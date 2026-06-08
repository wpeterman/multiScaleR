#' @title Estimate Memory Requirements for Multiscale Optimization
#'
#' @description
#' Estimate the memory footprint of a \code{\link{kernel_prep}} object and
#' provide a conservative recommendation for parallel worker counts. The helper
#' inspects the stored kernel inputs, optionally builds the same optimization
#' context used by \code{\link{multiScale_optim}}, detects the number of
#' physical CPU cores, and estimates a maximum \code{n_cores} value while
#' leaving a user-defined number of cores free.
#'
#' @param kernel_inputs A \code{"multiScaleR_data"} object created by
#'   \code{\link{kernel_prep}}.
#' @param fitted_mod Optional fitted model object. When supplied, the helper
#'   also measures the optimization context that would be serialized to
#'   parallel workers during \code{\link{multiScale_optim}}.
#' @param join_by Optional data frame passed through to
#'   \code{\link{multiScale_optim}} for \code{unmarked} workflows.
#' @param refit_fn Optional custom refit function passed through to
#'   \code{\link{multiScale_optim}}.
#' @param n_cores Optional positive integer. When supplied, the helper also
#'   reports the estimated peak memory for that specific worker count. When
#'   omitted, the reported peak estimate uses the recommended maximum worker
#'   count.
#' @param PSOCK Logical. Should parallel runs be budgeted as PSOCK workers.
#'   Default: \code{TRUE} on non-Unix platforms and \code{FALSE} on Unix-like
#'   platforms. Even on systems that support forking, the recommendation uses a
#'   conservative per-worker duplication assumption.
#' @param safety_factor Numeric scalar \eqn{\ge 1}. Multiplier applied to the
#'   per-worker payload to allow for temporary allocations during optimization.
#'   Default: \code{1.5}.
#' @param ram_fraction Numeric scalar in \code{(0, 1]}. Fraction of detected
#'   total system RAM considered usable for this run. The remainder is reserved
#'   as headroom for the OS and other R objects. Default: \code{0.75}.
#' @param leave_free Integer scalar \eqn{\ge 0}. Number of physical CPU cores
#'   to leave unused when computing the recommended maximum worker count.
#'   Default: \code{2}.
#'
#' @return A list of class \code{"multiScaleR_ram"} with:
#' \describe{
#'   \item{\code{requested_n_cores}}{Worker count used for the
#'     \code{peak_parallel_estimate}.}
#'   \item{\code{recommended_n_cores}}{Conservative recommended maximum worker
#'     count after applying both CPU and RAM limits.}
#'   \item{\code{max_cores_by_cpu}}{Maximum worker count allowed by detected
#'     physical CPU cores after leaving \code{leave_free} cores unused.}
#'   \item{\code{max_cores_by_ram}}{Maximum worker count allowed by the RAM
#'     budget, or \code{NA} if total system RAM could not be detected.}
#'   \item{\code{backend}}{Character string identifying the assumed parallel
#'     backend budget: \code{"PSOCK"} or \code{"fork"}.}
#'   \item{\code{system}}{Named list describing detected physical cores,
#'     logical cores, total RAM, usable RAM, and the reservation settings.}
#'   \item{\code{point_summary}}{One-row data frame summarizing the number of
#'     points, extracted cells per point, and source raster layer count.}
#'   \item{\code{component_bytes}}{Data frame reporting estimated byte sizes for
#'     the major stored objects and worker bundles.}
#'   \item{\code{notes}}{Character vector with interpretation notes and any
#'     detection warnings.}
#' }
#'
#' @details
#' The estimate is advisory rather than exact. It is based on
#' \code{\link[utils:object.size]{utils::object.size()}} and therefore reflects
#' the size of stored R objects, not all transient allocations made by
#' \code{exactextractr}, model refits, or BLAS/LAPACK.
#'
#' The recommendation is intentionally conservative. Even on Unix-like systems
#' where forking can share memory copy-on-write, worker counts are budgeted as
#' though each worker may require its own copy of the optimization payload.
#'
#' If \code{fitted_mod} is omitted, the RAM recommendation is based on
#' \code{kernel_inputs} alone and may understate the true optimization payload.
#'
#' @examples
#' library(terra)
#'
#' pts <- vect(cbind(c(3, 5, 7),
#'                   c(7, 5, 3)))
#' r <- rast(matrix(rnorm(100), nrow = 10))
#' names(r) <- "hab"
#'
#' kernel_inputs <- kernel_prep(
#'   pts = pts,
#'   raster_stack = r,
#'   max_D = 2,
#'   kernel = "gaussian",
#'   verbose = FALSE
#' )
#'
#' estimate_multiscale_ram(kernel_inputs)
#' @export
estimate_multiscale_ram <- function(kernel_inputs,
                                    fitted_mod = NULL,
                                    join_by = NULL,
                                    refit_fn = NULL,
                                    n_cores = NULL,
                                    PSOCK = (.Platform$OS.type != "unix"),
                                    safety_factor = 1.5,
                                    ram_fraction = 0.75,
                                    leave_free = 2L) {
  if (!inherits(kernel_inputs, "multiScaleR_data")) {
    stop("`kernel_inputs` must be a `multiScaleR_data` object from `kernel_prep()`.",
         call. = FALSE)
  }

  if (!is.null(join_by) && !is.data.frame(join_by)) {
    stop("`join_by` must be a data frame if provided.", call. = FALSE)
  }
  if (!is.null(refit_fn) && !is.function(refit_fn)) {
    stop("`refit_fn` must be a function if provided.", call. = FALSE)
  }
  validate_scalar_logical(PSOCK, "PSOCK")
  validate_scalar_numeric(safety_factor, "safety_factor", lower = 1)
  validate_scalar_numeric(ram_fraction, "ram_fraction", lower = 0, upper = 1)
  validate_scalar_numeric(leave_free, "leave_free", lower = 0, integerish = TRUE)
  if (!is.null(n_cores)) {
    validate_scalar_numeric(n_cores, "n_cores", lower = 1, integerish = TRUE)
  }

  d_list <- kernel_inputs$d_list
  raw_cov <- kernel_inputs$raw_cov
  kernel_dat <- kernel_inputs$kernel_dat

  point_summary <- data.frame(
    n_points = length(d_list),
    mean_cells_per_point = if (length(d_list) > 0) mean(lengths(d_list)) else NA_real_,
    median_cells_per_point = if (length(d_list) > 0) stats::median(lengths(d_list)) else NA_real_,
    max_cells_per_point = if (length(d_list) > 0) max(lengths(d_list)) else NA_real_,
    n_source_layers = if (length(raw_cov) > 0) ncol(raw_cov[[1]]) else 0L,
    stringsAsFactors = FALSE
  )

  opt_context <- NULL
  opt_context_error <- NULL
  if (!is.null(fitted_mod)) {
    opt_context <- tryCatch(
      build_opt_context(
        fitted_mod = fitted_mod,
        cov_df = raw_cov,
        join_by = join_by,
        refit_fn = refit_fn,
        scale_vars = kernel_inputs$scale_vars,
        unit_conv = kernel_inputs$unit_conv,
        resolution = kernel_inputs$resolution,
        n_cols = kernel_inputs$n_cols
      ),
      error = function(e) {
        opt_context_error <<- conditionMessage(e)
        NULL
      }
    )
  }

  component_sizes <- list(
    kernel_inputs_total = utils::object.size(kernel_inputs),
    kernel_dat = utils::object.size(kernel_dat),
    d_list = utils::object.size(d_list),
    raw_cov = utils::object.size(raw_cov),
    join_by = utils::object.size(join_by),
    fitted_mod = utils::object.size(fitted_mod),
    opt_context = utils::object.size(opt_context)
  )

  worker_bundle <- list(
    d_list = d_list,
    raw_cov = raw_cov,
    kernel = kernel_inputs$kernel,
    join_by = join_by,
    opt_context = opt_context,
    fitted_mod = fitted_mod
  )

  main_bundle <- list(
    kernel_inputs = kernel_inputs,
    join_by = join_by,
    opt_context = opt_context,
    fitted_mod = fitted_mod
  )

  worker_bundle_bytes <- .msr_as_numeric_size(utils::object.size(worker_bundle))
  main_bundle_bytes <- .msr_as_numeric_size(utils::object.size(main_bundle))

  physical_cores <- .msr_detect_cores(logical = FALSE)
  logical_cores <- .msr_detect_cores(logical = TRUE)
  total_ram_bytes <- .msr_total_ram_bytes()
  usable_ram_bytes <- if (is.finite(total_ram_bytes)) total_ram_bytes * ram_fraction else NA_real_

  max_cores_by_cpu <- .msr_max_cores_by_cpu(physical_cores, leave_free)
  max_cores_by_ram <- .msr_max_cores_by_ram(
    main_bundle_bytes = main_bundle_bytes,
    worker_bundle_bytes = worker_bundle_bytes,
    usable_ram_bytes = usable_ram_bytes,
    safety_factor = safety_factor
  )

  recommended_n_cores <- .msr_recommended_cores(
    max_cores_by_cpu = max_cores_by_cpu,
    max_cores_by_ram = max_cores_by_ram
  )

  requested_n_cores <- if (is.null(n_cores)) recommended_n_cores else as.integer(n_cores)
  peak_parallel_bytes <- .msr_peak_parallel_bytes(
    n_cores = requested_n_cores,
    main_bundle_bytes = main_bundle_bytes,
    worker_bundle_bytes = worker_bundle_bytes,
    safety_factor = safety_factor,
    PSOCK = PSOCK
  )

  component_bytes <- data.frame(
    component = c(
      "kernel_inputs_total",
      "kernel_dat",
      "d_list",
      "raw_cov",
      "join_by",
      "fitted_mod",
      "opt_context",
      "worker_bundle",
      "main_bundle",
      "peak_parallel_estimate"
    ),
    bytes = c(
      .msr_as_numeric_size(component_sizes$kernel_inputs_total),
      .msr_as_numeric_size(component_sizes$kernel_dat),
      .msr_as_numeric_size(component_sizes$d_list),
      .msr_as_numeric_size(component_sizes$raw_cov),
      .msr_as_numeric_size(component_sizes$join_by),
      .msr_as_numeric_size(component_sizes$fitted_mod),
      .msr_as_numeric_size(component_sizes$opt_context),
      worker_bundle_bytes,
      main_bundle_bytes,
      peak_parallel_bytes
    ),
    pretty = c(
      .msr_format_bytes(component_sizes$kernel_inputs_total),
      .msr_format_bytes(component_sizes$kernel_dat),
      .msr_format_bytes(component_sizes$d_list),
      .msr_format_bytes(component_sizes$raw_cov),
      .msr_format_bytes(component_sizes$join_by),
      .msr_format_bytes(component_sizes$fitted_mod),
      .msr_format_bytes(component_sizes$opt_context),
      .msr_format_bytes(worker_bundle_bytes),
      .msr_format_bytes(main_bundle_bytes),
      .msr_format_bytes(peak_parallel_bytes)
    ),
    stringsAsFactors = FALSE
  )

  notes <- c(
    if (is.null(fitted_mod)) {
      "No fitted model was supplied, so the recommendation is based on `kernel_inputs` alone."
    } else {
      NA_character_
    },
    if (!is.null(opt_context_error)) {
      paste0("Optimization context could not be built: ", opt_context_error)
    } else {
      NA_character_
    },
    if (!is.finite(total_ram_bytes)) {
      "Total system RAM could not be detected, so the CPU-based limit was used."
    } else {
      NA_character_
    },
    if (!is.finite(physical_cores)) {
      "Physical core detection failed; the helper fell back to at least one worker."
    } else {
      NA_character_
    },
    if (isTRUE(PSOCK)) {
      "Parallel RAM budgeting assumes PSOCK-style worker duplication."
    } else {
      "Parallel RAM budgeting remains conservative even when forked workers may share memory."
    }
  )
  notes <- notes[!is.na(notes)]

  out <- list(
    requested_n_cores = requested_n_cores,
    recommended_n_cores = recommended_n_cores,
    max_cores_by_cpu = max_cores_by_cpu,
    max_cores_by_ram = max_cores_by_ram,
    backend = if (isTRUE(PSOCK)) "PSOCK" else "fork",
    system = list(
      physical_cores = physical_cores,
      logical_cores = logical_cores,
      leave_free = as.integer(leave_free),
      total_ram_bytes = total_ram_bytes,
      total_ram_pretty = .msr_format_bytes(total_ram_bytes),
      usable_ram_fraction = ram_fraction,
      usable_ram_bytes = usable_ram_bytes,
      usable_ram_pretty = .msr_format_bytes(usable_ram_bytes),
      safety_factor = safety_factor
    ),
    point_summary = point_summary,
    component_bytes = component_bytes,
    notes = notes
  )

  class(out) <- "multiScaleR_ram"
  out
}


.msr_warn_unsafe_parallel_ram <- function(kernel_inputs,
                                          fitted_mod,
                                          join_by,
                                          refit_fn,
                                          n_cores,
                                          PSOCK) {
  if (is.null(n_cores) || as.integer(n_cores) <= 1L) {
    return(invisible(NULL))
  }
  if (!inherits(kernel_inputs, "multiScaleR_data")) {
    return(invisible(NULL))
  }

  ram <- tryCatch(
    estimate_multiscale_ram(
      kernel_inputs = kernel_inputs,
      fitted_mod = fitted_mod,
      join_by = join_by,
      refit_fn = refit_fn,
      n_cores = n_cores,
      PSOCK = PSOCK
    ),
    error = function(e) NULL
  )
  if (is.null(ram) || is.na(ram$max_cores_by_ram) ||
      as.integer(n_cores) <= ram$max_cores_by_ram) {
    return(invisible(NULL))
  }

  peak <- ram$component_bytes$pretty[
    ram$component_bytes$component == "peak_parallel_estimate"
  ]
  warning(
    "Requested `n_cores = ", as.integer(n_cores),
    "` may exceed the conservative RAM budget for this `kernel_prep()` payload. ",
    "`estimate_multiscale_ram()` recommends at most ",
    ram$max_cores_by_ram, " worker(s) by RAM for the ",
    ram$backend, " backend; estimated peak use is ", peak,
    ". Consider lowering `n_cores` or running `estimate_multiscale_ram()` first.",
    call. = FALSE
  )

  invisible(ram)
}


.msr_detect_cores <- function(logical = FALSE) {
  cores <- suppressWarnings(parallel::detectCores(logical = logical))
  if (length(cores) != 1 || is.na(cores) || !is.finite(cores) || cores < 1) {
    return(NA_real_)
  }

  as.integer(cores)
}


.msr_total_ram_bytes <- function() {
  sysname <- Sys.info()[["sysname"]]

  if (identical(.Platform$OS.type, "windows")) {
    out <- .msr_system_stdout(
      "powershell",
      c("-NoProfile", "-Command",
        "[int64](Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory")
    )
    ram <- .msr_parse_first_numeric(out)
    if (is.finite(ram)) {
      return(ram)
    }
  }

  if (identical(sysname, "Linux") && file.exists("/proc/meminfo")) {
    lines <- tryCatch(readLines("/proc/meminfo", n = 1L, warn = FALSE),
                      error = function(e) character())
    ram <- .msr_parse_first_numeric(lines)
    if (is.finite(ram)) {
      return(ram * 1024)
    }
  }

  if (identical(sysname, "Darwin")) {
    out <- .msr_system_stdout("sysctl", c("-n", "hw.memsize"))
    ram <- .msr_parse_first_numeric(out)
    if (is.finite(ram)) {
      return(ram)
    }
  }

  NA_real_
}


.msr_system_stdout <- function(command, args = character()) {
  tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = FALSE),
    warning = function(e) character(),
    error = function(e) character()
  )
}


.msr_parse_first_numeric <- function(x) {
  if (length(x) == 0) {
    return(NA_real_)
  }

  x <- trimws(x)
  x <- x[nzchar(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }

  matches <- regmatches(x, regexpr("[0-9]+(?:\\.[0-9]+)?", x, perl = TRUE))
  matches <- matches[lengths(matches) > 0]
  if (length(matches) == 0) {
    return(NA_real_)
  }

  out <- suppressWarnings(as.numeric(matches[[1]]))
  if (length(out) != 1 || is.na(out) || !is.finite(out)) {
    return(NA_real_)
  }

  out
}


.msr_max_cores_by_cpu <- function(physical_cores, leave_free) {
  if (!is.finite(physical_cores)) {
    return(1L)
  }

  max(1L, as.integer(physical_cores) - as.integer(leave_free))
}


.msr_max_cores_by_ram <- function(main_bundle_bytes,
                                  worker_bundle_bytes,
                                  usable_ram_bytes,
                                  safety_factor) {
  if (!is.finite(usable_ram_bytes)) {
    return(NA_integer_)
  }

  worker_bytes <- worker_bundle_bytes * safety_factor
  if (!is.finite(worker_bytes) || worker_bytes <= 0) {
    return(1L)
  }

  remaining <- usable_ram_bytes - main_bundle_bytes
  if (!is.finite(remaining) || remaining <= 0) {
    return(1L)
  }

  max(1L, as.integer(floor(remaining / worker_bytes)))
}


.msr_recommended_cores <- function(max_cores_by_cpu, max_cores_by_ram) {
  if (is.na(max_cores_by_ram)) {
    return(as.integer(max_cores_by_cpu))
  }

  as.integer(min(max_cores_by_cpu, max_cores_by_ram))
}


.msr_peak_parallel_bytes <- function(n_cores,
                                     main_bundle_bytes,
                                     worker_bundle_bytes,
                                     safety_factor,
                                     PSOCK) {
  n_cores <- max(1L, as.integer(n_cores))
  worker_bytes <- worker_bundle_bytes * safety_factor

  if (isTRUE(PSOCK)) {
    return(main_bundle_bytes + (n_cores * worker_bytes))
  }

  # Use the same linear worker budget on Unix-like systems to stay
  # conservative even when copy-on-write may reduce realized RAM use.
  main_bundle_bytes + (n_cores * worker_bytes)
}


.msr_as_numeric_size <- function(x) {
  if (is.null(x)) {
    return(0)
  }

  as.numeric(x)
}


.msr_format_bytes <- function(x) {
  if (length(x) != 1 || is.na(x) || !is.finite(x)) {
    return(NA_character_)
  }

  format(structure(x, class = "object_size"), units = "auto")
}


#' @export
#' @method print multiScaleR_ram
print.multiScaleR_ram <- function(x, ...) {
  cat("\nEstimated parallel RAM budget for `multiScaleR`\n")
  cat("Backend assumption: ", x$backend, "\n", sep = "")
  cat("Recommended max workers: ", x$recommended_n_cores, "\n", sep = "")
  cat("CPU limit after reserve: ", x$max_cores_by_cpu, "\n", sep = "")
  cat("RAM limit: ", x$max_cores_by_ram, "\n", sep = "")
  cat("Physical cores detected: ", x$system$physical_cores, "\n", sep = "")
  cat("Logical cores detected: ", x$system$logical_cores, "\n", sep = "")
  cat("Total RAM detected: ", x$system$total_ram_pretty, "\n", sep = "")
  cat("Usable RAM budget: ", x$system$usable_ram_pretty, "\n", sep = "")
  cat("Peak estimate at ", x$requested_n_cores, " worker(s): ",
      x$component_bytes$pretty[
        x$component_bytes$component == "peak_parallel_estimate"
      ],
      "\n",
      sep = "")

  if (length(x$notes) > 0) {
    cat("\nNotes:\n")
    cat(paste0("- ", x$notes, collapse = "\n"))
    cat("\n")
  }

  invisible(x)
}

#' @title Estimate Memory Requirements for Multiscale Optimization
#'
#' @description
#' Estimates how much memory a [kernel_prep()] object (and optionally a fitted
#' starting model) needs during [multiScale_optim()], and recommends a
#' conservative maximum number of parallel workers (`n_cores`). Use it before
#' running a large or slow model in parallel, because each worker holds its
#' own copy of the inputs.
#'
#' @param kernel_inputs A `multiScaleR_data` object created by
#'   [kernel_prep()] or [kernel_prep_by_group()].
#' @param fitted_mod Optional fitted starting model, as you would pass to
#'   [multiScale_optim()]. When supplied, the estimate also includes the
#'   optimization context sent to each worker. Without it, the estimate can
#'   understate the real payload. Default: `NULL`.
#' @param join_by Optional data frame passed through to [multiScale_optim()]
#'   for `unmarked` workflows. Default: `NULL`.
#' @param refit_fn Optional custom refit function passed through to
#'   [multiScale_optim()]. Default: `NULL`.
#' @param n_cores Optional positive integer. Worker count at which to report
#'   the peak memory estimate. When omitted, the recommended maximum is used.
#' @param PSOCK Logical. Budget workers as PSOCK (socket) workers. Default:
#'   `TRUE` on Windows and `FALSE` on Unix-like systems. The estimate is the
#'   same either way, because forked workers are also budgeted as though each
#'   holds its own copy.
#' @param safety_factor Numeric scalar, at least 1. Multiplier on each
#'   worker's payload to allow for temporary allocations during optimization.
#'   Larger values give fewer recommended workers. Default: `1.5`.
#' @param ram_fraction Numeric scalar between 0 and 1. Share of total system
#'   RAM treated as usable. The rest is left for the operating system and
#'   other R objects. Default: `0.75`.
#' @param leave_free Non-negative integer. Number of physical CPU cores to
#'   leave unused when computing the CPU limit. Default: `2`.
#' @param x A `multiScaleR_ram` object returned by `estimate_multiscale_ram()`.
#' @param ... Not used.
#'
#' @return `estimate_multiscale_ram()` returns a list of class
#' `multiScaleR_ram` with these elements:
#' * `requested_n_cores`: the worker count used for the peak estimate
#'   (`n_cores`, or the recommendation when `n_cores` is omitted).
#' * `recommended_n_cores`: the smaller of the CPU and RAM limits, or the CPU
#'   limit when RAM could not be detected. Never less than 1.
#' * `max_cores_by_cpu`: physical cores minus `leave_free`, at least 1.
#' * `max_cores_by_ram`: the number of workers that fit in the usable RAM, at
#'   least 1, or `NA` when total RAM could not be detected.
#' * `backend`: `"PSOCK"` or `"fork"`, the assumed parallel backend.
#' * `system`: list with the detected `physical_cores` and `logical_cores`,
#'   `leave_free`, `total_ram_bytes`, `usable_ram_fraction`,
#'   `usable_ram_bytes`, `safety_factor`, and readable versions of the RAM
#'   values (`total_ram_pretty`, `usable_ram_pretty`).
#' * `point_summary`: one-row data frame with the number of points, the
#'   mean, median, and maximum number of extracted cells per point (`NA` in
#'   lean mode), the number of source layers, whether cell-level data are
#'   stored (`cell_data_stored`), and the number of distance bins (`nbins`).
#' * `component_bytes`: data frame with columns `component`, `bytes`, and
#'   `pretty` (a readable size) for each stored object (`kernel_inputs_total`,
#'   `kernel_dat`, `d_list`, `raw_cov`, `binned`, `join_by`, `fitted_mod`,
#'   `opt_context`), the two bundles described in Details (`worker_bundle`
#'   and `main_bundle`), and `peak_parallel_estimate`.
#' * `notes`: character vector of interpretation notes and detection
#'   warnings.
#'
#' `print()` displays the main limits, the detected hardware, the peak
#' estimate, and the notes, and returns `x` invisibly.
#'
#' @details
#' # How the estimate is built
#'
#' Sizes come from [utils::object.size()], so they describe stored R objects,
#' not the transient memory used by `exactextractr`, model refits, or linear
#' algebra libraries. The *main bundle* is what the main R session holds (the
#' kernel inputs, `join_by`, the optimization context, and the model). The
#' *worker bundle* is what each parallel worker receives (the per-cell
#' distances and values, `join_by`, the optimization context, and the model).
#' The limits and the peak estimate follow from these sizes.
#'
#' * Usable RAM = total RAM x `ram_fraction`.
#' * RAM limit = (usable RAM - main bundle) / (worker bundle x
#'   `safety_factor`), rounded down.
#' * Peak estimate = main bundle + `n_cores` x worker bundle x
#'   `safety_factor`.
#'
#' Even where forking lets workers share memory, each worker is budgeted as
#' if it held its own copy. Hardware measurements that succeed are kept for
#' the rest of the R session. A failed RAM lookup is retried on the next call.
#'
#' # Read the report
#'
#' Compare the worker count you plan to use with `max_cores_by_ram`, and read
#' the peak estimate at that count. These are conservative estimates, not
#' measured memory use. The recommendation is never below one worker, but
#' that does not guarantee a serial run fits in memory. A run at or below the
#' limit can still run out of memory during a model refit. If RAM detection
#' or the estimate fails, no limit can be enforced.
#'
#' # Reduce the memory needed
#'
#' When every covariate is a kernel mean, rebuild the original
#' [kernel_prep()] call with `bin = TRUE, store_cell_data = FALSE`, keeping
#' its other settings, and check the estimate again. This lean mode keeps only
#' the distance-binned summaries that optimization uses. It is not available
#' for landscape or surface metrics, which need the cell-level data.
#' Otherwise, use fewer workers or `n_cores = 1`. Do not delete components of
#' the object by hand, and do not shrink a scientifically chosen `max_D` only
#' to pass the check.
#'
#' # How multiScale_optim() uses the estimate
#'
#' With `n_cores > 1`, [multiScale_optim()] calls this function with its
#' default settings before starting workers. It stops with an error when
#' `n_cores` exceeds `max_cores_by_ram`. It does not lower `n_cores` for you.
#' When RAM cannot be detected or the estimate fails, the check does not
#' block the run.
#'
#' @seealso [multiScale_optim()], [kernel_prep()]. The section "Estimate the
#'   scales of effect" in
#'   `vignette("v02-user-guide", package = "multiScaleR")` shows the check on
#'   a worked example.
#'
#' @examples
#' data("pts")
#' data("count_data")
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' kernel_inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                              kernel = "gaussian", verbose = FALSE)
#' mod <- glm(y ~ hab, family = poisson,
#'            data = data.frame(count_data["y"], kernel_inputs$kernel_dat))
#'
#' # Recommended worker limits and the peak memory estimate for two workers
#' ram <- estimate_multiscale_ram(kernel_inputs, fitted_mod = mod, n_cores = 2)
#' ram
#'
#' # Estimated size of each component, largest first
#' ram$component_bytes[order(-ram$component_bytes$bytes), ]
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
  binned <- kernel_inputs$binned
  kernel_dat <- kernel_inputs$kernel_dat

  # Cell-level data may have been dropped (`store_cell_data = FALSE`); fall back
  # to the binned summaries / kernel data to describe the point set.
  cells_per_point <- if (!is.null(d_list)) lengths(d_list) else integer(0)
  n_points <- if (!is.null(d_list)) {
    length(d_list)
  } else if (!is.null(binned)) {
    binned$n_points
  } else if (!is.null(kernel_dat)) {
    nrow(kernel_dat)
  } else {
    0L
  }
  n_source_layers <- if (!is.null(raw_cov) && length(raw_cov) > 0) {
    ncol(raw_cov[[1]])
  } else if (!is.null(binned)) {
    length(binned$sources)
  } else if (!is.null(kernel_dat)) {
    ncol(kernel_dat)
  } else {
    0L
  }

  point_summary <- data.frame(
    n_points = n_points,
    mean_cells_per_point = if (length(cells_per_point) > 0) mean(cells_per_point) else NA_real_,
    median_cells_per_point = if (length(cells_per_point) > 0) stats::median(cells_per_point) else NA_real_,
    max_cells_per_point = if (length(cells_per_point) > 0) max(cells_per_point) else NA_real_,
    n_source_layers = n_source_layers,
    cell_data_stored = !is.null(d_list),
    nbins = if (!is.null(binned)) binned$nbins else NA_integer_,
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
        n_cols = kernel_inputs$n_cols,
        binned = binned,
        d_list = d_list
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
    binned = utils::object.size(binned),
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
      "binned",
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
      .msr_as_numeric_size(component_sizes$binned),
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
      .msr_format_bytes(component_sizes$binned),
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


.msr_stop_unsafe_parallel_ram <- function(kernel_inputs,
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
  stop(
    "Parallel optimization was not started: requested `n_cores = ", as.integer(n_cores),
    "` exceeds the conservative RAM budget for this `kernel_prep()` payload.\n",
    "Estimated peak use: ", peak, " (", ram$backend,
    " backend); RAM-based worker limit: ", ram$max_cores_by_ram, ".\n",
    "Use `n_cores = 1` for serial execution, or rebuild `kernel_prep()` inputs ",
    "with `bin = TRUE, store_cell_data = FALSE` for lean mode when all ",
    "covariates are kernel-weighted (not landscape or surface metrics).\n",
    "Inspect `estimate_multiscale_ram(kernel_inputs, fitted_mod = ...)` ",
    "before retrying parallel execution.",
    call. = FALSE
  )
}


# Hardware (physical/logical cores, total RAM) does not change during an R
# session, but detecting it can be costly -- on Windows `.msr_total_ram_bytes()`
# spawns a PowerShell process. `multiScale_optim()` queries these on every call
# to populate its `next_run` recommendation, so the results are memoized to
# avoid repeated process spawns / system calls.
.msr_hw_cache <- new.env(parent = emptyenv())

.msr_detect_cores <- function(logical = FALSE) {
  key <- if (isTRUE(logical)) "logical_cores" else "physical_cores"
  cached <- .msr_hw_cache[[key]]
  if (!is.null(cached)) {
    return(cached)
  }

  cores <- suppressWarnings(parallel::detectCores(logical = logical))
  out <- if (length(cores) != 1 || is.na(cores) || !is.finite(cores) || cores < 1) {
    NA_real_
  } else {
    as.integer(cores)
  }
  .msr_hw_cache[[key]] <- out
  out
}


.msr_total_ram_bytes <- function() {
  cached <- .msr_hw_cache[["total_ram_bytes"]]
  if (is.numeric(cached) && length(cached) == 1L && is.finite(cached) && cached > 0) {
    return(cached)
  }
  out <- .msr_total_ram_bytes_uncached()
  # A transient discovery failure must not become a session-long hardware value.
  if (is.numeric(out) && length(out) == 1L && is.finite(out) && out > 0)
    .msr_hw_cache[["total_ram_bytes"]] <- out
  out
}


.msr_total_ram_bytes_uncached <- function() {
  sysname <- Sys.info()[["sysname"]]

  if (identical(.Platform$OS.type, "windows")) {
    ram <- .msr_windows_total_ram_bytes()
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


# The direct .NET query avoids depending on a responsive CIM provider. Both
# queries are bounded by .msr_system_stdout(); neither supplies a guessed total.
.msr_windows_total_ram_bytes <- function() {
  queries <- c(
    paste("Add-Type -AssemblyName Microsoft.VisualBasic;",
      "[Microsoft.VisualBasic.Devices.ComputerInfo]::new().TotalPhysicalMemory"),
    "[int64](Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory")
  for (query in queries) {
    out <- .msr_system_stdout("powershell", c("-NoProfile", "-Command", query))
    # Accept one positive numeric result, never digits embedded in diagnostics.
    out <- trimws(out)
    out <- out[grepl("^[0-9]+$", out)]
    if (length(out) != 1L) next
    ram <- .msr_parse_first_numeric(out)
    if (is.finite(ram) && ram > 0) return(ram)
  }
  NA_real_
}


.msr_system_stdout <- function(command, args = character()) {
  # Bound the external call (PowerShell on Windows can be slow to spawn under
  # load); on timeout/failure fall back to no output so RAM detection degrades
  # to NA rather than blocking. `timeout` is available since R 3.5.0.
  tryCatch(
    system2(command, args = args, stdout = TRUE, stderr = FALSE, timeout = 5),
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


#' @rdname estimate_multiscale_ram
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

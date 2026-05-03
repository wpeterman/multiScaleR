#' Define multiScaleR covariate transformations
#'
#' @description
#' These helpers define named model covariates as transformations of one or more
#' source raster layers. They are optional; if omitted, `kernel_prep()` preserves
#' the historical behavior where each raster layer becomes one optimized
#' kernel-weighted covariate with the same name as the raster layer.
#'
#' @param ... Named `kernel_var()` or `landscape_var()` specifications.
#' @param source Character. Name of the source raster layer.
#' @param metric Character. Landscape metric to calculate for `landscape_var()`.
#' @param radius Optional fixed buffer radius. If omitted, the radius is optimized.
#' @param base Logarithm base for diversity metrics.
#' @param classes_max Optional maximum number of classes for relative patch
#' richness (`rpr`).
#'
#' @return
#' `msr_vars()` returns a data frame of class `"multiScaleR_vars"` describing
#' the requested covariate transformations.
#'
#' @examples
#' vars <- msr_vars(
#'   forest_prop = kernel_var("forest"),
#'   forest_ed = landscape_var("forest", metric = "ed"),
#'   cover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500)
#' )
#'
#' @rdname msr_vars
#' @export
msr_vars <- function(...) {
  specs <- list(...)
  if (length(specs) == 0) {
    stop("At least one variable specification must be provided.", call. = FALSE)
  }

  spec_names <- names(specs)
  if (is.null(spec_names) || any(!nzchar(spec_names))) {
    stop("All `msr_vars()` specifications must be named.", call. = FALSE)
  }
  if (anyDuplicated(spec_names)) {
    stop("All `msr_vars()` names must be unique.", call. = FALSE)
  }

  for (i in seq_along(specs)) {
    if (!inherits(specs[[i]], "multiScaleR_var")) {
      stop(
        "All `msr_vars()` inputs must be created with `kernel_var()` or `landscape_var()`.",
        call. = FALSE
      )
    }
  }

  out <- data.frame(
    covariate = spec_names,
    source = vapply(specs, `[[`, character(1), "source"),
    type = vapply(specs, `[[`, character(1), "type"),
    metric = vapply(specs, function(x) .msr_chr_or_na(x$metric), character(1)),
    radius = vapply(specs, function(x) .msr_num_or_na(x$radius), numeric(1)),
    optimize = vapply(specs, `[[`, logical(1), "optimize"),
    base = vapply(specs, `[[`, numeric(1), "base"),
    classes_max = vapply(specs, function(x) .msr_num_or_na(x$classes_max), numeric(1)),
    stringsAsFactors = FALSE
  )

  class(out) <- c("multiScaleR_vars", "data.frame")
  out
}


#' @rdname msr_vars
#' @export
kernel_var <- function(source) {
  validate_character_scalar(source, "source")

  structure(
    list(
      type = "kernel",
      source = source,
      metric = NA_character_,
      radius = NA_real_,
      optimize = TRUE,
      base = exp(1),
      classes_max = NA_real_
    ),
    class = "multiScaleR_var"
  )
}


#' @rdname msr_vars
#' @export
landscape_var <- function(source,
                          metric,
                          radius = NULL,
                          base = exp(1),
                          classes_max = NULL) {
  validate_character_scalar(source, "source")
  metric <- match.arg(tolower(metric), .msr_landscape_metrics())
  validate_scalar_numeric(base, "base", positive = TRUE)
  if (base == 1) {
    stop("`base` must not equal 1.", call. = FALSE)
  }

  if (!is.null(radius)) {
    validate_scalar_numeric(radius, "radius", positive = TRUE)
  }
  if (!is.null(classes_max)) {
    validate_scalar_numeric(classes_max, "classes_max", positive = TRUE)
  }

  structure(
    list(
      type = "landscape",
      source = source,
      metric = metric,
      radius = if (is.null(radius)) NA_real_ else radius,
      optimize = is.null(radius),
      base = base,
      classes_max = if (is.null(classes_max)) NA_real_ else classes_max
    ),
    class = "multiScaleR_var"
  )
}


#' @export
print.multiScaleR_vars <- function(x, ...) {
  cat("multiScaleR variable specifications:\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}


.msr_chr_or_na <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    NA_character_
  } else {
    as.character(x[[1]])
  }
}


.msr_num_or_na <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    NA_real_
  } else {
    as.numeric(x[[1]])
  }
}


.msr_landscape_metrics <- function() {
  c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr", "prd", "rpr",
    "ta", "ed", "te", "lsi", "pladj", "contag")
}


.msr_composition_metrics <- function() {
  c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr", "prd", "rpr",
    "ta")
}


.msr_edge_metrics <- function() {
  c("ed", "te", "lsi")
}


.msr_adjacency_metrics <- function() {
  c("pladj", "contag")
}


.msr_default_scale_vars <- function(raster_stack) {
  vars <- data.frame(
    covariate = names(raster_stack),
    source = names(raster_stack),
    type = "kernel",
    metric = NA_character_,
    radius = NA_real_,
    optimize = TRUE,
    base = exp(1),
    classes_max = NA_real_,
    stringsAsFactors = FALSE
  )
  class(vars) <- c("multiScaleR_vars", "data.frame")
  vars
}


.msr_validate_scale_vars <- function(scale_vars, raster_stack, kernel = "gaussian") {
  if (is.null(scale_vars)) {
    scale_vars <- .msr_default_scale_vars(raster_stack)
  }

  if (!inherits(scale_vars, "multiScaleR_vars")) {
    stop("`scale_vars` must be created with `msr_vars()`.", call. = FALSE)
  }

  scale_vars <- as.data.frame(scale_vars, stringsAsFactors = FALSE)
  missing_sources <- setdiff(scale_vars$source, names(raster_stack))
  if (length(missing_sources) > 0) {
    stop(
      paste0(
        "The following `scale_vars` sources or optimized covariate layers are not raster layers: ",
        paste(missing_sources, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(scale_vars$covariate)) {
    stop("`scale_vars` covariate names must be unique.", call. = FALSE)
  }

  if (identical(kernel, "expow") &&
      any(scale_vars$type == "landscape" & scale_vars$optimize)) {
    stop(
      "Optimized `landscape_var()` specifications are not currently supported with `kernel = 'expow'`.",
      call. = FALSE
    )
  }

  class(scale_vars) <- c("multiScaleR_vars", "data.frame")
  scale_vars
}


.msr_optimized_scale_vars <- function(scale_vars) {
  if (is.null(scale_vars)) {
    return(NULL)
  }
  scale_vars[scale_vars$optimize, , drop = FALSE]
}


.msr_optimized_covariates <- function(scale_vars, cov_df = NULL) {
  if (is.null(scale_vars)) {
    return(colnames(cov_df[[1]]))
  }
  .msr_optimized_scale_vars(scale_vars)$covariate
}


.msr_needs_cells <- function(scale_vars) {
  any(scale_vars$type == "landscape" &
        scale_vars$metric %in% c(.msr_edge_metrics(), .msr_adjacency_metrics()))
}


.msr_extract_cells <- function(cov_df) {
  cells <- attr(cov_df, "cell", exact = TRUE)
  if (is.null(cells)) {
    stop(
      "Raster cell IDs are required for this landscape metric but were not cached.",
      call. = FALSE
    )
  }
  cells
}


.msr_eval_scale_vars <- function(d,
                                 cov_df,
                                 scale_vars,
                                 sigma,
                                 shape,
                                 kernel,
                                 unit_conv,
                                 resolution,
                                 n_cols,
                                 covariates = scale_vars$covariate) {
  out <- numeric(length(covariates))
  names(out) <- covariates

  opt_vars <- .msr_optimized_scale_vars(scale_vars)
  param_covariates <- covariates[covariates %in% opt_vars$covariate]

  for (j in seq_along(covariates)) {
    covariate <- covariates[[j]]
    spec <- scale_vars[scale_vars$covariate == covariate, , drop = FALSE]
    if (nrow(spec) != 1) {
      stop("Could not find a unique `scale_vars` specification for `", covariate, "`.",
           call. = FALSE)
    }

    values <- cov_df[, spec$source, drop = FALSE]

    if (identical(spec$type, "kernel")) {
      param_idx <- match(covariate, param_covariates)
      out[[j]] <- scale_type(
        d = d,
        kernel = kernel,
        sigma = sigma[[param_idx]],
        shape = if (is.null(shape)) NULL else shape[[param_idx]],
        r_stack.df = values
      )[[1]]
      next
    }

    radius <- if (isTRUE(spec$optimize)) {
      sigma[[match(covariate, param_covariates)]]
    } else {
      spec$radius / unit_conv
    }

    metric <- spec$metric
    if (metric %in% .msr_composition_metrics()) {
      out[[j]] <- .landscape_composition_by_buffer(
        d = d,
        r_stack.df = values,
        radius = radius,
        metric = metric,
        base = spec$base,
        resolution = resolution,
        classes_max = if (is.na(spec$classes_max)) NULL else spec$classes_max
      )[[1]]
    } else if (metric %in% .msr_edge_metrics()) {
      out[[j]] <- .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = .msr_extract_cells(cov_df),
        radius = radius,
        resolution = resolution,
        n_cols = n_cols,
        metric = metric
      )[[1]]
    } else if (metric %in% .msr_adjacency_metrics()) {
      out[[j]] <- .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = .msr_extract_cells(cov_df),
        radius = radius,
        n_cols = n_cols,
        metric = metric
      )[[1]]
    } else {
      stop("Unsupported landscape metric '", metric, "'.", call. = FALSE)
    }
  }

  out
}


.msr_merge_scl_params <- function(primary, fallback) {
  if (is.null(fallback)) {
    return(primary)
  }
  if (is.null(primary)) {
    return(fallback)
  }

  out <- fallback
  out$mean[names(primary$mean)] <- primary$mean
  out$sd[names(primary$sd)] <- primary$sd
  out
}


.msr_scale_vars_raster <- function(raster_stack,
                                   scale_vars,
                                   sigma,
                                   shape,
                                   kernel,
                                   pct_wt,
                                   fft,
                                   na.rm,
                                   verbose) {
  opt_vars <- .msr_optimized_scale_vars(scale_vars)
  out <- vector("list", nrow(scale_vars))
  names(out) <- scale_vars$covariate

  for (i in seq_len(nrow(scale_vars))) {
    spec <- scale_vars[i, , drop = FALSE]
    source_raster <- raster_stack[[spec$source]]

    if (identical(spec$type, "kernel")) {
      param_idx <- match(spec$covariate, opt_vars$covariate)
      if (isTRUE(verbose)) {
        cat(paste0(
          "\nSmoothing spatRaster variable '", spec$covariate,
          "' from layer '", spec$source,
          "' at sigma = ", floor(sigma[[param_idx]]), "\n"
        ))
      }
      out[[i]] <- .msr_kernel_raster_one(
        raster = source_raster,
        sigma = sigma[[param_idx]],
        shape = if (is.null(shape)) NULL else shape[[param_idx]],
        kernel = kernel,
        pct_wt = pct_wt,
        fft = fft,
        na.rm = na.rm
      )
      names(out[[i]]) <- spec$covariate
      next
    }

    radius <- if (isTRUE(spec$optimize)) {
      sigma[[match(spec$covariate, opt_vars$covariate)]]
    } else {
      spec$radius
    }

    if (isFALSE(fft)) {
      stop("Landscape metric raster projection currently requires `fft = TRUE`.",
           call. = FALSE)
    }
    if (isTRUE(verbose)) {
      cat(paste0(
        "\nCalculating landscape metric '", spec$metric,
        "' for variable '", spec$covariate,
        "' at radius = ", floor(radius), "\n"
      ))
    }

    if (spec$metric %in% .msr_composition_metrics()) {
      out[[i]] <- .landscape_composition_raster_fft(
        raster = source_raster,
        radius = radius,
        metric = spec$metric,
        base = spec$base,
        classes_max = if (is.na(spec$classes_max)) NULL else spec$classes_max,
        na.rm = na.rm
      )
    } else if (spec$metric %in% .msr_edge_metrics()) {
      out[[i]] <- .landscape_edge_raster_fft(
        raster = source_raster,
        radius = radius,
        metric = spec$metric,
        na.rm = na.rm
      )
    } else {
      out[[i]] <- .landscape_adjacency_raster_fft(
        raster = source_raster,
        radius = radius,
        metric = spec$metric,
        na.rm = na.rm
      )
    }
    names(out[[i]]) <- spec$covariate
  }

  terra::rast(out)
}


.msr_kernel_raster_one <- function(raster,
                                   sigma,
                                   shape,
                                   kernel,
                                   pct_wt,
                                   fft,
                                   na.rm) {
  mx <- kernel_dist(kernel = kernel,
                    sigma = sigma,
                    beta = shape,
                    prob = pct_wt)

  r_res <- terra::res(raster)[1]
  focal_d <- ceiling(mx / r_res) * 2
  if ((focal_d %% 2) == 0) {
    focal_d <- focal_d + 1
  }

  mat <- matrix(stats::rnorm(focal_d^2), focal_d, focal_d)
  r_wt <- terra::rast(mat)
  terra::crs(r_wt) <- terra::crs(raster)
  cntr_crd <- terra::xyFromCell(r_wt, ceiling(length(mat) / 2))
  cell_crds <- terra::crds(r_wt)
  r_wt[] <- fields::rdist(cntr_crd, cell_crds)[1, ] * r_res
  r_wt[] <- scale_type_r(d = as.vector(r_wt),
                         kernel = kernel,
                         sigma = sigma,
                         shape = shape,
                         output = "wts")

  wt_mat <- terra::as.matrix(r_wt, wide = TRUE)
  if (isTRUE(fft)) {
    mat <- terra::as.matrix(raster, wide = TRUE)
    out_mat <- fft_convolution(mat,
                               wt_mat,
                               fun = "mean",
                               na.rm = na.rm)
    out <- terra::rast(raster)
    terra::values(out) <- as.vector(out_mat)
    out
  } else {
    terra::focal(raster,
                 w = wt_mat,
                 fun = "mean",
                 na.rm = na.rm,
                 expand = FALSE)
  }
}

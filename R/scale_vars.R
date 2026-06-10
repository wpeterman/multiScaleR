#' Define multiScaleR covariate transformations
#'
#' @description
#' These helpers define named model covariates as transformations of one or
#' more source raster layers. They are optional; if omitted,
#' \code{\link{kernel_prep}} preserves the default behavior where each raster
#' layer becomes one optimized kernel-weighted covariate with the same name.
#'
#' Use \code{msr_vars()} to collect one or more \code{kernel_var()},
#' \code{landscape_var()}, or \code{surface_var()} specifications into a single
#' object that is passed to the \code{scale_vars} argument of
#' \code{\link{kernel_prep}} and \code{\link{kernel_scale.raster}}.
#'
#' @param ... Named \code{kernel_var()}, \code{landscape_var()}, or
#'   \code{surface_var()} specifications. Each argument must be named; the name
#'   becomes the covariate column name in the model data frame. All names must
#'   be unique. At least one specification is required.
#' @param source Character scalar. Name of the source raster layer in
#'   \code{raster_stack} from which the covariate is derived. Must exactly match
#'   a layer name in the raster stack provided to \code{\link{kernel_prep}}.
#' @param metric Character scalar. The metric to compute within the circular
#'   buffer. For \code{landscape_var()} this is a landscape ecology metric
#'   computed on a categorical raster; for \code{surface_var()} this is a
#'   continuous surface texture metric (\code{"sa"} or \code{"sq"}) computed on
#'   a continuous raster. Must be one of the supported metrics (see Details).
#'   Matched case-insensitively.
#' @param radius Optional positive numeric. Fixed buffer radius in the same
#'   units as the projection. When \code{NULL} (default), the radius is treated
#'   as a free parameter and optimized alongside the model. When supplied, the
#'   metric is computed at this fixed radius and no scale optimization is
#'   performed for this covariate.
#' @param weighted Logical, for \code{surface_var()} only. When \code{FALSE}
#'   (default), the surface metric is computed over a hard-edged circular
#'   neighborhood and the optimized scale parameter is that neighborhood radius.
#'   When \code{TRUE}, each cell is weighted by the model's distance kernel
#'   (the \code{kernel} passed to \code{\link{kernel_prep}}) and the optimized
#'   parameter is the kernel scale (sigma), exactly as for a kernel-weighted
#'   mean. The weighted form estimates the scale of effect of heterogeneity with
#'   a likelihood that is smooth in sigma. It is available only for the
#'   distribution metrics (\code{"sa"}, \code{"sq"}, \code{"ssk"}, \code{"sku"}),
#'   not the slope metric (\code{"sdq"}), and the scale is always optimized, so
#'   \code{radius} must be \code{NULL}.
#' @param base Positive numeric (not equal to 1). Logarithm base used when
#'   computing diversity metrics (\code{"shdi"}, \code{"shei"}, \code{"msidi"},
#'   \code{"msiei"}). Default: \code{exp(1)} (natural log). Use \code{2} for
#'   bits or \code{10} for Hartley units.
#' @param classes_max Optional positive numeric. Maximum number of possible
#'   patch types in the landscape, used only for relative patch richness
#'   (\code{"rpr"}). If \code{NULL} (default), the observed number of classes
#'   within each buffer is used as the denominator.
#'
#' @return
#' \describe{
#'   \item{\code{msr_vars()}}{A data frame of class \code{"multiScaleR_vars"}
#'     with one row per covariate. Columns are \code{covariate} (the name
#'     assigned in \code{...}), \code{source} (source raster layer name),
#'     \code{type} (\code{"kernel"}, \code{"landscape"}, or \code{"surface"}),
#'     \code{metric} (landscape or surface metric, or \code{NA}), \code{radius}
#'     (fixed radius or \code{NA} when optimized), \code{optimize} (logical
#'     indicating whether the scale is estimated), \code{base}, and
#'     \code{classes_max}.}
#'   \item{\code{kernel_var()}}{An internal list of class
#'     \code{"multiScaleR_var"} representing a kernel-weighted mean covariate.
#'     Pass one or more of these inside \code{msr_vars()}.}
#'   \item{\code{landscape_var()}}{An internal list of class
#'     \code{"multiScaleR_var"} representing a landscape metric covariate.
#'     Pass one or more of these inside \code{msr_vars()}.}
#'   \item{\code{surface_var()}}{An internal list of class
#'     \code{"multiScaleR_var"} representing a continuous surface texture metric
#'     covariate. Pass one or more of these inside \code{msr_vars()}.}
#' }
#'
#' @details
#' \strong{Kernel vs. landscape covariates}
#'
#' \code{kernel_var(source)} defines a kernel-weighted mean of the continuous
#' raster values within a circular neighborhood. The neighborhood radius (sigma)
#' is optimized by \code{\link{multiScale_optim}}.
#'
#' \code{landscape_var(source, metric)} derives a landscape ecology metric from
#' a categorical (or thresholded) raster layer within a circular buffer. Raster
#' classes must be encoded as finite integer-like values. The buffer radius can
#' be fixed or optimized.
#'
#' \code{surface_var(source, metric)} derives a continuous surface texture
#' metric from a continuous raster layer (elevation, NDVI, canopy height,
#' temperature, moisture, and similar surfaces) within a circular buffer. These
#' metrics summarize within-neighborhood heterogeneity rather than an average
#' value, so they complement \code{kernel_var()} (which summarizes the mean).
#' The buffer radius can be fixed or optimized. Note that these covariates
#' measure variability, not central tendency: changing the radius changes both
#' the ecological scale and the statistical quantity being summarized.
#'
#' By default a surface metric uses a hard-edged neighborhood. Setting
#' \code{weighted = TRUE} instead weights each cell by the model's distance
#' kernel and optimizes the kernel scale (sigma), which estimates the scale of
#' effect of heterogeneity with a smooth likelihood (see the \code{weighted}
#' argument). This applies to the distribution metrics (\code{"sa"},
#' \code{"sq"}, \code{"ssk"}, \code{"sku"}).
#'
#' \strong{Supported landscape metrics}
#'
#' Composition metrics (require a categorical raster):
#' \describe{
#'   \item{\code{"shdi"}}{Shannon diversity index.}
#'   \item{\code{"shei"}}{Shannon evenness index.}
#'   \item{\code{"sidi"}}{Simpson diversity index.}
#'   \item{\code{"siei"}}{Simpson evenness index.}
#'   \item{\code{"msidi"}}{Modified Simpson diversity index.}
#'   \item{\code{"msiei"}}{Modified Simpson evenness index.}
#'   \item{\code{"pr"}}{Patch richness (number of distinct classes).}
#'   \item{\code{"prd"}}{Patch richness density (pr per 100 ha).}
#'   \item{\code{"rpr"}}{Relative patch richness (pr / classes_max * 100).}
#'   \item{\code{"ta"}}{Total area of the buffer (ha).}
#' }
#'
#' Edge metrics (require adjacency information; cell IDs are cached internally):
#' \describe{
#'   \item{\code{"ed"}}{Edge density (m/ha).}
#'   \item{\code{"te"}}{Total edge length (m).}
#'   \item{\code{"lsi"}}{Landscape shape index.}
#' }
#'
#' Adjacency metrics:
#' \describe{
#'   \item{\code{"ai"}}{Aggregation index (\%).}
#'   \item{\code{"pladj"}}{Proportion of like adjacencies (\%).}
#'   \item{\code{"contag"}}{Contagion index (\%).}
#' }
#' Adjacency projections are class-count limited to avoid accidental use of
#' continuous rasters and to keep FFT projection times bounded. The current
#' ceilings are 200 classes for \code{"ai"} and \code{"pladj"}, and 50 classes
#' for \code{"contag"}, which scales with the square of the class count.
#'
#' \strong{Supported surface texture metrics}
#'
#' Surface metrics (require a continuous raster):
#' \describe{
#'   \item{\code{"sa"}}{Average roughness: the mean absolute deviation of the
#'     neighborhood values from their mean.}
#'   \item{\code{"sq"}}{Root mean square (RMS) roughness: the sample standard
#'     deviation of the neighborhood values. More sensitive to large deviations
#'     than \code{"sa"}.}
#'   \item{\code{"ssk"}}{Skewness: the asymmetry of the neighborhood value
#'     distribution. Positive when occasional high values dominate, negative
#'     when occasional low values do.}
#'   \item{\code{"sku"}}{Kurtosis (excess): the peakedness and tail weight of
#'     the neighborhood value distribution relative to a normal distribution
#'     (which has excess kurtosis 0).}
#'   \item{\code{"sdq"}}{Root mean square slope: the RMS of the local gradient
#'     magnitude, a measure of how steeply the surface changes (terrain or
#'     structural ruggedness).}
#' }
#' The dispersion and shape metrics (\code{"sa"}, \code{"sq"}, \code{"ssk"},
#' \code{"sku"}) match the \code{geodiv} package (\code{sq} uses an N - 1
#' denominator as \code{stats::sd()} does; \code{ssk} uses the bias-adjusted
#' skewness and \code{sku} the excess kurtosis, matching \code{geodiv}'s
#' defaults). \code{"sdq"} is reported as a true slope (the gradient is divided
#' by the cell resolution); \code{geodiv::sdq()} instead uses a cell spacing of
#' one, so the two agree after multiplying by the resolution.
#'
#' Raster projection of \code{"sq"}, \code{"ssk"}, \code{"sku"}, and \code{"sdq"}
#' uses fast Fourier transform (FFT) convolution; projection of \code{"sa"} uses
#' an exact compiled masked focal pass, because average roughness has no
#' closed-form FFT decomposition. \code{"sdq"} uses each cell's gradient, so
#' (like the landscape edge metrics) it caches raster cell IDs internally, and
#' its projection agrees with point extraction to within a small boundary
#' tolerance rather than bit-for-bit.
#'
#' \strong{Mixing covariate types}
#'
#' Multiple covariate types can be combined in one \code{msr_vars()} call. For
#' example, a kernel-weighted mean of forest cover and a fixed-radius edge
#' density metric can both be defined and passed together to
#' \code{\link{kernel_prep}}:
#'
#' \preformatted{
#' vars <- msr_vars(
#'   forest_mean  = kernel_var("forest"),
#'   forest_ed500 = landscape_var("forest", metric = "ed", radius = 500)
#' )
#' kernel_inputs <- kernel_prep(pts, rasters, max_D = 1000,
#'                              scale_vars = vars)
#' }
#'
#' Covariates with a fixed \code{radius} are computed once and not re-evaluated
#' during optimization, which can meaningfully reduce computation time.
#'
#' @seealso \code{\link{kernel_prep}}, \code{\link{multiScale_optim}},
#'   \code{\link{kernel_scale.raster}}
#'
#' @examples
#' ## Kernel-weighted mean only (equivalent to default behavior)
#' vars <- msr_vars(
#'   forest_prop = kernel_var("forest")
#' )
#'
#' ## Combining kernel, landscape, and surface covariates
#' vars <- msr_vars(
#'   forest_prop = kernel_var("forest"),
#'   forest_ed   = landscape_var("forest", metric = "ed"),
#'   cover_shdi_500 = landscape_var("landcover", metric = "shdi", radius = 500),
#'   elev_sd     = surface_var("elevation", metric = "sq")
#' )
#'
#' ## landscape_var with natural-log Shannon diversity (the default)
#' landscape_var("landcover", metric = "shdi")
#'
#' ## landscape_var with log2 Shannon diversity and a fixed 250 m radius
#' landscape_var("landcover", metric = "shdi", radius = 250, base = 2)
#'
#' ## surface_var: optimized RMS roughness, and fixed-radius average roughness
#' surface_var("elevation", metric = "sq")
#' surface_var("ndvi", metric = "sa", radius = 300)
#'
#' ## Kernel-weighted roughness: optimize the scale of effect of heterogeneity
#' surface_var("elevation", metric = "sq", weighted = TRUE)
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
    weighted = vapply(specs, function(x) isTRUE(x$weighted), logical(1)),
    base = vapply(specs, `[[`, numeric(1), "base"),
    classes_max = vapply(specs, function(x) .msr_num_or_na(x$classes_max), numeric(1)),
    stringsAsFactors = FALSE
  )

  class(out) <- c("multiScaleR_vars", "data.frame")
  out
}


#' @rdname msr_vars
#' @export
#' @examples
#' ## Define a single optimized kernel-weighted mean covariate
#' kv <- kernel_var("forest")
kernel_var <- function(source) {
  validate_character_scalar(source, "source")

  structure(
    list(
      type = "kernel",
      source = source,
      metric = NA_character_,
      radius = NA_real_,
      optimize = TRUE,
      weighted = FALSE,
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
      weighted = FALSE,
      base = base,
      classes_max = if (is.null(classes_max)) NA_real_ else classes_max
    ),
    class = "multiScaleR_var"
  )
}


#' @rdname msr_vars
#' @export
surface_var <- function(source,
                        metric,
                        radius = NULL,
                        weighted = FALSE) {
  validate_character_scalar(source, "source")
  metric <- match.arg(tolower(metric), .msr_surface_metrics())
  validate_scalar_logical(weighted, "weighted")

  if (isTRUE(weighted)) {
    if (metric %in% .msr_surface_neighbor_metrics()) {
      stop(
        sprintf(
          "`weighted = TRUE` is not supported for the slope metric `%s`; kernel weighting applies to the distribution metrics (%s).",
          metric,
          paste(setdiff(.msr_surface_metrics(), .msr_surface_neighbor_metrics()),
                collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(radius)) {
      stop(
        "`weighted = TRUE` surface metrics optimize the kernel scale, so `radius` must be NULL.",
        call. = FALSE
      )
    }
  }

  if (!is.null(radius)) {
    validate_scalar_numeric(radius, "radius", positive = TRUE)
  }

  structure(
    list(
      type = "surface",
      source = source,
      metric = metric,
      radius = if (is.null(radius)) NA_real_ else radius,
      optimize = is.null(radius),
      weighted = weighted,
      base = exp(1),
      classes_max = NA_real_
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
    "ta", "ed", "te", "lsi", "ai", "pladj", "contag")
}


.msr_composition_metrics <- function() {
  c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr", "prd", "rpr",
    "ta")
}


.msr_surface_metrics <- function() {
  c("sa", "sq", "ssk", "sku", "sdq")
}


# Surface metrics that require per-cell neighbor information (cell IDs), as
# opposed to the pointwise distribution metrics. Currently only `sdq` (RMS
# slope), which uses forward differences to each cell's neighbors.
.msr_surface_neighbor_metrics <- function() {
  c("sdq")
}


.msr_edge_metrics <- function() {
  c("ed", "te", "lsi")
}


.msr_adjacency_metrics <- function() {
  c("ai", "pladj", "contag")
}


.msr_default_scale_vars <- function(raster_stack) {
  vars <- data.frame(
    covariate = names(raster_stack),
    source = names(raster_stack),
    type = "kernel",
    metric = NA_character_,
    radius = NA_real_,
    optimize = TRUE,
    weighted = FALSE,
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

  # Optimized hard-neighborhood metrics (landscape metrics and unweighted
  # surface metrics) have no kernel, so the exponential-power shape parameter is
  # undefined for them. Weighted surface metrics do use the kernel (including
  # its shape), so they are permitted with `kernel = "expow"`.
  if (identical(kernel, "expow") &&
      any((scale_vars$type == "landscape" |
             (scale_vars$type == "surface" & !scale_vars$weighted)) &
            scale_vars$optimize)) {
    stop(
      "Optimized `landscape_var()` and unweighted `surface_var()` specifications are not currently supported with `kernel = 'expow'`.",
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


.msr_scale_vars_in_model <- function(scale_vars, multiScaleR) {
  if (is.null(scale_vars) ||
      is.null(multiScaleR) ||
      !inherits(multiScaleR, "multiScaleR")) {
    return(scale_vars)
  }

  used_covariates <- tryCatch(
    intersect(scale_vars$covariate,
              .model_predictors(.analysis_model(multiScaleR$opt_mod))),
    error = function(e) character(0)
  )

  if (length(used_covariates) == 0 && !is.null(multiScaleR$scale_est)) {
    used_covariates <- intersect(scale_vars$covariate,
                                 row.names(multiScaleR$scale_est))
  }

  if (length(used_covariates) == 0) {
    return(scale_vars)
  }

  scale_vars[scale_vars$covariate %in% used_covariates, , drop = FALSE]
}


.msr_optimized_covariates <- function(scale_vars, cov_df = NULL) {
  if (is.null(scale_vars)) {
    return(colnames(cov_df[[1]]))
  }
  .msr_optimized_scale_vars(scale_vars)$covariate
}


.msr_needs_cells <- function(scale_vars) {
  any(scale_vars$type == "landscape" &
        scale_vars$metric %in% c(.msr_edge_metrics(), .msr_adjacency_metrics())) ||
    any(scale_vars$type == "surface" &
          scale_vars$metric %in% .msr_surface_neighbor_metrics())
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
                                 covariates = scale_vars$covariate,
                                 validate = TRUE) {
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
    if (identical(spec$type, "surface")) {
      if (isTRUE(spec$weighted)) {
        param_idx <- match(covariate, param_covariates)
        out[[j]] <- .surface_weighted_by_buffer(
          d = d,
          r_stack.df = values,
          kernel = kernel,
          sigma = sigma[[param_idx]],
          shape = if (is.null(shape)) NULL else shape[[param_idx]],
          metric = metric,
          validate = validate
        )[[1]]
      } else if (metric %in% .msr_surface_neighbor_metrics()) {
        out[[j]] <- .surface_slope_by_buffer(
          d = d,
          r_stack.df = values,
          cells = .msr_extract_cells(cov_df),
          radius = radius,
          resolution = resolution,
          n_cols = n_cols,
          metric = metric,
          validate = validate
        )[[1]]
      } else {
        out[[j]] <- .surface_metric_by_buffer(
          d = d,
          r_stack.df = values,
          radius = radius,
          metric = metric,
          validate = validate
        )[[1]]
      }
    } else if (metric %in% .msr_composition_metrics()) {
      out[[j]] <- .landscape_composition_by_buffer(
        d = d,
        r_stack.df = values,
        radius = radius,
        metric = metric,
        base = spec$base,
        resolution = resolution,
        classes_max = if (is.na(spec$classes_max)) NULL else spec$classes_max,
        validate = validate
      )[[1]]
    } else if (metric %in% .msr_edge_metrics()) {
      out[[j]] <- .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = .msr_extract_cells(cov_df),
        radius = radius,
        resolution = resolution,
        n_cols = n_cols,
        metric = metric,
        validate = validate
      )[[1]]
    } else if (metric %in% .msr_adjacency_metrics()) {
      out[[j]] <- .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = .msr_extract_cells(cov_df),
        radius = radius,
        n_cols = n_cols,
        metric = metric,
        validate = validate
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
      stop("Landscape and surface metric raster projection currently requires `fft = TRUE`.",
           call. = FALSE)
    }
    if (isTRUE(verbose)) {
      if (identical(spec$type, "surface") && isTRUE(spec$weighted)) {
        cat(paste0(
          "\nCalculating weighted surface texture metric '", spec$metric,
          "' for variable '", spec$covariate,
          "' at sigma = ", floor(radius), "\n"
        ))
      } else {
        metric_label <- if (identical(spec$type, "surface")) {
          "surface texture metric"
        } else {
          "landscape metric"
        }
        cat(paste0(
          "\nCalculating ", metric_label, " '", spec$metric,
          "' for variable '", spec$covariate,
          "' at radius = ", floor(radius), "\n"
        ))
      }
    }

    if (identical(spec$type, "surface")) {
      if (isTRUE(spec$weighted)) {
        param_idx <- match(spec$covariate, opt_vars$covariate)
        out[[i]] <- .surface_weighted_raster_fft(
          raster = source_raster,
          kernel = kernel,
          sigma = sigma[[param_idx]],
          shape = if (is.null(shape)) NULL else shape[[param_idx]],
          metric = spec$metric,
          pct_wt = pct_wt,
          na.rm = na.rm
        )
      } else {
        out[[i]] <- .surface_metric_raster_fft(
          raster = source_raster,
          radius = radius,
          metric = spec$metric,
          na.rm = na.rm
        )
      }
    } else if (spec$metric %in% .msr_composition_metrics()) {
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

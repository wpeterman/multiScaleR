#' @title Apply Kernels to Raster Layers for Mapping
#'
#' @description Computes the model covariates (kernel-weighted means, or
#' declared landscape or surface metrics) for every cell of one or more raster
#' layers, at the scales estimated by [multiScale_optim()] or at scales you
#' supply. The result is a stack of covariate layers that a fitted model can
#' predict from, so this is the usual first step in mapping a fitted model.
#' `vignette("v03-map-predictions", package = "multiScaleR")` works through
#' projection, extrapolation checks, and clamping.
#'
#' @param raster_stack A `SpatRaster` with the source layer(s). Layer names
#'   must match the covariate names used in fitting, or the `source` names in
#'   `scale_vars`. Use the same map units and cell size as in fitting.
#' @param sigma Numeric vector of kernel scale parameters (or hard radii), in
#'   map units. Supply one value per optimized covariate in declaration order
#'   when `scale_vars` is used, or one per raster layer otherwise. Without
#'   `scale_vars`, a length that does not match the number of layers gives a
#'   warning and the first value is used for every layer. Ignored when
#'   `multiScaleR` is supplied.
#' @param multiScaleR A fitted `multiScaleR` object from [multiScale_optim()],
#'   or a prepared `multiScaleR_data` object from [kernel_prep()]. When
#'   supplied, `sigma`, `shape`, `kernel`, and `scale_vars` are taken from it.
#'   A fitted object supplies the estimated scales. A prepared object supplies
#'   the starting scales of [kernel_prep()]. Default: `NULL`.
#' @param shape Numeric vector of shape parameters for the exponential power
#'   kernel (`kernel = "expow"`), one per optimized covariate. A hard-radius
#'   covariate ignores its shape position. Ignored when `multiScaleR` is
#'   supplied. Default: `NULL`.
#' @param kernel Character. Kernel used for smoothing: `"gaussian"` (default),
#'   `"exp"`, `"expow"`, or `"fixed"`. See [multiScaleR-kernels]. Ignored when
#'   `multiScaleR` is supplied.
#' @param scale_vars Optional covariate specifications from [msr_vars()]. Use
#'   it to project landscape or surface metrics, or derived covariates,
#'   without a fitted object. When `multiScaleR` is supplied, its
#'   specifications are used, and for a fitted model only the covariates in
#'   the model formula are projected.
#' @param pct_wt Numeric between 0 and 1 (exclusive). Sets the size of the
#'   default smoothing window for kernel means and kernel-weighted surface
#'   metrics. The window extends to the line-based (1D) distance that holds
#'   `pct_wt` of the kernel weight, rounded up to whole cells. Larger values
#'   (for example `0.99`) include more of the kernel tail and take longer. It
#'   is not the share of planar weight within the window. Ignored when
#'   `match_buffer = TRUE` and by hard-radius metrics. Default: `0.975`.
#' @param fft Logical. If `TRUE` (default), smoothing uses fast Fourier
#'   transform (FFT) convolution, which is much faster for large rasters and
#'   wide kernels. Set `FALSE` for direct moving-window (focal) calculations.
#'   Both engines support kernel means and the declared surface and landscape
#'   metrics, with the same support and missing-value rules.
#' @param scale_center Logical. If `TRUE`, each smoothed layer is centered and
#'   scaled with the mean and standard deviation stored in
#'   `multiScaleR$scl_params`, so its values match the standardized
#'   covariates used in fitting. Required for prediction from a model fitted
#'   on standardized covariates, which is the default. Needs `multiScaleR`.
#'   Default: `FALSE`.
#' @param clamp Logical. If `TRUE`, standardized values are limited to the
#'   range observed at the survey points (moved by `pct_mx`), which prevents
#'   extrapolation beyond the fitted range. Active only when
#'   `scale_center = TRUE`. Default: `FALSE`.
#' @param pct_mx Numeric greater than -0.5 and no larger than 0.99. When
#'   `clamp = TRUE`, each clamping limit moves outward by `pct_mx` times the
#'   range of that covariate observed at the survey points:
#'   `lower = min - pct_mx * (max - min)` and
#'   `upper = max + pct_mx * (max - min)`, on the standardized scale. For
#'   example, a covariate observed from -1 to 3 is clamped to -1.4 and 3.4
#'   with `pct_mx = 0.1`, and to -0.6 and 2.6 with `pct_mx = -0.1`. Negative
#'   values move the limits inward. Default: `0` (the exact training range).
#' @param na.rm Logical. If `TRUE` (default), missing cells are left out and
#'   the weights of the remaining cells renormalize. If `FALSE`, any window
#'   with a missing cell, including cells beyond the raster edge, returns
#'   `NA`.
#' @param verbose Logical. Print a progress message for each layer. Default:
#'   `TRUE`.
#' @param offset_values Named numeric vector of raw input values for variables
#'   inside model offsets, such as `c(n_weeks = 4)` for
#'   `offset(log(n_weeks))`. When omitted, a log offset uses 1 (zero on the
#'   log scale) and a direct offset uses 0. Set it for a meaningful projection
#'   scenario. See Details.
#' @param include_fixed Logical. Set `TRUE`, with `scale_center = TRUE` and a
#'   fitted native model, to standardize modeled fixed-radius metrics with
#'   their verified original training preparation. The default `FALSE` uses
#'   only the scaling parameters stored on the supplied object. See Details.
#' @param match_buffer Logical. With the default `FALSE`, kernel means use the
#'   `pct_wt` window. Set `TRUE` to use the same circular footprint as
#'   [kernel_prep()]: every cell intersecting a circle of radius `max_D`.
#'   Requires a fitted or prepared `multiScaleR` object, a projected CRS with
#'   square cells at the prepared resolution, and `na.rm = TRUE`. Works with
#'   either value of `fft`. See Details.
#' @param site_values Optional named list of scalar values for nonspatial
#'   predictors, in their native units, for example
#'   `list(season = "spring", effort = 2, surveyed = TRUE)`. Requires a
#'   fitted object. Factor and character values must be fitted levels.
#'   Omitted numeric predictors default to zero, and categorical predictors
#'   need an explicit value. Unknown names and attempts to overwrite a
#'   spatial layer are errors. Use `offset_values` for variables that occur
#'   only inside offsets.
#' @param ... Reserved for deprecated arguments. Only `scale_opt`, a
#'   deprecated alias for `multiScaleR`, is handled (with a warning).
#'
#' @return A `SpatRaster` on the grid of `raster_stack`, with one layer per
#' covariate (as declared in `scale_vars` or stored in `multiScaleR`), or one
#' per raster layer when you supply `sigma` without `scale_vars`. Layers are
#' named as in the model. Read the values as follows:
#' * **Covariate layers.** With `scale_center = FALSE`, values are in the
#'   covariate's own units. A kernel-weighted mean is in the units of the
#'   source layer, and a metric in the units given in [msr_vars()]. With
#'   `scale_center = TRUE`, values are in standard deviations of the fitting
#'   data, so 0 is the average survey point.
#' * **Site-level layers.** With a fitted object, predictors not derived from
#'   rasters (for example survey effort) are added as constant layers, at 0
#'   unless set by `site_values`, so that the stack works with
#'   `terra::predict()`. These do not represent spatial variation.
#' * **Offset layers.** Variables inside offsets are added under their input
#'   names, at the values described for `offset_values`.
#'
#' With a fitted object, the returned R object also carries a
#' `"multiScaleR_projection"` attribute that records how the layers were
#' made. [predict.multiScaleR()] checks it before raster prediction. Raster
#' files written to disk do not store this attribute.
#'
#' @details
#' # Choose where the scales come from
#'
#' * **A fitted model** (`multiScaleR = fit`) is the usual choice for mapping.
#'   Each covariate is smoothed at its estimated sigma or radius with the
#'   fitted kernel.
#' * **A prepared object** (`multiScaleR = kernel_inputs`) uses the starting
#'   scales of [kernel_prep()]. Use it to check how a preparation maps before
#'   fitting.
#' * **Your own values** (`sigma`, `kernel`, and `shape`, with optional
#'   `scale_vars`) let you explore smoothing at any scale without a model.
#'
#' # Set the smoothing window
#'
#' By default (`match_buffer = FALSE`), each kernel mean and kernel-weighted
#' surface metric uses a square window that extends to the line-based
#' distance holding `pct_wt` of the kernel weight. For a Gaussian kernel and
#' the default `pct_wt = 0.975`, that is about 2.24 sigma. Hard-radius metrics
#' use the cells whose centers lie within their radius. This window need not
#' include the same cells that [kernel_prep()] used at the survey points,
#' which were every cell intersecting a circle of radius `max_D`, including
#' partly covered boundary cells.
#'
#' With `match_buffer = TRUE`, kernel means and kernel-weighted surface
#' metrics use that preparation footprint instead. The same projected
#' `sf::st_buffer()` circle and polygon-intersection cell selection are
#' applied around each raster cell center. Partly covered boundary cells get
#' their usual full-cell kernel weights, as in fitting, rather than
#' coverage-fraction weights. Weighted SQ remains a weighted population spread
#' (zero for constant values), and a hard-radius SQ remains a sample standard
#' deviation. Requirements and limits:
#'
#' * Supply a fitted or prepared object, a raster with a projected (not
#'   longitude and latitude) CRS, square cells at the prepared resolution, and
#'   `na.rm = TRUE`. Use the intended source layers in matching map units.
#' * Both smoothing engines support matching. Direct focal calculation
#'   (`fft = FALSE`) can be much slower with a large `max_D`, so check the
#'   cost on a small area first.
#' * A metric radius larger than the prepared `max_D` raises an error,
#'   because its truncated footprint cannot be reproduced by capping the
#'   radius. Prepare again with an adequate `max_D` before fitting such a
#'   model.
#' * Matching aligns the footprint at raster cell centers. It does not undo
#'   distance binning, reconstruct missing CRS metadata, replace missing
#'   source coverage, or remove the difference between a sampled point and a
#'   cell center.
#'
#' # Choose a smoothing engine
#'
#' FFT and direct focal calculations evaluate the same neighborhood with the
#' same finite-cell normalization, so they agree up to small floating-point
#' differences. Weighted means renormalize over the available mapped cells,
#' including at raster edges. Cells that are missing in the source layer stay
#' missing in the output, so smoothing never fills unmapped locations. With
#' `na.rm = FALSE`, incomplete windows, including those that extend beyond the
#' map, return `NA`. Weighted SQ over constant values returns zero, while the
#' standardized higher moments stay undefined when the spread is zero. Use
#' `fft = FALSE` mainly to confirm a result, because it can be much slower.
#'
#' # Standardize and clamp for prediction
#'
#' By default, [multiScale_optim()] fits the model to centered and scaled
#' covariates, so set `scale_center = TRUE` before predicting from the fitted
#' model. The scaling uses the means and standard deviations stored in the
#' supplied object (for a fitted model, those at the estimated scales). Map
#' cells can take covariate values outside the range seen at the survey
#' points. `clamp = TRUE` holds such values at the surveyed limits, and
#' `pct_mx` moves those limits in or out.
#' `vignette("v03-map-predictions", package = "multiScaleR")` shows how to
#' find extrapolated cells and choose a setting.
#'
#' # Fill site-level predictors and offsets
#'
#' With a fitted object, the function finds predictors that no raster layer
#' supplies, such as survey effort or a habitat type recorded in the field.
#' They cannot be mapped, so they are added as constant layers. A numeric
#' predictor defaults to zero, which is its mean when it was centered and
#' scaled. Set `site_values` for any other scenario. Categorical predictors
#' are skipped with a warning unless you supply a fitted level, and their
#' raster category labels keep the fitted levels. Checked prediction also
#' restores logical types. Native `unmarked` models require explicit
#' nonspatial prediction scenarios. When every predictor comes from a raster,
#' no constant layers are added. You can also supply remaining variables
#' through native prediction data or `terra::predict(..., const = scenario)`.
#'
#' Offset terms use their underlying input variables as layer names. For
#' `offset(log(n_weeks))`, the default layer is `n_weeks = 1`, one week of
#' effort. Set `offset_values = c(n_weeks = 4)` to predict for four weeks.
#' Direct offsets default to zero. Compound offsets, such as
#' `offset(log(effort * duration))`, need explicit raw values for each input.
#' Literal constant offsets need no layer. Ambiguous or conflicting default
#' scenarios raise an error.
#'
#' # Include verified fixed-radius metrics
#'
#' A fitted model can contain fixed-radius metrics whose scaling parameters
#' are kept only in its original preparation. Set `include_fixed = TRUE` with
#' `scale_center = TRUE` to standardize them with that preparation. The
#' function first checks that the original fixed transforms and training
#' normalization reproduce the native model covariates exactly. Changed
#' definitions or training inputs raise an error. The check needs:
#'
#' * The retained training buffers, matching fixed metric definitions, and a
#'   raster at the fitted resolution.
#' * A recoverable training frame (for `glm()` and `coxph()`, fit with
#'   `model = TRUE`), with the fixed predictors retained under their original
#'   names and row IDs matching the prepared sites.
#' * A fixed radius no larger than the fitted `max_D`.
#'
#' Transformations that do not retain a fixed predictor under its original
#' name are not reconstructed. This option does not alter your fitted
#' object, the metric definitions, or the interpretation of native
#' uncertainty. [prepare_unmarked_data()] describes the same verification for
#' `unmarked` models.
#'
#' # Predict from the result
#'
#' Pass the returned stack to [predict.multiScaleR()], which checks the
#' projection settings, or to `terra::predict()` with the fitted model
#' (`fit$opt_mod`). For a conditional-logistic step-selection model, native
#' `type = "risk"` with `reference = "zero"` gives a relative selection
#' score, not a choice or occurrence probability. Supply required stratum
#' variables through `const`. See [multiScaleR-prediction].
#'
#' @seealso [multiScaleR-kernels], [msr_vars()], [predict.multiScaleR()],
#'   [multiScaleR-prediction].
#'   `vignette("v03-map-predictions", package = "multiScaleR")` maps a fitted
#'   model, and `vignette("v04-landscape-metrics", package = "multiScaleR")`
#'   and `vignette("v05-surface-metrics", package = "multiScaleR")` map
#'   metrics.
#'
#' @examples
#' ## Smooth a raster at a chosen sigma, without a fitted model
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' hab_100 <- kernel_scale.raster(hab, sigma = 100, kernel = "gaussian",
#'                                verbose = FALSE)
#' # Each cell now holds the kernel-weighted proportion of habitat around it
#' plot(c(hab, hab_100))
#'
#' ## Map a fitted model
#' # A fit stored with the package (the Quick Start model of counts on
#' # land1 and land2, with a site-level covariate)
#' fit <- readRDS(system.file("extdata", "cache-quickstart.rds",
#'                            package = "multiScaleR"))$opt
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))
#' # Smooth each layer at its estimated sigma, standardize it as in fitting,
#' # and clamp it to the range seen at the survey points
#' r_scaled <- kernel_scale.raster(land, multiScaleR = fit,
#'                                 scale_center = TRUE, clamp = TRUE,
#'                                 verbose = FALSE)
#' # land1 and land2 in standard deviation units, plus a constant layer
#' # (0, the mean) for the site-level covariate
#' names(r_scaled)
#' # Expected count in every cell
#' pred <- terra::predict(r_scaled, fit$opt_mod, type = "response")
#' plot(pred)
#'
#' ## Use the preparation footprint of kernel_prep()
#' # match_buffer needs a projected CRS in meters. The bundled raster has
#' # none, so assign one for this example.
#' terra::crs(land) <- "EPSG:32617"
#' r_matched <- kernel_scale.raster(land, multiScaleR = fit,
#'                                  scale_center = TRUE, clamp = TRUE,
#'                                  match_buffer = TRUE, verbose = FALSE)
#' # Change in standardized land1 (standard deviation units) from using the
#' # wider preparation footprint instead of the pct_wt window
#' summary(abs(terra::values(r_matched$land1) -
#'               terra::values(r_scaled$land1))[, 1])
#'
#' @rdname kernel_scale.raster
#' @export
#' @importFrom terra crs rast subset cellFromRowCol crds focal res xyFromCell nlyr clamp nrow ncol values as.matrix

kernel_scale.raster <- function(raster_stack,
                                sigma = NULL,
                                multiScaleR = NULL,
                                shape = NULL,
                                kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                                scale_vars = NULL,
                                pct_wt = 0.975,
                                fft = TRUE,
                                scale_center = FALSE,
                                clamp = FALSE,
                                pct_mx = 0,
                                na.rm = TRUE,
                                verbose = TRUE,
                                offset_values = NULL,
                                include_fixed = FALSE,
                                match_buffer = FALSE,
                                site_values = NULL,
                                ...){

  args <- list(...)

  if ("scale_opt" %in% names(args)) {
    warning("Argument 'scale_opt' is deprecated. Use 'multiScaleR' instead.", call. = FALSE)
    if (is.null(multiScaleR) & !is.null(args$scale_opt)) {
      multiScaleR <- args$scale_opt
    }
  }

  if(!inherits(raster_stack, "SpatRaster")){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }
  validate_scalar_logical(fft, "fft")
  validate_scalar_logical(scale_center, "scale_center")
  validate_scalar_logical(clamp, "clamp")
  validate_scalar_logical(na.rm, "na.rm")
  validate_scalar_logical(verbose, "verbose")
  validate_scalar_logical(include_fixed, "include_fixed")
  validate_scalar_logical(match_buffer, "match_buffer")
  site_values <- .msr_site_values(site_values)
  if (length(site_values) && !inherits(multiScaleR, "multiScaleR"))
    stop("site_values requires a fitted multiScaleR object.", call. = FALSE)
  if (!is.null(offset_values)) {
    if (!is.numeric(offset_values) || is.null(names(offset_values)) ||
        any(!nzchar(names(offset_values))) || anyDuplicated(names(offset_values)) ||
        any(!is.finite(offset_values))) {
      stop("`offset_values` must be a named finite numeric vector.", call. = FALSE)
    }
  }
  validate_scalar_numeric(pct_wt,
                          "pct_wt",
                          lower = 0,
                          upper = 1,
                          inclusive_lower = FALSE,
                          inclusive_upper = FALSE)
  validate_pct_mx(pct_mx)

  if(!is.null(multiScaleR)){
    validate_multiScaleR_input(multiScaleR)
  }
  buffer_radius <- NULL
  if (match_buffer) {
    if (!(inherits(multiScaleR, "multiScaleR") ||
          inherits(multiScaleR, "multiScaleR_data")))
      stop("match_buffer requires a fitted or prepared multiScaleR object.", call. = FALSE)
    if (!isTRUE(na.rm))
      stop("match_buffer requires na.rm = TRUE.", call. = FALSE)
    prepared <- if (inherits(multiScaleR, "multiScaleR"))
      multiScaleR$kernel_inputs else multiScaleR
    buffer_radius <- multiScaleR$max_D
    validate_scalar_numeric(buffer_radius, "prepared max_D", positive = TRUE)
    validate_scalar_numeric(prepared$resolution, "prepared resolution", positive = TRUE)
    if (!nzchar(terra::crs(raster_stack)) || isTRUE(terra::is.lonlat(raster_stack)) ||
        !isTRUE(all.equal(terra::res(raster_stack), rep(prepared$resolution, 2L))))
      stop("match_buffer requires projected square cells at the prepared resolution.", call. = FALSE)
  }
  fixed_contract <- NULL
  if (include_fixed) {
    if (!isTRUE(scale_center) || !inherits(multiScaleR, "multiScaleR"))
      stop("include_fixed requires a supported fitted native model and scale_center = TRUE.", call. = FALSE)
    if (!isTRUE(all.equal(terra::res(raster_stack),
          rep(multiScaleR$kernel_inputs$resolution, 2L))))
      stop("Prediction raster resolution must match the fit.", call. = FALSE)
    fixed_contract <- .msr_fixed_training_contract(multiScaleR)
    fixed_specs <- multiScaleR$kernel_inputs$scale_vars
    fixed_specs <- fixed_specs[fixed_specs$covariate %in% fixed_contract$covariates, , drop = FALSE]
    if (any(fixed_specs$radius > multiScaleR$max_D))
      stop("A fixed radius exceeded the fitting buffer; reprepare with an adequate buffer or supply explicitly verified projection layers.", call. = FALSE)
    # This local list copy leaves the caller's fitted object unchanged.
    multiScaleR$scl_params$mean[fixed_contract$covariates] <- fixed_contract$center
    multiScaleR$scl_params$sd[fixed_contract$covariates] <- fixed_contract$spread
  }

  kernel <- match.arg(kernel)

  if (is.null(scale_vars) && !is.null(multiScaleR)) {
    if (inherits(multiScaleR, "multiScaleR") &&
        !is.null(multiScaleR$kernel_inputs$scale_vars)) {
      scale_vars <- multiScaleR$kernel_inputs$scale_vars
    } else if (inherits(multiScaleR, "multiScaleR_data") &&
               !is.null(multiScaleR$scale_vars)) {
      scale_vars <- multiScaleR$scale_vars
    }
  }

  if (!is.null(scale_vars)) {
    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      kernel <- multiScaleR$kernel_inputs$kernel
    } else if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")) {
      kernel <- multiScaleR$kernel
    }

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      scale_vars <- .msr_scale_vars_in_model(scale_vars = scale_vars,
                                             multiScaleR = multiScaleR)
    }

    scale_vars <- .msr_validate_scale_vars(scale_vars = scale_vars,
                                           raster_stack = raster_stack,
                                           kernel = kernel)
    if (include_fixed && length(fixed_contract$covariates)) {
      fixed <- fixed_contract$covariates
      expected <- multiScaleR$kernel_inputs$scale_vars
      expected <- expected[match(fixed, expected$covariate), , drop = FALSE]
      actual <- scale_vars[match(fixed, scale_vars$covariate), , drop = FALSE]
      if (anyNA(actual$covariate) ||
          !isTRUE(all.equal(as.data.frame(expected), as.data.frame(actual),
                            check.attributes = FALSE)))
        stop("Fixed-radius source and metric definitions must match the fit.", call. = FALSE)
    }
    opt_vars <- .msr_optimized_scale_vars(scale_vars)
    n_optimized <- nrow(opt_vars)

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      sigma <- multiScaleR$scale_est[opt_vars$covariate, "Mean"]
      shape <- if (!is.null(multiScaleR$shape_est)) {
        .msr_shape_values(multiScaleR)[opt_vars$covariate]
      } else {
        NULL
      }
    } else if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")) {
      sigma <- multiScaleR$sigma * multiScaleR$unit_conv
      shape <- multiScaleR$shape
    } else if (n_optimized == 0 && is.null(sigma)) {
      sigma <- numeric(0)
    }

    if (n_optimized > 0 && is.null(sigma)) {
      stop("sigma values must be specified for optimized `scale_vars`.\n",
           call. = FALSE)
    }
    if (n_optimized > 0) {
      validate_numeric_vector(sigma,
                              "sigma",
                              length_ = n_optimized,
                              positive = TRUE)
    }
    if (kernel == "expow" && n_optimized > 0) {
      validate_numeric_vector(shape,
                              "shape",
                              length_ = n_optimized,
                              positive = TRUE)
    }

    if(isTRUE(scale_center) && is.null(multiScaleR)){
      warning("`scale_center = TRUE` was ignored because no `multiScaleR` object was provided.",
              call. = FALSE)
    }

    if(isTRUE(clamp) && !isTRUE(scale_center)){
      warning("`clamp = TRUE` was ignored because `scale_center` is FALSE.",
              call. = FALSE)
    }

    smooth_stack <- .msr_scale_vars_raster(
      raster_stack = raster_stack,
      scale_vars = scale_vars,
      sigma = sigma,
      shape = shape,
      kernel = kernel,
      pct_wt = pct_wt,
      fft = fft,
      na.rm = na.rm,
      verbose = verbose,
      buffer_radius = buffer_radius
    )

    if (isTRUE(scale_center) &&
        (inherits(multiScaleR, "multiScaleR_data") || inherits(multiScaleR, "multiScaleR"))) {
      smooth_stack <- scale_center_raster(r = smooth_stack,
                                          multiScaleR = multiScaleR,
                                          clamp = clamp,
                                          pct_mx = pct_mx)
    }

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      smooth_stack <- .add_site_covariate_rasters(
        smooth_stack = smooth_stack,
        multiScaleR = multiScaleR,
        raster_covs = names(smooth_stack),
        offset_values = offset_values, site_values = site_values
      )
    }

    return(.msr_stamp_projection(smooth_stack, multiScaleR, scale_center,
      match_buffer, scale_vars, include_fixed, site_values,
      list(fft = fft, clamp = clamp, pct_mx = pct_mx, na.rm = na.rm,
           pct_wt = pct_wt, offset_values = offset_values)))
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")){
    covs <- rownames(multiScaleR$scale_est)
    sigma <- multiScaleR$scale_est[,1]
    shape <- .msr_shape_values(multiScaleR)
    kernel <- multiScaleR$kernel_inputs$kernel

    if(!any(covs %in% names(raster_stack))){
      stop('optimized covariate(s) not present in the provided SpatRaster!')
    } else {
      var <- intersect(covs, names(raster_stack))
      c <- which(covs %in% var)
      raster_stack <- subset(raster_stack, var)
      sigma <- sigma[c]
      shape <- shape[c]
    }
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")){
    covs <- colnames(multiScaleR$kernel_dat)
    var <- intersect(covs, names(raster_stack))
    c <- which(var %in% covs)

    sigma <- multiScaleR$sigma[c] * multiScaleR$unit_conv
    shape <- multiScaleR$shape[c]
    kernel <- multiScaleR$kernel

    if(!all(covs %in% names(raster_stack))){
      stop('optimized covariate is not present in the provided SpatRaster!')
    } else {
      raster_stack <- subset(raster_stack, covs)
    }
  } else {
    covs <- var <- names(raster_stack)
  }

  if(is.null(sigma)){
    stop("sigma values must be specified\n")
  }
  validate_numeric_vector(sigma,
                          "sigma",
                          positive = TRUE)

  if(length(sigma) != nlyr(raster_stack)){
    warning("Number of sigma values must equal the number of raster layers!!!  \n  All raster layers will be smoothed using the same sigma value")
    sigma <- rep(sigma[1], nlyr(raster_stack))
  }

  if(kernel == "expow"){
    validate_numeric_vector(shape,
                            "shape",
                            length_ = length(sigma),
                            positive = TRUE)
  }

  if(isTRUE(scale_center) && is.null(multiScaleR)){
    warning("`scale_center = TRUE` was ignored because no `multiScaleR` object was provided.",
            call. = FALSE)
  }

  if(isTRUE(clamp) && !isTRUE(scale_center)){
    warning("`clamp = TRUE` was ignored because `scale_center` is FALSE.",
            call. = FALSE)
  }


  smooth_list <- wt_list <-  vector('list', length(sigma))
  out <- rast(raster_stack[[1]])

  for(i in 1:length(sigma)){
    # lyr <- covs[i]
    lyr <- var[i]

    if (match_buffer) {
      smooth_list[[i]] <- .msr_kernel_raster_one(raster_stack[[lyr]], sigma[i],
        shape[i], kernel, pct_wt, fft, na.rm, buffer_radius = buffer_radius)
      next
    }

    mx <- kernel_dist(kernel = kernel,
                      dimension = "1d",
                      sigma = sigma[i],
                      shape = shape[i],
                      prob = pct_wt)

    r_res <- res(raster_stack)[1]
    focal_d <- ceiling(mx / r_res) * 2

    if((focal_d %% 2) == 0) {
      focal_d <- focal_d + 1
    }

    mat <- matrix(0, focal_d, focal_d)
    r_wt <- rast(mat)
    crs(r_wt) <- crs(raster_stack)
    cntr_crd <- cellFromRowCol(r_wt, focal_d/2, focal_d/2)
    cntr_crd <- xyFromCell(r_wt, ceiling(length(mat)/2))
    cell_crds <- crds(r_wt)
    r_wt[] <- rdist(cntr_crd, cell_crds)[1,] * r_res
    r_wt[] <- scale_type_r(d = as.vector(r_wt),
                           kernel = kernel,
                           sigma = sigma[i],
                           shape = shape[i],
                           output = 'wts')


    wt_mat <- as.matrix(r_wt, wide = T)

    if(verbose){
      cat(paste0("\nSmoothing spatRaster ",i, " of ", length(sigma), ": ",lyr," at sigma = ",floor(sigma[i]),"\n"))
    }
    if(isTRUE(fft)){
      mat <- as.matrix(raster_stack[[lyr]], wide = T)

      fft_mat <- fft_convolution(mat,
                                 wt_mat,
                                 fun = "mean",
                                 na.rm = na.rm)
      values(out) <- as.vector(fft_mat)

      smooth_list[[i]] <- out
    } else {
      smooth_list[[i]] <- focal(raster_stack[[lyr]],
                                w = wt_mat,
                                fun = "mean",
                                na.rm = na.rm,
                                na.policy = "omit",
                                expand = F)
    }

  }
  smooth_stack <- rast(smooth_list)
  names(smooth_stack) <- names(raster_stack)

  if (isTRUE(scale_center) &&
      (inherits(multiScaleR, "multiScaleR_data") || inherits(multiScaleR, "multiScaleR"))) {
    smooth_stack <- scale_center_raster(r = smooth_stack,
                                        multiScaleR = multiScaleR,
                                        clamp = clamp,
                                        pct_mx = pct_mx)
  }

  # add dummy rasters for site-level covariates in fitted model
  if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
    smooth_stack <- .add_site_covariate_rasters(
      smooth_stack = smooth_stack,
      multiScaleR = multiScaleR,
      raster_covs = var,
      offset_values = offset_values, site_values = site_values
    )
  }

  return(.msr_stamp_projection(smooth_stack, multiScaleR, scale_center,
    match_buffer, scale_vars, include_fixed, site_values,
    list(fft = fft, clamp = clamp, pct_mx = pct_mx, na.rm = na.rm,
         pct_wt = pct_wt, offset_values = offset_values)))
}

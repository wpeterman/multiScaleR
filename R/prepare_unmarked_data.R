#' Prepare spatial covariates for native unmarked prediction
#'
#' @description Calculate optimized spatial predictors at new sites and insert
#'   them into an explicitly keyed native prediction data frame. Training means
#'   and standard deviations are reused, so prediction sites are on the fitted
#'   model's covariate scale.
#' @param object A fitted \code{multiScaleR} object containing an
#'   \code{unmarkedFit}, fitted scales, training standardization, and retained
#'   \code{kernel_inputs}.
#' @param newdata A nonempty data frame in the row layout required by native
#'   unmarked prediction for your chosen component. Include site key columns
#'   and all nonspatial variables, including factors and offset inputs.
#'   Existing optimized spatial columns are replaced. Supply fixed-radius
#'   predictors on their native fitted scale explicitly unless
#'   \code{include_fixed = TRUE}. Other columns and row order remain unchanged.
#' @param kernel_inputs Spatial inputs created by \code{\link{kernel_prep}}
#'   for the prediction sites. Use the fitted kernel, source/metric definitions,
#'   buffer extent, and raster resolution. Cell-level inputs are required for
#'   landscape and surface metrics; lean binned inputs are supported for kernel
#'   means.
#' @param join_by A data frame with one unique key row per prepared prediction
#'   site, in \code{kernel_inputs} point order. Its column names must occur in
#'   \code{newdata}. Every prediction row must match a key. Repeated keys in
#'   \code{newdata} repeat that site's predictors, allowing explicit site/year
#'   or site/observation scenarios without changing their order.
#' @param include_fixed Logical. The default \code{FALSE} leaves fixed-radius
#'   predictors as explicit native inputs. Set \code{TRUE} to calculate modeled
#'   fixed-radius metrics too, after verifying that their retained preparation
#'   values and standardization match the fitted static site covariates. Missing
#'   or inconsistent training information raises an error.
#' @details
#' This function prepares data only. Pass its result to
#' \code{unmarked::predict(object$opt_mod, type = ..., newdata = ...)}, or supply
#' the same inputs and keys to \code{\link{prepare_prediction}} to retain a
#' compact prediction snapshot. Native prediction determines component meaning,
#' factor handling, offsets, missing-data behavior, and conditional intervals.
#' Neither function adds uncertainty in selected spatial scales.
#'
#' Spatial predictors are static site covariates, as required by the unmarked
#' fitting adapter. You must supply the native row layout for detection and
#' annual components explicitly. The function does not infer observation times,
#' expand sites into years, or modify an unmarked frame. It does not borrow
#' nonspatial values from the fitting data.
#'
#' Distances inside prepared inputs are converted using their internal
#' \code{unit_conv}; reported fitted scales remain in map units. Metric
#' definitions must match the fit. Fixed-radius predictors are recalculated only
#' with \code{include_fixed = TRUE}. Their normalization comes from the original
#' preparation object, and their fitted site values must match that object after
#' explicit training-site alignment. This rejects manually rescaled or stale
#' fixed predictors rather than guessing their fitted scale.
#' Binned versus exact
#' kernel evaluation must also match. Repreparing bins at new sites can change
#' their representative distances slightly, just as in other binned prediction
#' paths. Source names are checked, but raster contents are not fingerprinted.
#' Use the intended source rasters and coordinate system when preparing sites.
#' When internal fitted parameters are retained, they are reused directly to
#' avoid rounding from a conversion to reported units and back. Inconsistent
#' stored parameters and reported scales or shapes raise an error.
#'
#' Missing spatial values remain missing for the native prediction method.
#' Infinite values or invalid training standard deviations raise errors.
#' Prediction values are not clamped or standardized within the new sample.
#' @return A data frame with the same rows, row names, and nonspatial columns as
#'   \code{newdata}. Calculated spatial columns contain the fitted transforms in
#'   standard-deviation units relative to the training sample. These are native
#'   model inputs, not response predictions. Do not standardize them again.
#' @examples
#' \dontrun{
#' # Keys follow the order of points used in kernel_prep().
#' keys <- data.frame(site = prediction_sites$site)
#' inputs <- kernel_prep(prediction_sites, habitat,
#'   max_D = fit$max_D, kernel = fit$kernel_inputs$kernel,
#'   scale_vars = fit$kernel_inputs$scale_vars,
#'   bin = !is.null(fit$kernel_inputs$binned))
#' nd <- prepare_unmarked_data(fit, scenarios, inputs, join_by = keys)
#' unmarked::predict(fit$opt_mod, type = "state", newdata = nd)
#'
#' # Freeze the same native component and controls for repeated use.
#' prepared <- prepare_prediction(fit, scenarios, inputs, component = "state",
#'   join_by = keys)
#' predict(prepared)
#'
#' # Opt in to checked reconstruction of modeled fixed-radius metrics.
#' nd_fixed <- prepare_unmarked_data(fit, scenarios, inputs, join_by = keys,
#'   include_fixed = TRUE)
#' }
#' @seealso \code{\link{kernel_prep}}, \code{\link{prepare_prediction}},
#'   \code{\link{kernel_scale.raster}}
#' @export
prepare_unmarked_data <- function(object, newdata, kernel_inputs, join_by,
                                  include_fixed = FALSE) {
  validate_scalar_logical(include_fixed, "include_fixed")
  if (!inherits(object, "multiScaleR") ||
      !methods::is(object$opt_mod, "unmarkedFit"))
    stop("object must contain a fitted native unmarked model.", call. = FALSE)
  if (!is.data.frame(newdata) || !nrow(newdata) || anyDuplicated(names(newdata)))
    stop("newdata must be a nonempty data frame with unique column names.", call. = FALSE)
  if (missing(join_by) || is.null(join_by))
    stop("Supply explicit join_by keys in prepared spatial site order.", call. = FALSE)
  train <- object$kernel_inputs
  inputs <- kernel_inputs
  covs <- names(object$scl_params$mean)
  if (is.null(train) || !length(covs) || anyDuplicated(covs) ||
      !all(covs %in% .unmarked_model_predictors(object$opt_mod)))
    stop("The fit must retain spatial inputs and named training standardization.", call. = FALSE)
  if (!inherits(inputs, "multiScaleR_data") || !identical(inputs$kernel, train$kernel))
    stop("Supply kernel_prep() inputs with the fitted kernel.", call. = FALSE)
  validate_scalar_numeric(inputs$unit_conv, "kernel_inputs$unit_conv", positive = TRUE)
  if (!isTRUE(all.equal(inputs$max_D, object$max_D)) ||
      !isTRUE(all.equal(inputs$min_D, object$min_D)))
    stop("Prediction buffer extent must match the fitted min_D and max_D.", call. = FALSE)
  if (!isTRUE(all.equal(inputs$resolution, train$resolution)))
    stop("Prediction raster resolution must match the fit.", call. = FALSE)
  specs <- train$scale_vars
  actual <- inputs$scale_vars
  if (is.null(specs) != is.null(actual))
    stop("Spatial source and metric definitions must match the fit.", call. = FALSE)
  if (!is.null(specs)) {
    specs <- specs[match(covs, specs$covariate), , drop = FALSE]
    actual <- actual[match(covs, actual$covariate), , drop = FALSE]
    if (anyNA(specs$covariate) || anyNA(actual$covariate) ||
        any(!specs$optimize) ||
        !isTRUE(all.equal(as.data.frame(specs), as.data.frame(actual),
                          check.attributes = FALSE)))
      stop("Spatial source and metric definitions must match the fit.", call. = FALSE)
  }
  center <- object$scl_params$mean[covs]
  spread <- object$scl_params$sd[covs]
  if (any(!is.finite(center)) || length(spread) != length(covs) ||
      any(!is.finite(spread) | spread <= 0))
    stop("Training means and positive standard deviations must be finite.", call. = FALSE)
  n <- if (!is.null(inputs$d_list)) length(inputs$d_list) else inputs$binned$n_points
  if (length(n) != 1L || !is.finite(n) || n < 1L)
    stop("Prediction inputs must contain prepared spatial sites.", call. = FALSE)
  if (any(covs %in% names(join_by)))
    stop("join_by key columns cannot be spatial predictors.", call. = FALSE)
  rows <- .unmarked_site_map(newdata, join_by, n)
  # Keep the fit's parameter order; fixed-radius metrics have no parameter.
  optimized <- if (is.null(specs)) covs else specs$covariate[specs$optimize]
  sigma <- object$scale_est[optimized, "Mean"] / inputs$unit_conv
  if (length(sigma) != length(optimized) || any(!is.finite(sigma) | sigma <= 0))
    stop("Finite positive fitted spatial scales are required.", call. = FALSE)
  par <- sigma
  if (inputs$kernel == "expow") {
    shape <- object$shape_est[optimized, "Mean"]
    if (length(shape) != length(optimized) || any(!is.finite(shape) | shape <= 0))
      stop("Finite positive fitted kernel shapes are required.", call. = FALSE)
    par <- c(par, shape)
  }
  stored <- object$optim_results$par
  if (!is.null(stored)) {
    # Reuse the actual fitted parameters without a report-unit round trip.
    # Reported scales must still agree exactly with that retained parameter set.
    validate_scalar_numeric(train$unit_conv, "fitted unit_conv", positive = TRUE)
    parameter_names <- rownames(object$scale_est)
    p <- length(parameter_names)
    expected_length <- p * if (inputs$kernel == "expow") 2L else 1L
    if (length(stored) != expected_length || any(!is.finite(stored)) ||
        !isTRUE(all.equal(unname(stored[seq_len(p)] * train$unit_conv),
          unname(object$scale_est$Mean), tolerance = 0, check.attributes = FALSE)))
      stop("Stored fitted parameters disagree with reported spatial scales.", call. = FALSE)
    idx <- match(optimized, parameter_names)
    par <- stored[idx] * (train$unit_conv / inputs$unit_conv)
    if (inputs$kernel == "expow") {
      shape <- stored[p + idx]
      if (!isTRUE(all.equal(unname(shape), unname(object$shape_est[optimized, "Mean"]),
                            tolerance = 0, check.attributes = FALSE)))
        stop("Stored fitted parameters disagree with reported kernel shapes.", call. = FALSE)
      par <- c(par, shape)
    }
  }
  binned_covs <- function(x) intersect(covs, x$binned$covariates)
  if (!setequal(binned_covs(inputs), binned_covs(train)))
    stop("Prediction binned versus exact kernel evaluation must match the fit.", call. = FALSE)
  # Use the same transform evaluator as refitting, with no model-row filtering.
  # Each site is evaluated once before repeated native scenario rows are joined.
  context <- list(covs = covs, n_covs = length(optimized), scale_vars = specs,
    unit_conv = inputs$unit_conv, resolution = inputs$resolution,
    n_cols = inputs$n_cols, binned = inputs$binned)
  raw <- .msr_kernel_cov_w(par, inputs$d_list, inputs$raw_cov, inputs$kernel, context)
  if (any(is.infinite(raw))) stop("Spatial prediction values cannot be infinite.", call. = FALSE)
  scaled <- sweep(sweep(raw, 2L, center), 2L, spread, "/")
  newdata[covs] <- as.data.frame(scaled[rows, , drop = FALSE])
  if (include_fixed)
    newdata <- .msr_unmarked_fixed_prediction_data(object, newdata, inputs,
                                                  rows, context, par, join_by)
  newdata
}

# Fixed predictors are not refitted by the scale optimizer. Their original
# preparation normalization is reusable only when it reproduces native inputs.
.msr_unmarked_fixed_prediction_data <- function(object, newdata, inputs,
                                               rows, context, par, join_by) {
  train <- object$kernel_inputs
  specs <- train$scale_vars
  fixed <- specs$covariate[!specs$optimize &
    specs$covariate %in% .unmarked_model_predictors(object$opt_mod)]
  if (!length(fixed)) return(newdata)
  if (any(fixed %in% names(join_by)))
    stop("join_by key columns cannot be spatial predictors.", call. = FALSE)
  expected <- specs[match(fixed, specs$covariate), , drop = FALSE]
  actual <- inputs$scale_vars[match(fixed, inputs$scale_vars$covariate), , drop = FALSE]
  if (anyNA(actual$covariate) ||
      !isTRUE(all.equal(as.data.frame(expected), as.data.frame(actual),
                        check.attributes = FALSE)))
    stop("Fixed-radius source and metric definitions must match the fit.", call. = FALSE)
  contract <- .msr_fixed_training_contract(object)
  # Keep all definitions so the evaluator can resolve fixed radii without
  # assigning them optimization parameters. Each prepared site is evaluated once.
  context$scale_vars <- specs
  raw <- .msr_kernel_cov_w(par, inputs$d_list, inputs$raw_cov, inputs$kernel,
                          context, covariates = fixed)
  if (any(is.infinite(raw))) stop("Spatial prediction values cannot be infinite.", call. = FALSE)
  scaled <- sweep(sweep(raw, 2L, contract$center), 2L, contract$spread, "/")
  newdata[fixed] <- as.data.frame(scaled[rows, , drop = FALSE])
  newdata
}

# Shared by point preparation and opt-in raster standardization. A retained
# mean/SD alone is not proof that the fitted native predictor used that scale.
.msr_fixed_training_contract <- function(object) {
  train <- object$kernel_inputs
  specs <- train$scale_vars
  model <- .analysis_model(object$opt_mod)
  unmarked <- methods::is(model, "unmarkedFit")
  variables <- if (unmarked) .unmarked_model_predictors(model) else
    .model_predictors(model)
  fixed <- specs$covariate[!specs$optimize &
    specs$covariate %in% variables]
  if (!length(fixed)) return(list(covariates = character(),
                                  center = numeric(), spread = numeric()))
  if (unmarked) site_covs <- model@data@siteCovs else {
    if (!inherits(model, c("glm", "coxph")) || !is.data.frame(model$model))
      stop("Automatic fixed predictors require a supported native model with a retained model frame.", call. = FALSE)
    site_covs <- model$model
  }
  if (!all(fixed %in% names(site_covs)) ||
      !all(fixed %in% names(train$kernel_dat)))
    stop("Automatic fixed predictors require retained prepared and native numeric covariates.", call. = FALSE)
  center <- train$scl_params$mean[fixed]
  spread <- train$scl_params$sd[fixed]
  if (length(center) != length(fixed) || any(!is.finite(center)) ||
      length(spread) != length(fixed) || any(!is.finite(spread) | spread <= 0))
    stop("Fixed predictors require finite retained training means and positive standard deviations.", call. = FALSE)
  if (unmarked) {
    training_rows <- .unmarked_site_map(site_covs, object$join_by,
      nrow(train$kernel_dat), rownames(train$kernel_dat))
  } else {
    ids <- rownames(train$kernel_dat)
    training_rows <- match(rownames(site_covs), ids)
    if (is.null(ids) || anyDuplicated(ids) || anyNA(training_rows))
      stop("Native model row IDs must match retained prepared training sites.", call. = FALSE)
  }
  training_context <- list(covs = character(), n_covs = 0L,
    scale_vars = specs, unit_conv = train$unit_conv,
    resolution = train$resolution, n_cols = train$n_cols, binned = train$binned)
  # These transforms have fixed radii and need no optimization parameters.
  training_raw <- .msr_kernel_cov_w(numeric(), train$d_list, train$raw_cov,
    train$kernel, training_context, covariates = fixed)
  training_scaled <- sweep(sweep(training_raw, 2L, center), 2L, spread, "/")
  for (nm in fixed) {
    if (!is.numeric(site_covs[[nm]]) || !is.numeric(train$kernel_dat[[nm]]) ||
        !isTRUE(all.equal(as.numeric(site_covs[[nm]]),
          as.numeric(train$kernel_dat[[nm]][training_rows]),
          tolerance = 0, check.attributes = FALSE)) ||
        !isTRUE(all.equal(as.numeric(train$kernel_dat[[nm]]),
          as.numeric(training_scaled[, nm]), tolerance = 0,
          check.attributes = FALSE)))
      stop("Fixed predictor '", nm,
           "' does not match its retained training preparation; supply it explicitly.",
           call. = FALSE)
  }
  list(covariates = fixed, center = center, spread = spread)
}

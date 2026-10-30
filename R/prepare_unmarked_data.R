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
  covs <- intersect(rownames(object$scale_est), .unmarked_model_predictors(object$opt_mod))
  spatial <- .msr_refit_spatial_values(object, kernel_inputs, covs)
  if (any(covs %in% names(join_by)))
    stop("join_by key columns cannot be spatial predictors.", call. = FALSE)
  rows <- .unmarked_site_map(newdata, join_by, spatial$n_sites)
  newdata[covs] <- as.data.frame(spatial$value[rows, , drop = FALSE])
  if (include_fixed)
    newdata <- .msr_fixed_prediction_data(object, newdata, kernel_inputs,
      rows, spatial$context, spatial$par, join_by)
  newdata
}

# Fixed predictors are not refitted by the scale optimizer. Their original
# preparation normalization is reusable only when it reproduces native inputs.
.msr_fixed_prediction_data <- function(object, newdata, inputs,
                                               rows, context, par, join_by) {
  train <- object$kernel_inputs
  specs <- train$scale_vars
  if (is.null(specs)) return(newdata)
  model <- .analysis_model(object$opt_mod)
  variables <- if (methods::is(model, "unmarkedFit"))
    .unmarked_model_predictors(model) else .model_predictors(model)
  fixed <- specs$covariate[!specs$optimize & specs$covariate %in% variables]
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
  if (is.null(specs)) return(list(covariates = character(),
                                  center = numeric(), spread = numeric()))
  # Native wrappers retain bounds outside kernel_inputs to avoid duplication.
  if (is.null(train$max_D)) train$max_D <- object$max_D
  model <- .analysis_model(object$opt_mod)
  unmarked <- methods::is(model, "unmarkedFit")
  variables <- if (unmarked) .unmarked_model_predictors(model) else
    .model_predictors(model)
  fixed <- specs$covariate[!specs$optimize &
    specs$covariate %in% variables]
  if (!length(fixed)) return(list(covariates = character(),
                                  center = numeric(), spread = numeric()))
  if (unmarked) site_covs <- model@data@siteCovs else {
    site_covs <- if (inherits(model, "merMod")) model@frame else if (is.list(model)) {
      if (is.data.frame(model$model)) model$model else if (is.data.frame(model$frame))
        model$frame else if (inherits(model, "lme") && is.data.frame(model$data)) {
          # nlme retains original data separately; restrict it to fitted rows.
          used <- rownames(model$fitted)
          rows <- match(used, rownames(model$data))
          if (length(used) && !anyNA(rows)) model$data[rows, , drop = FALSE] else NULL
        } else NULL
    } else NULL
    if (!is.data.frame(site_covs))
      stop("Automatic fixed predictors require a supported native model with a retained model frame.", call. = FALSE)
  }
  if (!all(fixed %in% names(site_covs)) ||
      !all(fixed %in% names(train$kernel_dat)))
    stop("Automatic fixed predictors require retained prepared and native numeric covariates.", call. = FALSE)
  if (unmarked) {
    training_rows <- .unmarked_site_map(site_covs, object$join_by,
      nrow(train$kernel_dat), rownames(train$kernel_dat))
  } else {
    ids <- rownames(train$kernel_dat)
    training_rows <- match(rownames(site_covs), ids)
    if (is.null(ids) || anyDuplicated(ids) || anyNA(training_rows))
      stop("Native model row IDs must match retained prepared training sites.", call. = FALSE)
  }
  .msr_fixed_reference_contract(train, site_covs, training_rows, fixed)
}

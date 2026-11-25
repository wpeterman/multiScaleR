#' Prepare spatial covariates for native unmarked prediction
#'
#' @description
#' Calculate the optimized spatial predictors at new sites and insert them
#' into an explicitly keyed data frame for native `unmarked` prediction. The
#' training means and standard deviations are reused, so the new sites are on
#' the fitted model's covariate scale. Pass the result to
#' `unmarked::predict()`, or give the same inputs to [prepare_prediction()]
#' for a reusable snapshot.
#'
#' @param object A fitted `multiScaleR` object whose `opt_mod` is an
#'   `unmarkedFit`, with its fitted scales, training standardization, and
#'   retained `kernel_inputs`.
#' @param newdata A nonempty data frame of prediction rows, in the layout
#'   that native `unmarked` prediction needs for your chosen component.
#'   Include the site key columns and every non-spatial
#'   variable, including factors and offset inputs. Existing optimized spatial
#'   columns are replaced. Supply fixed-radius predictors on their fitted scale
#'   unless `include_fixed = TRUE`. Other columns and the row order are left
#'   unchanged.
#' @param kernel_inputs A `multiScaleR_data` object from [kernel_prep()] for
#'   the new sites. Use the fitted kernel, source and metric definitions,
#'   buffer extent (`max_D`), raster resolution, and the same binned or exact
#'   evaluation as the fit. Landscape and surface metrics need cell-level
#'   inputs, while lean binned inputs are enough for kernel means.
#' @param join_by A data frame with one unique key row per prepared site, in
#'   `kernel_inputs` point order. Its column names must occur in `newdata`,
#'   and every `newdata` row must match a key. Repeated keys in `newdata`
#'   repeat that site's spatial values, which lets you write site-by-year or
#'   site-by-survey scenarios without changing their order. Required.
#' @param include_fixed Logical. The default `FALSE` leaves fixed-radius
#'   predictors as native inputs that you supply. Set `TRUE` to calculate
#'   modeled fixed-radius metrics at the new sites too, after verifying that
#'   their retained preparation values and standardization match the fitted
#'   static site covariates. Missing or inconsistent training information
#'   raises an error.
#'
#' @details
#' # What the function does and does not do
#'
#' This function prepares data only. Native prediction determines the
#' component's meaning, factor handling, offsets, missing-data behavior, and
#' confidence limits, and neither this function nor [prepare_prediction()]
#' adds uncertainty in the selected scales. [multiScaleR-prediction] lists
#' what each `unmarked` component predicts.
#'
#' Spatial predictors are static site covariates, as the `unmarked` fitting
#' adapter requires. Supply the native row layout for detection and annual
#' components yourself. The function does not infer observation times,
#' expand sites into years, modify an `unmarkedFrame`, or borrow non-spatial
#' values from the fitting data.
#'
#' # Match the spatial inputs to the fit
#'
#' Distances inside `kernel_inputs` are converted with its own `unit_conv`,
#' and the reported scales remain in map units. Metric definitions and the
#' choice of binned or exact kernel evaluation must match the fit. Preparing
#' bins at new sites can shift their representative distances slightly, as in
#' other binned prediction routes. Source names are checked, but raster
#' contents are not fingerprinted, so prepare the sites from the intended
#' source rasters and coordinate system. When the internal fitted parameters
#' are retained, they are reused directly, which avoids rounding from a
#' conversion to reported units and back. Stored parameters that disagree
#' with the reported scales or shapes raise an error.
#'
#' The calculated columns use the training sample's mean and standard
#' deviation. Prediction values are never clamped or standardized within the
#' new sample. Missing spatial values stay missing for the native prediction
#' method, while infinite values or invalid training standard deviations
#' raise errors.
#'
#' # Recompute fixed-radius predictors
#'
#' Fixed-radius predictors are recalculated only with `include_fixed = TRUE`.
#' Their normalization comes from the original preparation object, and their
#' fitted site values must match that object after explicit training-site
#' alignment. This rejects manually rescaled or stale fixed predictors, changed
#' radii, ambiguous training keys, and missing training buffers rather than
#' guessing their fitted scale. Keep the original inputs if you need this
#' verification.
#'
#' # Distance-removal models
#'
#' For `unmarked::gdistremoval()` fits, only active component formulas
#' identify spatial predictors. Availability is inactive with one primary
#' period, and the distance formula is inactive for a uniform key. Optimized
#' static site names must not also occur in the native yearly or observation
#' covariate tables, so use distinct names for distinct inputs. Dot formulas use
#' the original component's available training-column names. GDR models with
#' active formula offsets or TMB zero inflation are rejected, as [aic_tab()]
#' explains.
#'
#' @return A data frame with the same rows, row names, and non-spatial columns
#'   as `newdata`. The calculated spatial columns hold the fitted kernel
#'   transforms in standard-deviation units of the training sample, the scale
#'   the native model was fitted on. They are model inputs, not predictions.
#'   Do not standardize them again.
#'
#' @seealso [prepare_prediction()], [predict.multiScaleR()],
#'   [multiScaleR-prediction], [kernel_prep()], [kernel_scale.raster()], and
#'   `vignette("v02-user-guide", package = "multiScaleR")` for fitting and
#'   mapping `unmarked` models.
#'
#' @examples
#' ## A stored N-mixture fit used in vignette("v02-user-guide"). Abundance
#' ## depends on the bin1 and cont2 layers of a simulated landscape.
#' nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
#'                             package = "multiScaleR"))
#' layers <- subset(sim_rast(user_seed = 999, dim = 250), c("bin1", "cont2"))
#'
#' ## Three new sites, prepared with the fitted kernel, max_D, covariate
#' ## definitions, and binning
#' new_pts <- terra::vect(cbind(x = c(800, 1250, 1700),
#'                              y = c(900, 1250, 1600)))
#' new_inputs <- kernel_prep(new_pts, layers, max_D = nmix$max_D,
#'                           kernel = nmix$kernel_inputs$kernel,
#'                           scale_vars = nmix$kernel_inputs$scale_vars,
#'                           bin = !is.null(nmix$kernel_inputs$binned),
#'                           verbose = FALSE)
#'
#' ## Keys follow the prepared point order. Site B appears twice, for
#' ## example in two survey years, and reuses the same spatial values.
#' keys <- data.frame(site = c("A", "B", "C"))
#' scenarios <- data.frame(site = c("A", "B", "B", "C"))
#' nd <- prepare_unmarked_data(nmix, scenarios, new_inputs, join_by = keys)
#' nd  # bin1 and cont2 in training standard-deviation units
#'
#' ## Expected abundance at each row, with native 95% confidence limits
#' unmarked::predict(nmix$opt_mod, type = "state", newdata = nd)
#'
#' ## The same rows and component as a reusable snapshot
#' prepared <- prepare_prediction(nmix, scenarios, new_inputs,
#'                                component = "state", join_by = keys)
#' predict(prepared)
#' @export
prepare_unmarked_data <- function(object, newdata, kernel_inputs, join_by,
                                  include_fixed = FALSE) {
  validate_scalar_logical(include_fixed, "include_fixed")
  if (!inherits(object, "multiScaleR") ||
      !methods::is(object$opt_mod, "unmarkedFit"))
    stop("object must contain a fitted native unmarked model.", call. = FALSE)
  .msr_gdr_native_guard(object$opt_mod)
  if (!is.data.frame(newdata) || !nrow(newdata) || anyDuplicated(names(newdata)))
    stop("newdata must be a nonempty data frame with unique column names.", call. = FALSE)
  if (missing(join_by) || is.null(join_by))
    stop("Supply explicit join_by keys in prepared spatial site order.", call. = FALSE)
  covs <- intersect(rownames(object$scale_est), .unmarked_model_predictors(object$opt_mod))
  .msr_gdr_spatial_guard(object$opt_mod, covs)
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
    site_covs <- if (inherits(model, "merMod")) model@frame else
      if (inherits(model, "HLfit")) .msr_selection_frame(model) else if (is.list(model)) {
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
    # A fixed metric used only inside an expression (poly(ed, 2), log(ed))
    # appears in the model frame as that expression. Take its raw column
    # from the model data, restricted to the model frame's rows.
    if (!all(fixed %in% names(site_covs))) {
      raw <- .model_data(model)
      rows <- if (is.data.frame(raw)) match(rownames(site_covs), rownames(raw)) else NA
      if (is.data.frame(raw) && all(fixed %in% names(raw)) && !anyNA(rows))
        site_covs <- raw[rows, , drop = FALSE]
    }
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

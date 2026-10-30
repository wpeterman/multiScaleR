# Recognize documented zero-effect controls without comparing formula environments.
# This describes conditioning at zero, never integration over random effects.
.msr_native_zero_random_effects <- function(model, controls, kind = NULL) {
  if (inherits(model, "lme") || identical(kind, "lme")) {
    level <- controls$level
    return(is.numeric(level) && length(level) == 1L && !is.na(level) && level == 0)
  }
  if (!inherits(model, c("merMod", "glmmTMB")) &&
      !methods::is(model, "unmarkedFit") &&
      !isTRUE(kind %in% c("glmmTMB", "unmarked", "merMod", "lmerMod", "glmerMod")))
    return(FALSE)
  value <- controls$re.form
  (is.atomic(value) && length(value) == 1L && is.na(value)) ||
    (inherits(value, "formula") && length(value) == 2L &&
      is.numeric(value[[2L]]) && length(value[[2L]]) == 1L &&
      !is.na(value[[2L]]) && value[[2L]] == 0)
}

# Shared native-refit transformation. Keep physical scales, source definitions,
# fitted parameter order, and training normalization identical across models.
.msr_refit_spatial_values <- function(object, inputs, covs) {
  train <- object$kernel_inputs
  if (is.null(train) || !length(covs) || anyDuplicated(covs))
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
  # Keep the fit's parameter order; fixed-radius metrics have no parameter.
  optimized <- if (is.null(specs)) covs else specs$covariate[specs$optimize]
  sigma <- object$scale_est[optimized, "Mean"] / inputs$unit_conv
  if (length(sigma) != length(optimized) || any(!is.finite(sigma) | sigma <= 0))
    stop("Finite positive fitted spatial scales are required.", call. = FALSE)
  par <- sigma
  if (inputs$kernel == "expow") {
    shape <- .msr_shape_values(object)[optimized]
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
      active <- !optimized %in% .msr_hard_radius_covariates(specs)
      if (!isTRUE(all.equal(unname(shape[active]), unname(object$shape_est[optimized[active], "Mean"]),
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
  list(value = scaled, n_sites = n, context = context, par = par)
}

# Construct native prediction data without refitting or restandardizing the
# prediction sample. Explicit keys allow multiple scenario rows for one site.
.msr_refit_prediction_data <- function(object, newdata, kernel_inputs,
                                       join_by = NULL, include_fixed = FALSE) {
  validate_scalar_logical(include_fixed, "include_fixed")
  if (!is.data.frame(newdata) || !nrow(newdata) || anyDuplicated(names(newdata)))
    stop("Supply nonempty newdata and matching kernel_inputs together.", call. = FALSE)
  model <- .analysis_model(object$opt_mod)
  variables <- if (methods::is(model, "unmarkedFit"))
    .unmarked_model_predictors(model) else .model_predictors(model)
  covs <- intersect(rownames(object$scale_est), variables)
  spatial <- .msr_refit_spatial_values(object, kernel_inputs, covs)
  if (is.null(join_by)) {
    rows <- .msr_profile_input_rows(newdata, kernel_inputs)
  } else {
    if (any(covs %in% names(join_by)))
      stop("join_by key columns cannot be spatial predictors.", call. = FALSE)
    rows <- .unmarked_site_map(newdata, join_by, spatial$n_sites)
  }
  newdata[covs] <- as.data.frame(spatial$value[rows, , drop = FALSE])
  if (include_fixed)
    newdata <- .msr_fixed_prediction_data(object, newdata, kernel_inputs, rows,
      spatial$context, spatial$par, join_by)
  newdata
}

# Native model classes keep their own prediction targets. Their uncertainty is
# conditional on estimated spatial parameters; no joint covariance is invented.
.msr_predict_refit <- function(object, newdata, kernel_inputs, type, se.fit,
                               uncertainty, interval, ...) {
  validate_scalar_logical(se.fit, "se.fit")
  if (interval != "none" || (se.fit && uncertainty != "conditional"))
    stop("Refit predictions support native conditional SEs only; use uncertainty = 'conditional' and interval = 'none'.", call. = FALSE)
  model <- .analysis_model(object$opt_mod)
  if (is.null(model)) stop("The fit must retain its native opt_mod.", call. = FALSE)
  dots <- list(...)
  native_args <- dots$native_args
  dots$native_args <- NULL
  if (!is.null(native_args)) {
    if (!is.list(native_args) || (length(native_args) &&
        (is.null(names(native_args)) || anyNA(names(native_args)) ||
         any(!nzchar(names(native_args))) || anyDuplicated(names(native_args)) ||
         any(names(native_args) %in% c("object", "newdata", "type", "se.fit",
                                       "join_by", "include_fixed", "native_args")))))
      stop("native_args must uniquely name native controls, excluding wrapper controls.", call. = FALSE)
    if (any(names(native_args) %in% names(dots)))
      stop("Native prediction controls must not be supplied twice.", call. = FALSE)
    dots <- c(dots, native_args)
  }
  join_by <- dots$join_by
  include_fixed <- if (is.null(dots$include_fixed)) FALSE else dots$include_fixed
  dots[c("join_by", "include_fixed")] <- NULL
  validate_scalar_logical(include_fixed, "include_fixed")
  if (is.null(newdata) != is.null(kernel_inputs))
    stop("Supply newdata and matching kernel_inputs together, or omit both for fitted rows.", call. = FALSE)
  if (is.null(newdata) && (!is.null(join_by) || include_fixed))
    stop("join_by and include_fixed require explicit newdata and kernel_inputs.", call. = FALSE)
  if (!is.null(newdata)) {
    newdata <- .msr_refit_prediction_data(object, newdata, kernel_inputs,
                                          join_by, include_fixed)

  }
  .msr_native_prediction_rows(object, newdata, type, se.fit, dots)
}

# Shared native target for prepared point rows and checked raster chunks.
.msr_native_prediction_rows <- function(object, newdata, type, se.fit, dots = list()) {
  model <- .analysis_model(object$opt_mod)
  unmarked <- methods::is(model, "unmarkedFit")
  # A survival linear predictor or risk is not a within-choice-set probability.
  # Keep its native names, including lp/risk/survival where the model permits.
  if (is.null(type)) type <- if (unmarked) "state" else
    if (inherits(model, "coxph")) "lp" else "response"
  validate_character_scalar(type, "type")
  if (inherits(model, "coxph")) {
    if (type == "response")
      stop("Native survival refits have no generic response probability; use type = 'lp' or 'risk', or the conditional_logit backend for choice-set probabilities.", call. = FALSE)
    if (type == "link") type <- "lp"
  }
  if (inherits(model, "lm") && !inherits(model, "glm") && type == "link")
    type <- "response"
  if (!is.null(newdata)) {
    if (unmarked) {
      component <- .msr_unmarked_marginal_formula(model, type)
      if (identical(type, "det") && "formula" %in% methods::slotNames(model) &&
          length(model@formula) == 3L) component <- model@formula[[2L]]
      required <- .msr_marginal_variables(component, offsets = TRUE)
    } else {
      formula <- stats::formula(model)
      required <- union(.model_predictors(model),
        .msr_marginal_variables(formula[[length(formula)]], offsets = TRUE))
      if (.msr_native_zero_random_effects(model, dots)) {
        # Inspect fixed syntax directly. Formula environments do not decide
        # whether grouping columns or random-only slopes are needed.
        required <- .msr_marginal_variables(formula[[length(formula)]], offsets = TRUE)
        if (inherits(model, "glmmTMB")) for (f in
            model$modelInfo$allForm[c("ziformula", "dispformula")])
          if (!is.null(f)) required <- union(required,
            .msr_marginal_variables(f[[length(f)]], offsets = TRUE))
      } else required <- union(required, all.vars(formula[[length(formula)]]))
      if (inherits(model, "coxph") && type %in% c("expected", "survival"))
        required <- union(required, all.vars(formula[[2L]]))
    }
    absent <- setdiff(required, names(newdata))
    if (length(absent)) stop("newdata must explicitly contain formula variables: ",
      paste(absent, collapse = ", "), call. = FALSE)
  }
  args <- c(list(object = model, type = type), dots)
  if (!is.null(newdata)) args$newdata <- newdata
  predictor <- if (unmarked) unmarked::predict else stats::predict
  if (!unmarked) args$se.fit <- se.fit
  value <- do.call(predictor, args)
  if (unmarked && !se.fit) value <- value$Predicted
  if (se.fit && !unmarked &&
      (!is.list(value) || is.null(value$se.fit)))
    stop("The native prediction method did not supply se.fit; this model's uncertainty needs its native inference interface.", call. = FALSE)
  attr(value, "prediction_contract") <- .msr_native_prediction_contract(
    object, if (unmarked) "unmarked" else class(model)[1L], type, dots)
  value
}

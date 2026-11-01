# Bounded frequentist adapter; independent of the Bayesian likelihood/prior API.
.msr_fit_mixed_backend <- function(model, inputs, join_by, par, start_strategy,
                                    n_cores, refit_fn, control, verbose, call) {
  for (pkg in c("RTMB", "TMB", "glmmTMB"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("joint_mixed requires the optional package ", pkg, ".", call. = FALSE)
  if (!inherits(model, "glmmTMB") ||
      !inputs$kernel %in% c("gaussian", "exp"))
    stop("joint_mixed requires a glmmTMB model and Gaussian or exponential inputs.", call. = FALSE)
  if (!is.null(join_by) || !is.null(refit_fn) || start_strategy != "single" ||
      (!is.null(n_cores) && n_cores != 1L))
    stop("joint_mixed uses retained rows and serial independent starts; join_by and refit_fn are unsupported.", call. = FALSE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!names(control) %in% c("starts", "n_starts")) || anyDuplicated(names(control)))))
    stop("joint_mixed profile_control accepts only starts and n_starts.", call. = FALSE)
  if (!is.null(model$call$offset))
    stop("Use a formula offset() rather than the separate offset argument.", call. = FALSE)
  mf <- stats::model.frame(model)
  specs <- inputs$scale_vars
  surface_names <- if (is.null(specs)) character() else specs$covariate[
    specs$type == "surface" & specs$metric == "sq" & specs$weighted & specs$optimize]
  spatial_names <- c(names(inputs$binned$vsum), surface_names)
  fixed_terms <- model$modelInfo$reTrms$cond$terms$fixed
  active <- .msr_active_formula_variables(fixed_terms)
  covs <- intersect(active, spatial_names)
  other_spatial <- intersect(spatial_names, unique(c(
    all.vars(model$modelInfo$allForm$ziformula), all.vars(model$modelInfo$allForm$dispformula))))
  if (length(other_spatial))
    stop("Spatial predictors in zero-inflation or dispersion formulas are unsupported.", call. = FALSE)
  if (!length(covs)) stop("No supported spatial predictors occur in the formula.", call. = FALSE)
  if (!is.null(specs)) {
    used <- specs[specs$covariate %in% active, , drop = FALSE]
    optimized <- used[used$optimize, , drop = FALSE]
    if (any(!(optimized$type == "kernel" | optimized$covariate %in% surface_names)) || !all(optimized$covariate %in% covs))
      stop("joint_mixed requires optimized kernel means or weighted surface SQ.", call. = FALSE)
    if (any(!used$optimize)) .msr_fixed_training_contract(list(opt_mod = model,
      kernel_inputs = inputs))
  }
  rows <- .msr_profile_input_rows(mf, inputs)
  engine <- if (any(covs %in% surface_names)) .msr_prepare_surface_engine(inputs, covs, rows) else
    .msr_prepare_kernel_engine(inputs$binned, covs, rows, inputs$kernel)
  # Recover raw formula variables (including offset/transform inputs), then fix
  # the sample to the native fit's retained rows. Never resolve missing data in
  # a formula environment or silently change the sample during optimization.
  dat <- .model_data(model, effects = "all")
  required <- unique(unlist(lapply(model$modelInfo$allForm[c("formula", "ziformula", "dispformula")], all.vars)))
  if (!is.data.frame(dat) || !all(required %in% names(dat)) ||
      anyDuplicated(rownames(dat)) || anyNA(match(rownames(mf), rownames(dat))))
    stop("Cannot recover explicit formula variables and retained rows from the native model.", call. = FALSE)
  dat <- dat[match(rownames(mf), rownames(dat)), , drop = FALSE]
  # Reuse deterministic start generation, without adding a GLM dispersion coordinate.
  ctl <- .msr_profile_control(control)
  bounds <- .msr_profile_starts(list(objective = list(engine = engine, family = "poisson",
    parameter_names = engine$parameter_names)), inputs, par, ctl, model)
  # Reject unsupported features on the ORIGINAL fit before constructing a fresh
  # call that intentionally omits starts, maps, priors, weights and subsets.
  .msr_native_mixed_design(model, engine, covs, bounds$starts[1L, ],
                           rownames(mf), check_kernel = FALSE)
  native_refit <- function(kp) {
    nd <- dat
    kd <- .msr_eval_kernel_engine(engine, unname(kp))
    for (nm in covs) nd[[nm]] <- kd$value[, nm]
    glmmTMB::glmmTMB(formula = stats::formula(model), data = nd,
      family = stats::family(model), ziformula = model$modelInfo$allForm$ziformula,
      dispformula = model$modelInfo$allForm$dispformula,
      contrasts = model$modelInfo$contrasts, na.action = stats::na.fail)
  }
  fits <- .msr_native_joint(engine, covs, rownames(mf), bounds$starts,
                            bounds$lower, bounds$upper, native_refit)
  if (is.na(fits$selected)) {
    condition <- structure(list(message = "No joint_mixed start passed joint Hessian and cold native audits.",
      call = NULL, start_results = fits$outcomes), class = c("msr_mixed_fit_error", "error", "condition"))
    stop(condition)
  }
  best <- fits$outcomes[[fits$selected]]
  kp <- best$log_scale
  kd <- .msr_eval_kernel_engine(engine, kp)
  scale <- exp(kp) * inputs$unit_conv
  missing_hessian <- matrix(NA_real_, length(kp), length(kp))
  # Keep only modeled raw columns, with fitted kernel values, so transformed
  # effect curves do not depend on a cold refit's temporary data environment.
  prediction_data <- dat[, required, drop = FALSE]
  for (nm in covs) prediction_data[[nm]] <- kd$value[, nm]
  out <- list(backend = "joint_mixed",
    scale_est = data.frame(Mean = scale, SE = NA_real_, row.names = covs), shape_est = NULL,
    optim_results = list(par = exp(kp), par_unscale = scale, value = best$objective,
      convergence = best$optim$convergence, message = best$optim$message,
      hessian = missing_hessian, hessian_unscale = missing_hessian),
    opt_mod = fits$model, fitted_mod_original = model, min_D = inputs$min_D, max_D = inputs$max_D,
    kernel_inputs = inputs[setdiff(names(inputs), c("min_D", "max_D"))],
    scl_params = list(mean = stats::setNames(kd$center, covs), sd = stats::setNames(kd$spread, covs)),
    join_by = NULL, opt_context = NULL, profile_scale_est = NULL,
    mixed_fit = list(engine = engine, log_scale = kp, prediction_data = prediction_data,
      starts = bounds$starts,
      lower = bounds$lower, upper = bounds$upper, start_results = fits$outcomes,
      selected = fits$selected, required = setdiff(required, c(all.vars(stats::formula(model)[[2L]]), covs))),
    uncertainty = list(available = FALSE, reasons = "Mixed-model scale uncertainty is not calibrated."),
    diagnostics = .empty_diagnostics(), next_run = NULL, warn_message = 0, call = call)
  class(out) <- "multiScaleR"
  out <- .msr_attach_fixed_reference(out)
  out$diagnostics$joint_mixed <- list(triggered = TRUE, uncertainty_available = FALSE,
    reasons = out$uncertainty$reasons, attempted_starts = length(fits$outcomes),
    eligible_starts = sum(vapply(fits$outcomes, function(x) x$eligible, TRUE)),
    boundary = any(kp - bounds$lower < 1e-5 | bounds$upper - kp < 1e-5))
  # Scale intervals are deliberately unavailable; do not run an interval-based
  # search-extent recommendation on missing uncertainty estimates.
  if (verbose) message("Joint mixed fit: ", out$diagnostics$joint_mixed$eligible_starts,
    "/", length(fits$outcomes), " starts passed. Scale uncertainty is unavailable.")
  out
}

.msr_predict_mixed <- function(object, newdata, kernel_inputs, type, se.fit,
                                uncertainty, interval, include_fixed = FALSE, ...) {
  dots <- list(...)
  native_args <- dots$native_args
  dots$native_args <- NULL
  if (!is.null(native_args)) {
    reserved <- c("object", "newdata", "type", "se.fit", "kernel_inputs", "include_fixed", "native_args")
    if (!is.list(native_args) || (length(native_args) &&
        (is.null(names(native_args)) || anyNA(names(native_args)) ||
         any(!nzchar(names(native_args))) || anyDuplicated(names(native_args)) ||
         any(names(native_args) %in% c(reserved, names(dots))))))
      stop("native_args must uniquely name native controls, excluding wrapper controls and duplicates.", call. = FALSE)
    dots <- c(dots, native_args)
  }
  dots <- .msr_native_control_names(object$opt_mod, dots)
  native_covariance <- isTRUE(as.logical(dots[["cov.fit"]]))
  if (interval != "none" || ((se.fit || native_covariance) && uncertainty != "conditional"))
    stop("joint_mixed supports conditional prediction SEs only; use uncertainty = 'conditional' and interval = 'none'.", call. = FALSE)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit))
    stop("se.fit must be one logical value.", call. = FALSE)
  newdata <- .msr_fixed_point_data(object, newdata, kernel_inputs, include_fixed)
  newdata <- .msr_mixed_prediction_data(object, newdata, kernel_inputs)
  do.call(stats::predict, c(list(object = object$opt_mod, newdata = newdata,
    type = type, se.fit = se.fit), dots))
}

# Shared by live prediction and frozen native-output snapshots. Keep spatial
# normalization, explicit formula variables, and row/unit checks in one place.
.msr_mixed_prediction_data <- function(object, newdata, kernel_inputs) {
  if (!is.null(newdata) || !is.null(kernel_inputs)) {
    e <- object$mixed_fit$engine
    if (!is.data.frame(newdata) || (is.null(kernel_inputs$binned) && !inherits(e, "msr_surface_engine")) ||
        !identical(kernel_inputs$kernel, e$kernel))
      stop("Supply newdata and matching kernel_inputs together.", call. = FALSE)
    if (!all(object$mixed_fit$required %in% names(newdata)))
      stop("newdata must explicitly contain all fixed, offset, and grouping variables.", call. = FALSE)
    validate_scalar_numeric(kernel_inputs$unit_conv, "unit_conv", positive = TRUE)
    b <- kernel_inputs$binned
    rows <- .msr_profile_input_rows(newdata, kernel_inputs)
    ratio <- kernel_inputs$unit_conv / object$kernel_inputs$unit_conv
    if (!is.null(b)) b$bin_rep <- b$bin_rep * ratio
    if (inherits(e, "msr_surface_engine")) {
      expected <- e$surface_specs
      if (is.null(kernel_inputs$scale_vars)) stop("New surface specifications are required.", call. = FALSE)
      actual <- kernel_inputs$scale_vars[match(expected$covariate, kernel_inputs$scale_vars$covariate), , drop = FALSE]
      fields <- c("covariate", "source", "type", "metric", "weighted", "optimize")
      if (!isTRUE(all.equal(actual[, fields], expected[, fields], check.attributes = FALSE)))
        stop("New surface source and metric specifications must match the fitted model.", call. = FALSE)
      prepared <- kernel_inputs; prepared$binned <- b
      prepared$d_list <- lapply(prepared$d_list, function(d) d * ratio)
      en <- .msr_prepare_surface_engine(prepared, e$covariates, rows, standardize = FALSE)
    } else en <- .msr_prepare_kernel_engine(b, e$covariates, rows, e$kernel, standardize = FALSE)
    kd <- .msr_eval_kernel_engine(en, object$mixed_fit$log_scale)
    for (nm in e$covariates) newdata[[nm]] <-
      (kd$raw[, nm] - object$scl_params$mean[nm]) / object$scl_params$sd[nm]
  }
  newdata
}

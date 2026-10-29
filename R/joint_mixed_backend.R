# Bounded frequentist adapter; independent of the Bayesian likelihood/prior API.
.msr_fit_mixed_backend <- function(model, inputs, join_by, par, start_strategy,
                                    n_cores, refit_fn, control, verbose, call) {
  for (pkg in c("RTMB", "TMB", "glmmTMB"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("joint_mixed requires the optional package ", pkg, ".", call. = FALSE)
  if (!inherits(model, "glmmTMB") || is.null(inputs$binned) ||
      !inputs$kernel %in% c("gaussian", "exp"))
    stop("joint_mixed requires a glmmTMB model and Gaussian or exponential binned inputs.", call. = FALSE)
  if (!is.null(join_by) || !is.null(refit_fn) || start_strategy != "single" ||
      (!is.null(n_cores) && n_cores != 1L))
    stop("joint_mixed uses retained rows and serial independent starts; join_by and refit_fn are unsupported.", call. = FALSE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      any(!names(control) %in% c("starts", "n_starts")) || anyDuplicated(names(control)))))
    stop("joint_mixed profile_control accepts only starts and n_starts.", call. = FALSE)
  if (!is.null(model$call$offset))
    stop("Use a formula offset() rather than the separate offset argument.", call. = FALSE)
  mf <- stats::model.frame(model)
  covs <- intersect(all.vars(stats::formula(model)), names(inputs$binned$vsum))
  other_spatial <- intersect(names(inputs$binned$vsum), unique(c(
    all.vars(model$modelInfo$allForm$ziformula), all.vars(model$modelInfo$allForm$dispformula))))
  if (length(other_spatial))
    stop("Spatial predictors in zero-inflation or dispersion formulas are unsupported.", call. = FALSE)
  if (!length(covs)) stop("No binned kernel predictors occur in the formula.", call. = FALSE)
  specs <- inputs$scale_vars
  if (!is.null(specs)) {
    used <- specs[specs$covariate %in% all.vars(stats::formula(model)), , drop = FALSE]
    if (any(used$type != "kernel" | !used$optimize) || !all(used$covariate %in% covs))
      stop("joint_mixed requires optimized kernel-type spatial predictors only.", call. = FALSE)
  }
  rows <- .msr_profile_rows(mf, inputs$binned)
  engine <- .msr_prepare_kernel_engine(inputs$binned, covs, rows, inputs$kernel)
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
  out <- list(backend = "joint_mixed",
    scale_est = data.frame(Mean = scale, SE = NA_real_, row.names = covs), shape_est = NULL,
    optim_results = list(par = exp(kp), par_unscale = scale, value = best$objective,
      convergence = best$optim$convergence, message = best$optim$message,
      hessian = missing_hessian, hessian_unscale = missing_hessian),
    opt_mod = fits$model, fitted_mod_original = model, min_D = inputs$min_D, max_D = inputs$max_D,
    kernel_inputs = inputs[setdiff(names(inputs), c("min_D", "max_D"))],
    scl_params = list(mean = stats::setNames(kd$center, covs), sd = stats::setNames(kd$spread, covs)),
    join_by = NULL, opt_context = NULL, profile_scale_est = NULL,
    mixed_fit = list(engine = engine, log_scale = kp, starts = bounds$starts,
      lower = bounds$lower, upper = bounds$upper, start_results = fits$outcomes,
      selected = fits$selected, required = setdiff(required, c(all.vars(stats::formula(model)[[2L]]), covs))),
    uncertainty = list(available = FALSE, reasons = "Mixed-model scale uncertainty is not calibrated."),
    diagnostics = .empty_diagnostics(), next_run = NULL, warn_message = 0, call = call)
  class(out) <- "multiScaleR"
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
                                uncertainty, interval, ...) {
  if (interval != "none" || (se.fit && uncertainty != "conditional"))
    stop("joint_mixed supports conditional prediction SEs only; use uncertainty = 'conditional' and interval = 'none'.", call. = FALSE)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit))
    stop("se.fit must be one logical value.", call. = FALSE)
  newdata <- .msr_mixed_prediction_data(object, newdata, kernel_inputs)
  stats::predict(object$opt_mod, newdata = newdata, type = type, se.fit = se.fit, ...)
}

# Shared by live prediction and frozen native-output snapshots. Keep spatial
# normalization, explicit formula variables, and row/unit checks in one place.
.msr_mixed_prediction_data <- function(object, newdata, kernel_inputs) {
  if (!is.null(newdata) || !is.null(kernel_inputs)) {
    e <- object$mixed_fit$engine
    if (!is.data.frame(newdata) || is.null(kernel_inputs$binned) ||
        !identical(kernel_inputs$kernel, e$kernel))
      stop("Supply newdata and matching binned kernel_inputs together.", call. = FALSE)
    if (!all(object$mixed_fit$required %in% names(newdata)))
      stop("newdata must explicitly contain all fixed, offset, and grouping variables.", call. = FALSE)
    validate_scalar_numeric(kernel_inputs$unit_conv, "unit_conv", positive = TRUE)
    b <- kernel_inputs$binned
    rows <- .msr_profile_rows(newdata, b)
    b$bin_rep <- b$bin_rep * kernel_inputs$unit_conv / object$kernel_inputs$unit_conv
    en <- .msr_prepare_kernel_engine(b, e$covariates, rows, e$kernel, standardize = FALSE)
    kd <- .msr_eval_kernel_engine(en, object$mixed_fit$log_scale)
    for (nm in e$covariates) newdata[[nm]] <-
      (kd$raw[, nm] - object$scl_params$mean[nm]) / object$scl_params$sd[nm]
  }
  newdata
}

# Public backend adapter. Keep model-frame/design decisions outside the numeric
# objective so every likelihood evaluation uses exactly the same sample.
.msr_profile_rows <- function(data, binned) {
  ids <- binned$point_ids
  if (is.null(ids)) ids <- as.character(seq_len(nrow(binned$vsum[[1L]])))
  if (length(ids) != nrow(binned$vsum[[1L]]) || anyNA(ids) ||
      anyDuplicated(ids) || anyDuplicated(rownames(data)))
    stop("Prepared point IDs and model row names must be unique.", call. = FALSE)
  rows <- match(rownames(data), as.character(ids))
  if (anyNA(rows))
    stop("Model/data row names must match prepared point IDs. Preserve row names through subsetting.",
         call. = FALSE)
  rows
}

# Surface-only preparation has no binned means. Its explicit preparation rows
# supply the same row-ID contract without inventing binned summaries.
.msr_profile_input_rows <- function(data, inputs) {
  if (!is.null(inputs$binned)) return(.msr_profile_rows(data, inputs$binned))
  ids <- rownames(inputs$kernel_dat)
  if (is.null(ids) || length(ids) != length(inputs$d_list) || anyNA(ids) ||
      anyDuplicated(ids) || anyDuplicated(rownames(data)))
    stop("Surface preparation requires matching explicit point IDs and raw buffers.", call. = FALSE)
  rows <- match(rownames(data), ids)
  if (anyNA(rows)) stop("Model/data row names must match prepared point IDs.", call. = FALSE)
  rows
}

.msr_profile_contract <- function(model, inputs, kernel_engine = "reference") {
  if (!inherits(model, "glm") || !is.data.frame(model$model))
    stop("profile_glm requires a glm or MASS::glm.nb fit with model = TRUE.", call. = FALSE)
  if (inputs$kernel == "fixed")
    stop("profile_glm requires smooth kernels.", call. = FALSE)
  if (!is.null(model$call$offset))
    stop("Use an offset() term in the formula, not the separate offset argument.", call. = FALSE)
  fam <- model$family
  family <- if (inherits(model, "negbin") ||
                grepl("^Negative Binomial", fam$family)) "nbinom2" else fam$family
  links <- c(gaussian = "identity", poisson = "log", binomial = "logit", nbinom2 = "log")
  if (!family %in% names(links) || fam$link != links[[family]])
    stop("profile_glm supports Gaussian identity, Poisson/NB2 log, and binomial logit models.",
         call. = FALSE)
  mf <- model$model
  tt <- stats::terms(model)
  specs <- inputs$scale_vars
  surface_names <- if (is.null(specs)) character() else specs$covariate[
    specs$type == "surface" & specs$metric == "sq" & specs$weighted & specs$optimize]
  covs <- intersect(all.vars(stats::delete.response(tt)), c(names(inputs$binned$vsum), surface_names))
  if (!length(covs))
    stop("No supported spatial predictors occur in the model formula; kernel means require bin = TRUE.", call. = FALSE)
  if (!is.null(specs)) {
    used <- specs[specs$covariate %in% all.vars(tt), , drop = FALSE]
    optimized <- used[used$optimize, , drop = FALSE]
    if (any(!(optimized$type == "kernel" | optimized$covariate %in% surface_names)))
      stop("profile_glm supports optimized kernel means and weighted surface SQ only; use refit for other surface metrics or hard-radius predictors.",
           call. = FALSE)
    if (!all(optimized$covariate %in% covs))
      stop("Every modeled kernel mean must have binned summaries.", call. = FALSE)
    if (any(!used$optimize)) .msr_fixed_training_contract(list(opt_mod = model,
      kernel_inputs = inputs))
  }
  # A kernel covariate must correspond to one additive design column. Fixed
  # covariates may have interactions/transforms, using the fitted terms object.
  for (label in attr(tt, "term.labels")) {
    expr <- str2lang(label)
    involved <- intersect(all.vars(expr), covs)
    if (length(involved) && !(length(involved) == 1L &&
                             identical(expr, as.name(involved))))
      stop("Kernel predictors must be bare additive terms; interactions and transformations involving them are unsupported.",
           call. = FALSE)
  }
  for (idx in attr(tt, "offset")) {
    if (length(intersect(all.vars(attr(tt, "variables")[[idx + 1L]]), covs)))
      stop("Kernel predictors cannot occur in offsets.", call. = FALSE)
  }
  if (!all(covs %in% names(mf)) ||
      !all(vapply(mf[covs], function(z) is.numeric(z) && is.null(dim(z)), logical(1))))
    stop("Kernel predictors must be numeric model-frame columns.", call. = FALSE)
  X <- stats::model.matrix(tt, mf, contrasts.arg = model$contrasts)
  if (anyDuplicated(colnames(X)))
    stop("Design column names must be distinct; rename colliding factor levels or predictors.",
         call. = FALSE)
  if (!all(covs %in% colnames(X)) || any(!is.finite(X)))
    stop("The fitted model must have a finite, recoverable design matrix.", call. = FALSE)
  fixed_columns <- setdiff(colnames(X), covs)
  rows <- .msr_profile_input_rows(mf, inputs)
  e <- if (any(covs %in% surface_names))
    .msr_prepare_surface_engine(inputs, covs, rows, implementation = kernel_engine) else
      .msr_prepare_kernel_engine(inputs$binned, covs, rows, inputs$kernel,
                                implementation = kernel_engine)
  y <- stats::model.response(mf)
  trials <- rep(1, nrow(mf))
  if (family == "binomial") {
    if (is.matrix(y) && ncol(y) == 2L) {
      trials <- rowSums(y)
      if (any(!is.finite(y) | y < 0 | y != floor(y)))
        stop("Grouped binomial responses must be integer success/failure counts.", call. = FALSE)
      y <- y[, 1L]
    } else {
      if (is.factor(y) && nlevels(y) == 2L) y <- as.numeric(y) - 1
      if (!is.numeric(y) || any(!y %in% c(0, 1)))
        stop("Use binary responses or cbind(successes, failures), not weighted proportions.",
             call. = FALSE)
    }
  }
  weights <- stats::model.weights(mf)
  if (is.null(weights)) weights <- rep(1, nrow(mf))
  if (any(!is.finite(weights) | weights <= 0 | weights != floor(weights)) ||
      (family == "gaussian" && any(weights != 1)))
    stop("profile_glm requires positive integer frequency weights; Gaussian weights must all equal one.",
         call. = FALSE)
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nrow(mf))
  objective <- .msr_profile_glm_objective(e, as.numeric(y), family,
                                        X[, fixed_columns, drop = FALSE],
                                        offset, weights, trials)
  list(objective = objective, frame = mf, terms = tt, original_columns = colnames(X),
       fixed_columns = fixed_columns, contrasts = model$contrasts,
       xlevels = model$xlevels, rows = rows, covariates = covs)
}

.msr_profile_control <- function(control, joint_nb2 = FALSE) {
  defaults <- list(n_starts = 6L, starts = NULL, theta_bounds = c(1e-4, 1e5),
                   theta_start = NULL, gradient_tolerance = 1e-3,
                   iter.max = 500L, eval.max = 1000L, kernel_engine = "reference",
                   coefficient_start = "cold", worker_pool = NULL)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(!nzchar(names(control))) ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("Unknown or unnamed entries in profile_control.", call. = FALSE)
  defaults[names(control)] <- control
  if (!is.null(defaults$worker_pool) && (!joint_nb2 ||
      !inherits(defaults$worker_pool, "SOCKcluster") || !length(defaults$worker_pool)))
    stop("worker_pool requires a nonempty PSOCK cluster and backend = 'joint_nb2'.", call. = FALSE)
  if (!is.character(defaults$coefficient_start) || length(defaults$coefficient_start) != 1L ||
      is.na(defaults$coefficient_start) || !defaults$coefficient_start %in% c("cold", "warm"))
    stop("coefficient_start must be 'cold' or 'warm'.", call. = FALSE)
  if (!is.character(defaults$kernel_engine) || length(defaults$kernel_engine) != 1L ||
      is.na(defaults$kernel_engine) || !defaults$kernel_engine %in% c("reference", "stream"))
    stop("kernel_engine must be 'reference' or 'stream'.", call. = FALSE)
  for (nm in c("n_starts", "iter.max", "eval.max"))
    validate_scalar_numeric(defaults[[nm]], paste0("profile_control$", nm),
                            positive = TRUE, integerish = TRUE)
  validate_scalar_numeric(defaults$gradient_tolerance, "gradient_tolerance", positive = TRUE)
  tb <- defaults$theta_bounds
  if (joint_nb2 && is.null(control$theta_bounds)) {
    tb <- defaults$theta_bounds <- c(1e-4, Inf)
  }
  if (joint_nb2) {
    if (!is.numeric(tb) || length(tb) != 2L || !is.finite(tb[1]) || tb[1] < 1e-8 ||
        !identical(tb[2], Inf))
      stop("joint_nb2 theta_bounds must be c(positive lower bound >= 1e-8, Inf).", call. = FALSE)
  } else if (!is.numeric(tb) || length(tb) != 2L || any(!is.finite(tb)) ||
             tb[1] <= 0 || tb[2] <= tb[1] || tb[2] > 1e6)
    stop("theta_bounds must be increasing, positive, and no larger than 1e6.", call. = FALSE)
  if (!is.null(defaults$theta_start))
    validate_scalar_numeric(defaults$theta_start, "theta_start", positive = TRUE)
  defaults
}

.msr_profile_starts <- function(contract, inputs, par, control, model) {
  obj <- contract$objective
  e <- obj$engine
  p <- e$p
  lower <- rep(log(inputs$min_D / inputs$unit_conv), p)
  upper <- rep(log(inputs$max_D / inputs$unit_conv), p)
  if (e$kernel == "expow") {
    lower <- c(lower, rep(log(0.75), p))
    upper <- c(upper, rep(log(50), p))
  }
  if (any(lower >= upper)) stop("min_D must be smaller than max_D.", call. = FALSE)
  nk <- length(lower)
  if (!is.null(control$starts)) {
    starts <- control$starts
    if (!is.matrix(starts) || !is.numeric(starts) || !nrow(starts) ||
        ncol(starts) != nk || any(!is.finite(starts) | starts <= 0))
      stop("profile_control$starts must be a positive matrix of kernel parameters, one start per row.",
           call. = FALSE)
    starts <- log(starts)
  } else {
    fractions <- seq(0.15, 0.85, length.out = control$n_starts)
    starts <- t(matrix(vapply(seq_along(fractions), function(i) {
      # Deterministic, partly staggered scale starts avoid requiring an RNG seed.
      f <- if (i %% 2L) rep(fractions[i], p) else
        pmin(0.95, pmax(0.05, fractions[i] + seq(-0.12, 0.12, length.out = p)))
      c(lower[seq_len(p)] + f * (upper[seq_len(p)] - lower[seq_len(p)]),
        if (e$kernel == "expow") rep(log(c(1, 2, 3)[(i-1L) %% 3L + 1L]), p))
    }, numeric(nk)), nrow = nk))
  }
  if (!is.null(par)) {
    if (length(par) != nk || any(!is.finite(par) | par <= 0))
      stop("par must contain positive kernel parameters in model predictor order.", call. = FALSE)
    starts <- rbind(log(unname(par)), starts)
  }
  if (obj$family == "nbinom2") {
    initial <- control$theta_start
    if (is.null(initial)) initial <- model$theta
    if (is.null(initial) || !is.finite(initial) || initial <= 0) initial <- 2
    ts <- rep(c(initial, 0.3, 2, 20), length.out = nrow(starts))
    ts <- pmax(control$theta_bounds[1], pmin(control$theta_bounds[2], ts))
    starts <- cbind(starts, log(ts))
    lower <- c(lower, log(control$theta_bounds[1]))
    upper <- c(upper, log(control$theta_bounds[2]))
  }
  colnames(starts) <- obj$parameter_names
  list(starts = starts, lower = lower, upper = upper)
}

.msr_profile_conditional_model <- function(contract, fitted, model, kernel_par) {
  obj <- contract$objective
  kd <- .msr_eval_kernel_engine(obj$engine, unname(kernel_par))
  mf <- contract$frame
  for (nm in contract$covariates) mf[[nm]] <- kd$value[, nm]
  X <- fitted$design[, contract$original_columns, drop = FALSE]
  fit <- stats::glm.fit(X, if (obj$family == "binomial") obj$y / obj$trials else obj$y,
                       start = if (isTRUE(obj$boundary_nb2)) unname(fitted$beta[colnames(X)]) else NULL,
                       weights = obj$weights * if (obj$family == "binomial") obj$trials else 1,
                       offset = obj$offset, family = fitted$fit$family,
                       intercept = attr(contract$terms, "intercept") == 1L,
                       control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
  # Frequency-weighted grouped binomial data are repeated groups, not a single
  # group with multiplied trials. Keep the likelihood's combinatorial constants.
  fit$aic <- 2 * fitted$value + 2 * (ncol(X) + (obj$family %in% c("gaussian", "nbinom2")))
  fit$model <- mf
  fit$terms <- contract$terms
  fit$formula <- stats::formula(model)
  fit$contrasts <- contract$contrasts
  fit$xlevels <- contract$xlevels
  fit$call <- model$call
  fit$na.action <- model$na.action
  fit$x <- X
  class(fit) <- c("glm", "lm")
  if (obj$family == "nbinom2") {
    fit$theta <- fitted$theta
    fit$SE.theta <- NA_real_
    fit$twologlik <- -2 * fitted$value
    fit$aic <- 2 * fitted$value + 2 * (ncol(X) + 1L)
    class(fit) <- c("negbin", class(fit))
  }
  attr(fit, "scale_uncertainty") <- "conditional on optimized kernel parameters"
  # The original call refers to pre-optimization covariates. Prevent update()
  # from silently refitting that stale design.
  class(fit) <- c("msr_conditional_glm", class(fit))
  fit
}

#' @exportS3Method stats::update
#' @noRd
update.msr_conditional_glm <- function(object, ...) {
  stop("Refit the original model and call multiScale_optim() again; update() on opt_mod would reuse pre-optimization kernel data.",
       call. = FALSE)
}

# MASS's summary printer computes decimal places from theta and cannot format
# theta = Inf. The exact boundary has Poisson coefficient information; retain
# NB2 likelihood metadata on the fitted object and use the Poisson summary.
#' @exportS3Method base::summary
#' @noRd
summary.msr_conditional_glm <- function(object, ...) {
  if (inherits(object, "negbin") && is.infinite(object$theta))
    return(stats::summary.glm(object, ...))
  NextMethod("summary")
}

.msr_fit_profile_backend <- function(model, inputs, join_by, par, start_strategy,
                                      n_cores, refit_fn, control, verbose, call,
                                      backend = "profile_glm") {
  if (!is.null(join_by) || !is.null(refit_fn))
    stop("profile_glm uses retained model rows and does not accept join_by or refit_fn.",
         call. = FALSE)
  if ((!is.null(n_cores) && n_cores > 1L) || start_strategy != "single")
    stop("profile_glm runs serial full multi-start fits. Use n_cores = 1, start_strategy = 'single', and profile_control.",
         call. = FALSE)
  joint_nb2 <- identical(backend, "joint_nb2")
  control <- .msr_profile_control(control, joint_nb2)
  pool <- control$worker_pool
  # Connections belong to the caller, never to the saved fit or inference controls.
  control["worker_pool"] <- list(NULL)
  if (!is.null(pool))
    .msr_stop_unsafe_parallel_ram(inputs, model, NULL, NULL, length(pool), TRUE)
  contract <- .msr_profile_contract(model, inputs, control$kernel_engine)
  if (joint_nb2) {
    old <- contract$objective
    contract$objective <- .msr_joint_nb2_objective(old$engine, old$y, old$family,
      old$X_fixed, old$offset, old$weights, old$trials,
      coefficient_start = control$coefficient_start)
  } else if (identical(control$coefficient_start, "warm")) {
    contract$objective <- .msr_profile_nb2_continuation(contract$objective)
  }
  obj <- contract$objective
  bounds <- if (joint_nb2) .msr_joint_nb2_bounds(contract, inputs, par, control, model) else
    .msr_profile_starts(contract, inputs, par, control, model)
  fits <- if (joint_nb2 && !is.null(pool)) .msr_joint_nb2_pool(obj, bounds, control, pool) else
    if (joint_nb2) .msr_joint_nb2_multistart(obj, bounds, control) else
    .msr_multistart(obj, bounds$starts, bounds$lower, bounds$upper,
    gradient_tolerance = control$gradient_tolerance,
    control = list(iter.max = control$iter.max, eval.max = control$eval.max,
                   rel.tol = 1e-11, x.tol = 1e-9))
  if (is.null(fits$best))
    stop("No analytic GLM start passed convergence and score checks. Try additional starts or inspect the model/design.",
         call. = FALSE)
  best <- fits$best
  fitted <- if (is.function(obj$cold_evaluate)) obj$cold_evaluate(best$par) else
    obj$evaluate(best$par)
  nk <- length(obj$engine$parameter_names)
  kp <- unname(best$par[seq_len(nk)])
  kd <- .msr_eval_kernel_engine(obj$engine, kp)
  info <- .msr_profile_information(obj, best$par, bounds$lower, bounds$upper, fits)
  p <- obj$engine$p
  scale <- exp(kp[seq_len(p)]) * inputs$unit_conv
  scale_se <- if (info$available) scale * sqrt(diag(info$vcov)[seq_len(p)]) else rep(NA_real_, p)
  scale_est <- data.frame(Mean = scale, SE = scale_se, row.names = contract$covariates)
  shape_est <- NULL
  if (inputs$kernel == "expow") {
    shape <- exp(kp[p + seq_len(p)])
    shape_est <- data.frame(Mean = shape,
      SE = if (info$available) shape * sqrt(diag(info$vcov)[p + seq_len(p)]) else NA_real_,
      row.names = contract$covariates)
  }
  opt_mod <- .msr_profile_conditional_model(contract, fitted, model, kp)
  # Keep the ordinary fitted model's dispersion SE conditional on spatial scales.
  theta_interior <- obj$family == "nbinom2" &&
    utils::tail(best$par - bounds$lower, 1L) > 1e-5 &&
    utils::tail(bounds$upper - best$par, 1L) > 1e-5
  if (theta_interior && !info$poisson_compatible && !is.null(info$nuisance_vcov)) {
    opt_mod$SE.theta <- fitted$theta * sqrt(utils::tail(diag(info$nuisance_vcov), 1L))
  }
  profile_hessian <- if (!joint_nb2) stats::optimHess(unname(best$par), obj$fn, obj$gr) else
    matrix(NA_real_, length(best$par), length(best$par))
  profile_vcov <- if (info$available) info$vcov[info$profile_index, info$profile_index, drop = FALSE] else
    matrix(NA_real_, length(best$par), length(best$par))
  # Legacy optimizer slots retain natural kernel units only. Dispersion and
  # all log-parameter information live in the explicitly named profile_fit.
  natural_cov <- if (info$available) diag(exp(kp), nk) %*%
    profile_vcov[seq_len(nk), seq_len(nk), drop = FALSE] %*% diag(exp(kp), nk) else
    matrix(NA_real_, nk, nk)
  natural_hess <- if (info$available) solve(natural_cov) else natural_cov
  to_map <- c(rep(1 / inputs$unit_conv, p),
              if (inputs$kernel == "expow") rep(1, p))
  out <- list(
    backend = backend, scale_est = scale_est, shape_est = shape_est,
    optim_results = list(par = exp(kp), par_unscale = c(scale,
      if (!is.null(shape_est)) shape_est$Mean), value = fitted$value,
      hessian = natural_hess,
      hessian_unscale = natural_hess * outer(to_map, to_map),
      convergence = best$convergence, message = best$message),
    opt_mod = opt_mod, fitted_mod_original = model,
    min_D = inputs$min_D, max_D = inputs$max_D,
    kernel_inputs = inputs[setdiff(names(inputs), c("min_D", "max_D"))],
    scl_params = list(mean = stats::setNames(kd$center, contract$covariates),
                      sd = stats::setNames(kd$spread, contract$covariates)),
    join_by = NULL, opt_context = NULL, profile_scale_est = NULL,
    profile_fit = list(contract = contract, par = unname(best$par),
      parameter_names = obj$parameter_names, lower = bounds$lower, upper = bounds$upper,
      starts = fits$starts, start_results = fits$records, winner = fits$winner,
      hessian = profile_hessian, vcov = profile_vcov, control = control),
    uncertainty = info,
    dispersion_est = if (obj$family == "nbinom2")
      list(family = "NB2", theta = fitted$theta,
           theta_bounds = control$theta_bounds, poisson_limit = info$poisson_limit) else NULL,
    diagnostics = .empty_diagnostics(), next_run = NULL, warn_message = 0, call = call)
  class(out) <- "multiScaleR"
  out <- .msr_attach_fixed_reference(out)
  if (is.function(obj$coefficient_audit))
    out$profile_fit$coefficient_audit <- obj$coefficient_audit()
  if (!is.null(pool)) {
    out$profile_fit$worker_count <- length(pool)
    out$profile_fit$worker_coefficient_audits <- fits$worker_coefficient_audits
  }
  if (joint_nb2) {
    out$dispersion_est$alpha <- fitted$alpha
    out$dispersion_est$exact_poisson <- identical(fitted$alpha, 0)
    out$profile_fit$joint_par <- fitted$joint_par
    # Joint information uses log(theta) for positive alpha. Convert its marginal
    # dispersion covariance to rho for the explicitly named outer-parameter slot.
    if (info$available) {
      jac <- c(rep(1, nk), -fitted$alpha / (1 + fitted$alpha))
      out$profile_fit$vcov <- profile_vcov * outer(jac, jac)
      out$profile_fit$hessian <- solve(out$profile_fit$vcov)
    }
  }
  out$profile_scale_est <- .msr_profile_wald_table(out)
  out$diagnostics$profile_glm <- list(
    triggered = !info$available || info$poisson_limit,
    uncertainty_available = info$available, reasons = info$reasons,
    weak_identifiability = info$weak, boundary = info$boundary,
    poisson_limit = info$poisson_limit, condition_number = info$condition_number,
    poisson_compatible = info$poisson_compatible, poisson_deviance = info$poisson_deviance,
    log_scale_se = info$log_scale_se, multiple_modes = info$multiple_modes,
    eligible_starts = sum(fits$records$eligible), attempted_starts = nrow(fits$records))
  out$diagnostics$max_distance <- .max_distance_diagnostic(
    kernel_dist(out, dimension = "1d"), inputs$max_D)
  if (isTRUE(out$diagnostics$max_distance$triggered)) out$warn_message <- 1L
  prepared_n <- if (is.null(inputs$binned)) length(inputs$d_list) else nrow(inputs$binned$vsum[[1L]])
  out$diagnostics$sample_size <- list(code = "sample_size",
    triggered = obj$engine$n < prepared_n,
    fitted_n = obj$engine$n, prepared_n = prepared_n)
  out$next_run <- .next_run_recommendation(scale_est, shape_est, inputs$max_D, out$diagnostics)
  if (!info$available)
    warning("Joint Wald uncertainty unavailable: ", paste(info$reasons, collapse = "; "),
            ". Inspect diagnostics and likelihood profiles.", call. = FALSE)
  if (verbose) cat(if (joint_nb2) "Analytic joint NB2 fit: " else "Analytic-profile GLM fit: ", sum(fits$records$eligible), "/",
                   nrow(fits$records), " starts passed; family ", obj$family, ".\n", sep = "")
  out
}

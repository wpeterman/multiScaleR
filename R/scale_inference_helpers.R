# Inference refits reuse the prepared design and analytic objective. They do not
# build fitted-model objects or Hessians inside bootstrap replicates.
.msr_inference_fit <- function(obj, starts, lower, upper, control,
                               fixed_scale = NULL) {
  free <- if (is.null(fixed_scale)) seq_along(lower) else seq_along(lower)[-1L]
  expand <- function(z) {
    if (is.null(fixed_scale)) z else c(fixed_scale, z)
  }
  reduced <- list(fn = function(z) obj$fn(expand(z)),
                  gr = function(z) obj$gr(expand(z))[free])
  if (is.function(obj$reset_start)) reduced$reset_start <- obj$reset_start
  if (is.function(obj$cold_evaluate)) reduced$cold_evaluate <- function(par) {
    z <- obj$cold_evaluate(expand(par)); z$gradient <- z$gradient[free]; z
  }
  starts <- unique(starts[, free, drop = FALSE])
  fits <- .msr_multistart(reduced, starts, lower[free], upper[free],
    gradient_tolerance = control$gradient_tolerance,
    control = list(iter.max = control$iter.max, eval.max = control$eval.max,
                   rel.tol = 1e-11, x.tol = 1e-9))
  if (is.null(fits$best)) stop("No inference refit passed convergence and score checks.", call. = FALSE)
  list(par = expand(fits$best$par), value = fits$best$objective,
       records = fits$records)
}

.msr_inference_setup <- function(object, n_starts) {
  if (!inherits(object, "multiScaleR") || !.msr_analytic_backend(object))
    stop("scale_inference requires a profile_glm or joint_nb2 fit.", call. = FALSE)
  pf <- object$profile_fit
  obj <- pf$contract$objective
  if (obj$family != "nbinom2" || obj$engine$p != 1L)
    stop("scale_inference currently supports NB2 with one spatial covariate.", call. = FALSE)
  # Isolate the likelihood cache as well as the returned fit metadata.
  factory <- if (isTRUE(obj$boundary_nb2)) .msr_joint_nb2_objective else .msr_profile_glm_objective
  obj <- factory(obj$engine, obj$y, obj$family,
    obj$X_fixed, obj$offset, obj$weights, obj$trials,
    response_design = obj$response_design)
  if (isTRUE(obj$boundary_nb2) && identical(pf$control$coefficient_start, "warm"))
    obj <- .msr_joint_nb2_objective(obj$engine, obj$y, obj$family,
      obj$X_fixed, obj$offset, obj$weights, obj$trials, coefficient_start = "warm",
      response_design = obj$response_design)
  if (!isTRUE(obj$boundary_nb2) && identical(pf$control$coefficient_start, "warm"))
    obj <- .msr_profile_nb2_continuation(obj)
  validate_scalar_numeric(n_starts, "n_starts", lower = 2, integerish = TRUE)
  control <- pf$control
  control$n_starts <- n_starts
  control$starts <- NULL
  inputs <- object$kernel_inputs
  inputs$min_D <- exp(pf$lower[1L]) * inputs$unit_conv
  inputs$max_D <- exp(pf$upper[1L]) * inputs$unit_conv
  start_builder <- if (isTRUE(obj$boundary_nb2)) .msr_joint_nb2_bounds else .msr_profile_starts
  extra <- start_builder(pf$contract, inputs, NULL,
                               control, object$opt_mod)$starts
  # Include both bounds explicitly; weak effects can favor a boundary basin
  # that dispersed interior starts fail to reach.
  ends <- rbind(pf$par, pf$par)
  ends[, 1L] <- c(pf$lower[1L], pf$upper[1L])
  starts <- unique(rbind(pf$par, pf$starts, extra, ends))
  list(obj = obj, starts = starts, lower = pf$lower, upper = pf$upper,
       control = control, unit = object$kernel_inputs$unit_conv)
}

.msr_inference_refit <- function(setup, obj = setup$obj, starts = setup$starts,
                                 fixed_scale = NULL) {
  .msr_inference_fit(obj, starts, setup$lower, setup$upper, setup$control, fixed_scale)
}

# Remove every response column depending on the sole optimized kernel. The
# remaining X_fixed model makes scale and shape irrelevant, including products
# and factor interactions. It bounds profile deviance; it is not an effect test.
.msr_nb2_zero_effect <- function(setup) {
  obj <- setup$obj
  if (isTRUE(obj$boundary_nb2)) return(.msr_joint_nb2_zero_effect(setup))
  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  evaluate <- function(z) {
    if (identical(z, cache$par)) return(cache$result)
    theta <- exp(z)
    fit <- stats::glm.fit(obj$X_fixed, obj$y, weights = obj$weights,
      offset = obj$offset, family = MASS::negative.binomial(theta),
      control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
    # The reduced model can oscillate too, especially after removing a strong
    # kernel effect. Apply the same checked coefficient-likelihood fallback.
    if (!isTRUE(fit$converged) && fit$rank == ncol(obj$X_fixed))
      fit <- .msr_nb2_refine_glm(fit, obj$X_fixed, obj$y, theta,
                                obj$offset, obj$weights)
    if (!isTRUE(fit$converged) || fit$rank != ncol(obj$X_fixed))
      stop("The zero-effect NB2 fit failed.", call. = FALSE)
    cache$par <- z
    cache$result <- list(value = -sum(obj$weights * stats::dnbinom(
      obj$y, mu = fit$fitted.values, size = theta, log = TRUE)),
      gradient = .msr_nb2_theta_score(obj$y, fit$fitted.values, theta, obj$weights))
    cache$result
  }
  objective <- list(fn = function(z) evaluate(z)$value,
                    gr = function(z) evaluate(z)$gradient)
  fits <- .msr_multistart(objective, unique(setup$starts[, ncol(setup$starts), drop = FALSE]),
    utils::tail(setup$lower, 1L), utils::tail(setup$upper, 1L),
    gradient_tolerance = setup$control$gradient_tolerance,
    control = list(iter.max = setup$control$iter.max, eval.max = setup$control$eval.max))
  if (is.null(fits$best)) stop("No zero-effect fit passed convergence checks.", call. = FALSE)
  fits$best$objective
}

# Preserve separate accepted components. Grid-only sets retain outside brackets
# rather than pretending their outer accepted points are exact endpoints.
.msr_scale_components <- function(x, accepted, refine = NULL) {
  empty <- data.frame(lower = numeric(), upper = numeric(),
    lower_outside = numeric(), upper_outside = numeric(),
    lower_bound_limited = logical(), upper_bound_limited = logical())
  if (anyNA(accepted)) return(NULL)
  if (!any(accepted)) return(empty)
  runs <- rle(accepted)
  ends <- cumsum(runs$lengths)
  begins <- ends - runs$lengths + 1L
  rows <- lapply(which(runs$values), function(k) {
    a <- begins[k]; b <- ends[k]
    lo <- x[a]; hi <- x[b]
    outside_lo <- if (a > 1L) x[a - 1L] else NA_real_
    outside_hi <- if (b < length(x)) x[b + 1L] else NA_real_
    if (!is.null(refine)) {
      if (a > 1L) lo <- refine(outside_lo, lo)
      if (b < length(x)) hi <- refine(hi, outside_hi)
    }
    data.frame(lower = lo, upper = hi, lower_outside = outside_lo,
      upper_outside = outside_hi, lower_bound_limited = a == 1L,
      upper_bound_limited = b == length(x))
  })
  do.call(rbind, rows)
}

.msr_bootstrap_tail <- function(lr, observed, level) {
  if (!length(lr) || any(!is.finite(lr)))
    return(c(p = NA_real_, mc_lower = NA_real_, mc_upper = NA_real_))
  B <- length(lr)
  count <- sum(lr >= observed - 1e-8)
  ci <- stats::binom.test(count, B, conf.level = .95)$conf.int
  c(p = (1 + count) / (B + 1),
    mc_lower = (1 + B * ci[1L]) / (B + 1),
    mc_upper = (1 + B * ci[2L]) / (B + 1))
}

.msr_scale_bootstrap <- function(setup, constrained, reference, B) {
  obj <- setup$obj
  # Bootstrap generators must not depend on the preceding continuation path.
  fitted <- if (is.function(obj$cold_evaluate)) obj$cold_evaluate(constrained$par) else
    obj$evaluate(constrained$par)
  starts <- unique(rbind(constrained$par, reference$par, setup$starts))
  records <- vector("list", B)
  for (b in seq_len(B)) {
    warnings <- character()
    ans <- tryCatch(withCallingHandlers({
      y <- if (is.infinite(fitted$theta)) stats::rpois(obj$engine$n, fitted$fit$fitted.values) else
        stats::rnbinom(obj$engine$n, mu = fitted$fit$fitted.values, size = fitted$theta)
      factory <- if (isTRUE(obj$boundary_nb2)) .msr_joint_nb2_objective else .msr_profile_glm_objective
      boot_obj <- factory(obj$engine, y, "nbinom2",
        obj$X_fixed, obj$offset, obj$weights, obj$trials,
    response_design = obj$response_design)
      if (isTRUE(obj$boundary_nb2) && identical(obj$coefficient_start, "warm"))
        boot_obj <- .msr_joint_nb2_objective(obj$engine, y, "nbinom2",
          obj$X_fixed, obj$offset, obj$weights, obj$trials, coefficient_start = "warm",
      response_design = obj$response_design)
      if (!isTRUE(obj$boundary_nb2) && identical(obj$coefficient_start, "warm"))
        boot_obj <- .msr_profile_nb2_continuation(boot_obj)
      null <- .msr_inference_refit(setup, boot_obj, starts, constrained$par[1L])
      full <- .msr_inference_refit(setup, boot_obj, rbind(null$par, starts))
      lr <- 2 * (null$value - full$value)
      if (lr < -2e-5) stop("Unrestricted bootstrap fit is worse than its constrained fit.")
      data.frame(replicate = b, ok = TRUE, lr = max(0, lr),
        sigma = exp(full$par[1L]) * setup$unit,
        theta = if (isTRUE(obj$boundary_nb2)) 1 / expm1(utils::tail(full$par, 1L)) else
          exp(utils::tail(full$par, 1L)),
        boundary = any(full$par - setup$lower < 1e-5 | setup$upper - full$par < 1e-5),
        full_nll = full$value, constrained_nll = null$value,
        starts_attempted = nrow(full$records) + nrow(null$records),
        starts_eligible = sum(full$records$eligible) + sum(null$records$eligible), error = "")
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) data.frame(replicate = b, ok = FALSE, lr = NA_real_,
      sigma = NA_real_, theta = NA_real_, boundary = NA,
      full_nll = NA_real_, constrained_nll = NA_real_,
      starts_attempted = NA_integer_, starts_eligible = NA_integer_,
      error = conditionMessage(e)))
    ans$warnings <- paste(unique(warnings), collapse = "; ")
    records[[b]] <- ans
  }
  do.call(rbind, records)
}

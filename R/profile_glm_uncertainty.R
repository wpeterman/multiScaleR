# Reconstruct the full likelihood at arbitrary coefficients. Differentiating
# its analytic score gives observed information including scale, coefficient,
# and dispersion cross derivatives. No conditional model fitting occurs here.
.msr_profile_joint <- function(obj, profile_par) {
  nk <- length(obj$engine$parameter_names)
  fitted <- obj$evaluate(profile_par)
  nb <- length(fitted$beta)
  ib <- nk + seq_len(nb)
  disp <- obj$family %in% c("gaussian", "nbinom2")
  par <- c(profile_par[seq_len(nk)], unname(fitted$beta),
    if (disp) log(if (obj$family == "nbinom2") fitted$theta else fitted$residual_sd))
  names(par) <- c(obj$engine$parameter_names, paste0("beta:", colnames(fitted$design)),
    if (disp) if (obj$family == "nbinom2") "log_theta" else "log_residual_sd")
  if (isTRUE(obj$boundary_nb2)) {
    evaluate <- function(par) {
      alpha <- exp(-utils::tail(par, 1L))
      ans <- obj$joint$evaluate(c(par[-length(par)], log1p(alpha)))
      ans$gradient[length(par)] <- ans$gradient[length(par)] * (-alpha / (1 + alpha))
      ans
    }
    return(list(par = par, beta_index = ib,
      profile_index = c(seq_len(nk), length(par)),
      fn = function(par) evaluate(par)$value, gr = function(par) evaluate(par)$gradient))
  }
  evaluate <- function(par) {
    kd <- .msr_eval_kernel_engine(obj$engine, unname(par[seq_len(nk)]))
    X <- if (is.null(obj$response_design)) cbind(obj$X_fixed, kd$value) else
      .msr_response_design_matrix(obj$response_design, kd$value)
    beta <- par[ib]
    eta <- drop(X %*% beta) + obj$offset
    y <- obj$y
    w <- obj$weights
    mu <- switch(obj$family, gaussian = eta, binomial = stats::plogis(eta), exp(eta))
    d <- if (disp) exp(utils::tail(par, 1L)) else NA_real_
    if (obj$family == "gaussian") {
      ll <- stats::dnorm(y, mu, d, log = TRUE)
      score <- w * (mu - y) / d^2
      ds <- sum(w * (1 - (y - mu)^2 / d^2))
    } else if (obj$family == "nbinom2") {
      ll <- stats::dnbinom(y, mu = mu, size = d, log = TRUE)
      score <- w * (mu - y) / (1 + mu / d)
      ds <- .msr_nb2_theta_score(y, mu, d, w)
    } else if (obj$family == "poisson") {
      ll <- stats::dpois(y, mu, log = TRUE)
      score <- w * (mu - y)
    } else {
      ll <- stats::dbinom(y, obj$trials, mu, log = TRUE)
      score <- w * (obj$trials * mu - y)
    }
    if (is.null(obj$response_design)) {
      bk <- utils::tail(beta, obj$engine$p)
      gk <- colSums(kd$d_scale * score) * bk
      if (obj$engine$kernel == "expow")
        gk <- c(gk, colSums(kd$d_shape * score) * bk)
    } else {
      sensitivity <- .msr_response_kernel_sensitivity(obj$response_design, kd$value, beta)
      gk <- colSums(kd$d_scale * sensitivity * score)
      if (obj$engine$kernel == "expow")
        gk <- c(gk, colSums(kd$d_shape * sensitivity * score))
    }
    list(value = -sum(w * ll),
         gradient = c(gk, drop(crossprod(X, score)), if (disp) ds))
  }
  list(par = par, beta_index = ib,
       profile_index = c(seq_len(nk), if (obj$family == "nbinom2") length(par)),
       fn = function(par) evaluate(par)$value,
       gr = function(par) evaluate(par)$gradient)
}

.msr_information_inverse <- function(H) {
  if (any(!is.finite(H)) || any(diag(H) <= 0))
    return(list(vcov = NULL, condition = Inf))
  scaled <- H / sqrt(outer(diag(H), diag(H)))
  ev <- eigen(scaled, symmetric = TRUE, only.values = TRUE)$values
  condition <- if (min(ev) <= 0) Inf else max(ev) / min(ev)
  list(vcov = if (is.finite(condition) && condition < 1e8) solve(H) else NULL,
       condition = condition)
}

# Joint covariance may condition on kernel parameters at a search bound when
# the bound is the only failed check (see .msr_profile_bound_joint()).
.msr_profile_information <- function(obj, par, lower, upper, fits) {
  if (isTRUE(obj$boundary_nb2) && utils::tail(par, 1L) == 0)
    return(.msr_poisson_boundary_information(obj, par, fits))
  joint <- .msr_profile_joint(obj, par)
  H <- stats::optimHess(joint$par, joint$fn, joint$gr)
  dimnames(H) <- list(names(joint$par), names(joint$par))
  inv <- .msr_information_inverse(H)
  nk <- length(obj$engine$parameter_names)
  nuisance <- .msr_information_inverse(H[-seq_len(nk), -seq_len(nk), drop = FALSE])$vcov
  boundary <- any(par - lower < 1e-5 | upper - par < 1e-5)
  logse <- if (is.null(inv$vcov)) rep(NA_real_, nk) else sqrt(diag(inv$vcov)[seq_len(nk)])
  good <- which(fits$records$eligible &
                  fits$records$value <= fits$best$objective + stats::qchisq(.95, 1) / 2)
  multiple <- any(vapply(good, function(i)
    max(abs(fits$fits[[i]]$par[seq_len(nk)] - par[seq_len(nk)])) > 0.5, logical(1)))
  weak <- is.null(inv$vcov) || any(logse > 1, na.rm = TRUE) || multiple
  fitted <- obj$evaluate(par)
  poisson <- obj$family == "nbinom2" &&
    max(fitted$fit$fitted.values / fitted$theta) < 0.01
  poisson_deviance <- NA_real_
  poisson_compatible <- FALSE
  if (obj$family == "nbinom2") {
    # A Poisson fit at the selected scales is an attainable limiting model.
    # This is a one-sided diagnostic: optimizing its scales could improve it.
    pois <- stats::glm.fit(fitted$design, obj$y, weights = obj$weights,
      offset = obj$offset, family = stats::poisson(),
      control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
    if (isTRUE(pois$converged)) {
      poisson_deviance <- 2 * (-sum(obj$weights * stats::dpois(
        obj$y, pois$fitted.values, log = TRUE)) - fitted$value)
      poisson_compatible <- poisson_deviance <= stats::qchisq(.95, 1)
    }
  }
  weak <- weak || poisson_compatible
  scale <- .msr_profile_scale_se(H, nk, par, lower, upper, fits,
    global = poisson_compatible)
  reasons <- c(if (boundary) "parameter at an optimization bound",
    if (is.null(inv$vcov)) "nonpositive or ill-conditioned observed information",
    if (any(logse > 1, na.rm = TRUE)) "log-kernel standard error exceeds one",
    if (multiple) "competing scale modes within the 95% likelihood threshold",
    if (poisson_compatible) "Poisson limit within the diagnostic likelihood threshold")
  available <- !boundary && !weak
  V <- inv$vcov
  held <- rep(FALSE, nk)
  if (!available) {
    conditional <- .msr_profile_bound_joint(H, scale, boundary, multiple, poisson_compatible)
    if (!is.null(conditional)) {
      V <- conditional$vcov
      held <- conditional$held
    }
  }
  if (is.null(V)) V <- matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H))
  reporting_names <- sub("^log_shape:", "shape:", sub("^log_scale:", "sigma:",
    obj$engine$parameter_names))
  list(available = available, joint_available = available || any(held),
       joint_conditional_on = reporting_names[held], vcov = V, hessian = H,
       par = joint$par, beta_index = joint$beta_index,
       profile_index = joint$profile_index, nuisance_vcov = nuisance,
       conditional_vcov = if (!is.null(nuisance))
         nuisance[seq_along(joint$beta_index), seq_along(joint$beta_index), drop = FALSE],
       boundary = boundary, weak = weak, multiple_modes = multiple,
       poisson_limit = poisson, poisson_compatible = poisson_compatible,
       poisson_deviance = poisson_deviance, condition_number = inv$condition,
       log_scale_se = stats::setNames(logse, obj$engine$parameter_names),
       scale_log_se = stats::setNames(scale$log_se, obj$engine$parameter_names),
       scale_conditional_on = reporting_names[scale$conditional_on],
       reasons = reasons)
}

# Per-parameter Wald SEs for the reported scale table. The joint gate above
# withholds every SE when any check fails, and joint prediction uncertainty
# unless the bound is the only failure (see .msr_profile_bound_joint()).
# Here a kernel parameter loses only its own SE when it sits at a search bound,
# its log-scale SE exceeds one, or competing modes within the 95 percent
# likelihood threshold move it by more than 0.5 log units. Bounded kernel
# parameters are dropped from the observed information, so the remaining SEs
# hold them fixed at their bounds. A dispersion bound, a Poisson-compatible NB2
# fit, or unusable information still withholds every scale SE. `vcov` is the
# inverse of the reduced information over the rows in `keep` (free kernel
# parameters, then every coefficient and dispersion row), or NULL.
.msr_profile_scale_se <- function(H, nk, par, lower, upper, fits, global = FALSE) {
  at_bound <- par - lower < 1e-5 | upper - par < 1e-5
  kernel_bound <- at_bound[seq_len(nk)]
  # `conditional_on` flags the bounded parameters held fixed, only when some
  # other parameter keeps an SE.
  out <- list(log_se = rep(NA_real_, nk), conditional_on = rep(FALSE, nk),
              vcov = NULL, keep = integer())
  if (global || any(at_bound[-seq_len(nk)]) || all(kernel_bound)) return(out)
  keep <- c(which(!kernel_bound), seq.int(nk + 1L, length.out = nrow(H) - nk))
  vcov <- .msr_information_inverse(H[keep, keep, drop = FALSE])$vcov
  if (is.null(vcov)) return(out)
  out$vcov <- vcov
  out$keep <- keep
  free <- which(!kernel_bound)
  log_se <- rep(NA_real_, nk)
  log_se[free] <- sqrt(diag(vcov)[seq_along(free)])
  good <- which(fits$records$eligible &
                  fits$records$value <= fits$best$objective + stats::qchisq(.95, 1) / 2)
  moved <- vapply(seq_len(nk), function(k) any(vapply(good, function(i)
    abs(fits$fits[[i]]$par[k] - par[k]) > 0.5, logical(1))), logical(1))
  log_se[!is.finite(log_se) | log_se > 1 | moved] <- NA_real_
  out$log_se <- log_se
  if (any(is.finite(log_se))) out$conditional_on <- kernel_bound
  out
}

# Joint covariance conditional on kernel parameters at a search bound
# (profile_glm and joint_nb2). `scale` is the result of
# .msr_profile_scale_se(). It replaces the all-or-nothing joint gate only
# when the bound is the sole failure: the
# fit is at a bound, no dispersion parameter is bounded, at least one kernel
# parameter is free, the reduced information is usable, every free log-kernel
# SE is at most one, no competing modes lie within the 95 percent likelihood
# threshold, and no Poisson fit is compatible with the NB2 fit. Otherwise NULL.
# The full-size matrix keeps exact zeros in the bounded kernel rows and
# columns, so positional delta-method predictions hold those parameters fixed,
# as conditional predictions hold every kernel parameter fixed. Its free block
# is the same reduced inverse that gives the per-parameter scale SEs.
.msr_profile_bound_joint <- function(H, scale, boundary, multiple, poisson_compatible) {
  held <- scale$conditional_on
  if (!boundary || multiple || poisson_compatible || !any(held) ||
      is.null(scale$vcov) || !all(is.finite(scale$log_se[!held])))
    return(NULL)
  V <- matrix(0, nrow(H), ncol(H), dimnames = dimnames(H))
  V[scale$keep, scale$keep] <- scale$vcov
  list(vcov = V, held = held)
}

.msr_profile_wald_table <- function(object, shape = FALSE) {
  tab <- if (shape) object$shape_est else object$scale_est
  # Limits follow the reported per-parameter SE (SE / Mean is the log-scale
  # SE), so a withheld SE gives missing limits for that parameter only.
  se <- tab$SE / tab$Mean
  ok <- is.finite(se)
  tab[["2.5%"]] <- ifelse(ok, tab$Mean * exp(-stats::qnorm(.975) * se), NA_real_)
  tab[["97.5%"]] <- ifelse(ok, tab$Mean * exp(stats::qnorm(.975) * se), NA_real_)
  attr(tab, "interval_method") <- ifelse(ok, "log_wald", "unavailable")
  tab
}

#' Covariance matrices for an analytic multiscale fit
#'
#' @description
#' Return the covariance matrix of the estimated parameters from a
#' [multiScale_optim()] fit made with an analytic backend. The joint matrix
#' covers the kernel parameters and the regression coefficients together, and
#' it is what [predict.multiScaleR()] uses for joint uncertainty. The
#' conditional matrix covers the coefficients with the kernel parameters held
#' at their estimates. [multiScaleR-prediction] explains when each kind of
#' uncertainty applies.
#'
#' @param object A `multiScaleR` fit from `backend = "profile_glm"` or
#'   `backend = "joint_nb2"`. A `backend = "joint_mixed"` fit returns only
#'   native conditional covariance; its scale standard errors are in
#'   `object$scale_est`. Other backends raise an error. Use
#'   `vcov(object$opt_mod)` for their coefficient covariance conditional on the
#'   selected scales.
#' @param type `"joint"` (the default) or `"conditional"`.
#'   * `"joint"` returns the observed-information covariance of the log kernel
#'     parameters, the regression coefficients, and the log dispersion when the
#'     family has one.
#'   * `"conditional"` returns the coefficient covariance with the kernel
#'     parameters held fixed and the dispersion estimated. At the exact Poisson
#'     boundary of a `"joint_nb2"` fit it also holds the dispersion at that
#'     boundary, and joint covariance is unavailable there.
#' @param ... Passed to native `vcov()` for `"joint_mixed"`. Otherwise unused.
#'
#' @details
#' # Read the joint matrix
#'
#' Rows and columns are named by parameter:
#' * `log_scale:<covariate>`: `log(sigma / u)`, the natural log of that
#'   covariate's sigma after dividing it by the internal distance unit
#'   `u = object$kernel_inputs$unit_conv`. Dividing by `u` shifts the log by a
#'   constant, so the variance equals the variance of log(sigma) in map units.
#' * `log_shape:<covariate>`: the log shape parameter, for `kernel = "expow"`.
#' * `beta:<term>`: a regression coefficient on the link scale, named after
#'   its model-matrix column, for example `beta:(Intercept)`. A standardized
#'   kernel-weighted covariate's coefficient describes one standard deviation
#'   of that covariate.
#' * A final dispersion row when the family has one: `log_theta` for a
#'   negative binomial model and `log_residual_sd` for a Gaussian model.
#'
#' The square root of a diagonal element is a standard error on that
#' parameter's scale. For sigma, a 95% Wald interval in map units is
#' `sigma x exp(-1.96 x SE)` to `sigma x exp(1.96 x SE)`, where `SE` is the
#' standard error of `log_scale`. These are the limits that `summary(object)`
#' reports, and `object$scale_est$SE` is `sigma x SE`. Off-diagonal elements
#' show how the estimates trade off, and [stats::cov2cor()] turns them into
#' correlations.
#'
#' Joint covariance is unavailable for fits at a parameter bound or when the
#' identifiability checks fail. In that case `type = "joint"` raises an error
#' that lists `object$uncertainty$reasons`, rather than reporting a misleading
#' Wald result. The conditional matrix remains available.
#'
#' One exception applies to both backends. When kernel parameters
#' at a search bound are the only failed check, the joint matrix holds those
#' parameters fixed at their bounds. Their rows and columns are zero, and the
#' other rows and columns come from the observed information with the bounded
#' parameters removed. That reduced information must pass the remaining
#' checks. The matrix then carries a `conditional_on` attribute that names the
#' bounded parameters, as does `object$uncertainty$joint_conditional_on`. A
#' zero there means the parameter was held fixed, not that it is known
#' exactly, and its standard error is missing in `object$scale_est`. The other
#' standard errors are no wider than they would be if the bounded scale were
#' free to vary, so treat them as conditional on that scale. A dispersion
#' bound, every kernel parameter at a bound, or any other failed check still
#' withholds joint covariance.
#'
#' # Read the conditional matrix
#'
#' The conditional matrix holds only the `beta:` rows. It answers how
#' uncertain the coefficients are if the scales of effect are exactly their
#' estimates, which is also what `vcov(object$opt_mod)` and
#' `summary(object$opt_mod)` describe.
#'
#' @return A numeric covariance matrix with named rows and columns, as
#'   described in Details. For `"joint_mixed"` with `type = "conditional"`,
#'   the native glmmTMB covariance list, one matrix per model component,
#'   conditional on the selected scales.
#' @seealso [predict.multiScaleR()] for predictions that use these matrices,
#'   [multiScaleR-prediction], [multiScale_optim()], and
#'   `vignette("v09-estimates-and-intervals", package = "multiScaleR")`.
#'
#' @examples
#' ## Running example: counts at 100 survey points, two landscape layers
#' data("landscape_counts")
#' data("surv_pts")
#' pts <- terra::vect(surv_pts)
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))[[c("land1", "land2")]]
#' inputs <- kernel_prep(pts, land, max_D = 1700, kernel = "gaussian",
#'                       bin = TRUE, verbose = FALSE)
#' dat <- data.frame(landscape_counts, inputs$kernel_dat)
#' start <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
#' fit <- multiScale_optim(start, inputs, backend = "profile_glm",
#'                         n_cores = 1, verbose = FALSE)
#'
#' ## Joint covariance of the two log scales and four coefficients
#' V <- vcov(fit)
#' round(V, 4)
#'
#' ## 95% Wald limits for sigma, in meters, from the log-scale standard errors.
#' ## They match the limits in summary(fit).
#' se_log <- sqrt(diag(V))[c("log_scale:land1", "log_scale:land2")]
#' sigma <- stats::setNames(fit$scale_est$Mean, rownames(fit$scale_est))
#' round(cbind(sigma, lower = sigma * exp(-qnorm(0.975) * se_log),
#'             upper = sigma * exp(qnorm(0.975) * se_log)))
#'
#' ## Correlations show how the scale and coefficient estimates trade off
#' round(cov2cor(V), 2)
#'
#' ## Coefficient standard errors with and without scale uncertainty
#' round(cbind(joint = sqrt(diag(V))[-(1:2)],
#'             conditional = sqrt(diag(vcov(fit, type = "conditional")))), 3)
#' @exportS3Method stats::vcov
vcov.multiScaleR <- function(object, type = c("joint", "conditional"), ...) {
  type <- match.arg(type)
  if (identical(object$backend, "joint_mixed")) {
    if (type != "conditional")
      stop("Joint covariance is unavailable for joint_mixed. Scale SEs are in ",
           "object$scale_est; use type = 'conditional' for coefficient covariance.", call. = FALSE)
    return(stats::vcov(object$opt_mod, ...))
  }
  if (!.msr_analytic_backend(object))
    stop("This method requires backend = 'profile_glm' or 'joint_nb2'; use vcov(object$opt_mod) for conditional model covariance.",
         call. = FALSE)
  if (type == "joint") {
    u <- object$uncertainty
    # Older saved fits have no `joint_available`; their gate is `available`.
    # Without `available`, joint covariance needs named held parameters.
    if (!isTRUE(u$available) &&
        !(isTRUE(u$joint_available) && length(u$joint_conditional_on)))
      stop("Joint Wald uncertainty unavailable: ", paste(u$reasons, collapse = "; "),
           call. = FALSE)
    V <- u$vcov
    # Name the bounded parameters whose zero rows and columns hold them fixed.
    if (!isTRUE(u$available) && length(u$joint_conditional_on))
      attr(V, "conditional_on") <- u$joint_conditional_on
    return(V)
  }
  V <- object$uncertainty$conditional_vcov
  if (is.null(V)) stop("Conditional coefficient covariance is unavailable.", call. = FALSE)
  V
}

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
  reasons <- c(if (boundary) "parameter at an optimization bound",
    if (is.null(inv$vcov)) "nonpositive or ill-conditioned observed information",
    if (any(logse > 1, na.rm = TRUE)) "log-kernel standard error exceeds one",
    if (multiple) "competing scale modes within the 95% likelihood threshold",
    if (poisson_compatible) "Poisson limit within the diagnostic likelihood threshold")
  V <- inv$vcov
  if (is.null(V)) V <- matrix(NA_real_, nrow(H), ncol(H), dimnames = dimnames(H))
  list(available = !boundary && !weak, vcov = V, hessian = H,
       par = joint$par, beta_index = joint$beta_index,
       profile_index = joint$profile_index, nuisance_vcov = nuisance,
       conditional_vcov = if (!is.null(nuisance))
         nuisance[seq_along(joint$beta_index), seq_along(joint$beta_index), drop = FALSE],
       boundary = boundary, weak = weak, multiple_modes = multiple,
       poisson_limit = poisson, poisson_compatible = poisson_compatible,
       poisson_deviance = poisson_deviance, condition_number = inv$condition,
       log_scale_se = stats::setNames(logse, obj$engine$parameter_names),
       reasons = reasons)
}

.msr_profile_wald_table <- function(object, shape = FALSE) {
  p <- nrow(object$scale_est)
  idx <- if (shape) p + seq_len(p) else seq_len(p)
  tab <- if (shape) object$shape_est else object$scale_est
  se <- sqrt(diag(object$uncertainty$vcov)[idx])
  tab[["2.5%"]] <- if (object$uncertainty$available)
    tab$Mean * exp(-stats::qnorm(.975) * se) else NA_real_
  tab[["97.5%"]] <- if (object$uncertainty$available)
    tab$Mean * exp(stats::qnorm(.975) * se) else NA_real_
  attr(tab, "interval_method") <- rep(if (object$uncertainty$available)
    "log_wald" else "unavailable", p)
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
#'   native conditional covariance. Other backends raise an error. Use
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
      stop("joint_mixed scale covariance is unavailable; use type = 'conditional'.", call. = FALSE)
    return(stats::vcov(object$opt_mod, ...))
  }
  if (!.msr_analytic_backend(object))
    stop("This method requires backend = 'profile_glm' or 'joint_nb2'; use vcov(object$opt_mod) for conditional model covariance.",
         call. = FALSE)
  if (type == "joint") {
    if (!object$uncertainty$available)
      stop("Joint Wald uncertainty unavailable: ", paste(object$uncertainty$reasons, collapse = "; "),
           call. = FALSE)
    return(object$uncertainty$vcov)
  }
  V <- object$uncertainty$conditional_vcov
  if (is.null(V)) stop("Conditional coefficient covariance is unavailable.", call. = FALSE)
  V
}

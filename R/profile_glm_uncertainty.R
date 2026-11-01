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

#' Covariance for an analytic multiscale GLM
#'
#' @param object A \code{multiScaleR} fit with \code{backend = "profile_glm"}
#'   or \code{backend = "joint_nb2"}. For \code{backend = "joint_mixed"},
#'   only conditional native covariance is available.
#' @param type \code{"joint"} returns observed-information covariance of log
#'   kernel parameters, regression coefficients, and log dispersion (when
#'   present). \code{"conditional"} returns coefficient covariance with kernel
#'   parameters held fixed and dispersion estimated. At an exact Poisson
#'   boundary it also holds dispersion at that boundary; joint covariance
#'   remains unavailable.
#' @param ... Passed to native \code{vcov()} for \code{joint_mixed}; otherwise reserved.
#' @return A named covariance matrix. Joint covariance is unavailable at
#'   parameter bounds or when identifiability diagnostics fail; in that case
#'   this method raises an error rather than reporting a misleading Wald result.
#'   For \code{joint_mixed} with \code{type = "conditional"}, returns the native
#'   glmmTMB covariance list by model component, conditional on selected scales.
#' @seealso \code{\link{multiScale_optim}}, \code{\link{predict.multiScaleR}}
#' @export
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

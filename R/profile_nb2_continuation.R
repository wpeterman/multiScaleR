# Finite-theta adapter: reuse coefficient continuation without changing the
# profile backend's dispersion coordinates or exact-boundary contract.
.msr_profile_nb2_continuation <- function(reference) {
  if (reference$family != "nbinom2")
    stop("Warm coefficient starts require an NB2 model.", call. = FALSE)
  warm <- .msr_joint_nb2_objective(reference$engine, reference$y, reference$family,
    reference$X_fixed, reference$offset, reference$weights, reference$trials,
    coefficient_start = "warm")
  evaluate <- function(par) {
    if (!is.numeric(par) || length(par) != length(reference$parameter_names) ||
        any(!is.finite(par)) ||
        (!is.null(names(par)) && !identical(names(par), reference$parameter_names)))
      stop("Invalid profile parameter vector.", call. = FALSE)
    theta <- exp(utils::tail(par, 1L))
    if (!is.finite(theta) || theta <= 0)
      stop("NB2 theta must be positive and finite.", call. = FALSE)
    alpha <- 1 / theta
    z <- warm$evaluate(c(unname(par[-length(par)]), log1p(alpha)))
    z$gradient[length(par)] <- z$gradient[length(par)] * (-alpha / (1 + alpha))
    z$gradient <- unname(z$gradient)
    z$theta <- theta
    z$residual_sd <- NA_real_
    z$joint_par <- NULL
    z
  }
  obj <- reference
  obj$evaluate <- evaluate
  obj$fn <- function(par) evaluate(par)$value
  obj$gr <- function(par) evaluate(par)$gradient
  obj$reset_start <- warm$reset_start
  # This reference has independent coefficient initialization and caches.
  # Return its complete GLM metadata for selected-fit reconstruction.
  obj$cold_evaluate <- function(par) {
    warm$reset_start()
    reference$evaluate(par)
  }
  obj$coefficient_audit <- warm$coefficient_audit
  obj$coefficient_start <- "warm"
  # None of the evaluators needs the constructor's output list itself.
  on.exit(rm("obj", envir = environment()), add = TRUE)
  obj
}

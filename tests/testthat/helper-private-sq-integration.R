# These tests integrate one prepared group at fixed population parameters.
# They never call an unsupported fitting interface or bypass its capability guard.
sq_group_native_response <- function(g, p) {
  eta <- private_sq_reference(g, p, response = FALSE)$eta
  if (length(g$response$random)) {
    index <- g$response$random[[1L]]$index
    covariance <- diag(exp(2 * p$log_dispersion), length(eta)) +
      outer(index, index, "==") * exp(2 * p$log_response_sd)
    ch <- chol(covariance)
    return(.5 * (length(eta) * log(2 * pi) + 2 * sum(log(diag(ch))) +
      sum(backsolve(ch, g$y - eta, transpose = TRUE)^2)))
  }
  -sum(switch(g$family,
    gaussian = dnorm(g$y, eta, exp(p$log_dispersion), log = TRUE),
    poisson = dpois(g$y, exp(eta), log = TRUE),
    binomial = dbinom(g$y, g$trials, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(g$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE)))
}

# Independent cold integration uses a single unbounded interval and a density
# value at zero for centering. No production response/integrand routine is used.
sq_group_native_integral <- function(g, p) {
  log_density <- function(z) {
    q <- p; q$scale_z <- z
    -sq_group_native_response(g, q) + dnorm(z, log = TRUE)
  }
  shift <- log_density(0)
  bound <- if (g$family == "gaussian")
    -length(g$y) * (p$log_dispersion + .5 * log(2 * pi)) else 0
  value <- integrate(function(z) vapply(z, function(v) {
    # Skip only exact floating-point zero from a valid density upper bound.
    # This avoids 0*Inf arithmetic; neither normal tail is truncated.
    if (exp(bound + dnorm(v, log = TRUE) - shift) == 0) return(0)
    exp(log_density(v) - shift)
  }, 0), -Inf, Inf, rel.tol = 1e-11, abs.tol = 0,
    subdivisions = 1000L, stop.on.error = TRUE)$value
  -log(value) - shift
}


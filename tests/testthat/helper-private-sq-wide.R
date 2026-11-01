# Independent log-weight reference for existing tiny fixtures. It never
# constructs sigma or sigma^2, and never calls a production response leaf.
wide_sq_raw <- function(buffers, log_scale, kernel) {
  vapply(seq_along(buffers), function(i) {
    b <- buffers[[i]]; nearest <- min(b$d)
    gap <- if (kernel == "gaussian")
      (b$d - nearest) * (b$d + nearest) / 2 else b$d - nearest
    positive <- gap > 0
    weight <- rep(1, length(gap))
    power <- if (kernel == "gaussian") 2 else 1
    weight[positive] <- exp(-exp(log(gap[positive]) - power * log_scale[i]))
    average <- sum(weight * b$v) / sum(weight)
    sqrt(sum(weight * (b$v - average)^2) / sum(weight))
  }, 0)
}

wide_sq_reference <- function(s, p) {
  sm <- s$scales[[1L]]
  log_scale <- drop(sm$X %*% p$scale_beta[sm$beta_index])
  for (term in sm$random) log_scale <- log_scale +
    exp(p$log_scale_sd[term$sd_index]) * p$scale_z[term$z_index][term$index]
  raw <- wide_sq_raw(s$surface$buffers[[1L]], log_scale, s$kernel)
  kernel <- (raw - s$center[[1L]]) / s$spread[[1L]]
  eta <- drop(s$X_fixed %*% p$beta[seq_len(ncol(s$X_fixed))]) + s$offset +
    p$beta[ncol(s$X_fixed) + 1L] * kernel
  list(raw = raw, kernel = kernel, eta = eta)
}

# Each value is a normalized GROUP density, including integrated Gaussian
# response covariance. Population priors are outside this likelihood target.
wide_sq_nll <- function(g, p) {
  eta <- wide_sq_reference(g, p)$eta
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

wide_sq_integral <- function(g, p) {
  log_density <- function(z) {
    q <- p; q$scale_z <- z
    -wide_sq_nll(g, q) + dnorm(z, log = TRUE)
  }
  shift <- log_density(0)
  bound <- if (g$family == "gaussian")
    -length(g$y) * (p$log_dispersion + .5 * log(2 * pi)) else 0
  result <- integrate(function(z) vapply(z, function(t) {
    if (exp(bound + dnorm(t, log = TRUE) - shift) == 0) return(0)
    exp(log_density(t) - shift)
  }, 0), -Inf, Inf, rel.tol = 1e-11, abs.tol = 0,
    subdivisions = 1000L, stop.on.error = TRUE)
  -log(result$value) - shift
}


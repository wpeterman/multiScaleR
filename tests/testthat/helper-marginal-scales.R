# Independent adaptive integration of the annual effects at fixed outer
# parameters. Shared response effects would couple years and invalidate this.
# Evaluate kernel weights in log coordinates so integrator tail probes cannot
# form 0 * Inf at a bin with zero shifted distance. This is an independent
# ordinary-double reference, not a change to either production objective.
marginal_reference_nll <- function(spec, parameters) {
  sm <- spec$scales[[1L]]
  log_s <- drop(sm$X %*% parameters$scale_beta[sm$beta_index])
  for (term in sm$random)
    log_s <- log_s + exp(parameters$log_scale_sd[term$sd_index]) *
      parameters$scale_z[term$z_index][term$index]
  power <- if (spec$kernel == "gaussian") 2 else 1
  log_exponent <- log(spec$delta[[1L]]) - power * log_s -
    if (power == 2) log(2) else 0
  weights <- exp(-exp(pmin(log_exponent, 700)))
  kernel <- rowSums(weights * spec$V[[1L]]) / rowSums(weights * spec$C[[1L]])
  kernel <- (kernel - spec$center[[1L]]) / spec$spread[[1L]]
  eta <- drop(spec$X_fixed %*% parameters$beta[seq_len(ncol(spec$X_fixed))]) +
    spec$offset + parameters$beta[ncol(spec$X_fixed) + 1L] * kernel
  -switch(spec$family,
    gaussian = dnorm(spec$y, eta, exp(parameters$log_dispersion), log = TRUE),
    poisson = dpois(spec$y, exp(eta), log = TRUE),
    binomial = dbinom(spec$y, spec$trials, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(spec$y, mu = exp(eta), size = exp(parameters$log_dispersion), log = TRUE))
}

marginal_annual_integral <- function(spec, parameters, rel.tol = 1e-9) {
  stopifnot(!length(spec$response$random), length(spec$scales[[1L]]$random) == 1L)
  term <- spec$scales[[1L]]$random[[1L]]
  values <- vapply(seq_along(term$levels), function(group) {
    rows <- which(term$index == group)
    log_density <- function(z) {
      p <- parameters; p$scale_z[term$z_index[group]] <- z
      -sum(marginal_reference_nll(spec, p)[rows]) + dnorm(z, log = TRUE)
    }
    modes <- lapply(c(-2, 0, 2), function(z) stats::nlminb(z,
      function(x) -log_density(x), control = list(rel.tol = 1e-12)))
    shift <- -min(vapply(modes, `[[`, 0, "objective"))
    cuts <- c(-Inf, sort(unique(vapply(modes, function(x) as.numeric(x$par), 0))), Inf)
    integrand <- function(z) vapply(z, function(x) exp(log_density(x) - shift), 0)
    integral <- sum(vapply(seq_len(length(cuts) - 1L), function(i)
      stats::integrate(integrand, cuts[i], cuts[i + 1L], rel.tol = rel.tol,
        abs.tol = 1e-12, subdivisions = 1000L)$value, 0))
    shift + log(integral)
  }, 0)
  -sum(values)
}

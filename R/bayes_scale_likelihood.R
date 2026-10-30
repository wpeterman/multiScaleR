# This calculation works with ordinary doubles and with RTMB AD vectors.
# Scales remain explicit parameters: full Bayesian sampling integrates over
# their uncertainty instead of replacing it with a conditional optimum.
.msr_bayes_components <- function(spec, par, include_response = TRUE) {
  sum_rows <- if (inherits(par$beta, "advector")) RTMB::rowSums else base::rowSums
  scales <- kernel <- vector("list", length(spec$covariates))
  names(scales) <- names(kernel) <- spec$covariates
  eta <- drop(spec$X_fixed %*% par$beta[seq_len(ncol(spec$X_fixed))]) + spec$offset
  for (j in seq_along(spec$covariates)) {
    nm <- spec$covariates[j]
    sm <- spec$scales[[nm]]
    log_s <- drop(sm$X %*% par$scale_beta[sm$beta_index])
    for (term in sm$random)
      if (is.null(term$gaussian_block))
        log_s <- log_s + exp(par$log_scale_sd[term$sd_index]) *
          par$scale_z[term$z_index][term$index] else
        log_s <- log_s + .msr_gaussian_effect(term$gaussian_block,
          par$scale_z[term$z_index], par$log_scale_sd[term$sd_index])
    scales[[nm]] <- exp(log_s) * spec$unit_conv
    power <- if (spec$kernel == "gaussian") 2 else 1
    if (nm %in% names(spec$surface$buffers)) {
      x <- .msr_scale_surface_raw(spec$surface$buffers[[nm]], log_s,
        spec$kernel, inherits(par$beta, "advector"))
    } else {
      w <- exp(-spec$delta[[nm]] * exp(-power * log_s) /
                 if (power == 2) 2 else 1)
      x <- sum_rows(w * spec$V[[nm]]) / sum_rows(w * spec$C[[nm]])
    }
    kernel[[nm]] <- (x - spec$center[[nm]]) / spec$spread[[nm]]
    eta <- eta + kernel[[nm]] * par$beta[ncol(spec$X_fixed) + j]
  }
  if (include_response) for (term in spec$response$random)
    if (is.null(term$gaussian_block))
      eta <- eta + exp(par$log_response_sd[term$sd_index]) *
        par$response_z[term$z_index][term$index] else
      eta <- eta + .msr_gaussian_effect(term$gaussian_block,
        par$response_z[term$z_index], par$log_response_sd[term$sd_index])
  list(eta = eta, scales = scales, kernel = kernel)
}

.msr_bayes_objective <- function(spec, ad = FALSE) {
  if (ad) .msr_scale_surface_ram(spec)
  density <- if (ad) getExportedValue("RTMB", "dnorm") else stats::dnorm
  function(par) {
    p <- spec$priors
    eta <- .msr_bayes_components(spec, par)$eta
    ll <- .msr_scale_response_loglik(spec$family, spec$y, eta,
      par$log_dispersion, spec$trials, ad)
    nlp <- -sum(ll) - sum(density(par$beta, 0, p$beta_sd, log = TRUE)) -
      sum(density(par$scale_z, 0, 1, log = TRUE)) -
      sum(density(par$response_z, 0, 1, log = TRUE))
    for (nm in spec$covariates) {
      sm <- spec$scales[[nm]]
      idx <- sm$beta_index
      # Preserve evaluation of saved pre-cell-means specifications.
      means <- if (is.null(sm$prior_mean))
        c(log(p$scale_median[[nm]] / spec$unit_conv), rep(0, length(idx) - 1L)) else
        sm$prior_mean
      sds <- if (is.null(sm$prior_sd))
        c(p$scale_intercept_sd, rep(p$scale_coef_sd, length(idx) - 1L)) else sm$prior_sd
      nlp <- nlp - sum(density(par$scale_beta[idx], means, sds, log = TRUE))
    }
    # Half-normal priors on SDs, with the Jacobian for log-SD coordinates.
    nlp <- nlp - sum(density(exp(par$log_scale_sd), 0,
      p$scale_group_sd, log = TRUE) + log(2) + par$log_scale_sd)
    nlp <- nlp - sum(density(exp(par$log_response_sd), 0,
      p$response_group_sd, log = TRUE) + log(2) + par$log_response_sd)
    nlp - sum(density(par$log_dispersion, p$dispersion_log_mean,
                       p$dispersion_log_sd, log = TRUE))
  }
}

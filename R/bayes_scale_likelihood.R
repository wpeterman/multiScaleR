# This calculation works with ordinary doubles and with RTMB AD vectors.
# Scales remain explicit parameters: full Bayesian sampling integrates over
# their uncertainty instead of replacing it with a conditional optimum.
.msr_bayes_components <- function(spec, par, include_response = TRUE) {
  ad <- inherits(par$beta, "advector")
  scales <- kernel <- vector("list", length(spec$covariates))
  names(scales) <- names(kernel) <- spec$covariates
  expanded <- !is.null(spec$response_design)
  eta <- if (expanded) spec$offset else
    drop(spec$X_fixed %*% par$beta[seq_len(ncol(spec$X_fixed))]) + spec$offset
  for (j in seq_along(spec$covariates)) {
    nm <- spec$covariates[j]
    log_s <- .msr_bayes_log_scale(spec$scales[[nm]], par)
    scales[[nm]] <- exp(log_s) * spec$unit_conv
    kernel[[nm]] <- .msr_bayes_kernel_value(spec, nm, log_s, ad)
    if (!expanded) eta <- eta + kernel[[nm]] * par$beta[ncol(spec$X_fixed) + j]
  }
  if (expanded) eta <- eta + .msr_response_design_eta(spec$response_design, kernel, par$beta)
  if (include_response) for (term in spec$response$random)
    if (is.null(term$gaussian_block))
      eta <- eta + exp(par$log_response_sd[term$sd_index]) *
        par$response_z[term$z_index][term$index] else
      eta <- eta + .msr_gaussian_effect(term$gaussian_block,
        par$response_z[term$z_index], par$log_response_sd[term$sd_index])
  list(eta = eta, scales = scales, kernel = kernel)
}

.msr_bayes_log_scale <- function(sm, par) {
  log_s <- drop(sm$X %*% par$scale_beta[sm$beta_index])
  for (term in sm$random)
    if (is.null(term$gaussian_block))
      log_s <- log_s + exp(par$log_scale_sd[term$sd_index]) *
        par$scale_z[term$z_index][term$index] else
      log_s <- log_s + .msr_gaussian_effect(term$gaussian_block,
        par$scale_z[term$z_index], par$log_scale_sd[term$sd_index])
  log_s
}

.msr_bayes_kernel_value <- function(spec, nm, log_s, ad = FALSE) {
  if (nm %in% names(spec$surface$buffers)) {
    x <- .msr_scale_surface_raw(spec$surface$buffers[[nm]], log_s, spec$kernel, ad)
  } else {
    sum_rows <- if (ad) RTMB::rowSums else base::rowSums
    power <- if (spec$kernel == "gaussian") 2 else 1
    multiplier <- exp(-power * log_s)
    delta <- spec$delta[[nm]]
    w <- exp(-delta * multiplier / if (power == 2) 2 else 1)
    # In ordinary arithmetic, an underflowed distance scale has the exact
    # nearest-bin limit. Keep zeros from being multiplied by infinity.
    if (!ad && any(is.infinite(multiplier))) {
      rows <- which(is.infinite(multiplier))
      w[rows, ] <- delta[rows, , drop = FALSE] == 0
    }
    x <- sum_rows(w * spec$V[[nm]]) / sum_rows(w * spec$C[[nm]])
  }
  (x - spec$center[[nm]]) / spec$spread[[nm]]
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

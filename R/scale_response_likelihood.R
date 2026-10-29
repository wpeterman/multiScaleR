# Conditional response densities shared by posterior and likelihood targets.
# No population prior, latent density, integration or kernel arithmetic belongs
# here. Keep the normalized per-observation log density for ordinary and AD use.
.msr_scale_response_loglik <- function(family, y, eta, log_dispersion,
                                       trials = NULL, ad = FALSE) {
  if (ad) return(switch(family,
    gaussian = RTMB::dnorm(y, eta, exp(log_dispersion), log = TRUE),
    poisson = RTMB::dpois(y, exp(eta), log = TRUE),
    binomial = RTMB::dbinom_robust(y, trials, eta, log = TRUE),
    nbinom2 = RTMB::dnbinom_robust(y, eta, 2 * eta - log_dispersion, log = TRUE),
    stop("Unsupported scale response family.", call. = FALSE)))
  switch(family,
    gaussian = stats::dnorm(y, eta, exp(log_dispersion), log = TRUE),
    poisson = stats::dpois(y, exp(eta), log = TRUE),
    # The logit expression stays finite even when plogis(eta) rounds to 0 or 1.
    binomial = lchoose(trials, y) - y * pmax(-eta, 0) -
      (trials - y) * pmax(eta, 0) - trials * log1p(exp(-abs(eta))),
    # NB2 size is exp(log_dispersion). Avoid reconstructing it by subtracting
    # two predictors, which can lose the size at very large predictor values.
    nbinom2 = stats::dnbinom(y, mu = exp(eta), size = exp(log_dispersion), log = TRUE),
    stop("Unsupported scale response family.", call. = FALSE))
}

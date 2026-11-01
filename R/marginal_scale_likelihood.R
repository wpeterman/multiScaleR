# Separate unpenalized likelihood: only latent standard-normal densities remain.
# Kernel evaluation is shared read-only with the Bayesian implementation.
.msr_marginal_response_nll <- function(spec, parameters, ad = FALSE) {
  eta <- .msr_bayes_components(spec, parameters)$eta
  -.msr_scale_response_loglik(spec$family, spec$y, eta,
    parameters$log_dispersion, spec$trials, ad)
}

.msr_marginal_objective <- function(spec, ad = FALSE) {
  if (ad) .msr_scale_surface_ram(spec)
  normal <- if (ad) RTMB::dnorm else stats::dnorm
  function(parameters) {
    sum(.msr_marginal_response_nll(spec, parameters, ad)) -
      sum(normal(parameters$scale_z, 0, 1, log = TRUE)) -
      sum(normal(parameters$response_z, 0, 1, log = TRUE))
  }
}

.msr_marginal_adfun <- function(spec, parameters = spec$parameters) {
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Marginal likelihood requires the optional package 'RTMB'.", call. = FALSE)
  .msr_scale_surface_ram(spec)
  random <- intersect(c("scale_z", "response_z"), names(parameters)[lengths(parameters) > 0L])
  RTMB::MakeADFun(.msr_marginal_objective(spec, ad = TRUE), parameters,
    random = if (length(random)) random else NULL,
    inner.control = list(maxit = 1000L, tol = 1e-10, grad.tol = 1e-10, step.tol = 1e-12),
    silent = TRUE)
}

# Conditional fixed-parameter reference only. This preserves the native finite-K
# SUM, not a prior renormalized after truncation and not a Gaussian integral.
.msr_gdr_finite_score_bins <- function(expected_abundance, pd, pr, y_distance,
    y_removal, K, dispersion = Inf, probability_method = "caller_supplied_unqualified") {
  if (any(vapply(list(expected_abundance, pd, pr, y_distance, y_removal, K, dispersion),
      is.object, logical(1))))
    stop("This is a numeric-only conditional reference; AD and fitting admission remain guarded.", call. = FALSE)
  if (!is.numeric(K) || !is.null(dim(K)) || length(K) != 1L || !is.finite(K) || K < 0 || K != floor(K) || K > .Machine$integer.max)
    stop("Supply one finite nonnegative integer-valued cutoff K; its storage type is retained.", call. = FALSE)
  exact <- .msr_gdr_exact_score_bins(expected_abundance, pd, pr, y_distance,
    y_removal, dispersion, probability_method)
  probability <- rowSums(pd) * rowSums(pr)
  if (any(!is.finite(probability)) || any(probability < 0 | probability > 1))
    stop("Finite-K correction requires recorded probabilities within [0, 1].", call. = FALSE)
  observed <- rowSums(y_distance)
  if (any(observed > K)) stop("The finite cutoff must cover every observed total.", call. = FALSE)
  if (is.infinite(dispersion)) {
    missed_mean <- expected_abundance * (1 - probability)
    missed_size <- rep(Inf, length(observed))
    # Retain the already-qualified Poisson primitive and its original guards.
    correction <- .msr_gdr_poisson_cutoff_correction(expected_abundance, probability, observed, K)
  } else {
    # Gamma-Poisson thinning: M | Y is NB(size + Y), with its updated rate.
    # Work on the log scale so a large prior mean/size does not overflow a ratio.
    add <- function(a, b) pmax(a, b) + log1p(exp(pmin(a, b) - pmax(a, b)))
    missed_size <- dispersion + observed
    missed_mean <- exp(log(expected_abundance) + log1p(-probability) +
      add(log(dispersion), log(observed)) -
      add(log(dispersion), log(exact$per_event$expected_observed)))
    if (any(!is.finite(missed_size)) || any(!is.finite(missed_mean)))
      stop("Conditional missed-count moments are not representable; no target approximation is substituted.", call. = FALSE)
    correction <- stats::pnbinom(K - observed, size = missed_size, mu = missed_mean, log.p = TRUE)
  }
  if (any(!is.finite(correction)) || any(correction > 0))
    stop("The finite-K conditional CDF correction is invalid.", call. = FALSE)
  finite <- exact$per_event$logLik + correction
  list(logLik = sum(finite), infinite_logLik = exact$logLik,
    per_event = data.frame(observed_total = observed, expected_circle_abundance = expected_abundance,
      recorded_probability = probability, infinite_logLik = exact$per_event$logLik,
      log_CDF_correction = correction, finite_logLik = finite,
      infinite_parent_missed_mean_given_observed = missed_mean,
      infinite_parent_missed_size_given_observed = missed_size),
    cutoff = K, cutoff_storage = typeof(K), probability_method = probability_method,
    target = "finite_K_sum_via_conditional_missed_count_CDF",
    native_target_equivalence_requires_identical_bins = TRUE,
    prior_truncation_renormalized = FALSE, family = if (is.infinite(dispersion)) "P" else "NB",
    qualification = "private_conditional_fixed_parameter_reference_not_Gaussian_integrated",
    fit_admission_established = FALSE, uncertainty = "none_no_fitting_or_Gaussian_integration")
}

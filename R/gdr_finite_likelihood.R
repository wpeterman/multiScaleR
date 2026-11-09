# Private single-primary Poisson prerequisite. Preserve the native finite sum,
# which is computational truncation, not a renormalized abundance distribution.
.msr_gdr_finite_poisson_score_bins <- function(expected_abundance, pd, pr,
    y_distance, y_removal, cutoff, probability_method = "caller_supplied_unqualified") {
  if (!is.numeric(cutoff) || !is.null(dim(cutoff)) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff < 0 || cutoff != floor(cutoff) ||
      cutoff > .Machine$integer.max)
    stop("Supply one finite nonnegative integer abundance cutoff.")
  out <- .msr_gdr_exact_score_bins(expected_abundance, pd, pr,
    y_distance, y_removal, dispersion = Inf, probability_method = probability_method)
  correction <- .msr_gdr_poisson_cutoff_correction(expected_abundance,
    rowSums(pd) * rowSums(pr), rowSums(y_distance), cutoff)
  out$per_event$total <- out$per_event$total + correction
  out$per_event$logLik <- out$per_event$logLik + correction
  out$per_event$log_abundance_cdf <- correction
  out$logLik <- sum(out$per_event$logLik)
  out$abundance_cutoff <- as.integer(cutoff)
  out$finite_count_target <- "native_finite_K_not_renormalized"
  out$qualification <- "private_single_primary_fixed_point_P_finite_K_bins_only"
  out$uncertainty <- "none_no_fitting_or_backend_admission"
  out
}

# Reuse a validated invariant point without duplicating all event predictions
# across a cutoff ladder. This is the same per-event correction, in row order.
.msr_gdr_poisson_cutoff_correction <- function(expected_abundance, detection, observed, cutoff) {
  if (!is.numeric(cutoff) || !is.null(dim(cutoff)) || length(cutoff) != 1L ||
      !is.finite(cutoff) || cutoff < 0 || cutoff != floor(cutoff) || cutoff > .Machine$integer.max)
    stop("Supply one finite nonnegative integer abundance cutoff.")
  vectors <- list(expected_abundance, detection, observed)
  if (!all(vapply(vectors, function(x) is.numeric(x) && is.null(dim(x)) &&
      length(x) > 0L && all(is.finite(x)), logical(1))) ||
      length(unique(lengths(vectors))) != 1L || any(expected_abundance < 0) ||
      any(detection < 0) || any(observed < 0 | observed != floor(observed)))
    stop("Supply aligned finite abundance, detection and observed-count vectors.")
  if (any(observed > cutoff)) stop("The finite cutoff must cover every observed total.")
  if (any(detection > 1)) stop("Finite-K detection probabilities must not exceed one.")
  correction <- stats::ppois(cutoff - observed, expected_abundance * (1 - detection), log.p = TRUE)
  if (any(!is.finite(correction))) stop("Finite-K log probabilities are outside the representable private domain.")
  correction
}

# Ephemeral serial R callbacks, never a saved objective or backend admission.
# For x=log(mu), d log P(Pois(mu)<=m)/dx = -mu*Pois(m;mu)/P(Pois(mu)<=m).
# Using the stable log-CDF in the adjoint retains differentiability through all
# orders without taking log() of an underflowed nonlog probability.
.msr_gdr_poisson_logcdf_program <- function(upper) {
  if (!is.numeric(upper) || !is.null(dim(upper)) || !length(upper) ||
      any(!is.finite(upper)) || any(upper < 0 | upper != floor(upper)) ||
      any(upper > .Machine$integer.max))
    stop("Use a nonempty vector of finite nonnegative integer Poisson cutoffs.")
  if (!requireNamespace("RTMB", quietly = TRUE)) stop("This private AD primitive requires RTMB.")
  upper <- as.integer(upper)
  # Avoid capturing the factory caller or a later objective in either callback.
  value_env <- list2env(list(m = upper), parent = baseenv())
  value <- eval(quote(function(x) {
    if (length(x) != length(m) || any(!is.finite(x))) stop("Invalid finite log-rate vector.")
    result <- stats::ppois(m, exp(x), log.p = TRUE)
    if (any(!is.finite(result))) stop("Log-CDF outside the representable private domain.")
    result
  }), value_env)
  derivative_env <- list2env(list(m = upper), parent = baseenv())
  derivative <- eval(quote(function(x, y, dy)
    -exp((m + 1) * x - exp(x) - lgamma(m + 1) - y) * dy), derivative_env)
  list(evaluate = RTMB::ADjoint(value, derivative), upper = upper,
    requires_serial_R_callbacks = TRUE, portable_objective = FALSE,
    fit_admission_established = FALSE,
    qualification = "private_log_CDF_derivative_prerequisite_not_Gaussian_integration")
}

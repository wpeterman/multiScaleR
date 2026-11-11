# Evaluate the native closed-population finite sum at fixed conditional
# parameters. Each site has one abundance shared by its primary surveys.
.msr_gdr_finite_primary_score <- function(expected_abundance, availability, pd, pr,
    y_distance, y_removal, K, dispersion = Inf, probability_method = "caller_supplied_unqualified") {
  inputs <- list(expected_abundance, availability, pd, pr, y_distance, y_removal, K, dispersion)
  if (any(vapply(inputs, is.object, logical(1))))
    stop("This finite-primary reference requires ordinary numeric inputs.", call. = FALSE)
  if (!is.numeric(K) || !is.null(dim(K)) || length(K) != 1L || !is.finite(K) ||
      K < 0 || K != floor(K) || K > .Machine$integer.max)
    stop("Supply one finite nonnegative integer-valued cutoff K.", call. = FALSE)
  if (!is.matrix(availability) || !is.numeric(availability) || !length(availability) ||
      any(!is.finite(availability)) || any(availability < 0 | availability > 1) ||
      !is.numeric(expected_abundance) || !is.null(dim(expected_abundance)) ||
      length(expected_abundance) != nrow(availability) ||
      any(!is.finite(expected_abundance)) || any(expected_abundance < 0))
    stop("Supply site abundance and site-by-primary availability within [0, 1].", call. = FALSE)
  if (!is.matrix(y_distance) || nrow(y_distance) != length(availability))
    stop("Supply one paired margin per site and primary in native site-major order.", call. = FALSE)
  if (length(availability) * (as.double(K) + 1) > 1e7)
    stop("The requested finite-primary reference exceeds the diagnostic work bound.", call. = FALSE)

  primary <- ncol(availability)
  margin <- .msr_gdr_exact_score_bins(rep(expected_abundance, each = primary), pd, pr,
    y_distance, y_removal, dispersion, probability_method)$per_event
  observed <- margin$observed_total
  if (any(observed > K)) stop("The finite cutoff must cover every observed total.", call. = FALSE)
  phi <- as.vector(t(availability))
  recorded <- phi * rowSums(pd) * rowSums(pr)
  if (any(recorded < 0 | recorded > 1))
    stop("Combined recorded probabilities must lie within [0, 1].", call. = FALSE)
  poisson <- is.infinite(dispersion)
  site_logLik <- vapply(seq_along(expected_abundance), function(site) {
    events <- (site - 1L) * primary + seq_len(primary)
    candidates <- seq.int(max(observed[events]), K)
    mu <- expected_abundance[site]
    weight <- if (poisson) stats::dpois(candidates, mu, log = TRUE) else
      stats::dnbinom(candidates, size = dispersion, mu = mu, log = TRUE)
    for (event in events)
      weight <- weight + stats::dbinom(observed[event], candidates, recorded[event], log = TRUE)
    peak <- max(weight)
    if (peak == -Inf) return(-Inf)
    peak + log(sum(exp(weight - peak))) + sum(margin$distance[events] + margin$removal[events])
  }, numeric(1))
  list(logLik = sum(site_logLik), site_logLik = site_logLik, cutoff = K,
    cutoff_storage = typeof(K), family = if (poisson) "P" else "NB",
    primary_periods = primary, probability_method = probability_method,
    recorded_probability = matrix(recorded, nrow(availability), primary, byrow = TRUE),
    target = "conditional_shared_abundance_finite_K", prior_truncation_renormalized = FALSE,
    fit_admission_established = FALSE, uncertainty = "none_fixed_conditional_reference")
}

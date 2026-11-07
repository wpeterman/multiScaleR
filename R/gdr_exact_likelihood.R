# Private numerical prerequisite. No fitted-model adapter or backend admission.
# Probabilities are integrated over point-transect annuli, normalized by the
# complete sampled circle. Work in radius units to avoid squared-scale overflow.
.msr_gdr_exact_distance_bins <- function(distance_scale, dist_breaks,
    keyfun = c("halfnorm", "exp", "uniform")) {
  keyfun <- match.arg(keyfun)
  if (!is.numeric(dist_breaks) || !is.null(dim(dist_breaks)) || length(dist_breaks) < 2L ||
      any(!is.finite(dist_breaks)) || dist_breaks[1L] != 0 || any(diff(dist_breaks) <= 0))
    stop("Use finite increasing point-transect distance breaks beginning at zero.")
  if (!is.numeric(distance_scale) || !is.null(dim(distance_scale)) || !length(distance_scale) || anyNA(distance_scale) ||
      any(distance_scale <= 0)) stop("Distance scales must be positive; Inf requests the uniform limit.")
  radius <- utils::tail(dist_breaks, 1L)
  lower <- utils::head(dist_breaks, -1L) / radius
  upper <- utils::tail(dist_breaks, -1L) / radius
  out <- matrix(0, length(distance_scale), length(lower))
  beta <- radius / distance_scale
  if (any(!is.finite(beta))) stop("Distance ratios are outside the representable private domain.")
  for (j in seq_along(lower)) {
    width <- (upper[j] - lower[j]) * (upper[j] + lower[j])
    if (keyfun == "uniform") {
      out[, j] <- width
    } else if (keyfun == "halfnorm") {
      a <- 0.5 * beta^2
      if (any(!is.finite(a))) stop("Squared distance ratios are outside the private domain.")
      x <- a * width
      relative <- rep(1, length(x))
      positive <- x > 0
      relative[positive] <- -expm1(-x[positive]) / x[positive]
      out[, j] <- exp(-a * lower[j]^2) * width * relative
    } else {
      d <- upper[j] - lower[j]
      x <- beta * d
      e0 <- rep(1, length(x))
      positive <- x > 0
      e0[positive] <- -expm1(-x[positive]) / x[positive]
      # E1(x) = integral_0^1 t exp(-x*t) dt. Its closed expression subtracts
      # nearly equal numbers for thin bins or very large distance scales.
      small <- x <= 0.1
      e1 <- numeric(length(x))
      if (any(small)) {
        term <- rep(1, sum(small))
        for (k in 0:12) {
          if (k > 0L) term <- term * (-x[small]) / k
          e1[small] <- e1[small] + term / (k + 2)
        }
      }
      if (any(!small)) {
        z <- x[!small]
        e1[!small] <- (-expm1(-z) - z * exp(-z)) / z^2
      }
      out[, j] <- 2 * exp(-beta * lower[j]) * (lower[j] * d * e0 + d^2 * e1)
    }
  }
  if (any(!is.finite(out)) || any(out < 0) || any(rowSums(out) > 1 + 1e-12))
    stop("The private distance integration did not produce valid probabilities.")
  out
}

.msr_gdr_exact_removal_bins <- function(probability, interval_lengths) {
  if (!is.matrix(probability) || !is.numeric(probability) || any(!is.finite(probability)) ||
      any(probability < 0 | probability > 1)) stop("Supply a site-by-interval probability matrix.")
  if (!is.numeric(interval_lengths) || !is.null(dim(interval_lengths)) || length(interval_lengths) != ncol(probability) ||
      any(!is.finite(interval_lengths)) || any(interval_lengths <= 0) ||
      any(interval_lengths != floor(interval_lengths))) stop("Supply positive integer interval lengths.")
  out <- probability
  survival <- rep(1, nrow(out))
  for (j in seq_len(ncol(out))) {
    interval <- -expm1(interval_lengths[j] * log1p(-probability[, j]))
    out[, j] <- survival * interval
    survival <- survival * (1 - interval)
  }
  out
}

.msr_gdr_exact_observed <- function(expected_abundance, distance_scale, removal,
    y_distance, y_removal, dist_breaks, interval_lengths,
    keyfun = c("halfnorm", "exp", "uniform"), dispersion = Inf) {
  keyfun <- match.arg(keyfun)
  counts <- function(y) is.matrix(y) && is.numeric(y) && nrow(y) > 0L && ncol(y) > 0L &&
    all(is.finite(y)) && all(y >= 0) && all(y == floor(y))
  if (!counts(y_distance) || !counts(y_removal) || nrow(y_distance) != nrow(y_removal) ||
      any(rowSums(y_distance) != rowSums(y_removal)))
    stop("The private identity requires complete nonnegative integer paired margins with equal totals.")
  n <- nrow(y_distance)
  if (!is.numeric(expected_abundance) || !is.null(dim(expected_abundance)) || length(expected_abundance) != n ||
      any(!is.finite(expected_abundance)) || any(expected_abundance < 0))
    stop("Supply one finite nonnegative expected circle abundance per event.")
  if (!is.numeric(distance_scale) || !is.null(dim(distance_scale)) || !length(distance_scale) ||
      !length(distance_scale) %in% c(1L, n)) stop("Supply one distance scale or one per event.")
  distance_scale <- rep_len(distance_scale, n)
  if (!is.matrix(removal)) {
    if (!is.numeric(removal) || !length(removal) %in% c(1L, n))
      stop("Supply a removal probability per event or a site-by-interval matrix.")
    removal <- matrix(rep_len(removal, n), n, length(interval_lengths))
  }
  if (nrow(removal) != n) stop("Removal rows must match survey events.")
  pd <- .msr_gdr_exact_distance_bins(distance_scale, dist_breaks, keyfun)
  pr <- .msr_gdr_exact_removal_bins(removal, interval_lengths)
  if (ncol(pd) != ncol(y_distance) || ncol(pr) != ncol(y_removal))
    stop("Distance and removal bins must match the observed margins.")
  .msr_gdr_exact_score_bins(expected_abundance, pd, pr, y_distance, y_removal,
    dispersion, probability_method = "stable_analytic_point_annuli")
}

# Separate the thinning identity from its distance integration. Native finite-K
# and native-bin approximation errors can then be measured independently.
.msr_gdr_exact_score_bins <- function(expected_abundance, pd, pr, y_distance,
    y_removal, dispersion = Inf, probability_method = "caller_supplied_unqualified") {
  counts <- function(y) is.matrix(y) && is.numeric(y) && nrow(y) > 0L && ncol(y) > 0L &&
    all(is.finite(y)) && all(y >= 0) && all(y == floor(y))
  if (!counts(y_distance) || !counts(y_removal) || nrow(y_distance) != nrow(y_removal) ||
      any(rowSums(y_distance) != rowSums(y_removal)))
    stop("The private identity requires complete nonnegative integer paired margins with equal totals.")
  n <- nrow(y_distance)
  if (!is.numeric(expected_abundance) || !is.null(dim(expected_abundance)) || length(expected_abundance) != n ||
      any(!is.finite(expected_abundance)) || any(expected_abundance < 0))
    stop("Supply one finite nonnegative expected circle abundance per event.")
  valid_bins <- function(p, y) is.matrix(p) && is.numeric(p) && identical(dim(p), dim(y)) &&
    all(is.finite(p)) && all(p >= 0) && all(rowSums(p) <= 1 + 1e-12)
  if (!valid_bins(pd, y_distance) || !valid_bins(pr, y_removal))
    stop("Supply finite nonnegative bin probabilities aligned with their response margins.")
  if (!is.character(probability_method) || length(probability_method) != 1L ||
      is.na(probability_method) || !nzchar(probability_method))
    stop("Identify the supplied probability method explicitly.")
  if (!is.numeric(dispersion) || !is.null(dim(dispersion)) || length(dispersion) != 1L || is.na(dispersion) || dispersion <= 0)
    stop("Dispersion is one positive NB size or Inf for the exact Poisson boundary.")
  pd_total <- rowSums(pd)
  pr_total <- rowSums(pr)
  distance_conditional <- pd / pd_total
  removal_conditional <- pr / pr_total
  observed_mean <- expected_abundance * pd_total * pr_total
  observed <- rowSums(y_distance)
  total <- if (is.infinite(dispersion)) stats::dpois(observed, observed_mean, log = TRUE) else
    stats::dnbinom(observed, size = dispersion, mu = observed_mean, log = TRUE)
  multinomial <- function(y, probability) {
    # A zero total contributes exactly zero, even when detection is impossible.
    terms <- matrix(0, nrow(y), ncol(y))
    positive <- y > 0
    p <- probability[positive]
    p[is.nan(p)] <- 0
    terms[positive] <- y[positive] * log(p)
    lgamma(rowSums(y) + 1) - rowSums(lgamma(y + 1)) + rowSums(terms)
  }
  distance <- multinomial(y_distance, distance_conditional)
  removal_part <- multinomial(y_removal, removal_conditional)
  expected_distance <- distance_conditional * observed_mean
  expected_removal <- removal_conditional * observed_mean
  expected_distance[observed_mean == 0, ] <- 0
  expected_removal[observed_mean == 0, ] <- 0
  conditional_distance <- distance_conditional * observed
  conditional_removal <- removal_conditional * observed
  conditional_distance[observed == 0, ] <- 0
  conditional_removal[observed == 0, ] <- 0
  list(logLik = sum(total + distance + removal_part),
    per_event = data.frame(total = total, distance = distance, removal = removal_part,
      logLik = total + distance + removal_part, expected_observed = observed_mean,
      expected_abundance = expected_abundance, observed_total = observed),
    expected_distance = expected_distance, expected_removal = expected_removal,
    conditional_distance = conditional_distance, conditional_removal = conditional_removal,
    distance_probability = pd, removal_probability = pr,
    distance_conditional = distance_conditional, removal_conditional = removal_conditional,
    zero_detection_events = which(pd_total == 0 | pr_total == 0),
    probability_method = probability_method,
    quantities = list(expected_bands = "unconditional_expected_observed_band_counts",
      conditional_bands = "band_counts_given_the_single_observed_total",
      expected_abundance = "expected_circle_abundance_not_density_or_realized_N"),
    qualification = "private_single_primary_fixed_point_P_NB_complete_point_transect_only",
    uncertainty = "none_no_fitting_or_backend_admission")
}

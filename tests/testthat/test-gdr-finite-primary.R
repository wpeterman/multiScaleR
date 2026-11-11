test_that("shared abundance across primary periods agrees with the native finite sum", {
  skip_if_not_installed("unmarked")
  skip_if(as.character(utils::packageVersion("unmarked")) != "1.5.2", "Native GDR reference needs qualification")
  # Distinct totals and covariates identify the site-major primary layout.
  yd <- rbind(c(1, 0, 0, 2), c(0, 0, 1, 0), c(1, 1, 0, 0))
  yr <- rbind(c(0, 1, 0, 1, 0, 1), c(0, 0, 0, 0, 1, 0), c(1, 0, 1, 0, 0, 0))
  breaks <- c(0, 25, 75)
  durations <- c(1, 2, 1)
  temperature <- c(-.2, .5, -.5, .2, .7, -.4)
  habitat <- c(-.5, .5, .2)
  frame <- unmarked::unmarkedFrameGDR(yd, yr, numPrimary = 2L, dist.breaks = breaks,
    unitsIn = "m", period.lengths = durations, siteCovs = data.frame(habitat = habitat),
    yearlySiteCovs = data.frame(temperature = temperature))
  abundance <- exp(log(8) + .2 * habitat)
  phi <- matrix(plogis(qlogis(.6) + .7 * temperature), 3L, 2L, byrow = TRUE)
  # Reproduce the native bins so this check isolates the count summation.
  annuli <- pi * diff(breaks^2)
  area <- matrix(rep(annuli, each = length(phi)), nrow = length(phi))
  pd <- .msr_gdr_audit_native_bins(rep(40, length(phi)), breaks, "halfnorm", area, area / sum(annuli))
  pr <- .msr_gdr_exact_removal_bins(matrix(.3, length(phi), length(durations)), durations)
  distance <- matrix(as.vector(t(yd)), ncol = 2L, byrow = TRUE)
  removal <- matrix(as.vector(t(yr)), ncol = 3L, byrow = TRUE)
  for (K in c(3L, 40L)) for (size in c(Inf, 2)) {
    mixture <- if (is.infinite(size)) "P" else "NB"
    start <- c(log(8), .2, if (mixture == "NB") log(size), qlogis(.6), .7, log(40), qlogis(.3))
    fit <- unmarked::gdistremoval(lambdaformula = ~habitat, phiformula = ~temperature,
      data = frame, keyfun = "halfnorm", output = "abund", K = K, starts = start,
      mixture = mixture, engine = "C", se = FALSE, control = list(maxit = 0L), threads = 1L)
    value <- .msr_gdr_finite_primary_score(abundance, phi, pd, pr, distance, removal, K, size, "native_bins")
    expect_equal(value$logLik, -fit@nllFun(start), tolerance = 1e-11)
    expect_identical(value$cutoff, K)
    expect_false(value$prior_truncation_renormalized)
    expect_false(value$fit_admission_established)
  }
})

test_that("exact zero availability recovers the retained-primary parent", {
  yd <- rbind(c(1, 1), c(0, 0))
  yr <- rbind(c(0, 2), c(0, 0))
  pd <- rbind(c(.2, .3), c(.3, .2))
  pr <- rbind(c(.3, .4), c(.4, .3))
  for (size in c(Inf, 2)) {
    parent <- .msr_gdr_finite_score_bins(8, .65 * pd[1L, , drop = FALSE], pr[1L, , drop = FALSE],
      yd[1L, , drop = FALSE], yr[1L, , drop = FALSE], 40L, size)
    exact <- .msr_gdr_finite_primary_score(8, matrix(c(.65, 0), 1L), pd, pr, yd, yr, 40L, size)
    nearby <- .msr_gdr_finite_primary_score(8, matrix(c(.65, 1e-10), 1L), pd, pr, yd, yr, 40L, size)
    expect_equal(exact$logLik, parent$logLik, tolerance = 1e-12)
    expect_lt(abs(nearby$logLik - exact$logLik), 1e-8)
    expect_lt(nearby$logLik, exact$logLik)
    impossible <- .msr_gdr_finite_primary_score(8, matrix(c(0, .5), 1L), pd, pr, yd, yr, 40L, size)
    expect_identical(impossible$logLik, -Inf)
    both_zero <- .msr_gdr_finite_primary_score(8, matrix(0, 1L, 2L), pd, pr,
      0 * yd, 0 * yr, 4L, size)
    prior_mass <- if (is.infinite(size)) ppois(4, 8, log.p = TRUE) else pnbinom(4, size, mu = 8, log.p = TRUE)
    expect_equal(both_zero$logLik, prior_mass, tolerance = 1e-12)
  }
  phi <- matrix(c(.65, .2), 1L)
  value <- .msr_gdr_finite_primary_score(8, phi, pd, pr, yd, yr, 40L)
  near_poisson <- .msr_gdr_finite_primary_score(8, phi, pd, pr, yd, yr, 40L, 1e8)
  expect_equal(near_poisson$logLik, value$logLik, tolerance = 1e-6)
  expect_error(.msr_gdr_finite_primary_score(8, phi, pd, pr, yd, yr, 1L), "cover every observed total")
  expect_error(.msr_gdr_finite_primary_score(8, matrix(c(.65, -.1), 1L), pd, pr, yd, yr, 40L), "within")
})

test_that("rare availability separates a small cutoff from the shared-count limit", {
  counts <- matrix(c(2, 1), 2L)
  bins <- matrix(1, 2L, 1L)
  mean <- 10000
  rates <- c(1.2, .8)
  phi <- matrix(rates / mean, 1L)
  for (size in c(Inf, 2)) {
    # In the rare-detection limit the Poisson parent gives independent counts.
    # The NB parent keeps one shared gamma intensity across both periods.
    limit <- if (is.infinite(size)) sum(dpois(counts, rates, log = TRUE)) else {
      denominator <- size + sum(rates)
      lgamma(size + sum(counts)) - lgamma(size) - sum(lgamma(counts + 1)) +
        size * log(size / denominator) + sum(counts * log(rates / denominator))
    }
    reference <- .msr_gdr_finite_primary_score(mean, phi, bins, bins, counts, counts, 100000L, size)
    capped <- .msr_gdr_finite_primary_score(mean, phi, bins, bins, counts, counts, 40L, size)
    expect_lt(abs(reference$logLik - limit), 2e-3)
    expect_gt(limit - capped$logLik, 10)
  }
})

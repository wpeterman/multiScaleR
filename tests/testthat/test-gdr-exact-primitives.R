testthat::test_that("stable distance bins agree with independent point quadrature and uniform limits", {
  breaks <- c(0, 20, 50, 75)
  scales <- 75 * c(0.1, 0.5, 1, 10, 1000, 1e6, Inf)
  for (key in c("halfnorm", "exp")) {
    actual <- .msr_gdr_exact_distance_bins(scales, breaks, key)
    reference <- t(vapply(scales, function(sigma) vapply(seq_len(3), function(j) {
      f <- if (key == "halfnorm") function(r) 2 * r / 75^2 * exp(-0.5 * (r / sigma)^2) else
        function(r) 2 * r / 75^2 * exp(-r / sigma)
      stats::integrate(f, breaks[j], breaks[j + 1L], rel.tol = 1e-12,
        abs.tol = 1e-13, subdivisions = 200L)$value
    }, numeric(1)), numeric(3)))
    testthat::expect_lte(max(abs(actual - reference)), 1e-12)
    testthat::expect_equal(actual[nrow(actual), ], diff((breaks / 75)^2), tolerance = 1e-14)
    testthat::expect_equal(actual, .msr_gdr_exact_distance_bins(scales / 1000, breaks / 1000, key),
      tolerance = 1e-14)
    testthat::expect_true(all(actual >= 0) && all(rowSums(actual) <= 1 + 1e-12))
  }
  uniform <- .msr_gdr_exact_distance_bins(1, breaks, "uniform")
  testthat::expect_equal(as.numeric(uniform), diff((breaks / 75)^2), tolerance = 1e-14)
})

testthat::test_that("exponential thin-bin and E1 cutoff cases retain probability mass", {
  for (delta in c(0.099999, 0.1, 0.100001)) {
    sigma <- 0.1 / delta
    actual <- .msr_gdr_exact_distance_bins(sigma, c(0, 0.1, 1), "exp")[1, 1]
    reference <- stats::integrate(function(r) 2 * r * exp(-r / sigma),
      0, 0.1, rel.tol = 1e-12, abs.tol = 1e-13, subdivisions = 200L)$value
    testthat::expect_lte(abs(actual - reference), 1e-14)
  }
  actual <- .msr_gdr_exact_distance_bins(1e-5, c(0, 1e-6, 1), "exp")[1, 1]
  reference <- stats::integrate(function(r) 2 * r * exp(-r / 1e-5),
    0, 1e-6, rel.tol = 1e-12, abs.tol = 1e-13, subdivisions = 200L)$value
  testthat::expect_lte(abs(actual - reference), 1e-24)
})

testthat::test_that("first-removal probabilities preserve unequal durations and endpoint limits", {
  p <- matrix(c(.2, .4, .6, .1, .3, .5), 2, 3, byrow = TRUE)
  lengths <- c(1, 2, 1)
  actual <- .msr_gdr_exact_removal_bins(p, lengths)
  interval <- 1 - (1 - p)^rep(lengths, each = 2)
  reference <- cbind(interval[, 1], (1 - interval[, 1]) * interval[, 2],
    (1 - interval[, 1]) * (1 - interval[, 2]) * interval[, 3])
  testthat::expect_equal(actual, reference, tolerance = 1e-14)
  testthat::expect_equal(rowSums(actual), 1 - apply(1 - interval, 1, prod), tolerance = 1e-14)
  testthat::expect_identical(.msr_gdr_exact_removal_bins(matrix(0, 1, 3), lengths), matrix(0, 1, 3))
  testthat::expect_equal(.msr_gdr_exact_removal_bins(matrix(1, 1, 3), lengths), matrix(c(1, 0, 0), 1))
})

testthat::test_that("one observed total and paired conditional margins reproduce P and NB probabilities", {
  yd <- rbind(c(0, 0, 0), c(1, 0, 1), c(2, 1, 0), c(0, 2, 1), c(1, 1, 1))
  yr <- rbind(c(0, 0, 0), c(1, 1, 0), c(1, 1, 1), c(2, 0, 1), c(0, 1, 2))
  abundance <- c(8, 10, 6, 9, 4)
  for (size in c(2, 20, Inf)) {
    actual <- .msr_gdr_exact_observed(abundance, 40, .3, yd, yr,
      c(0, 20, 50, 75), c(1, 2, 1), "halfnorm", size)
    reference <- vapply(seq_len(5), function(i) {
      total <- if (is.infinite(size)) stats::dpois(sum(yd[i, ]),
        actual$per_event$expected_observed[i], log = TRUE) else
          stats::dnbinom(sum(yd[i, ]), size = size,
            mu = actual$per_event$expected_observed[i], log = TRUE)
      total + stats::dmultinom(yd[i, ], prob = actual$distance_conditional[i, ], log = TRUE) +
        stats::dmultinom(yr[i, ], prob = actual$removal_conditional[i, ], log = TRUE)
    }, numeric(1))
    testthat::expect_equal(actual$per_event$logLik, reference, tolerance = 1e-12)
    testthat::expect_equal(actual$logLik, sum(reference), tolerance = 1e-12)
    testthat::expect_equal(rowSums(actual$expected_distance), actual$per_event$expected_observed)
    testthat::expect_equal(rowSums(actual$expected_removal), actual$per_event$expected_observed)
    testthat::expect_equal(rowSums(actual$conditional_distance), rowSums(yd))
    testthat::expect_equal(rowSums(actual$conditional_removal), rowSums(yr))
    testthat::expect_identical(actual$uncertainty, "none_no_fitting_or_backend_admission")
  }
  zero <- .msr_gdr_exact_observed(8, Inf, 0, matrix(0, 1, 3), matrix(0, 1, 3),
    c(0, 20, 50, 75), c(1, 2, 1), "exp")
  testthat::expect_equal(zero$logLik, 0)
  testthat::expect_equal(zero$zero_detection_events, 1L)
  testthat::expect_true(all(zero$expected_removal == 0))
  impossible <- .msr_gdr_exact_observed(0, 40, .3, matrix(c(1, 0, 0), 1),
    matrix(c(1, 0, 0), 1), c(0, 20, 50, 75), c(1, 2, 1), "halfnorm")
  testthat::expect_identical(impossible$logLik, -Inf)
})

testthat::test_that("unsupported private inputs cannot silently alter the response target", {
  testthat::expect_error(.msr_gdr_exact_distance_bins(0, c(0, 75)), "positive")
  testthat::expect_error(.msr_gdr_exact_distance_bins(40, c(0, Inf)), "finite increasing")
  testthat::expect_error(.msr_gdr_exact_distance_bins(40, c(5, 75)), "beginning at zero")
  testthat::expect_error(.msr_gdr_exact_removal_bins(matrix(.3, 1, 3), c(1, 0, 1)), "positive integer")
  testthat::expect_error(.msr_gdr_exact_observed(8, 40, .3, matrix(c(1, NA, 0), 1),
    matrix(c(1, 0, 0), 1), c(0, 20, 50, 75), c(1, 2, 1)), "complete")
  testthat::expect_error(.msr_gdr_exact_observed(8, 40, .3, matrix(c(1, 0, 0), 1),
    matrix(c(0, 0, 0), 1), c(0, 20, 50, 75), c(1, 2, 1)), "equal totals")
  testthat::expect_error(.msr_gdr_exact_observed(matrix(8, 1), 40, .3, matrix(c(1, 0, 0), 1),
    matrix(c(1, 0, 0), 1), c(0, 20, 50, 75), c(1, 2, 1)), "one finite")
  testthat::expect_error(.msr_gdr_exact_distance_bins(matrix(40, 1), c(0, 75)), "positive")
})

testthat::test_that("density and distance units yield the same explicit circle abundance target", {
  yd <- matrix(c(1, 1, 0), 1)
  yr <- matrix(c(1, 0, 1), 1)
  density_kmsq <- 800
  abundance_m <- density_kmsq * pi * (75 / 1000)^2
  density_ha <- density_kmsq / 100
  abundance_km <- density_ha * pi * 0.075^2 * 100
  meters <- .msr_gdr_exact_observed(abundance_m, 40, .3, yd, yr,
    c(0, 20, 50, 75), c(1, 2, 1), "halfnorm")
  kilometers <- .msr_gdr_exact_observed(abundance_km, .04, .3, yd, yr,
    c(0, .02, .05, .075), c(1, 2, 1), "halfnorm")
  testthat::expect_equal(abundance_m, abundance_km, tolerance = 1e-14)
  testthat::expect_equal(meters$logLik, kilometers$logLik, tolerance = 1e-12)
  testthat::expect_equal(meters$per_event$expected_observed,
    kilometers$per_event$expected_observed, tolerance = 1e-14)
})

testthat::test_that("supplied native-bin targets retain their own observed-process identity", {
  pd <- matrix(c(.1, .2, .1), 1)
  pr <- matrix(c(.3, .3, .2), 1)
  yd <- matrix(c(1, 1, 0), 1)
  yr <- matrix(c(1, 0, 1), 1)
  for (size in c(2, Inf)) {
    actual <- .msr_gdr_exact_score_bins(8, pd, pr, yd, yr, size,
      probability_method = "supplied_native_bin_fixture")
    mean <- 8 * sum(pd) * sum(pr)
    total <- if (is.infinite(size)) stats::dpois(2, mean, log = TRUE) else
      stats::dnbinom(2, size = size, mu = mean, log = TRUE)
    reference <- total + stats::dmultinom(as.numeric(yd), prob = as.numeric(pd), log = TRUE) +
      stats::dmultinom(as.numeric(yr), prob = as.numeric(pr), log = TRUE)
    testthat::expect_equal(actual$logLik, reference, tolerance = 1e-12)
    testthat::expect_identical(actual$probability_method, "supplied_native_bin_fixture")
  }
  testthat::expect_error(.msr_gdr_exact_score_bins(8, pd * 4, pr, yd, yr), "bin probabilities")
  testthat::expect_error(.msr_gdr_exact_score_bins(8, pd, pr, yd, yr,
    probability_method = NA_character_), "method explicitly")
})

testthat::test_that("distance and removal margins can have different numbers of bins", {
  actual <- .msr_gdr_exact_observed(8, 40, .3, matrix(c(1, 1), 1),
    matrix(c(1, 0, 1), 1), c(0, 25, 75), c(1, 2, 1), "halfnorm")
  testthat::expect_identical(dim(actual$expected_distance), c(1L, 2L))
  testthat::expect_identical(dim(actual$expected_removal), c(1L, 3L))
  testthat::expect_true(is.finite(actual$logLik))
})

test_that("finite-K P and NB identities match independent latent-count sums", {
  # Same five fixed paired margins as the existing exact-observation fixture.
  yd <- rbind(c(0, 0, 0), c(1, 0, 1), c(2, 1, 0), c(0, 2, 1), c(1, 1, 1))
  yr <- rbind(c(0, 0, 0), c(1, 1, 0), c(1, 1, 1), c(2, 0, 1), c(0, 1, 2))
  abundance <- c(8, 10, 6, 9, 4)
  pd <- .msr_gdr_exact_distance_bins(rep(40, 5), c(0, 20, 50, 75), "halfnorm")
  pr <- .msr_gdr_exact_removal_bins(matrix(.3, 5, 3), c(1, 2, 1))
  for (size in c(2, 20, Inf)) for (K in c(3, 10, 40)) {
    actual <- .msr_gdr_finite_score_bins(abundance, pd, pr, yd, yr, K, size)
    reference <- vapply(seq_len(nrow(yd)), function(i) {
      y <- sum(yd[i, ]); n <- y:K; p <- sum(pd[i, ]) * sum(pr[i, ])
      prior <- if (is.infinite(size)) dpois(n, abundance[i], log = TRUE) else
        dnbinom(n, size = size, mu = abundance[i], log = TRUE)
      terms <- prior + dbinom(y, n, p, log = TRUE)
      top <- max(terms)
      top + log(sum(exp(terms - top))) +
        dmultinom(yd[i, ], prob = pd[i, ], log = TRUE) +
        dmultinom(yr[i, ], prob = pr[i, ], log = TRUE)
    }, 0)
    expect_equal(unname(actual$per_event$finite_logLik), reference, tolerance = 1e-12)
    expect_equal(actual$logLik, sum(reference), tolerance = 1e-12)
    expect_true(all(actual$per_event$log_CDF_correction <= 0))
    expect_false(actual$prior_truncation_renormalized)
    expect_false(actual$fit_admission_established)
    expect_identical(actual$uncertainty, "none_no_fitting_or_Gaussian_integration")
    if (is.infinite(size)) expect_equal(actual$logLik,
      .msr_gdr_finite_poisson_score_bins(abundance, pd, pr, yd, yr, K)$logLik, tolerance = 0)
  }
  integer_cutoff <- .msr_gdr_finite_score_bins(abundance, pd, pr, yd, yr, 10L)
  numeric_cutoff <- .msr_gdr_finite_score_bins(abundance, pd, pr, yd, yr, 10)
  expect_identical(integer_cutoff$cutoff_storage, "integer")
  expect_identical(numeric_cutoff$cutoff_storage, "double")
  expect_identical(integer_cutoff$cutoff, 10L)
  expect_identical(numeric_cutoff$cutoff, 10)
  expect_equal(integer_cutoff$logLik, numeric_cutoff$logLik, tolerance = 0)
  cat("GDR_FINITE_THINNING P_NB=2 cutoffs=3 conditional_CDF=exact prior_renormalization=0 fitting_admission=0\n")
})

test_that("cutoff and exact Poisson parents retain impossible and zero-detection cases", {
  yd <- yr <- matrix(c(1, 0, 0), 1)
  pd <- matrix(c(.1, .2, .1), 1); pr <- matrix(c(.3, .3, .2), 1)
  expect_error(.msr_gdr_finite_score_bins(8, pd, pr, yd, yr, 0), "cover every observed")
  impossible <- .msr_gdr_finite_score_bins(0, pd, pr, yd, yr, 10)
  expect_identical(impossible$logLik, -Inf)
  expect_equal(impossible$per_event$log_CDF_correction, 0)
  for (size in c(2, Inf)) {
    zero <- .msr_gdr_finite_score_bins(8, pd, matrix(0, 1, 3),
      matrix(0, 1, 3), matrix(0, 1, 3), 10, size)
    reference <- if (is.infinite(size)) ppois(10, 8, log.p = TRUE) else pnbinom(10, size, mu = 8, log.p = TRUE)
    expect_equal(zero$logLik, reference, tolerance = 1e-12)
    expect_equal(zero$infinite_logLik, 0)
  }
  poisson <- .msr_gdr_finite_score_bins(8, pd, pr, yd, yr, 10, Inf)
  nb <- .msr_gdr_finite_score_bins(8, pd, pr, yd, yr, 10, 1e8)
  expect_equal(nb$logLik, poisson$logLik, tolerance = 1e-6)
  expect_identical(poisson$family, "P")
  expect_identical(nb$family, "NB")
  for (K in list(-1, 2.5, Inf, NA_real_, c(1, 2), .Machine$integer.max + 1))
    expect_error(.msr_gdr_finite_score_bins(8, pd, pr, yd, yr, K), "cutoff K")
  expect_error(.msr_gdr_finite_score_bins(structure(8, class = "unqualified_numeric"),
    pd, pr, yd, yr, 10), "numeric-only conditional")
})

test_that("rare-availability sections expose finite-cutoff artifacts without fitting", {
  yd <- rbind(c(0, 0, 0), c(1, 0, 1), c(2, 1, 0), c(0, 2, 1), c(1, 1, 1))
  yr <- rbind(c(0, 0, 0), c(1, 1, 0), c(1, 1, 1), c(2, 0, 1), c(0, 1, 2))
  pd <- .msr_gdr_exact_distance_bins(rep(40, 5), c(0, 20, 50, 75), "halfnorm")
  durations <- c(1, 2, 1); intercept <- qlogis(.3); base <- 10
  for (size in c(2, Inf)) {
    mean_limit <- base * exp(intercept) * sum(durations) * rowSums(pd)
    expected <- sum(vapply(1:5, function(i) {
      total <- if (is.infinite(size)) dpois(sum(yd[i, ]), mean_limit[i], log = TRUE) else
        dnbinom(sum(yd[i, ]), size, mu = mean_limit[i], log = TRUE)
      total + dmultinom(yd[i, ], prob = pd[i, ], log = TRUE) +
        dmultinom(yr[i, ], prob = durations, log = TRUE)
    }, 0))
    sections <- lapply(c(0, 2, 4, 8, 12, 16), function(shift) {
      pr <- .msr_gdr_exact_removal_bins(matrix(plogis(intercept - shift), 5, 3), durations)
      .msr_gdr_finite_score_bins(rep(base * exp(shift), 5), pd, pr, yd, yr, 40, size)
    })
    last <- tail(sections, 1L)[[1L]]
    expect_lt(abs(last$infinite_logLik - expected), 1e-5)
    expect_gt(last$infinite_logLik - last$logLik, 100)
    expect_true(all(vapply(sections, function(x) !x$fit_admission_established, logical(1))))
    expect_true(all(vapply(sections, function(x) identical(x$qualification,
      "private_conditional_fixed_parameter_reference_not_Gaussian_integrated"), logical(1))))
  }
})

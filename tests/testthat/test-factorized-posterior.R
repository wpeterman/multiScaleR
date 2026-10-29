test_that("conditional integration preserves known normal probabilities and both tails", {
  d <- multiScaleR:::.msr_factorized_distribution(function(z) dnorm(z, log = TRUE), 0)
  q <- c(-8, -2, 0, 2, 8)
  expect_equal(d$cdf(q), pnorm(q), tolerance = 1e-10)
  expect_equal(d$cdf(q, FALSE), pnorm(q, lower.tail = FALSE), tolerance = 1e-10)
  expect_equal(d$log_normalizer, 0, tolerance = 1e-10)
  expect_lte(d$relative_error, 1e-10)
  p <- c(1e-15, 1e-10, .025, .5, .975, 1 - 1e-10, 1 - 1e-15)
  expect_equal(d$quantile(p), qnorm(p), tolerance = 1e-7)
  expect_equal(log(d$cdf(8, FALSE)), pnorm(8, lower.tail = FALSE, log.p = TRUE), tolerance = 1e-9)
  expect_identical(d$quantile(c(0, 1)), c(-Inf, Inf))
  expect_identical(d$cdf(c(-Inf, Inf)), c(0, 1))
  expect_error(d$cdf(NA_real_), "nonmissing")
  expect_error(d$quantile(c(-.1, 1.1)), "between")
  expect_error(multiScaleR:::.msr_factorized_distribution(function(z) -z^2, 0, step = 0), "positive")
})

test_that("conditional CDFs retain separated modes and asymmetric probability mass", {
  logdensity <- function(z) {
    a <- log(.3) + dnorm(z, -3, .4, log = TRUE)
    b <- log(.7) + dnorm(z, 2, .8, log = TRUE)
    m <- max(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  d <- multiScaleR:::.msr_factorized_distribution(logdensity, 10, knots = c(-3, 2))
  q <- c(-5, -3, 0, 2, 5)
  truth <- .3 * pnorm(q, -3, .4) + .7 * pnorm(q, 2, .8)
  expect_equal(d$cdf(q), truth, tolerance = 1e-9)
  expect_equal(d$log_normalizer, 0, tolerance = 1e-9)
  p <- c(.025, .25, .5, .975)
  expect_equal(d$cdf(d$quantile(p)), p, tolerance = 1e-8)
})

test_that("annual reconstruction preserves the prior when landscape response is absent", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  for (unit in c(1, 1000)) {
    s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
      list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
    s$unit_conv <- unit
    p <- s$parameters; p$beta[length(p$beta)] <- 0; p$scale_z <- numeric()
    p$log_scale_sd <- log(.4)
    set.seed(1712); before <- .Random.seed
    a <- multiScaleR:::.msr_factorized_reconstruct(s, p, c(.2, .8))
    expect_identical(.Random.seed, before)
    expect_equal(a$latent_draw, qnorm(c(.2, .8)), tolerance = 1e-7)
    fixed <- p$scale_beta + log(unit)
    expect_equal(a$log_quantiles[, 1], fixed + .4 * qnorm(a$probs), tolerance = 1e-7)
    expect_equal(a$scale_draw, exp(a$log_draw), tolerance = 1e-12)
    expect_error(multiScaleR:::.msr_factorized_reconstruct(s, p, c(0, .5)), "open-unit")
    expect_error(multiScaleR:::.msr_factorized_annual_conditional(s, p, 3), "existing")
    bad <- p; bad$log_scale_sd <- -Inf
    expect_error(multiScaleR:::.msr_factorized_annual_conditional(s, bad, 1), "positive finite")
  }
})

test_that("shared scale reconstruction is deterministic and rejects spurious uniforms", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1), "gaussian", "key", 200)
  p <- s$parameters; p$scale_z <- numeric()
  a <- multiScaleR:::.msr_factorized_reconstruct(s, p)
  expected <- multiScaleR:::.msr_marginal_scale_predictions(s, p)$log_estimate
  expect_equal(a$log_quantiles, matrix(rep(expected, each = 3), nrow = 3))
  expect_null(a$log_draw)
  expect_length(a$log_group_normalizers, 0)
  expect_error(multiScaleR:::.msr_factorized_reconstruct(s, p, .5), "per annual")
})

test_that("shifted references agree with direct integration across response families", {
  skip_if_not_installed("RTMB")
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
      list(habitat = ~1 + (1 | year)), family, "key", 200)
    p <- s$parameters; p$beta <- c(.2, -.1, .4); p$scale_z <- numeric()
    shifted <- multiScaleR:::.msr_factorized_shifted_reference(s, p)
    direct <- multiScaleR:::.msr_factorized_objective(s, FALSE, 1e-10, .5)(p)
    expect_equal(shifted$value, direct, tolerance = 1e-8)
    expect_true(all(shifted$relative_error <= 1e-10))
    expect_true(shifted$infinite_support)
  }
})

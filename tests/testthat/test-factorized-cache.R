test_that("cached joint densities and conditional probabilities match ordinary evaluation", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 3, n_year = 2, kernel = kernel)
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
      list(habitat = ~1 + (1 | year)), family, "key", 200)
    p <- s$parameters; p$scale_z <- numeric(); p$beta <- c(.2, -.1, .4)
    before <- .Random.seed; config <- RTMB::TapeConfig()
    cache <- multiScaleR:::.msr_factorized_conditional_cache(s)
    expect_identical(.Random.seed, before)
    expect_identical(RTMB::TapeConfig(), config)
    for (group in 1:2) {
      d <- cache$log_density(p, group)
      g <- multiScaleR:::.msr_factorized_group(s, group)
      z <- c(-8, -2, 0, 2, 8)
      ordinary <- vapply(z, function(t) {
        q <- p; q$scale_z <- t
        -sum(multiScaleR:::.msr_factorized_response_nll(g, q)) + dnorm(t, log = TRUE)
      }, 0)
      expect_equal(vapply(z, d, 0), ordinary, tolerance = 1e-8)
      a <- multiScaleR:::.msr_factorized_annual_conditional(s, p, group)
      b <- multiScaleR:::.msr_factorized_annual_conditional(s, p, group, cache = cache)
      expect_equal(b$cdf(c(-3, 0, 3)), a$cdf(c(-3, 0, 3)), tolerance = 1e-9)
      expect_equal(b$cdf(c(-3, 0, 3), FALSE), a$cdf(c(-3, 0, 3), FALSE), tolerance = 1e-9)
      expect_equal(b$quantile(c(.025, .5, .975)), a$quantile(c(.025, .5, .975)), tolerance = 1e-7)
      expect_equal(b$log_normalizer, a$log_normalizer, tolerance = 1e-8)
      expect_lte(b$relative_error, 1e-10)
    }
    expect_identical(.Random.seed, before)
    expect_identical(RTMB::TapeConfig(), config)
  }
})

test_that("cached reconstruction preserves collapsed Gaussian groups and SD extremes", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  # Response groups are nested within years, so their Gaussian integral does
  # not connect the independent annual scale integrals.
  f$data$nested <- factor(rep(1:2, each = 3))
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | nested),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  p <- s$parameters; p$scale_z <- p$response_z <- numeric()
  p$beta <- c(.2, -.1, .4)
  cache <- multiScaleR:::.msr_factorized_conditional_cache(s)
  before <- .Random.seed
  for (tau in c(1e-7, .4, 20)) {
    p$log_scale_sd[] <- log(tau)
    a <- multiScaleR:::.msr_factorized_reconstruct(s, p, c(.2, .8))
    b <- multiScaleR:::.msr_factorized_reconstruct(s, p, c(.2, .8), cache = cache)
    expect_equal(b$log_quantiles, a$log_quantiles, tolerance = 1e-6)
    expect_equal(b$log_draw, a$log_draw, tolerance = 1e-6)
    expect_equal(b$log_group_normalizers, a$log_group_normalizers, tolerance = 1e-8)
    expect_identical(b$conditional_only, TRUE)
  }
  expect_identical(.Random.seed, before)
})

test_that("conditional caches reject stale specifications and incompatible inputs", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  p <- s$parameters; p$scale_z <- numeric()
  cache <- multiScaleR:::.msr_factorized_conditional_cache(s)
  altered <- s; altered$y[1] <- altered$y[1] + 1
  expect_error(multiScaleR:::.msr_factorized_annual_conditional(altered, p, 1, cache = cache), "match")
  expect_error(multiScaleR:::.msr_factorized_annual_conditional(s, p, 1, cache = list()), "match")
  for (group in list(0, 3, 1.5, 1e30, NA_real_, "1", 1:2))
    expect_error(cache$log_density(p, group), "compatible")
  missing <- p; missing$beta <- NULL
  for (bad in list(missing, unname(p), c(p, p[1])))
    expect_error(cache$log_density(bad, 1), "compatible")
  bad <- p; bad$beta <- c(p$beta, 0)
  expect_error(cache$log_density(bad, 1), "compatible")
  bad <- p; bad$beta[1] <- Inf
  expect_error(cache$log_density(bad, 1), "finite")
  bad <- p; bad$scale_z <- 0
  expect_error(cache$log_density(bad, 1), "compatible")
  density <- cache$log_density(p, 1)
  for (z in list(numeric(), c(0, 1), NA_real_, Inf, "0"))
    expect_error(density(z), "finite latent")
  shared <- multiScaleR:::.msr_marginal_reduce(s, "scale")
  expect_error(multiScaleR:::.msr_factorized_conditional_cache(shared), "one annual")
  # The independent likelihood audit must continue to use ordinary evaluation,
  # even when a cache has been built for this specification.
  expect_false("cache" %in% names(formals(multiScaleR:::.msr_factorized_shifted_reference)))
  expect_equal(multiScaleR:::.msr_factorized_shifted_reference(s, p)$value,
    multiScaleR:::.msr_factorized_objective(s, FALSE, 1e-10, .5)(p), tolerance = 1e-8)
})

test_that("prior simulation uses the prepared design without conditioning on responses", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 3)
  p <- prepare_bayesian_fixture(f)
  original <- serialize(p, NULL); rng <- .Random.seed
  a <- simulate_bayes_prior(p, nsim = 5, seed = 629)
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(p, NULL), original)
  expect_identical(a, simulate_bayes_prior(p, nsim = 5, seed = 629))
  invisible(simulate_bayes_prior(p, nsim = 1))
  expect_false(identical(.Random.seed, rng))
  expect_equal(dim(a$responses), c(nrow(f$data), 5))
  expect_equal(dim(a$means), dim(a$responses))
  expect_identical(a$observation_ids, p$specification$observation_ids)
  expect_identical(a$priors, p$specification$priors)
  expect_identical(rownames(a$log_scales$habitat), p$specification$scales$habitat$labels)
  p$specification$y <- p$specification$y + 100
  p$draws <- list(irrelevant = "posterior draws are ignored")
  expect_identical(a, simulate_bayes_prior(p, nsim = 5, seed = 629))
  expect_error(simulate_bayes_prior(list()), "object")
  expect_error(simulate_bayes_prior(p, 0), "nsim")
  expect_error(simulate_bayes_prior(p, 1.1), "nsim")
  expect_error(simulate_bayes_prior(p, seed = 0), "seed")
})

test_that("prior log distances retain map units under internal distance conversion", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 3)
  original <- prepare_bayesian_fixture(f)
  f$inputs$unit_conv <- f$inputs$unit_conv * 2
  f$inputs$binned$bin_rep <- f$inputs$binned$bin_rep / 2
  converted <- prepare_bayesian_fixture(f)
  a <- simulate_bayes_prior(original, nsim = 7, seed = 636)
  b <- simulate_bayes_prior(converted, nsim = 7, seed = 636)
  expect_equal(a$log_scales, b$log_scales, tolerance = 1e-12)
  expect_equal(a$means, b$means, tolerance = 1e-12)
  expect_equal(a$responses, b$responses)
})

test_that("prior scale draws have the declared shared and annual dependence", {
  f <- bayesian_scale_fixture(n_site = 2, n_year = 3)
  pr <- list(scale_median = 200, scale_intercept_sd = .3, scale_group_sd = .4,
    beta_sd = .1, dispersion_log_sd = .1)
  shared <- prepare_bayesian_fixture(f, scales = list(habitat = ~1), priors = pr)
  x <- simulate_bayes_prior(shared, nsim = 4000, seed = 630)$log_scales$habitat
  expect_equal(dim(x), c(1L, 4000L))
  expect_lt(abs(mean(x) - log(200)), .03)
  expect_lt(abs(stats::sd(as.numeric(x)) - .3), .02)
  pooled <- prepare_bayesian_fixture(f, priors = pr)
  x <- simulate_bayes_prior(pooled, nsim = 4000, seed = 631)$log_scales$habitat
  expect_lt(max(abs(rowMeans(x) - log(200))), .04)
  # Var(alpha + tau*z) = .3^2 + .4^2; covariance between years = Var(alpha).
  covariance <- stats::cov(t(x))
  expect_lt(max(abs(diag(covariance) - .25)), .03)
  expect_lt(max(abs(covariance[lower.tri(covariance)] - .09)), .03)
  independent <- prepare_bayesian_fixture(f, scales = list(habitat = ~0 + year), priors = pr)
  x <- simulate_bayes_prior(independent, nsim = 4000, seed = 632)$log_scales$habitat
  covariance <- stats::cov(t(x))
  expect_lt(max(abs(diag(covariance) - .09)), .02)
  expect_lt(max(abs(covariance[lower.tri(covariance)])), .02)
})

test_that("prior response means preserve group effects and likelihood units", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 3)
  f$data$effort <- 1
  p <- prepare_bayesian_fixture(f, family = "gaussian", priors = list(beta_sd = 1e-8))
  a <- simulate_bayes_prior(p, nsim = 6, seed = 633)
  same_site <- which(f$data$site == f$data$site[1])
  expect_lt(max(abs(as.matrix(a$means[same_site, ]) -
    matrix(as.numeric(unlist(a$means[1, ])), length(same_site), 6, byrow = TRUE))), 1e-6)
  expect_gt(max(abs(unlist(a$means[1, ]) - unlist(a$means[2, ]))), .01)
  local_mocked_bindings(.msr_bayes_components = function(spec, par)
    list(eta = rep(.4, length(spec$y))), .package = "multiScaleR")
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    p <- prepare_bayesian_fixture(f, family = "gaussian")
    p$specification$family <- family
    if (family == "binomial") p$specification$trials <- rep(7, nrow(f$data))
    if (family %in% c("poisson", "binomial")) p$specification$parameters$log_dispersion <- numeric()
    a <- simulate_bayes_prior(p, nsim = 4, seed = 634)
    expected <- switch(family, gaussian = .4, binomial = plogis(.4), exp(.4))
    expect_equal(unname(as.matrix(a$means)), matrix(expected, nrow(f$data), 4))
    expect_true(all(is.finite(as.matrix(a$responses))))
    if (family != "gaussian") {
      expect_true(all(as.matrix(a$responses) >= 0))
      expect_true(all(as.matrix(a$responses) == round(as.matrix(a$responses))))
    }
    if (family == "binomial") expect_true(all(as.matrix(a$responses) <= 7))
  }
})

test_that("prior simulation restores RNG state when a draw cannot be evaluated", {
  p <- prepare_bayesian_fixture(bayesian_scale_fixture())
  rng <- .Random.seed
  local_mocked_bindings(.msr_bayes_components = function(spec, par)
    list(eta = rep(Inf, length(spec$y))), .package = "multiScaleR")
  expect_error(simulate_bayes_prior(p, nsim = 2, seed = 635), "replicate 1.*No draw was replaced")
  expect_identical(.Random.seed, rng)
  rm(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", rng, envir = .GlobalEnv), add = TRUE)
  expect_error(simulate_bayes_prior(p, seed = 635), "Nonfinite")
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

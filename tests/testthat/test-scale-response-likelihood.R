test_that("shared response densities retain normalization and parameter meaning", {
  density <- multiScaleR:::.msr_scale_response_loglik
  eta <- c(-2, 0, 1.5)
  y <- c(0, 2, 4)
  ld <- log(1.7)
  expect_equal(density("gaussian", y, eta, ld), dnorm(y, eta, 1.7, log = TRUE))
  expect_equal(density("poisson", y, eta, numeric()), dpois(y, exp(eta), log = TRUE))
  expect_equal(density("binomial", y, eta, numeric(), c(0, 5, 7)),
    dbinom(y, c(0, 5, 7), plogis(eta), log = TRUE))
  expect_equal(density("nbinom2", y, eta, ld),
    dnbinom(y, mu = exp(eta), size = 1.7, log = TRUE))
  expect_equal(sum(exp(density("binomial", 0:7, .8, numeric(), 7))), 1)
  expect_equal(sum(exp(density("poisson", 0:100, .8, numeric()))), 1)
  expect_equal(sum(exp(density("nbinom2", 0:200, .8, ld))), 1, tolerance = 1e-10)
  expect_error(density("unknown", y, eta, ld), "Unsupported")
  for (e in c(-1000, 1000)) {
    expect_true(all(is.finite(density("binomial", c(0, 2, 7), e, numeric(), 7))))
    expect_equal(density("binomial", 0, e, numeric(), 0), 0)
  }
})

test_that("shared AD densities have independent values, scores and curvature", {
  skip_if_not_installed("RTMB")
  density <- multiScaleR:::.msr_scale_response_loglik
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    y <- c(0, 2, 4); trials <- c(0, 5, 7)
    p <- list(eta = c(-2, .2, 1.5), log_dispersion = log(1.7))
    ad <- RTMB::MakeADFun(function(q)
      -sum(density(family, y, q$eta, q$log_dispersion, trials, TRUE)), p, silent = TRUE)
    reference <- function(x) {
      eta <- x[1:3]; ld <- x[4]
      -sum(switch(family,
        gaussian = dnorm(y, eta, exp(ld), log = TRUE),
        poisson = dpois(y, exp(eta), log = TRUE),
        binomial = dbinom(y, trials, plogis(eta), log = TRUE),
        nbinom2 = dnbinom(y, mu = exp(eta), size = exp(ld), log = TRUE)))
    }
    expect_equal(unname(ad$fn(ad$par)), reference(ad$par), tolerance = 1e-11)
    h <- 1e-5
    fd <- vapply(seq_along(ad$par), function(j) {
      hi <- lo <- ad$par; hi[j] <- hi[j] + h; lo[j] <- lo[j] - h
      (reference(hi) - reference(lo)) / (2 * h)
    }, 0)
    expect_equal(as.numeric(ad$gr(ad$par)), fd, tolerance = 1e-7)
    expect_equal(unname(ad$he(ad$par)), unname(optimHess(ad$par, reference)), tolerance = 1e-5)
  }
})

test_that("posterior and likelihood adapters retain their different targets", {
  skip_if_not_installed("RTMB")
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    fit <- prepare_bayesian_fixture(f, family = family)
    s <- fit$specification; p <- s$parameters
    p$beta <- c(.3, -.1, .7)
    p$scale_z[] <- c(-.2, .1, .3)
    p$response_z[] <- seq(-.3, .3, length.out = length(p$response_z))
    eta <- multiScaleR:::.msr_bayes_components(s, p)$eta
    response <- -switch(family,
      gaussian = dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE),
      poisson = dpois(s$y, exp(eta), log = TRUE),
      binomial = dbinom(s$y, s$trials, plogis(eta), log = TRUE),
      nbinom2 = dnbinom(s$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE))
    expect_equal(unname(multiScaleR:::.msr_marginal_response_nll(s, p)),
      unname(response), tolerance = 1e-11)
    latent <- -sum(dnorm(p$scale_z, log = TRUE)) - sum(dnorm(p$response_z, log = TRUE))
    expect_equal(multiScaleR:::.msr_marginal_objective(s)(p), sum(response) + latent)
    prior <- s$priors
    sm <- s$scales$habitat
    penalty <- -sum(dnorm(p$beta, 0, prior$beta_sd, log = TRUE)) -
      sum(dnorm(p$scale_beta, sm$prior_mean, sm$prior_sd, log = TRUE)) -
      sum(dnorm(exp(p$log_scale_sd), 0, prior$scale_group_sd, log = TRUE) + log(2) + p$log_scale_sd) -
      sum(dnorm(exp(p$log_response_sd), 0, prior$response_group_sd, log = TRUE) + log(2) + p$log_response_sd) -
      sum(dnorm(p$log_dispersion, prior$dispersion_log_mean, prior$dispersion_log_sd, log = TRUE))
    expect_equal(multiScaleR:::.msr_bayes_objective(s)(p), sum(response) + latent + penalty,
      tolerance = 1e-10)
    expect_gt(abs(penalty), .01)
  }
})

test_that("marginal preparation freezes the design without statistical priors", {
  f <- bayesian_scale_fixture(n_site = 5, n_year = 3)
  spec <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | site) + offset(log(effort)),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), row_id = "key", reference_scale = 200)
  b <- multiScale_bayes(count ~ habitat + temp + (1 | site) + offset(log(effort)), f$inputs, f$data,
    scales = list(habitat = ~1 + (1 | year)), row_id = "key", method = "prepare",
    priors = list(scale_median = 200))$specification
  expect_null(spec$priors)
  expect_null(spec$scales[[1]]$prior_mean)
  expect_null(spec$scales[[1]]$prior_sd)
  for (nm in c("rows", "V", "C", "delta", "center", "spread", "X_fixed", "y", "trials", "offset"))
    expect_identical(spec[[nm]], b[[nm]])
  expect_identical(spec$scales[[1]]$X, b$scales[[1]]$X)
  expect_equal(spec$reference_scale, 200)
  expect_equal(spec$parameters$log_scale_sd, log(.2))
  reverse <- nrow(f$data):1
  reversed <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | site) + offset(log(effort)),
    f$inputs, f$data[reverse, ], list(habitat = ~1 + (1 | year)), row_id = "key")
  expect_identical(reversed$rows, rev(spec$rows))
  expect_equal(unname(spec$offset), log(f$data$effort))
  expect_error(multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    row_id = "key", reference_scale = 0), "reference_scale")
  expect_error(multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    row_id = "key", priors = list()), "unused argument")
})

test_that("shared scales, annual cells and fixed trends retain their formula meaning", {
  f <- bayesian_scale_fixture(n_site = 5, n_year = 3)
  f$data$trend <- as.numeric(f$data$year) - 2
  shared <- multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    list(habitat = ~1), "gaussian", "key", 200)
  cells <- multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    list(habitat = ~0 + year), "gaussian", "key", 200)
  trend <- multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    list(habitat = ~1 + trend + (1 | year)), "gaussian", "key", 200)
  expect_length(shared$parameters$scale_z, 0)
  expect_length(cells$parameters$scale_beta, 3)
  expect_length(cells$parameters$scale_z, 0)
  expect_length(trend$parameters$scale_beta, 2)
  expect_length(trend$parameters$scale_z, 3)
  trend$parameters$scale_beta[2] <- .1
  reduced <- multiScaleR:::.msr_marginal_reduce(trend, "scale")
  estimates <- multiScaleR:::.msr_marginal_scale_predictions(reduced, reduced$parameters)
  expect_equal(diff(estimates$log_estimate), c(.1, .1))
  expect_identical(multiScaleR:::.msr_marginal_reduce(reduced, "scale"), reduced)
})

test_that("the joint marginal objective retains only the required random densities", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 4, n_year = 3, kernel = kernel)
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    spec <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | site),
      f$inputs, f$data, list(habitat = ~1 + (1 | year)), family, "key", 200)
    p <- spec$parameters; p$beta <- c(.2, -.1, .4)
    p$scale_z[] <- c(-.2, 0, .3); p$response_z[] <- seq_along(p$response_z) / 10
    p$log_scale_sd[] <- log(.3); p$log_response_sd[] <- log(.4); p$log_dispersion[] <- log(.7)
    eta <- multiScaleR:::.msr_bayes_components(spec, p)$eta
    ll <- switch(family, gaussian = dnorm(spec$y, eta, .7, log = TRUE),
      poisson = dpois(spec$y, exp(eta), log = TRUE),
      binomial = dbinom(spec$y, spec$trials, plogis(eta), log = TRUE),
      nbinom2 = dnbinom(spec$y, mu = exp(eta), size = .7, log = TRUE))
    expected <- -sum(ll) - sum(dnorm(p$scale_z, log = TRUE)) - sum(dnorm(p$response_z, log = TRUE))
    fn <- multiScaleR:::.msr_marginal_objective(spec)
    expect_equal(fn(p), expected, tolerance = 1e-10)
    reference_spec <- multiScaleR:::.msr_marginal_reduce(spec, "response")
    no_response <- p; no_response$response_z <- no_response$log_response_sd <- numeric()
    expect_equal(unname(marginal_reference_nll(reference_spec, no_response)),
      unname(multiScaleR:::.msr_marginal_response_nll(reference_spec, no_response)), tolerance = 1e-10)
    ad <- RTMB::MakeADFun(multiScaleR:::.msr_marginal_objective(spec, TRUE), p, silent = TRUE)
    expect_equal(as.numeric(ad$fn(ad$par)), expected, tolerance = 1e-9)
    h <- 1e-5
    fd <- vapply(seq_along(ad$par), function(j) {
      hi <- lo <- ad$par; hi[j] <- hi[j] + h; lo[j] <- lo[j] - h
      (ad$fn(hi) - ad$fn(lo)) / (2 * h)
    }, 0)
    expect_equal(as.numeric(ad$gr(ad$par)), fd, tolerance = 1e-5)

    b <- multiScale_bayes(count ~ habitat + temp + (1 | site), f$inputs, f$data,
      scales = list(habitat = ~1 + (1 | year)), row_id = "key", family = family,
      method = "prepare", priors = list(scale_median = 200))$specification
    pr <- b$priors; sm <- b$scales[[1]]
    penalty <- -sum(dnorm(p$beta, 0, pr$beta_sd, log = TRUE)) -
      sum(dnorm(p$scale_beta, sm$prior_mean, sm$prior_sd, log = TRUE)) -
      sum(dnorm(exp(p$log_scale_sd), 0, pr$scale_group_sd, log = TRUE) + log(2) + p$log_scale_sd) -
      sum(dnorm(exp(p$log_response_sd), 0, pr$response_group_sd, log = TRUE) + log(2) + p$log_response_sd) -
      sum(dnorm(p$log_dispersion, pr$dispersion_log_mean, pr$dispersion_log_sd, log = TRUE))
    expect_equal(multiScaleR:::.msr_bayes_objective(b)(p) - fn(p), penalty, tolerance = 1e-9)
  }
})

test_that("Gaussian response-effect integration matches the analytic likelihood", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | site),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  s <- multiScaleR:::.msr_marginal_reduce(s, "scale")
  p <- s$parameters; p$beta <- c(.5, .2, 1.4)
  p$log_response_sd[] <- log(.4); p$log_dispersion[] <- log(.7)
  s$parameters <- p
  obj <- multiScaleR:::.msr_marginal_adfun(s)
  value <- as.numeric(obj$fn(obj$par))
  p$response_z[] <- 0
  eta <- multiScaleR:::.msr_bayes_components(s, p)$eta
  group <- s$response$random[[1]]$index
  Sigma <- diag(.7^2, length(group)) + .4^2 * outer(group, group, `==`)
  residual <- s$y - eta
  exact <- length(group) / 2 * log(2 * pi) +
    as.numeric(determinant(Sigma, logarithm = TRUE)$modulus) / 2 +
    drop(crossprod(residual, solve(Sigma, residual))) / 2
  expect_equal(value, exact, tolerance = 1e-8)

  centered_fn <- function(par) {
    converted <- par
    converted$response_z <- par$response_z / exp(par$log_response_sd)
    sum(multiScaleR:::.msr_marginal_response_nll(s, converted, TRUE)) -
      sum(RTMB::dnorm(par$response_z, 0, exp(par$log_response_sd), log = TRUE))
  }
  centered <- RTMB::MakeADFun(centered_fn, s$parameters, random = "response_z", silent = TRUE)
  expect_equal(as.numeric(centered$fn(centered$par)), exact, tolerance = 1e-8)
  expect_equal(as.numeric(centered$gr(centered$par)), as.numeric(obj$gr(obj$par)), tolerance = 1e-7)
})

test_that("annual integration has the right zero-SD limit and a quadrature reference", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  s$parameters$beta <- c(.5, .1, 1)
  s$parameters$log_dispersion[] <- log(.8)
  s$parameters$log_scale_sd[] <- log(1e-7)
  a <- multiScaleR:::.msr_marginal_adfun(s)
  reduced <- multiScaleR:::.msr_marginal_reduce(s, "scale")
  b <- multiScaleR:::.msr_marginal_adfun(reduced)
  expect_equal(as.numeric(a$fn(a$par)), as.numeric(b$fn(b$par)), tolerance = 1e-8)
  expect_length(reduced$parameters$scale_z, 0)
  expect_length(reduced$parameters$log_scale_sd, 0)
  expect_true(reduced$zero_sd[["scale"]])
  expect_equal(length(unique(multiScaleR:::.msr_marginal_scale_predictions(reduced,
    reduced$parameters)$estimate)), 1L)
  s$parameters$log_scale_sd[] <- log(.02)
  a <- multiScaleR:::.msr_marginal_adfun(s)
  q1 <- marginal_annual_integral(s, s$parameters)
  q2 <- marginal_annual_integral(s, s$parameters, rel.tol = 1e-10)
  expect_equal(q1, q2, tolerance = 1e-7)
  expect_equal(as.numeric(a$fn(a$par)), q1, tolerance = 1e-4)
  coupled <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + (1 | site),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), "gaussian", "key")
  expect_error(marginal_annual_integral(coupled, coupled$parameters))
  # Tail probes over the real line must preserve the nearest supported bin.
  s$parameters$scale_z[] <- -1e6
  expect_true(all(is.finite(marginal_reference_nll(s, s$parameters))))
})

test_that("independent conditional refinement verifies a unique Gaussian mode", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 10, n_year = 4)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | site),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  s <- multiScaleR:::.msr_marginal_reduce(s, "scale")
  p <- s$parameters; p$beta <- c(.5, .2, 1.4)
  p$log_response_sd[] <- log(.4); p$log_dispersion[] <- log(.07)
  p$response_z[] <- 0
  eta <- multiScaleR:::.msr_bayes_components(s, p)$eta
  term <- s$response$random[[1L]]
  residual_sum <- as.numeric(rowsum(s$y - eta, term$index, reorder = TRUE))
  sizes <- as.numeric(table(term$index))
  p$response_z[] <- .4 * residual_sum / (.07^2 + sizes * .4^2)
  audit <- multiScaleR:::.msr_marginal_inner_audit(s, p)
  expect_true(audit$passed)
  expect_lt(audit$gradient, 1e-7)
  expect_true(all(vapply(audit$starts, function(x) isTRUE(x$refinement$passed), TRUE)))
  expect_true(all(vapply(audit$starts, function(x) max(abs(x$par - p$response_z)) < 1e-7, TRUE)))
})

test_that("the private fitter recovers a positive year SD with shared site effects", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 20, n_year = 5, seed = 8321)
  response <- count ~ habitat + temp + (1 | site)
  s <- multiScaleR:::.msr_marginal_prepare(response, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  p <- s$parameters; p$beta <- c(.5, .2, 2)
  p$scale_z[] <- c(-1.3, -.6, .1, .7, 1.1); p$log_scale_sd[] <- log(.4)
  p$response_z[] <- cos(seq_along(p$response_z)); p$log_response_sd[] <- log(.25)
  truth <- multiScaleR:::.msr_marginal_scale_predictions(s, p)$log_estimate
  f$data$count <- multiScaleR:::.msr_bayes_components(s, p)$eta + .15 * sin(seq_len(nrow(f$data)))
  s <- multiScaleR:::.msr_marginal_prepare(response, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  fit <- multiScaleR:::.msr_marginal_fit(s)
  expect_true(fit$diagnostics$passed)
  expect_true(all(!fit$zero_sd))
  expect_gt(unname(fit$group_sd["scale"]), .2)
  expect_lt(sqrt(mean((fit$annual$log_estimate - truth)^2)), .1)
  expect_false(fit$uncertainty$available)
})

test_that("the private fitter retains zero-SD fits without claiming calibrated intervals", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 12, n_year = 4, annual = FALSE)
  f$data$temp <- rep(seq(-1, 1, length.out = 12), 4)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  p <- s$parameters; p$beta <- c(.5, .2, 2)
  f$data$count <- multiScaleR:::.msr_bayes_components(s, p)$eta +
    rep(.15 * sin(1:12), 4)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  before <- .Random.seed
  fit <- multiScaleR:::.msr_marginal_fit(s)
  expect_identical(.Random.seed, before)
  expect_s3_class(fit, "msr_marginal_prototype")
  expect_true(fit$diagnostics$passed)
  expect_true(fit$zero_sd[["scale"]])
  expect_equal(unname(fit$group_sd["scale"]), 0)
  expect_length(fit$attempts, 6)
  expect_false(fit$uncertainty$available)
  expect_true(all(!fit$annual$interval_available))
  expect_equal(length(unique(fit$annual$estimate)), 1L)
  repeat_fit <- multiScaleR:::.msr_marginal_fit(s, list(starts = 2L))
  expect_true(repeat_fit$diagnostics$passed)
  expect_equal(fit$log_likelihood, repeat_fit$log_likelihood, tolerance = 1e-5)
  expect_equal(unserialize(serialize(fit, NULL))$annual, fit$annual)
  expect_error(multiScaleR:::.msr_marginal_fit(s, list(prior_sd = 2)), "control")
})

test_that("unbounded and unsupported likelihood cases remain explicit", {
  f <- bayesian_scale_fixture(n_site = 5, n_year = 3)
  f$data$count[] <- 0
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "poisson", "key")
  fit <- multiScaleR:::.msr_marginal_fit(s)
  expect_false(fit$diagnostics$passed)
  expect_match(fit$diagnostics$reason, "No finite")
  expect_length(fit$attempts, 0)
  expect_false(multiScaleR:::.msr_marginal_curvature(diag(c(1, -1)))$passed)
  expect_false(multiScaleR:::.msr_marginal_curvature(diag(c(1, 1e-14)))$passed)
  expect_error(multiScaleR:::.msr_marginal_reduce(s, "unknown"), "zero")
  expect_error(multiScaleR:::.msr_marginal_reduce(s, "response"), "absent")
})

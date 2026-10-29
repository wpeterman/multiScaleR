test_that("marginal Laplace reproduces an exactly integrable random-effect model", {
  skip_if_not_installed("RTMB")
  y <- .7; se <- .3; prior_sd <- .5
  f <- function(p) {
    tau <- exp(p$log_scale_sd)
    -RTMB::dnorm(y, tau * p$z, se, log = TRUE) -
      RTMB::dnorm(p$z, 0, 1, log = TRUE) -
      RTMB::dnorm(tau, 0, prior_sd, log = TRUE) - log(2) - p$log_scale_sd
  }
  obj <- RTMB::MakeADFun(f, list(z = 0, log_scale_sd = log(.5)), silent = TRUE)
  hyper <- which(names(obj$par) == "log_scale_sd")
  ctl <- list(starts = 3L, iter.max = 1000L, max_grid = 129L,
              cdf_tol = .005, quantile_tol = .02, logz_tol = .01)
  evaluate <- function(x, warm) {
    tau <- if (x == 0) prior_sd * 1e-8 else if (x == 1) prior_sd * 8 else
      prior_sd * qnorm((1 + x) / 2)
    z <- multiScaleR:::.msr_bayes_marginal_node(obj, hyper, log(tau), prior_sd, warm, ctl)
    expect_equal(z$log_density, dnorm(y, 0, sqrt(se^2 + tau^2), log = TRUE), tolerance = 1e-8)
    z
  }
  grid <- multiScaleR:::.msr_bayes_marginal_grid(evaluate, prior_sd, ctl)
  expect_true(grid$converged)
  expect_equal(sum(grid$weights), 1, tolerance = 1e-12)
  expect_equal(multiScaleR:::.msr_bayes_marginal_cdf(grid, c(0, 1)), c(0, 1), tolerance = 1e-12)
  density <- function(tau) dnorm(y, 0, sqrt(se^2 + tau^2)) * 2 * dnorm(tau, 0, prior_sd)
  normalizer <- integrate(density, 0, Inf, rel.tol = 1e-10)$value
  quantiles <- multiScaleR:::.msr_bayes_marginal_quantile(grid, c(.025, .5, .975))
  truth <- vapply(quantiles, function(q) integrate(density, 0, q)$value / normalizer, 0)
  expect_equal(truth, c(.025, .5, .975), tolerance = .004)
  set.seed(42); before <- .Random.seed
  draws <- multiScaleR:::.msr_bayes_marginal_draws(grid, hyper, 30000L, 112)
  expect_identical(.Random.seed, before)
  tau <- exp(draws[, hyper])
  expect_gt(length(unique(tau)), 29000)
  expect_equal(mean(tau < .5), integrate(density, 0, .5)$value / normalizer, tolerance = .012)
  true_z <- integrate(function(t) y * t / (se^2 + t^2) * density(t), 0, Inf)$value / normalizer
  expect_equal(mean(draws[, -hyper]), true_z, tolerance = .025)
})

test_that("marginal SD quadrature retains the half-normal prior when data are uninformative", {
  ctl <- list(max_grid = 65L, cdf_tol = .005, quantile_tol = .02, logz_tol = .01)
  evaluate <- function(x, warm) list(log_density = 0, mode = c(z = 0, log_scale_sd = 0),
    root = matrix(1), competing = FALSE)
  grid <- multiScaleR:::.msr_bayes_marginal_grid(evaluate, .5, ctl)
  expect_true(grid$converged)
  expect_equal(multiScaleR:::.msr_bayes_marginal_cdf(grid, seq(0, 1, .1)), seq(0, 1, .1), tolerance = 1e-12)
  expect_equal(multiScaleR:::.msr_bayes_marginal_quantile(grid, c(.025, .5, .975)),
               .5 * qnorm((1 + c(.025, .5, .975)) / 2), tolerance = 1e-10)
})

test_that("marginal public inference rejects unsupported quadrature dimensions", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture()
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", scales = list(habitat = ~1), method = "marginal"), "exactly one")
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", scales = list(habitat = ~1 + (1 | year)),
    method = "marginal", control = list(max_grid = 20)), "max_grid")
})

test_that("marginal public fits retain continuous pooling uncertainty and fitted methods", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 10, n_year = 3)
  prepared <- multiScale_bayes(count ~ habitat + temp, f$inputs, f$data,
    family = "gaussian", row_id = "key", scales = list(habitat = ~1 + (1 | year)),
    method = "prepare")
  p <- prepared$parameters
  p$beta <- c(.5, .3, 1.5)
  p$scale_z <- c(-.4, 0, .5)
  p$log_scale_sd <- log(.3)
  f$data$count <- multiScaleR:::.msr_bayes_components(prepared$specification, p)$eta +
    .3 * sin(seq_len(nrow(f$data)))
  set.seed(713); before <- .Random.seed
  fit <- suppressWarnings(multiScale_bayes(count ~ habitat + temp, f$inputs, f$data,
    family = "gaussian", row_id = "key", scales = list(habitat = ~1 + (1 | year)),
    method = "marginal", seed = 913, laplace_draws = 300))
  expect_identical(.Random.seed, before)
  expect_true(fit$diagnostics$grid_converged)
  expect_identical(fit$inference_metadata$method, "marginal_laplace")
  expect_identical(tail(names(fit), 1), "inference_metadata")
  expect_identical(fit$inference_metadata$target, "posterior")
  expect_null(fit$sampler)
  expect_gt(length(unique(fit$draws$log_scale_sd[, 1])), 290)
  expect_equal(nrow(posterior_scales(fit)), 3)
  expect_equal(nrow(predict(fit, draws = 3)), nrow(f$data))
  expect_equal(dim(simulate(fit, nsim = 3, seed = 1)), c(nrow(f$data), 3))
  expect_output(print(fit), "Marginal Laplace numerical checks")
  expect_equal(posterior_scales(unserialize(serialize(fit, NULL))), posterior_scales(fit))
})

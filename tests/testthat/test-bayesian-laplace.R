test_that("Laplace covariance and draws retain joint Gaussian dependence", {
  H <- matrix(c(4, 1, 1, 2), 2)
  g <- multiScaleR:::.msr_bayes_laplace_geometry(H)
  expect_equal(g$covariance, solve(H), tolerance = 1e-12)
  set.seed(917)
  before <- .Random.seed
  x <- multiScaleR:::.msr_bayes_laplace_draws(c(2, -1), g$root, 30000, 61)
  expect_identical(.Random.seed, before)
  expect_equal(colMeans(x), c(2, -1), tolerance = .02)
  expect_equal(stats::cov(x), solve(H), tolerance = .015)
  expect_identical(x, multiScaleR:::.msr_bayes_laplace_draws(c(2, -1), g$root, 30000, 61))
  expect_error(multiScaleR:::.msr_bayes_laplace_geometry(diag(c(1, -1))),
               "not positive definite")
  expect_error(multiScaleR:::.msr_bayes_laplace_geometry(diag(c(1, 1e-14))),
               "ill-conditioned")
  failed <- tryCatch(multiScaleR:::.msr_bayes_laplace_failure("Failed curvature",
    list(list(par = c(1, 2), convergence = 1L))), msr_laplace_error = function(e) e)
  expect_s3_class(failed, "error")
  expect_identical(failed$starts[[1]]$par, c(1, 2))
})

test_that("Laplace inference uses the sampler posterior and existing fitted methods", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 12, n_year = 3)
  # A strong continuous response gives a stable small integration example.
  prepared <- multiScale_bayes(count ~ habitat + temp, f$inputs, f$data,
    family = "gaussian", scales = list(habitat = ~0 + year), row_id = "key",
    method = "prepare", priors = list(scale_median = 200))
  p <- prepared$parameters
  p$beta <- c(.5, .3, 2)
  f$data$count <- multiScaleR:::.msr_bayes_components(prepared$specification, p)$eta +
    .1 * sin(seq_len(nrow(f$data)))
  set.seed(312)
  previous <- .Random.seed
  fit <- multiScale_bayes(count ~ habitat + temp, f$inputs, f$data,
    family = "gaussian", scales = list(habitat = ~0 + year), row_id = "key",
    method = "laplace", laplace_draws = 300, seed = 31,
    priors = list(scale_median = 200))
  expect_identical(.Random.seed, previous)
  expect_true(fit$diagnostics$passed)
  expect_null(fit$sampler)
  expect_identical(fit$inference_metadata$method, "joint_laplace")
  expect_identical(tail(names(fit), 1), "inference_metadata")
  expect_identical(fit$inference_metadata$uncertainty$intervals, "marginal_credible")
  expect_equal(nrow(fit$draws$beta), 300)
  expect_equal(ncol(fit$draws$scale_z), 0)
  expect_equal(fit$laplace$covariance %*% fit$laplace$hessian,
               diag(length(fit$laplace$mode)), tolerance = 1e-8)
  obj <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(fit$specification, TRUE),
                         fit$parameters, silent = TRUE)
  expect_lt(max(abs(obj$gr(fit$laplace$mode))), 1e-4)
  expect_true(any(vapply(fit$laplace$starts, function(z) z$newton_steps > 0L, logical(1))))
  expect_equal(obj$he(fit$laplace$mode), fit$laplace$hessian, tolerance = 1e-8)
  expect_equal(nrow(posterior_scales(fit)), 3)
  expect_equal(nrow(predict(fit, draws = 3)), nrow(f$data))
  expect_equal(dim(simulate(fit, nsim = 3, seed = 1)), c(nrow(f$data), 3))
  expect_identical(summary(fit)$inference, "joint Laplace Gaussian approximation")
  expect_output(print(fit), "Laplace numerical checks")
  expect_equal(predict(unserialize(serialize(fit, NULL)), draws = 3),
               predict(fit, draws = 3))
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key",
    method = "laplace", family = "gaussian", control = list(adapt_delta = .9)), "Laplace control")
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key",
    method = "laplace", family = "gaussian", laplace_draws = 1), "laplace_draws")
})

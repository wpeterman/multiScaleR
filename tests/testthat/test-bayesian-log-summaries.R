test_that("log-distance summaries preserve default distances and define their units", {
  f <- bayesian_scale_fixture()
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f), n = 9)
  location <- fit$parameters$scale_beta[1]
  fit$draws$scale_beta[, 1] <- location + seq(-.4, .4, length.out = 9)
  distance <- posterior_scales(fit, level = .5)
  expect_identical(distance, posterior_scales(fit, level = .5, scale = "distance"))
  logs <- posterior_scales(fit, level = .5, scale = "log")
  expected <- fit$draws$scale_beta[, 1] + log(fit$specification$unit_conv)
  expect_equal(logs$mean, rep(mean(expected), nrow(logs)))
  expect_equal(logs$sd, rep(sd(expected), nrow(logs)))
  # With nine draws these interval probabilities land on order statistics,
  # avoiding differences caused by quantile interpolation after transformation.
  expect_equal(exp(logs$median), distance$median)
  expect_equal(exp(logs$lower), distance$lower)
  expect_equal(exp(logs$upper), distance$upper)
  expect_equal(logs$p_below_min_D, distance$p_below_min_D)
  expect_equal(logs$p_above_max_D, distance$p_above_max_D)
  expect_error(posterior_scales(fit, scale = "unknown"), "arg")
})

test_that("log summaries remain finite beyond natural-distance floating-point range", {
  f <- bayesian_scale_fixture()
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f), n = 9)
  fit$draws$scale_beta[, 1] <- 800 + seq(-.4, .4, length.out = 9)
  logs <- posterior_scales(fit, scale = "log")
  expect_true(all(is.finite(as.matrix(logs[c("mean", "median", "sd", "lower", "upper")]))))
  expect_true(all(logs$p_above_max_D == 1))
  expect_true(all(logs$p_below_min_D == 0))
  expect_gt(min(logs$median), 800)
})

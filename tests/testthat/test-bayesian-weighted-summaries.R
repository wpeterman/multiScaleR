test_that("weighted inverse CDF respects exact cumulative boundaries", {
  expect_equal(multiScaleR:::.msr_bayes_weighted_quantile(c(10, 20, 30),
    c(.1, .4, .5), c(0, .1, .5, 1)), c(10, 10, 20, 30))
})

test_that("direct weighted annual summaries use joint candidates and preserve defaults", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f), n = 4)
  fit$draws$scale_beta[, 1] <- log(c(.2, .4, .8, 1.6))
  fit$draws$scale_z[, ] <- matrix(c(-.3, 0, .2, .6, .4, -.2, 0, -.1), 4)
  candidates <- do.call(cbind, fit$draws)
  w <- c(.1, .2, .3, .4)
  fit$importance <- list(candidates = candidates, weights = w)
  original <- posterior_scales(fit)
  expect_identical(original, posterior_scales(fit, estimator = "draws"))
  weighted <- posterior_scales(fit, estimator = "weighted", scale = "log", level = .8)
  sm <- fit$specification$scales$habitat
  logs <- fit$draws$scale_beta[, 1] + fit$draws$scale_z *
    exp(fit$draws$log_scale_sd[, 1]) + log(fit$specification$unit_conv)
  expected <- apply(logs, 2, multiScaleR:::.msr_bayes_weighted_quantile,
                    weights = w, probs = c(.1, .5, .9))
  expect_equal(weighted$median, unname(expected[2, ]))
  expect_equal(weighted$lower, unname(expected[1, ]))
  expect_equal(weighted$upper, unname(expected[3, ]))
  expect_equal(weighted$mean, unname(colSums(logs * w)))
  expect_equal(weighted$sd, unname(sqrt(colSums(sweep(logs, 2, colSums(logs*w), "-")^2*w))))
  distance <- posterior_scales(fit, estimator = "weighted", level = .8)
  expect_equal(exp(weighted$median), distance$median)
  expect_equal(exp(weighted$lower), distance$lower)
  expect_equal(exp(weighted$upper), distance$upper)
  expect_equal(weighted$p_above_max_D, colSums((logs > log(f$inputs$max_D))*w))
  expect_equal(weighted$p_below_min_D, colSums((logs < log(f$inputs$min_D))*w))
  weighted_before <- posterior_scales(fit, estimator = "weighted", scale = "log")
  fit$draws$scale_beta[, 1] <- 800
  expect_identical(posterior_scales(fit, estimator = "weighted", scale = "log"), weighted_before)
  expect_true(all(is.finite(weighted_before$mean)))
  expect_error(posterior_scales(fit, estimator = "weighted", draws = 2), "draws = NULL")
  fit$importance$weights[1] <- -1
  expect_error(posterior_scales(fit, estimator = "weighted"), "Invalid retained")
  fit$importance <- NULL
  expect_error(posterior_scales(fit, estimator = "weighted"), "importance-corrected")
})

test_that("zero target mass does not contaminate weighted distance summaries", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f), n = 3)
  fit$draws$scale_beta[, 1] <- c(800, log(.2), log(.4))
  fit$importance <- list(candidates = do.call(cbind, fit$draws), weights = c(0, .4, .6))
  out <- posterior_scales(fit, estimator = "weighted")
  expect_true(all(is.finite(as.matrix(out[c("mean", "median", "sd", "lower", "upper")]))))
  expect_equal(out$mean, rep(.32 * fit$specification$unit_conv, 2))
  expect_equal(out$median, rep(.4 * fit$specification$unit_conv, 2))
})

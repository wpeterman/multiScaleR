test_that("Bayesian effect curves retain links, offsets, and fixed reference units", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  for (kernel in c("gaussian", "exp")) {
    f <- refit_spatial_fixture(kernel)
    for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
      dat <- f$data
      if (family == "binomial") dat$y <- dat$y %% 2
      fit <- multiScale_bayes(y ~ rough + z + group + (1 | site) + offset(log(exposure)),
        f$inputs, dat, scales = list(rough = ~1), family = family, method = "prepare")
      fit <- fake_bayesian_draws(fit, n = 7)
      columns <- fit$specification$beta_names
      fit$draws$beta[,] <- 0
      fit$draws$beta[, match("(Intercept)", columns)] <- seq(-.3, .3, length.out = 7)
      fit$draws$beta[, match("rough", columns)] <- .4
      fit$draws$beta[, match("z", columns)] <- .2
      fit$diagnostics <- list(passed = FALSE)
      before <- fit
      plots <- expect_no_warning(plot_marginal_effects(fit, length.out = 5))
      expect_setequal(names(plots), c("rough", "z", "group"))
      expect_identical(fit, before)
      expect_identical(as.character(plots$group$data$x), c("a", "b", "c"))
      x <- seq(min(fit$specification$response$frame$rough),
        max(fit$specification$response$frame$rough), length.out = 5)
      eta <- outer(.4 * x + .2 * mean(dat$z) + log(mean(dat$exposure)),
                    fit$draws$beta[, match("(Intercept)", columns)], "+")
      mu <- switch(family, gaussian = eta, binomial = plogis(eta), exp(eta))
      expected <- t(apply(mu, 1, quantile, probs = c(.025, .5, .975), names = FALSE))
      expect_equal(plots$rough$data$fit, expected[, 2], tolerance = 1e-12)
      expect_equal(plots$rough$data$lwr, expected[, 1], tolerance = 1e-12)
      expect_equal(plots$rough$data$upr, expected[, 3], tolerance = 1e-12)
      expect_equal(plots$rough$data$x, x * fit$specification$spread[["rough"]] +
        fit$specification$center[["rough"]], tolerance = 1e-12)
      expect_match(plots$rough$labels$caption, "Fit checks need attention")
      contract <- attr(plots, "prediction_contract")
      expect_identical(contract$inference_target, "posterior")
      expect_identical(contract$uncertainty$kernel_parameters,
                        "covariate_values_fixed_not_reevaluated")
      expect_false(contract$uncertainty$observation_noise)
      expect_identical(contract$diagnostics_status, "failed")
      changed <- fit
      changed$draws$scale_beta[,] <- changed$draws$scale_beta + 2
      other <- plot_marginal_effects(changed, length.out = 5)
      expect_equal(plots$rough$data, other$rough$data)
    }
  }
})

test_that("Bayesian curve draw selection is bounded and explicit", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- prepare_bayesian_fixture(f, family = "poisson")
  expect_error(plot_marginal_effects(fit), "no posterior draws")
  fit <- fake_bayesian_draws(fit, n = 1002)
  plots <- plot_marginal_effects(fit, length.out = 3)
  expect_length(attr(plots, "draw_ids"), 1000)
  expect_match(plots$habitat$labels$caption, "1000 stored coefficient draws")
  selected <- plot_marginal_effects(fit, length.out = 3, draws = 5)
  expect_length(attr(selected, "draw_ids"), 5)
  expect_error(plot_marginal_effects(fit, draws = 0), "draws")
  expect_error(plot_marginal_effects(fit, draws = 1003), "draws")
  fit$draws$beta[1, 1] <- Inf
  expect_error(plot_marginal_effects(fit, draws = 1002), "Stored finite response coefficients")
})

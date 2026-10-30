test_that("likelihood selection checks observed responses, weights, and measures", {
  f <- refit_spatial_fixture()
  a <- f$object$opt_mod
  same <- glm(y ~ mean + rough, data = f$data, family = poisson())
  expect_s3_class(aic_tab(list(a, same), AICc = FALSE), "aictab")
  reversed <- f$data[nrow(f$data):1, ]
  b <- glm(y ~ mean, data = reversed, family = poisson())
  expect_s3_class(bic_tab(list(a, b)), "bictab")
  changed <- f$data
  changed$y[3] <- changed$y[3] + 1
  other <- glm(y ~ mean, data = changed, family = poisson())
  expect_error(aic_tab(list(a, other), AICc = FALSE), "response values")
  expect_error(bic_tab(list(a, other)), "response values")
  weighted <- glm(y ~ mean, data = f$data, weights = rep(2, nrow(f$data)),
                   family = poisson())
  expect_error(aic_tab(list(a, weighted), AICc = FALSE), "likelihood weights")
  continuous <- lm(y ~ mean, data = f$data)
  expect_error(aic_tab(list(a, continuous), AICc = FALSE), "likelihood measures")
  fake <- structure(list(), class = "multiScaleR_bayes")
  expect_error(aic_tab(list(fake)), "Bayesian posterior fits")
  expect_error(bic_tab(list(fake)), "Bayesian posterior fits")
})

test_that("equivalent binomial count and proportion inputs share a selection contract", {
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$success <- seq_len(nrow(dat)) %% 4
  dat$trials <- 5 + seq_len(nrow(dat)) %% 3
  a <- glm(cbind(success, trials - success) ~ mean + rough, dat, family = binomial())
  b <- glm(I(success / trials) ~ mean, dat, weights = trials, family = binomial())
  c <- glm(cbind(success, trials - success) ~ mean, dat, family = binomial())
  expect_equal(as.numeric(logLik(b)), as.numeric(logLik(c)), tolerance = 1e-10)
  expect_s3_class(aic_tab(list(a, b), AICc = FALSE), "aictab")
  dat$trials <- dat$trials + 1
  other <- glm(cbind(success, trials - success) ~ mean, dat, family = binomial())
  expect_error(aic_tab(list(a, other), AICc = FALSE), "trial counts")
})

test_that("zero-weight and missing observations preserve the actual retained cohort", {
  f <- refit_spatial_fixture()
  n <- nrow(f$data)
  dat <- f$data
  weights <- rep(1, n)
  weights[4] <- 0
  a <- glm(y ~ mean, dat, family = poisson(), weights = weights)
  b <- glm(y ~ mean + rough, dat[-4, ], family = poisson())
  expect_equal(.msr_model_nobs(a), n - 1L)
  expect_s3_class(aic_tab(list(a, b), AICc = FALSE), "aictab")
  dat$mean[6] <- NA_real_
  dat$rough[7] <- NA_real_
  a <- glm(y ~ mean, dat, family = poisson())
  b <- glm(y ~ rough, dat, family = poisson())
  expect_error(aic_tab(list(a, b), AICc = FALSE), "different observation sets")
})

test_that("restricted likelihoods require matching objectives and fixed designs", {
  skip_if_not_installed("nlme")
  f <- refit_spatial_fixture()
  a <- nlme::gls(y ~ mean + rough, f$data, method = "REML")
  b <- nlme::gls(y ~ mean, f$data, method = "REML")
  ml <- nlme::gls(y ~ mean + rough, f$data, method = "ML")
  expect_error(aic_tab(list(a, b), AICc = FALSE), "same recoverable fixed-effect design")
  expect_error(bic_tab(list(a, ml)), "ML and REML")
  expect_warning(result <- aic_tab(list(a, a), AICc = FALSE), "redundant")
  expect_s3_class(result, "aictab")
})

test_that("partial likelihoods retain their method and conditioning sets", {
  skip_if_not_installed("survival")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$time <- seq_len(nrow(dat)) %% 4 + 1
  dat$status <- rep(c(1, 0, 1), length.out = nrow(dat))
  a <- survival::coxph(survival::Surv(time, status) ~ mean + strata(group),
    dat, model = TRUE, ties = "efron")
  b <- survival::coxph(survival::Surv(time, status) ~ rough + strata(group),
    dat, model = TRUE, ties = "breslow")
  c <- survival::coxph(survival::Surv(time, status) ~ rough + strata(site),
    dat, model = TRUE, ties = "efron")
  expect_error(aic_tab(list(a, b), AICc = FALSE), "Partial-likelihood methods")
  expect_error(aic_tab(list(a, c), AICc = FALSE), "different strata")
  expect_warning(result <- aic_tab(list(a, a), AICc = FALSE), "redundant")
  expect_s3_class(result, "aictab")
})


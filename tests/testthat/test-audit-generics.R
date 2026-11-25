# Standard model generics on fitted objects (integration audit 2026-11-24,
# findings X-2 and B-2). Every fit here comes from the public API.

audit_glm_fit <- function() {
  data(pts, package = "multiScaleR", envir = environment())
  data(count_data, package = "multiScaleR", envir = environment())
  hab <- terra::rast(system.file("extdata", "hab.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                        kernel = "gaussian", verbose = FALSE)
  start <- glm(y ~ hab, family = poisson, data = count_data)
  fit <- multiScale_optim(start, inputs, verbose = FALSE)
  list(fit = fit, start = start, inputs = inputs)
}

test_that("generics on a Poisson glm fit read the final model, never the start", {
  a <- audit_glm_fit()
  fit <- a$fit
  final <- fit$opt_mod
  # The optimized model differs from the starting model.
  expect_false(isTRUE(all.equal(coef(final), coef(a$start))))

  expect_equal(coef(fit), coef(final))
  expect_equal(fitted(fit), fitted(final))
  expect_type(fitted(fit), "double")
  expect_false(isTRUE(all.equal(unname(fitted(fit)), unname(fitted(a$start)))))
  expect_equal(residuals(fit), residuals(final))
  expect_equal(resid(fit, type = "pearson"), residuals(final, type = "pearson"))
  expect_equal(deviance(fit), deviance(final))
  expect_equal(df.residual(fit), df.residual(final))
  expect_identical(nobs(fit), nobs(final))
  expect_equal(formula(fit), formula(final))
  expect_equal(family(fit)$family, "poisson")
})

test_that("logLik counts the scale, so AIC and BIC match aic_tab and bic_tab", {
  fit <- audit_glm_fit()$fit
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), as.numeric(logLik(fit$opt_mod)))
  # Two glm coefficients plus one estimated sigma.
  expect_identical(as.numeric(attr(ll, "df")), attr(logLik(fit$opt_mod), "df") + 1)
  expect_identical(attr(ll, "nobs"), nobs(fit$opt_mod))
  expect_equal(AIC(fit), aic_tab(list(fit), AICc = FALSE)$AIC, tolerance = 1e-10)
  expect_equal(BIC(fit), bic_tab(list(fit))$BIC, tolerance = 1e-10)
  expect_equal(AIC(fit), AIC(fit$opt_mod) + 2, tolerance = 1e-10)
})

test_that("an unmarked fit delegates to unmarked and matches aic_tab", {
  skip_if_not_installed("unmarked")
  nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
                              package = "multiScaleR"))
  expect_equal(coef(nmix), unmarked::coef(nmix$opt_mod))
  expect_equal(AIC(nmix), aic_tab(list(nmix), AICc = FALSE)$AIC, tolerance = 1e-10)
  expect_equal(BIC(nmix), bic_tab(list(nmix))$BIC, tolerance = 1e-10)
  # Native unmarked AIC leaves the two estimated scales out.
  expect_equal(AIC(nmix), nmix$opt_mod@AIC + 2 * nrow(nmix$scale_est), tolerance = 1e-8)
  # unmarked has no family(); the native error passes through.
  expect_error(family(nmix))
})

test_that("Bayesian fits refuse point-estimate generics and name the alternatives", {
  f <- bayesian_scale_fixture()
  fit <- prepare_bayesian_fixture(f)
  for (g in list(coef, fitted, residuals, vcov))
    expect_error(g(fit), "posterior_scales\\(\\)")
  expect_error(logLik(fit), "not a maximized likelihood")
  expect_error(AIC(fit), "not a maximized likelihood")
  expect_identical(nobs(fit), nrow(f$data))
  expect_error(diagnostics(fit), "diagnose_bayes_scales\\(\\)")
})

test_that("diagnostics() on an unsupported object gives a clear message", {
  expect_error(diagnostics(list(a = 1)), "needs a multiScaleR fit")
  expect_error(diagnostics(lm(dist ~ speed, cars)), "class `lm`")
})

# Raster maps from conditional_logit fits. kernel_scale.raster() must add the
# strata(choice_id) variable as a constant layer holding a fitted ID, so that
# terra::predict() on the native coxph model works without `const` and without
# a spurious categorical-covariate warning. ssf_fixture() is in
# helper-step-selection.R.

test_that("kernel_scale.raster() supplies a fitted choice-set ID for clogit maps", {
  skip_on_cran()
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  m <- survival::clogit(case_ ~ forest + sl_km + strata(step_id),
                        data = fx$df, model = TRUE)
  fit <- multiScale_optim(m, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE,
                          profile_control = list(starts = matrix(0.3, ncol = 1)))

  # Default: the ID layer holds the first fitted ID; movement terms are 0.
  expect_no_warning(r <- kernel_scale.raster(fx$forest, multiScaleR = fit,
                                             scale_center = TRUE, verbose = FALSE))
  expect_setequal(names(r), c("forest", "sl_km", "step_id"))
  expect_equal(unname(unlist(terra::global(r$step_id, "max"))), 1)
  expect_equal(unname(unlist(terra::global(r$sl_km, "max"))), 0)

  # The native model predicts directly and gives exp(b * forest).
  p <- terra::predict(r, fit$opt_mod, type = "risk", reference = "zero",
                      na.rm = TRUE)
  b <- stats::coef(fit$opt_mod)[["forest"]]
  expect_equal(terra::values(p)[, 1], exp(b * terra::values(r$forest)[, 1]),
               tolerance = 1e-10)

  # Any fitted ID may be chosen; an unfitted ID is rejected.
  r7 <- kernel_scale.raster(fx$forest, multiScaleR = fit, scale_center = TRUE,
                            site_values = list(step_id = 7), verbose = FALSE)
  expect_equal(unname(unlist(terra::global(r7$step_id, "max"))), 7)
  expect_error(kernel_scale.raster(fx$forest, multiScaleR = fit,
                                   scale_center = TRUE, verbose = FALSE,
                                   site_values = list(step_id = 9999)),
               "fitted choice-set ID")

  # Character IDs become a one-level categorical layer with a fitted label.
  df <- fx$df
  df$step_id <- paste0("s", df$step_id)
  mc <- survival::clogit(case_ ~ forest + sl_km + strata(step_id),
                         data = df, model = TRUE)
  fc <- multiScale_optim(mc, fx$inputs, backend = "conditional_logit",
                         n_cores = 1, verbose = FALSE,
                         profile_control = list(starts = matrix(0.3, ncol = 1)))
  rc <- kernel_scale.raster(fx$forest, multiScaleR = fc, scale_center = TRUE,
                            verbose = FALSE)
  expect_equal(terra::levels(rc$step_id)[[1]]$label, "s1")
  pc <- terra::predict(rc, fc$opt_mod, type = "lp", reference = "zero",
                       na.rm = TRUE)
  expect_equal(terra::values(pc)[, 1],
               stats::coef(fc$opt_mod)[["forest"]] * terra::values(rc$forest)[, 1],
               tolerance = 1e-10)
})

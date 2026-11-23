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

  # Character IDs become a categorical layer carrying every fitted label,
  # set to the first one.
  df <- fx$df
  df$step_id <- paste0("s", df$step_id)
  mc <- survival::clogit(case_ ~ forest + sl_km + strata(step_id),
                         data = df, model = TRUE)
  fc <- multiScale_optim(mc, fx$inputs, backend = "conditional_logit",
                         n_cores = 1, verbose = FALSE,
                         profile_control = list(starts = matrix(0.3, ncol = 1)))
  rc <- kernel_scale.raster(fx$forest, multiScaleR = fc, scale_center = TRUE,
                            verbose = FALSE)
  expect_identical(terra::levels(rc$step_id)[[1]]$label,
                   levels(factor(df$step_id)))
  expect_identical(unique(as.character(terra::as.data.frame(rc$step_id)$step_id)),
                   levels(factor(df$step_id))[1L])
  pc <- terra::predict(rc, fc$opt_mod, type = "lp", reference = "zero",
                       na.rm = TRUE)
  expect_equal(terra::values(pc)[, 1],
               stats::coef(fc$opt_mod)[["forest"]] * terra::values(rc$forest)[, 1],
               tolerance = 1e-10)
})

test_that("factor choice-set IDs keep every fitted level in clogit maps", {
  skip_on_cran()
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  # Numeric-looking IDs stored as a factor. A numeric ID layer would not match
  # the fitted factor coding, and predict.coxph() would report a new level.
  df <- fx$df
  df$step_id <- factor(df$step_id)
  m <- survival::clogit(case_ ~ forest + sl_km + strata(step_id),
                        data = df, model = TRUE)
  fit <- multiScale_optim(m, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE,
                          profile_control = list(starts = matrix(0.3, ncol = 1)))
  b <- stats::coef(fit$opt_mod)[["forest"]]

  # Checks one map: a categorical ID layer carrying all fitted levels, set to
  # `id` in every cell, that terra::predict() and predict.coxph() both accept.
  check_map <- function(r, model, id) {
    expect_true(terra::is.factor(r$step_id))
    expect_identical(terra::levels(r$step_id)[[1]]$label, levels(df$step_id))
    expect_identical(unique(as.character(terra::as.data.frame(r$step_id)$step_id)), id)
    p <- terra::predict(r, model, type = "risk", reference = "zero", na.rm = TRUE)
    expect_equal(terra::values(p)[, 1], exp(b * terra::values(r$forest)[, 1]),
                 tolerance = 1e-10)
    nd <- terra::as.data.frame(r, na.rm = TRUE)
    expect_equal(unname(stats::predict(model, nd, type = "risk", reference = "zero")),
                 exp(b * nd$forest), tolerance = 1e-10)
  }

  # Default: the first fitted level, with no categorical-covariate warning.
  expect_no_warning(r <- kernel_scale.raster(fx$forest, multiScaleR = fit,
                                             scale_center = TRUE, verbose = FALSE))
  check_map(r, fit$opt_mod, "1")

  # site_values selects another fitted level; an unfitted level is rejected.
  r7 <- kernel_scale.raster(fx$forest, multiScaleR = fit, scale_center = TRUE,
                            site_values = list(step_id = "7"), verbose = FALSE)
  check_map(r7, fit$opt_mod, "7")
  expect_error(kernel_scale.raster(fx$forest, multiScaleR = fit,
                                   scale_center = TRUE, verbose = FALSE,
                                   site_values = list(step_id = "9999")),
               "fitted choice-set ID")

  # The factor coding survives when the original training data can no longer
  # be found, so the ID layer is built from the retained model frame instead.
  expired <- fit
  mod <- fit$opt_mod
  mod$call$data <- quote(expired_clogit_training_data)
  attr(mod$terms, ".Environment") <- asNamespace("survival")
  attr(attr(mod$model, "terms"), ".Environment") <- asNamespace("survival")
  expired$opt_mod <- mod
  expect_no_warning(re <- kernel_scale.raster(fx$forest, multiScaleR = expired,
                                              scale_center = TRUE, verbose = FALSE))
  check_map(re, mod, "1")
})

test_that("numeric choice-set IDs stay numeric beside transformed fixed terms", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("amt")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  # A fixed, non-optimized predictor transformed in the formula, fitted with
  # amt::fit_clogit(). The fit keeps a retained model frame, whose recovered
  # data hold step_id as a factor, while survival's fitted strata labels
  # ("step_id=3") are numeric. IDs start at 3 so the first fitted ID is not 1.
  df <- fx$df
  df$step_id <- df$step_id + 2
  df$river_m <- abs(fx$steps$x - 600)
  m <- amt::fit_clogit(df, case_ ~ forest + scale(river_m) + sl_km + strata(step_id),
                       model = TRUE, na.action = stats::na.fail)
  fit <- multiScale_optim(m, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE,
                          profile_control = list(starts = matrix(0.3, ncol = 1)))
  expect_identical(fit$profile_fit$contract$native_data_mode, "retained_model_frame")
  expect_identical(.msr_strata_ids(fit$opt_mod$model)$ids[1:2], c(3, 4))

  # Effect curves hold the ID at a fitted value (previously "new level 3").
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot_marginal_effects(fit, link = TRUE)
  expect_setequal(names(curves), c("forest", "river_m", "sl_km"))
  expect_equal(range(curves$river_m$data[[1]]), range(df$river_m))

  # The map uses a numeric ID layer and the native model predicts from it.
  r <- kernel_scale.raster(fx$forest, multiScaleR = fit, scale_center = TRUE,
                           site_values = list(river_m = 100), verbose = FALSE)
  expect_false(terra::is.factor(r$step_id))
  expect_equal(unname(unlist(terra::global(r$step_id, "max"))), 3)
  native <- fit$opt_mod$model
  b <- stats::coef(native)
  p <- terra::predict(r, native, type = "lp", reference = "zero", na.rm = TRUE)
  z100 <- (100 - mean(df$river_m)) / stats::sd(df$river_m)
  expect_equal(terra::values(p)[, 1],
               b[["forest"]] * terra::values(r$forest)[, 1] +
                 b[["scale(river_m)"]] * z100, tolerance = 1e-10)
})

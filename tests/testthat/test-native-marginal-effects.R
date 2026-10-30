test_that("native GLM marginal effects preserve factors, contrasts, and spatial units", {
  f <- refit_spatial_fixture()
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(f$object, length.out = 7))
  expect_identical(as.character(plots$group$data$x), c("a", "b", "c"))
  vars <- c("mean", "rough", "fixed", "group", "z", "exposure")
  defaults <- lapply(f$data[vars], function(v)
    if (is.numeric(v)) mean(v) else v[1])
  nd <- as.data.frame(lapply(defaults, rep, times = 7))
  nd$rough <- seq(min(f$data$rough), max(f$data$rough), length.out = 7)
  p <- predict(f$object$opt_mod, nd, type = "link", se.fit = TRUE)
  expect_equal(plots$rough$data$fit, unname(exp(p$fit)), tolerance = 1e-12)
  expect_equal(plots$rough$data$lwr, unname(exp(p$fit + qnorm(.025) * p$se.fit)),
               tolerance = 1e-12)
  expect_equal(plots$rough$data$x, nd$rough * f$object$scl_params$sd[["rough"]] +
    f$object$scl_params$mean[["rough"]], tolerance = 1e-12)
  expect_match(plots$rough$labels$caption, "condition on.*spatial scales")
  nd <- as.data.frame(lapply(defaults, rep, times = 3))
  nd$group <- factor(c("a", "b", "c"), levels = levels(f$data$group))
  expect_equal(plots$group$data$fit,
    unname(predict(f$object$opt_mod, nd, type = "response")), tolerance = 1e-12)
})

test_that("glmmTMB curves use its full response target and zero random effects", {
  skip_if_not_installed("glmmTMB")
  f <- refit_spatial_fixture()
  model <- glmmTMB::glmmTMB(
    y ~ mean + rough + fixed + group + z + (1 | site), ziformula = ~1,
    data = f$data, family = poisson(),
    start = list(betazi = qlogis(.3), theta = log(.4)),
    map = list(betazi = factor(NA), theta = factor(NA)))
  f$object$opt_mod <- model
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- expect_no_warning(plot_marginal_effects(f$object, length.out = 7))
  vars <- c("mean", "rough", "fixed", "group", "z", "site")
  defaults <- lapply(f$data[vars], function(v)
    if (is.numeric(v)) mean(v) else v[1])
  nd <- as.data.frame(lapply(defaults, rep, times = 7))
  nd$rough <- seq(min(f$data$rough), max(f$data$rough), length.out = 7)
  expected <- predict(model, nd, type = "response", re.form = NA, se.fit = TRUE)
  expect_equal(plots$rough$data$fit, unname(expected$fit), tolerance = 1e-10)
  expect_equal(plots$rough$data$lwr,
    unname(expected$fit + qnorm(.025) * expected$se.fit), tolerance = 1e-10)
  expect_match(plots$rough$labels$caption, "zero, not integrated")
  expect_false("site" %in% names(plots))
  expect_identical(as.character(plots$group$data$x), c("a", "b", "c"))
  rows <- c(5, 23, 17)
  expected <- predict(model, f$data[rows, ], type = "response",
                       re.form = NA, se.fit = TRUE)
  actual <- predict(f$object, f$data[rows, ], f$inputs, re.form = NA,
    se.fit = TRUE, uncertainty = "conditional", include_fixed = TRUE)
  expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-10)
  expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-10)
  expect_identical(attr(actual, "prediction_contract")$integration$random_effects,
                   "set_to_zero_not_integrated")
})

test_that("native fixed-radius normalization also supports linear model rasters", {
  f <- refit_spatial_fixture()
  f$object$opt_mod <- lm(y ~ mean + rough + fixed + group, data = f$data)
  expect_warning(raster <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE),
    "group.*categorical")
  actual <- terra::extract(raster, f$points, ID = FALSE)
  expect_true(all(c("mean", "rough", "fixed") %in% names(raster)))
  expect_equal(actual$fixed, f$data$fixed, tolerance = 1e-9)
  expect_equal(actual$rough, f$data$rough, tolerance = 1e-8)
})

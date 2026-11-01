test_that("native lm abbreviations retain conditional interval guards and levels", {
  f <- refit_spatial_fixture()
  model <- f$object$opt_mod <- lm(y ~ mean + rough + fixed + z, data = f$data)
  nd <- f$data[c(5, 23, 17), ]
  controls <- list(int = "confidence", lev = .9)
  expect_error(predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = controls), "require uncertainty = 'conditional'")
  actual <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = controls, uncertainty = "conditional")
  native <- predict(model, nd, int = "confidence", lev = .9)
  expect_equal(as.numeric(actual), as.numeric(native), tolerance = 1e-12)
  expect_equal(attr(actual, "prediction_contract")$uncertainty$interval_level, .9)
  snapshot <- prepare_prediction(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = controls)
  expect_identical(names(snapshot$provenance$native_args), c("interval", "level"))
  expect_equal(predict(snapshot, uncertainty = "conditional"), actual)
  expect_error(predict(snapshot), "conditional")
  expect_error(predict(f$object, native_args = list(newd = f$data)), "cannot replace wrapper")
  # Longer unknown names are not abbreviations of shorter native formals.
  plain <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = list(intervals = "confidence"))
  expect_equal(as.numeric(plain), as.numeric(predict(model, nd)), tolerance = 1e-12)
  expect_identical(attr(plain, "prediction_contract")$uncertainty$intervals, "not_provided")
  limits <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
    uncertainty = "conditional", native_args = list(interval = "confidence", level.extra = .8))
  expect_equal(attr(limits, "prediction_contract")$uncertainty$interval_level, .95)
  expect_equal(as.numeric(limits), as.numeric(predict(model, nd, interval = "confidence")), tolerance = 1e-12)
})

test_that("native covariance abbreviations cannot escape point or raster policies", {
  skip_if_not_installed("glmmTMB")
  f <- refit_spatial_fixture()
  model <- f$object$opt_mod <- glmmTMB::glmmTMB(y ~ mean + rough + fixed + z,
    data = f$data, family = gaussian())
  nd <- f$data[c(5, 23, 17), ]
  expect_error(predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = list(cov = TRUE)), "require uncertainty = 'conditional'")
  actual <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
    uncertainty = "conditional", native_args = list(cov = TRUE))
  native <- predict(model, nd, type = "response", cov = TRUE)
  for (nm in c("fit", "se.fit", "cov.fit")) expect_equal(actual[[nm]], native[[nm]], tolerance = 1e-10)
  expect_identical(attr(actual, "prediction_contract")$uncertainty$standard_errors, "provided")
  map <- kernel_scale.raster(f$rasters, multiScaleR = f$object, include_fixed = TRUE,
    scale_center = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  expect_error(predict(f$object, map, native_args = list(cov = TRUE)), "Cross-row native covariance")
  mixed <- f$object; mixed$backend <- "joint_mixed"
  expect_error(predict(mixed, native_args = list(cov = TRUE)), "supports conditional prediction SEs only")
  with_mocked_bindings({
    matched <- predict(mixed, uncertainty = "conditional", native_args = list(cov = TRUE))
    for (nm in c("fit", "se.fit", "cov.fit")) expect_equal(matched[[nm]], native[[nm]], tolerance = 1e-10)
  }, .msr_fixed_point_data = function(object, newdata, kernel_inputs, include_fixed) newdata,
     .msr_mixed_prediction_data = function(object, newdata, kernel_inputs) nd)
})

test_that("HLfit abbreviated variance controls retain native attributes and conditional gates", {
  skip_if_not_installed("spaMM")
  f <- native_gamma_fixture("HLfit")
  nd <- f$data[c(5, 23, 17), ]
  controls <- list(variance = list(respVar = TRUE))
  expect_error(predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = controls), "require uncertainty = 'conditional'")
  actual <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = controls, uncertainty = "conditional")
  native <- predict(f$object$opt_mod, nd, type = "response", variance = list(respVar = TRUE))
  for (nm in c("fixefVar", "predVar", "residVar", "respVar"))
    expect_equal(attr(actual, nm), attr(native, nm), tolerance = 1e-10)
  snapshot <- prepare_prediction(f$object, nd, f$inputs, component = "response",
    include_fixed = TRUE, native_args = controls)
  expect_identical(names(snapshot$provenance$native_args), "variances")
  expect_error(predict(snapshot), "conditional")
  expect_equal(predict(snapshot, uncertainty = "conditional"), actual, tolerance = 1e-10)
  map <- kernel_scale.raster(f$rasters, multiScaleR = f$object, include_fixed = TRUE,
    scale_center = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  expect_error(predict(f$object, map, native_args = controls), "require point predictions")
})

test_that("S4 lme4 prediction abbreviations use their inherited S3 method", {
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  f$data$y <- f$data$y + as.integer(f$data$site)^2
  model <- f$object$opt_mod <- lme4::lmer(y ~ mean + rough + fixed + z + (1 | site),
    data = f$data, REML = FALSE)
  nd <- f$data[c(5, 23, 17), ]; nd$site <- NULL
  actual <- predict(f$object, nd, f$inputs, include_fixed = TRUE, re.f = NA)
  native <- predict(model, nd, re.f = NA)
  expect_equal(as.numeric(actual), as.numeric(native), tolerance = 1e-10)
  expect_identical(attr(actual, "prediction_contract")$integration$random_effects,
    "set_to_zero_not_integrated")
})

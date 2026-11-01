test_that("native returned point and SE contracts describe actual uncertainty", {
  f <- refit_spatial_fixture()
  point <- predict(f$object, f$data, f$inputs, include_fixed = TRUE)
  info <- attr(point, "prediction_contract")$uncertainty
  expect_identical(info$standard_errors, "not_requested")
  expect_identical(info$intervals, "not_provided")
  expect_identical(info$population_parameters, "fixed_at_estimate")
  expect_true("population_parameter_uncertainty" %in% info$omitted)
  se <- predict(f$object, f$data, f$inputs, include_fixed = TRUE,
    se.fit = TRUE, uncertainty = "conditional")
  info <- attr(se, "prediction_contract")$uncertainty
  expected <- predict(f$object$opt_mod, f$data, type = "response", se.fit = TRUE)
  expect_equal(as.numeric(point), as.numeric(expected$fit), tolerance = 1e-12)
  expect_equal(as.numeric(se$se.fit), as.numeric(expected$se.fit), tolerance = 1e-12)
  expect_identical(info$standard_errors, "provided")
  expect_identical(info$intervals, "conditional_standard_errors_only")
  expect_identical(info$kernel_parameters, "fixed_at_estimate")
  expect_false(info$observation_noise)
  expect_false("population_parameter_uncertainty" %in% info$omitted)
  snapshot <- prepare_prediction(f$object, f$data, f$inputs, include_fixed = TRUE)
  expect_identical(snapshot$prediction_contract$response$uncertainty$intervals, "not_provided")
  expect_true(snapshot$prediction_contract$response$uncertainty$conditional_se_available)
  expect_equal(predict(snapshot), point)
  expect_equal(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"), se)
})

test_that("basis snapshots retain explicit unavailable SEs without claiming uncertainty", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture()
  f$object$opt_mod <- mgcv::gam(y ~ mean + rough + fixed + s(z, k = 4),
    data = f$data, method = "ML")
  snapshot <- prepare_prediction(f$object, f$data, f$inputs, component = "lpmatrix",
    include_fixed = TRUE)
  point <- predict(snapshot)
  expect_equal(as.numeric(point), as.numeric(predict(f$object$opt_mod,
    f$data, type = "lpmatrix")), tolerance = 1e-12)
  expect_identical(attr(point, "prediction_contract")$uncertainty$intervals, "not_provided")
  expect_identical(attr(point, "prediction_contract")$uncertainty$standard_errors, "not_requested")
  expect_false(snapshot$prediction_contract$lpmatrix$uncertainty$conditional_se_available)
  expect_error(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"), "did not supply se.fit")
})


test_that("native covariance controls report actual SEs and cannot lose cross-row covariance in maps", {
  skip_if_not_installed("glmmTMB")
  f <- refit_spatial_fixture()
  model <- f$object$opt_mod <- glmmTMB::glmmTMB(y ~ mean + rough + fixed,
    data = f$data, family = gaussian())
  nd <- f$data[c(5, 23, 17), ]
  expected <- predict(model, nd, type = "response", se.fit = FALSE, cov.fit = TRUE)
  expect_error(predict(f$object, nd, f$inputs, include_fixed = TRUE,
    native_args = list(cov.fit = TRUE)), "require uncertainty = 'conditional'")
  point <- predict(f$object, nd, f$inputs, include_fixed = TRUE, uncertainty = "conditional",
    native_args = list(cov.fit = TRUE))
  for (nm in c("fit", "se.fit", "cov.fit")) expect_equal(point[[nm]], expected[[nm]], tolerance = 1e-10)
  mixed <- f$object; mixed$backend <- "joint_mixed"
  with_mocked_bindings({
    actual <- predict(mixed, native_args = list(cov.fit = TRUE), uncertainty = "conditional")
    expect_equal(actual, expected, tolerance = 1e-10)
  }, .msr_fixed_point_data = function(object, newdata, kernel_inputs, include_fixed) newdata,
     .msr_mixed_prediction_data = function(object, newdata, kernel_inputs) nd)
  info <- attr(point, "prediction_contract")$uncertainty
  expect_false(info$se.fit_requested)
  expect_identical(info$standard_errors, "provided")
  expect_identical(info$covariance, "native_conditional_covariance")
  expect_identical(info$intervals, "conditional_standard_errors_only")
  expect_false("population_parameter_uncertainty" %in% info$omitted)
  map <- kernel_scale.raster(f$rasters, multiScaleR = f$object, include_fixed = TRUE,
    scale_center = TRUE, match_buffer = TRUE, verbose = FALSE)
  expect_error(predict(f$object, map, native_args = list(cov.fit = TRUE)), "Cross-row native covariance")
})


test_that("joint mixed covariance controls retain the conditional-uncertainty gate", {
  fake <- structure(list(backend = "joint_mixed"), class = "multiScaleR")
  expect_error(predict(fake, cov.fit = TRUE), "supports conditional prediction SEs only")
  expect_error(predict(fake, native_args = list(cov.fit = TRUE)), "supports conditional prediction SEs only")
  expect_error(predict(fake, cov.fit = TRUE, native_args = list(cov.fit = TRUE)), "uniquely name")
})


test_that("native residual variance is noise without manufactured coefficient uncertainty", {
  contract <- .msr_native_prediction_contract(NULL, "HLfit", "response", list())
  info <- .msr_returned_prediction_uncertainty(contract, FALSE,
    native_variances = "residVar")$uncertainty
  expect_true(info$observation_noise)
  expect_identical(info$population_parameters, "fixed_at_estimate")
  expect_true("population_parameter_uncertainty" %in% info$omitted)
})

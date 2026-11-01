# Provider parity checks hold selected scales and reference normalization fixed.
# Mean intervals and future-observation intervals retain different targets.
contains_snapshot_environment <- function(x) {
  if (is.environment(x) || is.function(x)) return(TRUE)
  values <- if (is.list(x)) vapply(x, contains_snapshot_environment, TRUE) else FALSE
  any(values) || any(vapply(attributes(x), contains_snapshot_environment, TRUE))
}

native_interval_fixture <- function() {
  f <- refit_spatial_fixture()
  f$object$opt_mod <- lm(y ~ mean + rough + fixed + z, data = f$data)
  f
}

test_that("native lm interval columns preserve values and identify the uncertainty target", {
  f <- native_interval_fixture()
  nd <- f$data[c(5, 23, 17), ]
  for (kind in c("confidence", "prediction")) for (se in c(FALSE, TRUE)) {
    args <- list(interval = kind, level = .9)
    expect_error(predict(f$object, nd, f$inputs, include_fixed = TRUE,
      native_args = args), "require uncertainty = 'conditional'")
    actual <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
      native_args = args, se.fit = se, uncertainty = "conditional")
    expected <- predict(f$object$opt_mod, nd, interval = kind, level = .9, se.fit = se)
    if (se) {
      expect_equal(actual$fit, expected$fit, tolerance = 1e-12)
      expect_equal(actual$se.fit, expected$se.fit, tolerance = 1e-12)
      expect_equal(actual$residual.scale, expected$residual.scale)
    } else expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
    info <- attr(actual, "prediction_contract")$uncertainty
    expect_identical(info$intervals, if (kind == "confidence")
      "native_mean_confidence" else "native_prediction_interval")
    expect_identical(info$native_intervals, kind)
    expect_identical(info$interval_level, .9)
    expect_identical(info$se.fit_requested, se)
    expect_identical(info$standard_errors, if (se) "provided" else "not_requested")
    expect_identical(info$observation_noise, kind == "prediction")
    expect_false("population_parameter_uncertainty" %in% info$omitted)
    expect_true("kernel_parameter_uncertainty" %in% info$omitted)
    expect_identical("observation_noise" %in% info$omitted, kind == "confidence")
  }
  expect_error(predict(f$object, nd, f$inputs, type = "terms", include_fixed = TRUE,
    uncertainty = "conditional", native_args = list(interval = "confidence")),
    "require response predictions")
})

test_that("lm interval snapshots keep native levels, columns and conditional gates", {
  f <- native_interval_fixture()
  nd <- f$data[c(5, 23, 17), ]
  for (kind in c("confidence", "prediction")) {
    snapshot <- prepare_prediction(f$object, nd, f$inputs, include_fixed = TRUE,
      native_args = list(interval = kind, level = .9))
    expect_equal(snapshot$n_rows, 3)
    expect_null(snapshot$opt_mod)
    expect_null(snapshot$kernel_inputs)
    expect_false(contains_snapshot_environment(snapshot))
    expect_error(predict(snapshot), "conditional")
    restored <- unserialize(serialize(snapshot, NULL))
    for (type in c("response", "link")) for (se in c(FALSE, TRUE)) {
      actual <- predict(restored, type = type, se.fit = se, uncertainty = "conditional")
      expected <- predict(f$object, nd, f$inputs, type = type, se.fit = se,
        include_fixed = TRUE, uncertainty = "conditional",
        native_args = list(interval = kind, level = .9))
      expect_equal(actual, expected, tolerance = 1e-12)
    }
    expect_error(prepare_prediction(f$object, nd, f$inputs, include_fixed = TRUE,
      native_args = list(interval = kind), max_bytes = 1), "max_bytes")
  }
})

test_that("native lm interval maps preserve columns and mean versus observation semantics", {
  f <- native_interval_fixture()
  map <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    include_fixed = TRUE, scale_center = TRUE, match_buffer = TRUE,
    site_values = list(z = .2), verbose = FALSE)
  reference <- f$data; reference$z <- .2
  for (kind in c("confidence", "prediction")) for (se in c(FALSE, TRUE)) {
    expect_error(predict(f$object, map, native_args = list(interval = kind)),
      "require uncertainty = 'conditional'")
    actual <- predict(f$object, map, se.fit = se, uncertainty = "conditional",
      native_args = list(interval = kind, level = .9))
    native <- predict(f$object$opt_mod, reference, se.fit = se, interval = kind, level = .9)
    expected <- if (se) cbind(native$fit, se.fit = native$se.fit) else native
    expect_equal(as.numeric(as.matrix(terra::extract(actual, f$points, ID = FALSE))),
      as.numeric(expected), tolerance = 1e-8)
    expect_equal(terra::nlyr(actual), if (se) 4 else 3)
    info <- attr(actual, "prediction_contract")$uncertainty
    expect_identical(info$intervals, if (kind == "confidence")
      "native_mean_confidence" else "native_prediction_interval")
    expect_identical(info$observation_noise, kind == "prediction")
    expect_false("population_parameter_uncertainty" %in% info$omitted)
    expect_identical("observation_noise" %in% info$omitted, kind == "confidence")
    expect_identical(info$interval_level, .9)
  }
})

test_that("HLfit snapshots preserve compact native variance and interval attributes", {
  skip_if_not_installed("spaMM")
  f <- refit_spatial_fixture()
  f$data$y <- f$data$y + as.integer(f$data$site)^2
  f$object$opt_mod <- spaMM::fitme(y ~ mean + rough + fixed + z + (1 | site),
    data = f$data, family = poisson(), method = "ML")
  nd <- f$data[c(5, 23, 17), ]; nd$site <- NULL
  for (controls in list(list(variances = list(respVar = TRUE)), list(intervals = "predVar"))) {
    controls$re.form <- NA
    snapshot <- prepare_prediction(f$object, nd, f$inputs, component = "response",
      include_fixed = TRUE, native_args = controls)
    expect_false(snapshot$prediction_contract$response$uncertainty$conditional_se_available)
    expect_error(predict(snapshot), "conditional")
    expect_error(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"), "did not supply se.fit")
    native <- predict(f$object, nd, f$inputs, include_fixed = TRUE,
      uncertainty = "conditional", native_args = controls)
    actual <- predict(unserialize(serialize(snapshot, NULL)), uncertainty = "conditional")
    expect_equal(actual, native, tolerance = 1e-10)
    expect_null(snapshot$opt_mod)
    expect_null(snapshot$kernel_inputs)
    expect_false(contains_snapshot_environment(snapshot))
    expect_lt(length(serialize(snapshot, NULL)), 1024^2)
    info <- attr(actual, "prediction_contract")$uncertainty
    if (!is.null(controls$variances)) {
      expect_true(info$observation_noise)
      expect_false("observation_noise" %in% info$omitted)
    } else expect_identical(info$intervals, "native_attribute_intervals")
  }
  expect_error(prepare_prediction(f$object, nd, f$inputs, component = "response",
    include_fixed = TRUE, native_args = list(variances = list(cov = TRUE))), "Cross-row native covariance")
  expect_error(prepare_prediction(f$object, nd, f$inputs, component = "response",
    include_fixed = TRUE, native_args = 1), "native_args must uniquely name")
})

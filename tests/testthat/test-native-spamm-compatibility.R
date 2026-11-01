# Real provider parity at declared spatial scales; these fixtures do not assess
# scale recovery, approximation accuracy, or interval calibration.
spamm_native_fixture <- local({
  cache <- list()
  function(family = "gaussian") {
    if (!is.null(cache[[family]])) return(cache[[family]])
    f <- refit_spatial_fixture()
    f$data$y <- f$data$y + as.integer(f$data$site)^2
    f$object$opt_mod <- spaMM::fitme(y ~ mean + rough + fixed + z + (1 | site),
      data = f$data, family = if (family == "gaussian") gaussian() else poisson(), method = "ML")
    cache[[family]] <<- f
    f
  }
})

test_that("real HLfit fixed features, zero effects and snapshots retain native values", {
  skip_if_not_installed("spaMM")
  for (family in c("gaussian", "poisson")) {
    f <- spamm_native_fixture(family)
    model <- f$object$opt_mod
    nd <- f$data[c(5, 23, 17), ]; nd$site <- NULL; nd$fixed <- NULL
    reference <- f$data[c(5, 23, 17), ]; reference$site <- NULL
    for (type in c("response", "link")) {
      expected <- predict(model, reference, type = type, re.form = NA)
      actual <- predict(f$object, nd, f$inputs, type = type, re.form = NA, include_fixed = TRUE)
      expect_s3_class(actual, "spaMM_predictions")
      expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
      expect_identical(attr(actual, "prediction_contract")$integration$random_effects,
        "set_to_zero_not_integrated")
    }
    expect_error(predict(f$object, nd, f$inputs, include_fixed = TRUE), "site")
    snapshot <- prepare_prediction(f$object, nd, f$inputs, component = "response",
      include_fixed = TRUE, native_args = list(re.form = NA))
    expect_equal(predict(snapshot), predict(f$object, nd, f$inputs,
      re.form = NA, include_fixed = TRUE))
    expect_false(snapshot$prediction_contract$response$uncertainty$conditional_se_available)
    expect_error(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"), "did not supply se.fit")
    expect_null(snapshot$opt_mod)
    expect_equal(predict(unserialize(serialize(snapshot, NULL))), predict(snapshot))
    expect_warning(map <- kernel_scale.raster(f$rasters, multiScaleR = f$object, include_fixed = TRUE,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE),
      "site.*categorical")
    actual <- predict(f$object, map, re.form = NA)
    reference <- f$data; reference$z <- .2; reference$site <- NULL
    expect_equal(terra::extract(actual, f$points, ID = FALSE)[[1]],
      as.numeric(predict(model, reference, type = "response", re.form = NA)), tolerance = 1e-8)
    expect_error(predict(f$object, map, re.form = NA,
      native_args = list(variances = list(predVar = TRUE))), "require point predictions")
  }
})

test_that("HLfit plots transform native link-scale uncertainty and preserve variance attributes", {
  skip_if_not_installed("spaMM")
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (family in c("gaussian", "poisson")) {
    f <- spamm_native_fixture(family)
    model <- f$object$opt_mod
    plots <- plot_marginal_effects(f$object, length.out = 7)
    nd <- as.data.frame(lapply(f$data[c("mean", "rough", "fixed", "z")],
      function(x) rep(mean(x), 7)))
    nd$rough <- seq(min(f$data$rough), max(f$data$rough), length.out = 7)
    native <- predict(model, nd, type = "link", re.form = NA,
      variances = list(fixefVar = TRUE))
    inverse <- stats::family(model)$linkinv
    expect_equal(plots$rough$data$fit, as.numeric(inverse(native)), tolerance = 1e-10)
    expect_equal(plots$rough$data$lwr,
      as.numeric(inverse(native + qnorm(.025) * sqrt(attr(native, "fixefVar")))), tolerance = 1e-10)
    expect_equal(plots$rough$data$upr,
      as.numeric(inverse(native + qnorm(.975) * sqrt(attr(native, "fixefVar")))), tolerance = 1e-10)
    expect_match(plots$rough$labels$caption, "zero, not integrated")
    point <- predict(f$object, f$data, f$inputs, re.form = NA, include_fixed = TRUE,
      uncertainty = "conditional", native_args = list(variances = list(respVar = TRUE)))
    expect_error(predict(f$object, f$data, f$inputs, re.form = NA, include_fixed = TRUE,
      native_args = list(variances = list(respVar = TRUE))), "require uncertainty = 'conditional'")
    native <- predict(model, f$data, type = "response", re.form = NA,
      variances = list(respVar = TRUE))
    for (nm in c("fixefVar", "predVar", "residVar", "respVar"))
      expect_equal(attr(point, nm), attr(native, nm), tolerance = 1e-10)
    info <- attr(point, "prediction_contract")$uncertainty
    expect_true(info$observation_noise)
    expect_true("predVar" %in% info$native_variances$linear_predictor_scale)
    expect_true("respVar" %in% info$native_variances$response_scale)
    expect_identical(info$standard_errors, "not_requested")
    limits <- predict(f$object, f$data, f$inputs, re.form = NA, include_fixed = TRUE,
      uncertainty = "conditional", native_args = list(intervals = "predVar"))
    native_limits <- predict(model, f$data, re.form = NA, type = "response", intervals = "predVar")
    expect_equal(attr(limits, "intervals"), attr(native_limits, "intervals"), tolerance = 1e-10)
    expect_identical(attr(limits, "prediction_contract")$uncertainty$intervals, "native_attribute_intervals")
    expect_true("predVar" %in% attr(limits, "prediction_contract")$uncertainty$native_intervals)
    expect_false(attr(limits, "prediction_contract")$uncertainty$observation_noise)
    link <- predict(f$object, f$data, f$inputs, type = "link", re.form = NA,
      include_fixed = TRUE, uncertainty = "conditional",
      native_args = list(variances = list(respVar = TRUE)))
    noise_present <- length(attr(link, "respVar")) > 0 || length(attr(link, "residVar")) > 0
    expect_identical(attr(link, "prediction_contract")$uncertainty$observation_noise, noise_present)
  }
})

test_that("HLfit penalties include every native DoF component and retain restricted comparisons", {
  skip_if_not_installed("spaMM")
  for (family in c("gaussian", "poisson")) {
    f <- spamm_native_fixture(family)
    model <- f$object$opt_mod
    actual <- aic_tab(list(f$object), AICc = FALSE)
    expect_equal(actual$K, if (family == "gaussian") 9 else 8)
    expect_equal(actual$AIC, unname(extractAIC(model)["AIC"]) + 4, tolerance = 1e-9)
    expect_equal(bic_tab(list(f$object))$BIC,
      -2 * as.numeric(logLik(model)) + log(nrow(f$data)) * actual$K, tolerance = 1e-9)
    expect_identical(rownames(.msr_selection_frame(model)), rownames(f$data))
  }
  f <- spamm_native_fixture()
  reml <- spaMM::fitme(y ~ mean + rough + fixed + z + (1 | site), data = f$data,
    family = gaussian(), method = "REML")
  other <- spaMM::fitme(y ~ mean + rough + fixed + I(z^2) + (1 | site), data = f$data,
    family = gaussian(), method = "REML")
  expect_error(aic_tab(list(f$object$opt_mod, reml), AICc = FALSE), "ML and REML")
  expect_error(bic_tab(list(reml, other)), "same recoverable fixed-effect design")
  signature <- .msr_selection_signature(reml)
  expect_equal(unname(as.matrix(signature$design)),
    unname(as.matrix(stats::model.matrix(reml))[order(rownames(f$data)), , drop = FALSE]))
})

test_that("HLfit partial-effect rows retain only selected grouping inputs", {
  skip_if_not_installed("spaMM")
  f <- spamm_native_fixture()
  dat <- f$data; dat$other <- factor(rep(1:4, 9)); dat$y <- dat$y + as.integer(dat$other)^2
  model <- spaMM::fitme(y ~ mean + rough + fixed + z + (1 | site) + (1 | other),
    data = dat, family = gaussian(), method = "ML")
  f$object$opt_mod <- model
  nd <- dat; nd$other <- NULL
  control <- ~ (1 | site)
  expected <- predict(model, nd, re.form = control, type = "response")
  actual <- predict(f$object, nd, f$inputs, re.form = control, include_fixed = TRUE)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
  expect_identical(attr(actual, "prediction_contract")$effects$response, "selected_conditional_effects")
})


test_that("real HLfit two-worker maps preserve zero-effect point predictions", {
  skip_on_cran()
  skip_if_not_installed("spaMM")
  f <- spamm_native_fixture()
  expect_warning(map <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    include_fixed = TRUE, scale_center = TRUE, match_buffer = TRUE,
    site_values = list(z = .2), verbose = FALSE), "site.*categorical")
  serial <- predict(f$object, map, re.form = NA)
  parallel <- predict(f$object, map, re.form = NA, raster_args = list(cores = 2))
  expect_equal(terra::values(parallel), terra::values(serial), tolerance = 1e-12)
  expect_identical(attr(parallel, "prediction_contract")$integration$random_effects,
    "set_to_zero_not_integrated")
})

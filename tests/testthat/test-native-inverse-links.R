# Independent native references use fixed selected scales, never a new scale fit.
check_native_gamma_compatibility <- function(kind) {
  f <- native_gamma_fixture(kind)
  model <- f$object$opt_mod
  nd <- f$data[c(5, 23, 17), ]
  for (type in c("response", "link")) {
    actual <- predict(f$object, nd, f$inputs, type = type, include_fixed = TRUE)
    native <- predict(model, nd, type = type)
    expect_equal(as.numeric(actual), as.numeric(native), tolerance = 1e-10)
    expect_false(attr(actual, "prediction_contract")$uncertainty$observation_noise)
    snapshot <- prepare_prediction(f$object, nd, f$inputs, component = type, include_fixed = TRUE)
    expect_equal(predict(unserialize(serialize(snapshot, NULL))), actual, tolerance = 1e-10)
    expect_null(snapshot$opt_mod)
    if (kind == "glm") {
      se <- predict(f$object, nd, f$inputs, type = type, se.fit = TRUE,
        uncertainty = "conditional", include_fixed = TRUE)
      native_se <- predict(model, nd, type = type, se.fit = TRUE)
      expect_equal(se$fit, native_se$fit, tolerance = 1e-10)
      expect_equal(se$se.fit, native_se$se.fit, tolerance = 1e-10)
      expect_equal(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"), se)
    } else expect_error(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"), "did not supply se.fit")
  }
  map <- kernel_scale.raster(f$rasters, multiScaleR = f$object, include_fixed = TRUE,
    scale_center = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  reference <- f$data; reference$z <- .2
  for (type in c("response", "link")) {
    actual <- predict(f$object, map, type = type)
    expect_equal(terra::extract(actual, f$points, ID = FALSE)[[1]],
      as.numeric(predict(model, reference, type = type)), tolerance = 1e-8)
    if (kind == "glm") {
      actual_se <- predict(f$object, map, type = type, se.fit = TRUE, uncertainty = "conditional")
      native_se <- predict(model, reference, type = type, se.fit = TRUE)
      expect_equal(as.numeric(as.matrix(terra::extract(actual_se, f$points, ID = FALSE))),
        as.numeric(cbind(native_se$fit, native_se$se.fit)), tolerance = 1e-8)
    }
  }
  native_k <- if (kind == "glm") attr(logLik(model), "df") else sum(spaMM::DoF(model))
  expect_equal(native_k, 6)
  expect_equal(aic_tab(list(f$object), AICc = FALSE)$K, native_k + 2)
  expect_equal(bic_tab(list(f$object))$BIC,
    -2 * as.numeric(logLik(model)) + log(nrow(f$data)) * (native_k + 2), tolerance = 1e-9)

  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(f$object, length.out = 7)
  grid <- as.data.frame(lapply(f$data[c("mean", "rough", "fixed", "z")], function(x) rep(mean(x), 7)))
  grid$rough <- seq(min(f$data$rough), max(f$data$rough), length.out = 7)
  native <- if (kind == "glm") predict(model, grid, type = "link", se.fit = TRUE) else
    predict(model, grid, type = "link", re.form = NA, variances = list(fixefVar = TRUE))
  eta <- if (kind == "glm") native$fit else as.numeric(native)
  se <- if (kind == "glm") native$se.fit else sqrt(attr(native, "fixefVar"))
  family <- stats::family(model)
  expect_true(all(eta + qnorm(.025) * se > 0))
  expect_equal(plots$rough$data$fit, as.numeric(family$linkinv(eta)), tolerance = 1e-10)
  # An inverse link reverses the endpoint order on this strictly positive branch.
  expect_equal(plots$rough$data$lwr, as.numeric(family$linkinv(eta + qnorm(.975) * se)), tolerance = 1e-10)
  expect_equal(plots$rough$data$upr, as.numeric(family$linkinv(eta + qnorm(.025) * se)), tolerance = 1e-10)
  expect_true(all(plots$rough$data$lwr <= plots$rough$data$fit))
  expect_true(all(plots$rough$data$fit <= plots$rough$data$upr))
  expect_true(all(plots$rough$data$lwr > 0))
}

test_that("Gamma inverse glm features, maps, snapshots, penalties and plots match native targets", {
  check_native_gamma_compatibility("glm")
})

test_that("Gamma inverse HLfit features and ordered mean limits match native targets", {
  skip_if_not_installed("spaMM")
  check_native_gamma_compatibility("HLfit")
})

test_that("native link confidence limits never bridge a pole or manufacture finite SEs", {
  f <- native_gamma_fixture("glm")
  expect_warning(bounds <- .msr_native_link_ci(f$object$opt_mod,
    c(1, 1, 1), c(1, Inf, NA_real_), stats::family(f$object$opt_mod)$linkinv),
    "not shown for 3 row")
  expect_true(all(is.na(bounds$lower)))
  expect_true(all(is.na(bounds$upper)))
  gaussian_inverse <- glm(y ~ mean + rough + fixed + z,
    data = f$data, family = gaussian(link = "inverse"))
  expect_warning(bounds <- .msr_native_link_ci(gaussian_inverse, 1, 1,
    stats::family(gaussian_inverse)$linkinv), "not shown")
  expect_true(is.na(bounds$lower))
  expect_true(is.na(bounds$upper))
})


test_that("plots retain native fitted values when inverse-link confidence limits are unavailable", {
  f <- native_gamma_fixture("glm")
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  warnings <- character()
  plots <- withCallingHandlers(with_mocked_bindings(
    plot_marginal_effects(f$object, length.out = 7),
    predict = function(object, newdata, ...) list(fit = rep(1, nrow(newdata)),
      se.fit = rep(1, nrow(newdata)))), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
  expect_length(warnings, length(plots))
  expect_true(all(grepl("Native mean confidence limits were not shown for 7 row(s)",
    warnings, fixed = TRUE)))
  expect_equal(plots$rough$data$fit, rep(1, 7))
  expect_true(all(is.na(plots$rough$data$lwr)))
  expect_true(all(is.na(plots$rough$data$upr)))
  expect_match(plots$rough$labels$caption, "Confidence limits unavailable")
})

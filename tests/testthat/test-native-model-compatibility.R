test_that("native LMM and GLMM point/raster adapters retain their prediction target", {
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  for (kind in c("lmer", "glmer")) {
    model <- if (kind == "lmer")
      lme4::lmer(y ~ mean + rough + fixed + z + (1 | site), dat) else
      lme4::glmer(y ~ mean + rough + fixed + z + (1 | site), dat, family = poisson())
    object <- f$object; object$opt_mod <- model
    nd <- dat[c(2, 27, 13), ]
    expected <- predict(model, nd, type = "response", re.form = NA)
    actual <- predict(object, nd, f$inputs, re.form = NA, include_fixed = TRUE)
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
    inputs <- kernel_scale.raster(f$rasters, multiScaleR = object,
      scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
      site_values = list(z = .2, site = levels(dat$site)[1]))
    output <- predict(object, inputs, re.form = NA)
    reference <- dat; reference$z <- .2
    expect_equal(terra::extract(output, f$points, ID = FALSE)[[1]],
      as.numeric(predict(model, reference, type = "response", re.form = NA)), tolerance = 1e-8)
  }
})

test_that("native lme controls and retained fixed normalization work together", {
  skip_if_not_installed("nlme")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  model <- nlme::lme(y ~ mean + rough + fixed + z, random = ~1 | site,
                      data = dat, method = "ML", keep.data = TRUE)
  object <- f$object; object$opt_mod <- model
  nd <- dat[c(2, 27, 13), ]
  actual <- predict(object, nd, f$inputs, include_fixed = TRUE,
                      native_args = list(level = 0))
  expect_equal(as.numeric(actual), as.numeric(predict(model, nd, level = 0)), tolerance = 1e-10)
  snap <- prepare_prediction(object, nd, f$inputs, include_fixed = TRUE,
                             native_args = list(level = 0))
  expect_equal(as.numeric(predict(snap)), as.numeric(actual))
})

test_that("NB2, hurdle, and zero-inflated refit adapters retain native mean predictions", {
  f <- refit_spatial_fixture("exp")
  dat <- f$data
  dat$y <- ifelse(seq_len(nrow(dat)) %% 4 == 0, 0, 4 + seq_len(nrow(dat)) %% 5)
  models <- list()
  if (requireNamespace("MASS", quietly = TRUE)) models$nbinom2 <-
    MASS::glm.nb(y ~ mean + rough + fixed + z, dat)
  if (requireNamespace("pscl", quietly = TRUE)) {
    models$zeroinfl <- pscl::zeroinfl(y ~ mean + rough + fixed + z | 1, dat)
    models$hurdle <- pscl::hurdle(y ~ mean + rough + fixed + z | 1, dat)
  }
  skip_if(!length(models), "No optional count-model providers installed")
  for (model in models) {
    object <- f$object; object$opt_mod <- model
    nd <- dat[c(4, 29, 18), ]
    expect_equal(as.numeric(predict(object, nd, f$inputs, include_fixed = TRUE)),
      as.numeric(predict(model, nd, type = "response")), tolerance = 1e-10)
    snap <- prepare_prediction(object, nd, f$inputs, component = "response", include_fixed = TRUE)
    expect_equal(as.numeric(predict(snap)), as.numeric(predict(model, nd, type = "response")),
                   tolerance = 1e-10)
  }
})

test_that("unmarked shared point prediction preserves its native component", {
  skip_if_not_installed("unmarked")
  f <- refit_unmarked_spatial_fixture()
  object <- f$object
  x <- f$inputs
  dat <- f$data
  actual <- predict(object, dat, x, type = "state")
  expected <- unmarked::predict(object$opt_mod, newdata = dat, type = "state")
  expect_equal(as.numeric(actual), as.numeric(expected$Predicted), tolerance = 1e-10)
  uncertain <- predict(object, dat, x, type = "state", se.fit = TRUE,
                         uncertainty = "conditional")
  expect_equal(uncertain$Predicted, expected$Predicted, tolerance = 1e-10)
  expect_equal(uncertain$SE, expected$SE, tolerance = 1e-10)
})

test_that("analytic SQ fits use checked rasters without claiming joint scale uncertainty", {
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(442)
  dat$y <- rpois(nrow(dat), exp(3 + 1.5 * dat$rough + .1 * dat$z))
  model <- glm(y ~ rough + z, dat, family = poisson())
  object <- suppressMessages(multiScale_optim(model, f$inputs, backend = "profile_glm",
                                              verbose = FALSE))
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = object,
    scale_center = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  nd <- dat; nd$z <- .2
  expected <- predict(object, nd, f$inputs, se.fit = TRUE, uncertainty = "conditional")
  output <- predict(object, inputs, se.fit = TRUE, uncertainty = "conditional")
  info <- attr(output, "prediction_contract")$uncertainty
  expect_identical(info$standard_errors, "provided")
  expect_identical(info$intervals, "conditional_standard_errors_only")
  point_map <- predict(object, inputs)
  info <- attr(point_map, "prediction_contract")$uncertainty
  expect_identical(info$standard_errors, "not_requested")
  expect_identical(info$intervals, "not_provided")
  expect_equal(terra::extract(point_map, f$points, ID = FALSE)[[1]],
    as.numeric(expected$fit), tolerance = 1e-8)
  actual <- terra::extract(output, f$points, ID = FALSE)
  expect_equal(actual$fit, as.numeric(expected$fit), tolerance = 1e-8)
  expect_equal(actual$se.fit, as.numeric(expected$se.fit), tolerance = 1e-8)
  expect_error(predict(object, inputs, se.fit = TRUE, uncertainty = "joint"), "conditional SEs")
})

test_that("joint NB2 weighted-SQ rasters use the analytic conditional covariance", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("RTMB")
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(443)
  dat$y <- rnbinom(nrow(dat), mu = exp(3 + 1.5 * dat$rough + .1 * dat$z), size = 3)
  model <- MASS::glm.nb(y ~ rough + z, data = dat)
  # This bounded software fixture reaches a scale bound. Its point and
  # conditional predictions remain available, while joint uncertainty must fail.
  expect_warning(object <- suppressMessages(multiScale_optim(model, f$inputs,
    backend = "joint_nb2", verbose = FALSE)), "Joint Wald uncertainty unavailable")
  expect_false(object$uncertainty$available)
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = object, scale_center = TRUE,
    match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  nd <- dat; nd$z <- .2
  expected <- predict(object, nd, f$inputs, se.fit = TRUE, uncertainty = "conditional")
  output <- predict(object, inputs, se.fit = TRUE, uncertainty = "conditional")
  info <- attr(output, "prediction_contract")$uncertainty
  expect_identical(info$standard_errors, "provided")
  expect_identical(info$intervals, "conditional_standard_errors_only")
  point_map <- predict(object, inputs)
  info <- attr(point_map, "prediction_contract")$uncertainty
  expect_identical(info$standard_errors, "not_requested")
  expect_identical(info$intervals, "not_provided")
  expect_equal(terra::extract(point_map, f$points, ID = FALSE)[[1]],
    as.numeric(expected$fit), tolerance = 1e-8)
  actual <- terra::extract(output, f$points, ID = FALSE)
  expect_equal(actual$fit, as.numeric(expected$fit), tolerance = 1e-8)
  expect_equal(actual$se.fit, as.numeric(expected$se.fit), tolerance = 1e-8)
  expect_error(predict(object, nd, f$inputs, se.fit = TRUE, uncertainty = "joint"),
    "Joint Wald uncertainty unavailable")
})

test_that("two-worker raster prediction uses the same package and native controls", {
  skip_on_cran()
  f <- refit_spatial_fixture()
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
    site_values = list(group = "b", z = .2), offset_values = c(exposure = 2))
  serial <- predict(f$object, inputs)
  parallel <- predict(f$object, inputs, raster_args = list(cores = 2))
  expect_equal(terra::values(parallel), terra::values(serial), tolerance = 1e-12)
})

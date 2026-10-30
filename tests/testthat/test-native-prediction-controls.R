test_that("zero random effects use syntax and retain fixed-variable checks", {
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(182)
  dat$only_random <- 1.5 + cos(seq_len(nrow(dat)) / 3)
  dat$y <- 3 + dat$mean + dat$rough + as.integer(dat$site) * dat$only_random +
    rnorm(nrow(dat), sd = .2)
  model <- lme4::lmer(y ~ mean + rough + fixed + z + (0 + only_random | site),
    data = dat)
  object <- f$object
  object$opt_mod <- model
  rows <- c(3, 22, 11)
  complete <- dat[rows, ]
  nd <- complete
  nd[c("site", "only_random", "fixed")] <- NULL
  controls <- list(NA, stats::as.formula("~0", env = new.env(parent = baseenv())))
  for (control in controls) {
    expected <- predict(model, complete, re.form = control)
    actual <- predict(object, nd, f$inputs, include_fixed = TRUE, re.form = control)
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
    contract <- attr(actual, "prediction_contract")
    expect_identical(contract$integration$random_effects, "set_to_zero_not_integrated")
    snap <- prepare_prediction(object, nd, f$inputs, include_fixed = TRUE,
      native_args = list(re.form = control))
    expect_equal(predict(snap), actual)
    expect_identical(snap$prediction_contract$response$effects$response,
      "zero_deviation_reference")
  }
  expect_error(predict(object, nd, f$inputs, include_fixed = TRUE),
    "formula variables")
  bad <- nd
  bad$z <- NULL
  expect_error(predict(object, bad, f$inputs, include_fixed = TRUE, re.form = NA),
    "formula variables")
  expect_false(.msr_native_zero_random_effects(model, list(re.form = ~1 | site)))
  expect_false(.msr_native_zero_random_effects(model, list(re.form = NULL)))

  rasters <- suppressWarnings(kernel_scale.raster(f$rasters, multiScaleR = object,
    include_fixed = TRUE, scale_center = TRUE, match_buffer = TRUE,
    site_values = list(z = .2), verbose = FALSE))
  expect_false("site" %in% names(rasters))
  map <- predict(object, rasters, re.form = controls[[2]])
  reference <- dat
  reference$z <- .2
  expect_equal(terra::extract(map, f$points, ID = FALSE)[[1]],
    as.numeric(predict(model, reference, re.form = controls[[2]])), tolerance = 1e-8)
})

test_that("nlme population controls accept integer and double zero levels", {
  skip_if_not_installed("nlme")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  model <- nlme::lme(y ~ mean + rough + fixed + z, random = ~1 | site,
    data = dat, method = "ML", keep.data = TRUE)
  object <- f$object
  object$opt_mod <- model
  nd <- dat[c(3, 22, 11), ]
  nd$site <- NULL
  for (level in list(0, 0L)) {
    actual <- predict(object, nd, f$inputs, include_fixed = TRUE,
      native_args = list(level = level))
    expect_equal(as.numeric(actual), as.numeric(predict(model, nd, level = level)),
      tolerance = 1e-10)
    expect_identical(attr(actual, "prediction_contract")$integration$random_effects,
      "set_to_zero_not_integrated")
    snap <- prepare_prediction(object, nd, f$inputs, include_fixed = TRUE,
      native_args = list(level = level))
    expect_equal(predict(snap), actual)
  }
  expect_error(predict(object, nd, f$inputs, include_fixed = TRUE),
    "group|formula variables")
})

test_that("GAM snapshots retain native components and covariance choices", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(185)
  dat$y <- 2 + dat$mean + dat$rough + sin(dat$z * 2) + rnorm(nrow(dat), sd = .2)
  model <- mgcv::gam(y ~ mean + rough + fixed + s(z, k = 4), data = dat, method = "ML")
  object <- f$object
  object$opt_mod <- model
  nd <- dat[c(2, 25, 13), ]
  for (type in c("response", "link", "terms", "iterms", "lpmatrix")) {
    controls <- list(exclude = "s(z)", unconditional = TRUE, block.size = 2)
    expected <- do.call(stats::predict, c(list(object = model, newdata = nd,
      type = type, se.fit = FALSE), controls))
    direct <- predict(object, nd, f$inputs, type = type, include_fixed = TRUE,
      native_args = controls)
    snap <- prepare_prediction(object, nd, f$inputs, component = type,
      include_fixed = TRUE, native_args = controls)
    expect_equal(as.numeric(direct), as.numeric(expected), tolerance = 1e-10)
    expect_equal(predict(snap), direct)
    expect_null(snap$opt_mod)
    expect_identical(snap$prediction_contract[[type]]$uncertainty$smoothing_parameters,
      "native_unconditional_requested")
    if (type == "lpmatrix") {
      expect_error(predict(snap, se.fit = TRUE, uncertainty = "conditional"),
        "did not supply se.fit")
    } else {
      expected_se <- do.call(stats::predict, c(list(object = model, newdata = nd,
        type = type, se.fit = TRUE), controls))
      got <- predict(snap, se.fit = TRUE, uncertainty = "conditional")
      expect_equal(as.numeric(got$fit), as.numeric(expected_se$fit), tolerance = 1e-10)
      expect_equal(as.numeric(got$se.fit), as.numeric(expected_se$se.fit), tolerance = 1e-10)
    }
  }
  snap <- prepare_prediction(object, nd, f$inputs, component = "iterms", include_fixed = TRUE,
    native_args = list(terms = "s(z)", iterms.type = 2))
  expect_identical(snap$prediction_contract$iterms$uncertainty$smoothing_parameters,
    "native_default_conditional")
  expected <- predict(model, nd, type = "iterms", terms = "s(z)",
    iterms.type = 2, se.fit = TRUE)
  actual <- predict(snap, se.fit = TRUE, uncertainty = "conditional")
  expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-10)
  expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-10)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(snap, path)
  expect_equal(predict(readRDS(path)), predict(snap))
  expect_error(prepare_prediction(f$object, nd, f$inputs,
    native_args = list(unconditional = TRUE)), "supported native")
  expect_error(prepare_prediction(object, nd, f$inputs,
    native_args = list(newdata.guaranteed = TRUE)), "supported native")
})

test_that("parallel zero-effect rasters use installed native row controls", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  object <- f$object
  object$opt_mod <- lme4::lmer(y ~ mean + rough + fixed + z + (1 | site), dat)
  rasters <- suppressWarnings(kernel_scale.raster(f$rasters, multiScaleR = object,
    include_fixed = TRUE, scale_center = TRUE, match_buffer = TRUE,
    site_values = list(z = .2), verbose = FALSE))
  expect_false("site" %in% names(rasters))
  control <- stats::as.formula("~0", env = new.env(parent = baseenv()))
  serial <- predict(object, rasters, re.form = control)
  parallel <- predict(object, rasters, re.form = control, raster_args = list(cores = 2))
  expect_equal(terra::values(parallel), terra::values(serial), tolerance = 1e-12)
})

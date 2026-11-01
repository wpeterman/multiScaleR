test_that("GAM raster components and native covariance match prepared points", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(1381)
  dat$y <- 2 + dat$mean + dat$rough + sin(2 * dat$z) + rnorm(nrow(dat), sd = .2)
  object <- f$object
  model <- object$opt_mod <- mgcv::gam(y ~ mean + rough + fixed + s(z, k = 4),
    data = dat, method = "ML")
  rasters <- kernel_scale.raster(f$rasters, multiScaleR = object,
    include_fixed = TRUE, scale_center = TRUE, match_buffer = TRUE,
    site_values = list(z = .2), verbose = FALSE)
  reference <- dat
  reference$z <- .2
  controls <- list(exclude = "s(z)", unconditional = TRUE)
  for (type in c("response", "terms", "iterms", "lpmatrix")) {
    se <- type != "lpmatrix"
    expected <- do.call(stats::predict, c(list(object = model, newdata = reference,
      type = type, se.fit = se), controls))
    point <- predict(object, reference, f$inputs, type = type, se.fit = se,
      uncertainty = "conditional", include_fixed = TRUE, native_args = controls)
    map <- predict(object, rasters, type = type, se.fit = se,
      uncertainty = "conditional", native_args = controls)
    actual <- as.matrix(terra::extract(map, f$points, ID = FALSE))
    target <- if (se) cbind(expected$fit, expected$se.fit) else as.matrix(expected)
    expect_equal(dim(actual), dim(target))
    expect_equal(as.numeric(actual), as.numeric(target), tolerance = 1e-8)
    if (se) {
      expect_equal(as.numeric(point$fit), as.numeric(expected$fit), tolerance = 1e-10)
      expect_equal(as.numeric(point$se.fit), as.numeric(expected$se.fit), tolerance = 1e-10)
    }
    contract <- attr(map, "prediction_contract")
    expect_identical(contract$uncertainty$standard_errors,
      if (se) "provided" else "not_requested")
    expect_identical(contract$uncertainty$intervals,
      if (se) "conditional_standard_errors_only" else "not_provided")
    expect_identical(contract$term_selection$exclude, "s(z)")
    expect_identical(contract$uncertainty$smoothing_parameters,
      "native_unconditional_requested")
    expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
    expect_false(contract$uncertainty$observation_noise)
    expect_true(all(nzchar(names(map))))
    expect_false(anyDuplicated(names(map)) > 0)
  }
  expect_error(predict(object, rasters, se.fit = TRUE, uncertainty = "joint"),
    "conditional SEs")
})

test_that("parallel partial-effect maps retain selected native conditioning", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$other <- factor(rep(1:4, 9))
  dat$y <- dat$y + as.integer(dat$site)^2 + as.integer(dat$other)^2
  object <- f$object
  model <- object$opt_mod <- lme4::lmer(y ~ mean + rough + fixed + z +
    (1 | site) + (1 | other), data = dat)
  controls <- list(re.form = stats::as.formula("~(1|site)",
    env = new.env(parent = baseenv())))
  rasters <- suppressWarnings(kernel_scale.raster(f$rasters, multiScaleR = object,
    include_fixed = TRUE, scale_center = TRUE, match_buffer = TRUE,
    site_values = list(z = .2, site = levels(dat$site)[1]), verbose = FALSE))
  expect_false("other" %in% names(rasters))
  serial <- predict(object, rasters, native_args = controls)
  parallel <- predict(object, rasters, native_args = controls, raster_args = list(cores = 2))
  expect_equal(terra::values(parallel), terra::values(serial), tolerance = 1e-12)
  reference <- dat
  reference$z <- .2
  reference$site <- factor(levels(dat$site)[1], levels = levels(dat$site))
  expect_equal(terra::extract(parallel, f$points, ID = FALSE)[[1]],
    as.numeric(predict(model, reference, re.form = controls$re.form)), tolerance = 1e-8)
  expect_identical(attr(parallel, "prediction_contract")$effects$response,
    "selected_conditional_effects")
})

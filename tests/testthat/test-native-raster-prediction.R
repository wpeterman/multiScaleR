test_that("checked raster prediction preserves native means, SEs, and fixed features", {
  f <- refit_spatial_fixture()
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
    site_values = list(group = "b", z = .2), offset_values = c(exposure = 2))
  nd <- f$data
  nd$group <- factor("b", levels = levels(nd$group))
  nd$z <- .2; nd$exposure <- 2
  expected <- predict(f$object$opt_mod, nd, type = "response", se.fit = TRUE)
  output <- predict(f$object, inputs, se.fit = TRUE, uncertainty = "conditional")
  actual <- terra::extract(output, f$points, ID = FALSE)
  expect_equal(actual$fit, as.numeric(expected$fit), tolerance = 1e-8)
  expect_equal(actual$se.fit, as.numeric(expected$se.fit), tolerance = 1e-8)
  expect_identical(attr(output, "prediction_contract")$uncertainty$kernel_parameters,
                    "fixed_at_estimate")
  expect_identical(attr(output, "projection_contract")$scenario_values,
                    list(group = "b", z = .2))
  expect_identical(attr(output, "projection_contract")$controls$offset_values,
                    c(exposure = 2))
  expect_false(attr(output, "projection_contract")$controls$clamp)
  expect_error(predict(f$object, inputs, raster_args = list(const = data.frame(rough = 0))),
    "const inputs must add variables")
  expect_error(predict(f$object, inputs, native_args = stats::setNames(list(TRUE), "")),
    "uniquely name native")
  expect_error(predict(f$object, inputs, se.fit = TRUE), "conditional SEs")
  expect_error(predict(f$object, inputs, interval = "confidence"), "conditional SEs")
  expect_error(predict(f$object, inputs, kernel_inputs = f$inputs), "do not take kernel_inputs")
  expect_error(predict(f$object, f$rasters), "kernel_scale.raster")
  wrong <- f$object
  wrong$scale_est$Mean[1] <- wrong$scale_est$Mean[1] + 1
  expect_error(predict(wrong, inputs), "this spatial fit")
  descriptor <- attr(inputs, "multiScaleR_projection")
  altered <- inputs
  descriptor$specifications$metric[descriptor$specifications$covariate == "rough"] <- "sa"
  attr(altered, "multiScaleR_projection") <- descriptor
  expect_error(predict(f$object, altered), "source and metric specifications")
})

test_that("raster file and constant controls retain native logical inputs", {
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$flag <- rep(c(TRUE, FALSE), length.out = nrow(dat))
  f$object$opt_mod <- glm(y ~ mean + rough + fixed + flag, dat, family = poisson())
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE,
    site_values = list(flag = TRUE), verbose = FALSE)
  path <- tempfile(fileext = ".tif")
  on.exit(unlink(path), add = TRUE)
  output <- predict(f$object, inputs, raster_args = list(filename = path))
  dat$flag <- TRUE
  expect_true(file.exists(path))
  expect_equal(terra::extract(output, f$points, ID = FALSE)[[1]],
    as.numeric(predict(f$object$opt_mod, dat, type = "response")), tolerance = 1e-7)
  expect_error(predict(f$object, inputs, raster_args = list(unrecognized = 1)),
                "supported terra")
})

test_that("additional categorical constants work without replacing computed raster predictors", {
  f <- refit_spatial_fixture()
  expect_warning(inputs <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, clamp = TRUE,
    site_values = list(z = .2), offset_values = c(exposure = 2), verbose = FALSE),
    "Site-level covariate 'group' is categorical")
  scenario <- data.frame(group = factor("b", levels = levels(f$data$group)))
  output <- predict(f$object, inputs, raster_args = list(const = scenario))
  nd <- f$data
  nd$group <- factor("b", levels = levels(nd$group))
  nd$z <- .2; nd$exposure <- 2
  expect_equal(terra::extract(output, f$points, ID = FALSE)[[1]],
    as.numeric(predict(f$object$opt_mod, nd, type = "response")), tolerance = 1e-8)
  expect_true(attr(output, "projection_contract")$controls$clamp)
})

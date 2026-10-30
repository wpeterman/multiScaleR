test_that("explicit numeric and categorical raster scenarios preserve native predictions", {
  f <- refit_spatial_fixture()
  rasters <- expect_no_warning(kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
    site_values = list(group = "b", z = .2), offset_values = c(exposure = 2)))
  expect_true(terra::is.factor(rasters[["group"]]))
  expect_identical(terra::levels(rasters[["group"]])[[1]]$label, c("a", "b", "c"))
  nd <- f$data
  nd$group <- factor("b", levels = levels(nd$group))
  nd$z <- .2
  nd$exposure <- 2
  expected <- predict(f$object$opt_mod, nd, type = "response")
  prediction <- predict(f$object, rasters, type = "response")
  actual <- terra::extract(prediction, f$points, ID = FALSE)[[1]]
  expect_equal(actual, as.numeric(expected), tolerance = 1e-8)
})

test_that("raster scenarios use raw transform and multi-input offset variables", {
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$area <- 2 + seq_len(nrow(dat)) / nrow(dat)
  f$object$opt_mod <- glm(y ~ mean + rough + fixed + group + I(z^2) +
    offset(log(exposure / area)), data = dat, family = poisson())
  args <- list(raster_stack = f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
    site_values = list(group = "c", z = 2), offset_values = c(exposure = 2, area = 4))
  rasters <- expect_no_warning(do.call(kernel_scale.raster, args))
  expect_true(all(c("z", "exposure", "area", "group") %in% names(rasters)))
  expect_false("I(z^2)" %in% names(rasters))
  dat$group <- factor("c", levels = levels(dat$group))
  dat$z <- 2
  dat$exposure <- 2
  dat$area <- 4
  prediction <- predict(f$object, rasters, type = "response")
  expect_equal(terra::extract(prediction, f$points, ID = FALSE)[[1]],
    as.numeric(predict(f$object$opt_mod, dat, type = "response")), tolerance = 1e-8)
  args$offset_values <- c(exposure = 2)
  expect_error(do.call(kernel_scale.raster, args), "explicit raw values")
})

test_that("raster scenario controls reject ambiguous or unsupported values", {
  f <- refit_spatial_fixture()
  args <- list(raster_stack = f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE)
  expect_error(do.call(kernel_scale.raster, c(args, list(site_values = list(group = "new")))),
                "fitted level")
  expect_error(do.call(kernel_scale.raster, c(args, list(site_values = list(rough = 0)))),
                "nonraster model predictors")
  expect_error(do.call(kernel_scale.raster, c(args, list(site_values = list(unknown = 0)))),
                "nonraster model predictors")
  expect_error(do.call(kernel_scale.raster, c(args, list(site_values = list(z = 1:2)))),
                "scalar scenario")
  expect_error(do.call(kernel_scale.raster, c(args, list(site_values = list(exposure = 2),
    offset_values = c(exposure = 2)))), "not both")
  expect_error(do.call(kernel_scale.raster, c(args, list(offset_values = c(typo = 2)))),
                "names not used as formula offsets")
  expect_error(kernel_scale.raster(f$rasters, sigma = c(19, 29),
    site_values = list(group = "b"), verbose = FALSE), "requires a fitted")
})

test_that("logical raster scenarios retain the native logical predictor meaning", {
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$flag <- rep(c(TRUE, FALSE), length.out = nrow(dat))
  f$object$opt_mod <- glm(y ~ mean + rough + fixed + flag, dat, family = poisson())
  rasters <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE,
    site_values = list(flag = TRUE), verbose = FALSE)
  prediction <- predict(f$object, rasters, type = "response")
  dat$flag <- TRUE
  expect_equal(terra::extract(prediction, f$points, ID = FALSE)[[1]],
    as.numeric(predict(f$object$opt_mod, dat, type = "response")), tolerance = 1e-8)
})

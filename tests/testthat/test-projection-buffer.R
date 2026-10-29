# Established raster fixture and its cell centers only. These are numerical
# transform checks, not a new ecological dataset or optimizer assessment.
projection_buffer_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- make_unmarked_fixture()
    source <- terra::subset(f$rs_two, "cont1")
    # The existing simulated grid is in arbitrary projected meter units but
    # has no CRS label. Label this test copy explicitly without resampling it.
    terra::crs(source) <- "EPSG:3857"
    cells <- c(1L, terra::ncol(source), terra::ncell(source),
               20L * terra::ncol(source) + 22L, 31L * terra::ncol(source) + 17L)
    pts <- terra::vect(terra::xyFromCell(source, cells), crs = terra::crs(source))
    cache <<- list(source = source, cells = cells, pts = pts)
    cache
  }
})

test_that("prepared circular buffers match exact point means and weighted SQ", {
  f <- projection_buffer_fixture()
  specs <- msr_vars(mean = kernel_var("cont1"),
                    rough = surface_var("cont1", "sq", weighted = TRUE))
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    inputs <- kernel_prep(f$pts, f$source, max_D = 150, sigma = c(60, 100),
      shape = if (kernel == "expow") c(1.5, 2.5) else NULL,
      kernel = kernel, scale_vars = specs, bin = FALSE, verbose = FALSE)
    projected <- kernel_scale.raster(f$source, multiScaleR = inputs,
      match_buffer = TRUE, verbose = FALSE)
    expected <- sweep(sweep(as.matrix(inputs$kernel_dat), 2L,
      inputs$scl_params$sd, "*"), 2L, inputs$scl_params$mean, "+")
    actual <- as.matrix(terra::extract(projected, f$pts)[, -1L])
    expect_equal(actual, expected, tolerance = 1e-10, ignore_attr = TRUE)
    # pct_wt does not change the explicit prepared support.
    same <- kernel_scale.raster(f$source, multiScaleR = inputs,
      match_buffer = TRUE, pct_wt = .9, verbose = FALSE)
    expect_identical(terra::values(same), terra::values(projected))
  }
})

test_that("matched higher weighted moments retain their native point definitions", {
  f <- projection_buffer_fixture()
  for (metric in c("sa", "ssk", "sku")) {
    specs <- msr_vars(texture = surface_var("cont1", metric, weighted = TRUE))
    inputs <- kernel_prep(f$pts, f$source, max_D = 150, sigma = 100,
      kernel = "gaussian", scale_vars = specs, bin = FALSE, verbose = FALSE)
    projected <- kernel_scale.raster(f$source, multiScaleR = inputs,
      match_buffer = TRUE, verbose = FALSE)
    expected <- inputs$kernel_dat$texture * inputs$scl_params$sd["texture"] +
      inputs$scl_params$mean["texture"]
    expect_equal(as.numeric(terra::extract(projected, f$pts)[, 2L]),
      as.numeric(expected), tolerance = 1e-10)
  }
})

test_that("matched projection normalizes finite neighbors but preserves source masks", {
  f <- projection_buffer_fixture()
  source <- f$source
  z <- terra::values(source)
  z[c(2L, 51L, f$cells[5L]), 1L] <- NA_real_
  terra::values(source) <- z
  specs <- msr_vars(mean = kernel_var("cont1"),
                    rough = surface_var("cont1", "sq", weighted = TRUE))
  inputs <- kernel_prep(f$pts, source, max_D = 150, sigma = c(80, 80),
    kernel = "gaussian", scale_vars = specs, bin = FALSE, verbose = FALSE)
  projected <- kernel_scale.raster(source, multiScaleR = inputs,
    match_buffer = TRUE, verbose = FALSE)
  expected <- sweep(sweep(as.matrix(inputs$kernel_dat), 2L,
    inputs$scl_params$sd, "*"), 2L, inputs$scl_params$mean, "+")
  actual <- as.matrix(terra::extract(projected, f$pts)[, -1L])
  expect_equal(actual[1:4, ], expected[1:4, ], tolerance = 1e-10, ignore_attr = TRUE)
  expect_true(all(is.na(actual[5L, ])))
  expect_true(all(is.na(terra::values(projected)[is.na(z[, 1L]), ])))
  # Point preparation permits neighbors around a missing center. The raster
  # intentionally keeps its existing domain mask, as documented by this opt-in.
  expect_true(all(is.finite(expected[5L, ])))
})

test_that("matched constant fields retain constant means and zero population SQ", {
  f <- projection_buffer_fixture()
  specs <- msr_vars(mean = kernel_var("cont1"),
                    rough = surface_var("cont1", "sq", weighted = TRUE))
  inputs <- kernel_prep(f$pts, f$source, max_D = 150, sigma = c(80, 80),
    kernel = "gaussian", scale_vars = specs, bin = FALSE, verbose = FALSE)
  constant <- f$source
  terra::values(constant) <- .5
  projected <- kernel_scale.raster(constant, multiScaleR = inputs,
    match_buffer = TRUE, verbose = FALSE)
  expect_equal(terra::values(projected)[, "mean"], rep(.5, terra::ncell(constant)),
               tolerance = 1e-12)
  expect_equal(terra::values(projected)[, "rough"], rep(0, terra::ncell(constant)),
               tolerance = 0)
})

test_that("matching preserves supported fixed radii and rejects truncated larger ones", {
  f <- projection_buffer_fixture()
  specs <- msr_vars(mean = kernel_var("cont1"),
                    rough = surface_var("cont1", "sq", radius = 100))
  inputs <- kernel_prep(f$pts, f$source, max_D = 150, sigma = 80,
    kernel = "gaussian", scale_vars = specs, bin = FALSE, verbose = FALSE)
  projected <- kernel_scale.raster(f$source, multiScaleR = inputs,
    match_buffer = TRUE, verbose = FALSE)
  expected <- sweep(sweep(as.matrix(inputs$kernel_dat), 2L,
    inputs$scl_params$sd, "*"), 2L, inputs$scl_params$mean, "+")
  expect_equal(as.matrix(terra::extract(projected, f$pts)[, -1L]), expected,
               tolerance = 1e-10, ignore_attr = TRUE)
  larger <- msr_vars(mean = kernel_var("cont1"),
                     rough = surface_var("cont1", "sq", radius = 250))
  truncated <- kernel_prep(f$pts, f$source, max_D = 150, sigma = 80,
    kernel = "gaussian", scale_vars = larger, bin = FALSE, verbose = FALSE)
  expect_error(kernel_scale.raster(f$source, multiScaleR = truncated,
    match_buffer = TRUE, verbose = FALSE), "radius exceeds the prepared buffer")
})

test_that("prepared matching refuses incompatible modes and leaves defaults intact", {
  f <- projection_buffer_fixture()
  inputs <- kernel_prep(f$pts, f$source, max_D = 150, sigma = 80,
    kernel = "gaussian", bin = FALSE, verbose = FALSE)
  expect_error(kernel_scale.raster(f$source, sigma = 80, match_buffer = TRUE,
    verbose = FALSE), "requires a fitted or prepared")
  expect_error(kernel_scale.raster(f$source, multiScaleR = inputs,
    match_buffer = TRUE, fft = FALSE, verbose = FALSE), "fft = TRUE")
  expect_error(kernel_scale.raster(f$source, multiScaleR = inputs,
    match_buffer = TRUE, na.rm = FALSE, verbose = FALSE), "na.rm = TRUE")
  incompatible <- inputs
  incompatible$resolution <- incompatible$resolution * 2
  expect_error(kernel_scale.raster(f$source, multiScaleR = incompatible,
    match_buffer = TRUE, verbose = FALSE), "prepared resolution")
  unknown_crs <- f$source
  terra::crs(unknown_crs) <- ""
  expect_error(kernel_scale.raster(unknown_crs, multiScaleR = inputs,
    match_buffer = TRUE, verbose = FALSE), "projected square cells")
  legacy <- kernel_scale.raster(f$source, multiScaleR = inputs, verbose = FALSE)
  explicit <- kernel_scale.raster(f$source, multiScaleR = inputs,
    match_buffer = FALSE, verbose = FALSE)
  expect_identical(terra::values(explicit), terra::values(legacy))
})

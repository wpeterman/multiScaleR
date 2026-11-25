# Regression tests for raster refusals that name the failing condition
# (audit findings X-1, M-8, and B-4 for kernel_scale.raster()).

# A small sim_rast() landscape. Like the shipped example rasters, it has an
# empty CRS, which is the usual cause of the match_buffer refusal.
audit_raster_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    rs <- terra::subset(sim_rast(dim = 40, resolution = 10, user_seed = 11), "cont1")
    cells <- c(410L, 620L, 830L, 1015L, 1225L)
    pts <- terra::vect(terra::xyFromCell(rs, cells))
    inputs <- suppressWarnings(kernel_prep(pts, rs, max_D = 60, sigma = 30,
      kernel = "gaussian", bin = FALSE, verbose = FALSE))
    cache <<- list(rs = rs, pts = pts, inputs = inputs)
    cache
  }
})

test_that("match_buffer names an empty CRS and how to assign one", {
  f <- audit_raster_fixture()
  expect_false(nzchar(terra::crs(f$rs)))
  err <- tryCatch(kernel_scale.raster(f$rs, multiScaleR = f$inputs,
    match_buffer = TRUE, verbose = FALSE), error = conditionMessage)
  expect_match(err, "`raster_stack` has none", fixed = TRUE)
  expect_match(err, "terra::crs(r) <- \"EPSG:", fixed = TRUE)
  expect_match(err, "\"local\"", fixed = TRUE)
  expect_no_match(err, "longitude|resolution")

  # Both suggested assignments let the same call run.
  for (definition in c("EPSG:32617", "local")) {
    labeled <- f$rs
    terra::crs(labeled) <- definition
    expect_s4_class(kernel_scale.raster(labeled, multiScaleR = f$inputs,
      match_buffer = TRUE, verbose = FALSE), "SpatRaster")
  }
})

test_that("match_buffer names longitude/latitude and resolution mismatches", {
  f <- audit_raster_fixture()
  lonlat <- terra::rast(nrows = 40, ncols = 40, xmin = -83, xmax = -82.6,
    ymin = 40, ymax = 40.4, crs = "EPSG:4326", vals = 1)
  expect_error(kernel_scale.raster(lonlat, multiScaleR = f$inputs,
    match_buffer = TRUE, verbose = FALSE), "uses longitude/latitude")

  coarse <- terra::aggregate(f$rs, 2)
  terra::crs(coarse) <- "EPSG:32617"
  err <- tryCatch(kernel_scale.raster(coarse, multiScaleR = f$inputs,
    match_buffer = TRUE, verbose = FALSE), error = conditionMessage)
  expect_match(err, "prepared resolution (10 map units)", fixed = TRUE)
  expect_match(err, "has 20 x 20 cells", fixed = TRUE)
  expect_no_match(err, "has none|longitude")
})

test_that("fft = FALSE names an oversized kernel window instead of the terra error", {
  f <- audit_raster_fixture()
  small <- f$rs[1:10, 1:10, drop = FALSE]
  # A large sigma on a 10 x 10 cell raster gives a window more than twice its size.
  err <- tryCatch(kernel_scale.raster(small, sigma = 300, fft = FALSE,
    verbose = FALSE), error = conditionMessage)
  expect_match(err, "more than twice the size of the raster (10 x 10 cells)", fixed = TRUE)
  expect_match(err, "fft = TRUE", fixed = TRUE)
  expect_no_match(err, "nrow(w)", fixed = TRUE)
  # The suggested route works for the same call.
  expect_s4_class(kernel_scale.raster(small, sigma = 300, fft = TRUE,
    verbose = FALSE), "SpatRaster")
  # A window within terra's limit still uses focal and agrees with the FFT.
  direct <- kernel_scale.raster(small, sigma = 10, fft = FALSE, verbose = FALSE)
  expect_s4_class(direct, "SpatRaster")
})

test_that("the metric focal path checks the window before calling terra", {
  x <- matrix(1, 5, 5)
  expect_error(.msr_focal_convolution(x, matrix(1, 13, 13)),
    "kernel window \\(13 x 13 cells\\) is more than twice")
  expect_error(.msr_focal_convolution(x, matrix(1, 3, 13)), "fft = TRUE")
  # Windows at the limit keep the existing result.
  expect_equal(.msr_focal_convolution(x, matrix(1, 3, 3))[3, 3], 9)
  expect_true(.msr_check_focal_window(matrix(1, 10, 10), x))
})

test_that("kernel_scale.raster points a Bayesian fit to the Bayesian tools", {
  f <- audit_raster_fixture()
  bayes <- structure(list(), class = "multiScaleR_bayes")
  err <- tryCatch(kernel_scale.raster(f$rs, multiScaleR = bayes, verbose = FALSE),
    error = conditionMessage)
  expect_match(err, "`kernel_scale.raster()` works on fits from `multiScale_optim()`",
    fixed = TRUE)
  expect_match(err, "posterior_scales()", fixed = TRUE)
  expect_match(err, "kernel_scale.raster(sigma = )", fixed = TRUE)
  expect_no_match(err, "must be a `multiScaleR` or", fixed = TRUE)
  # Other inputs keep the class refusal.
  expect_error(kernel_scale.raster(f$rs, multiScaleR = list(), verbose = FALSE),
    "must be a `multiScaleR` or `multiScaleR_data` object")
})

test_that("kernel means retain a constant surface at raster edges in either engine", {
  r <- terra::rast(nrows = 18, ncols = 23, xmin = 0, xmax = 230,
                   ymin = 0, ymax = 180, crs = "EPSG:3857")
  r[] <- 10
  names(r) <- "hab"
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    for (fft in c(TRUE, FALSE)) {
      result <- kernel_scale.raster(r, sigma = 12, kernel = kernel,
        shape = if (kernel == "expow") 1.5 else NULL, fft = fft, verbose = FALSE)
      expect_equal(as.numeric(terra::values(result)), rep(10, terra::ncell(r)),
        tolerance = 1e-10, info = paste(kernel, fft))
    }
  }
})

test_that("kernel mean source masks and missing-support controls agree across engines", {
  r <- terra::rast(nrows = 18, ncols = 23, xmin = 0, xmax = 230,
                   ymin = 0, ymax = 180, crs = "EPSG:3857")
  terra::values(r) <- seq_len(terra::ncell(r)) %% 13
  hole <- terra::cellFromRowCol(r, 8, 10)
  r[hole] <- NA_real_
  names(r) <- "hab"
  for (remove in c(TRUE, FALSE)) {
    a <- kernel_scale.raster(r, sigma = 12, fft = TRUE, na.rm = remove, verbose = FALSE)
    b <- kernel_scale.raster(r, sigma = 12, fft = FALSE, na.rm = remove, verbose = FALSE)
    expect_equal(terra::values(a), terra::values(b), tolerance = 1e-10)
    expect_true(is.na(terra::values(b)[hole, 1]))
    if (!remove) expect_true(is.na(terra::values(a)[1, 1]))
  }
})

test_that("surface and categorical na.rm FALSE retain missing neighborhood failures", {
  r <- terra::rast(nrows = 18, ncols = 23, xmin = 0, xmax = 230,
                   ymin = 0, ymax = 180, crs = "EPSG:3857")
  terra::values(r) <- seq_len(terra::ncell(r)) %% 13
  hole <- terra::cellFromRowCol(r, 8, 10)
  near <- terra::cellFromRowCol(r, 8, 11)
  far <- terra::cellFromRowCol(r, 13, 16)
  r[hole] <- NA_real_
  names(r) <- "hab"
  specs <- list(msr_vars(metric = surface_var("hab", "sq", radius = 15)),
    msr_vars(metric = surface_var("hab", "sa", radius = 15)),
    msr_vars(metric = landscape_var("hab", "shdi", radius = 15)))
  for (spec in specs) for (fft in c(TRUE, FALSE)) {
    value <- terra::values(kernel_scale.raster(r, scale_vars = spec, fft = fft,
      na.rm = FALSE, verbose = FALSE))[, 1]
    expect_true(is.na(value[hole]))
    expect_true(is.na(value[near]))
    expect_true(is.na(value[1]))
    expect_true(is.finite(value[far]))
  }
})

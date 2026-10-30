test_that("direct metric convolution preserves orientation, masks, and zero support", {
  x <- matrix(seq_len(35), 5, 7)
  x[3, 2] <- NA_real_
  k <- matrix(0, 3, 3)
  k[1, 2] <- 1
  k[2, 3] <- 2
  for (stencil in list(k, k * 0, matrix(1, 3, 5))) {
    expect_equal(.msr_focal_convolution(x, stencil),
      fft_convolution(x, stencil, fun = "sum"), tolerance = 1e-11)
  }
})

test_that("weighted surface raster metrics support focal smoothing under all kernels", {
  for (kernel in c("gaussian", "exp", "expow")) {
    f <- refit_spatial_fixture(kernel)
    raster <- f$rasters[["elev"]]
    raster[c(20, 44, 136)] <- NA_real_
    for (metric in c("sa", "sq", "ssk", "sku")) {
      specs <- msr_vars(texture = surface_var("elev", metric, weighted = TRUE))
      args <- list(raster_stack = raster, sigma = 29,
        shape = if (kernel == "expow") 1.7 else NULL,
        kernel = kernel, scale_vars = specs, verbose = FALSE)
      a <- do.call(kernel_scale.raster, c(args, list(fft = TRUE)))
      b <- do.call(kernel_scale.raster, c(args, list(fft = FALSE)))
      expect_equal(terra::values(a), terra::values(b), tolerance = 2e-8)
      expect_identical(names(b), "texture")
    }
  }
})

test_that("hard-radius texture and categorical metrics support focal smoothing", {
  f <- refit_spatial_fixture()
  for (metric in c("sa", "sq", "ssk", "sku", "sdq", "sdr")) {
    specs <- msr_vars(texture = surface_var("elev", metric, radius = 31))
    args <- list(raster_stack = f$rasters[["elev"]], scale_vars = specs, verbose = FALSE)
    a <- do.call(kernel_scale.raster, c(args, list(fft = TRUE)))
    b <- do.call(kernel_scale.raster, c(args, list(fft = FALSE)))
    expect_equal(terra::values(a), terra::values(b), tolerance = 1e-8)
  }
  r <- terra::rast(nrows = 22, ncols = 26, xmin = 0, xmax = 260,
                   ymin = 0, ymax = 220, crs = "EPSG:3857")
  ij <- terra::rowColFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- 1 + (ij[, 1] + ij[, 2] %/% 2) %% 4
  r[c(23, 99, 139)] <- NA_real_
  names(r) <- "cover"
  metrics <- c(.msr_composition_metrics(), "pland", "ca", .msr_edge_metrics(),
               .msr_adjacency_metrics())
  for (metric in metrics) {
    spec <- landscape_var("cover", metric, radius = 31,
      classes_max = if (metric == "rpr") 4 else NULL,
      class = if (metric %in% c("pland", "ca", "clumpy")) 1 else NULL)
    args <- list(raster_stack = r, scale_vars = msr_vars(land = spec), verbose = FALSE)
    a <- do.call(kernel_scale.raster, c(args, list(fft = TRUE)))
    b <- do.call(kernel_scale.raster, c(args, list(fft = FALSE)))
    expect_equal(terra::values(a), terra::values(b), tolerance = 1e-9,
                 info = metric)
  }
})

test_that("focal projection can retain the exact prepared buffer and aliases", {
  for (kernel in c("gaussian", "exp", "expow")) {
    f <- refit_spatial_fixture(kernel)
    x <- f$inputs
    x$sigma <- c(19, 29, 17) / x$unit_conv
    if (kernel == "expow") x$shape <- c(1.3, 1.7, 1.2)
    a <- kernel_scale.raster(f$rasters, multiScaleR = x, match_buffer = TRUE,
                             fft = TRUE, verbose = FALSE)
    b <- kernel_scale.raster(f$rasters, multiScaleR = x, match_buffer = TRUE,
                             fft = FALSE, verbose = FALSE)
    expect_equal(terra::values(a), terra::values(b), tolerance = 2e-8)
    expect_identical(names(b), c("mean", "rough", "fixed", "unused"))
  }
})

test_that("roughness is stable under large offsets and zero on constant support", {
  v <- c(1, 2, 4, 8, 16)
  w <- exp(-seq_along(v) / 2)
  for (metric in c("sa", "sq", "ssk", "sku"))
    expect_equal(.surface_weighted_metric(v + 2^40, w, metric),
                 .surface_weighted_metric(v, w, metric), tolerance = 1e-13)
  r <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 160,
                   ymin = 0, ymax = 160, crs = "EPSG:3857")
  terra::values(r) <- rep(seq_len(16), 16)
  names(r) <- "elev"
  shifted <- r + 2^40
  for (metric in c("sa", "sq")) {
    for (weighted in c(FALSE, TRUE)) {
      specs <- msr_vars(texture = surface_var("elev", metric,
        radius = if (!weighted) 31 else NULL, weighted = weighted))
      args <- list(scale_vars = specs, sigma = if (weighted) 19 else NULL,
                    verbose = FALSE)
      a <- do.call(kernel_scale.raster, c(list(raster_stack = r), args))
      b <- do.call(kernel_scale.raster, c(list(raster_stack = shifted), args))
      expect_equal(terra::values(a), terra::values(b), tolerance = 1e-10)
    }
  }
  r[] <- 7
  specs <- msr_vars(texture = surface_var("elev", "sq", weighted = TRUE))
  for (fft in c(TRUE, FALSE)) {
    result <- kernel_scale.raster(r, sigma = 19, scale_vars = specs,
                                  fft = fft, verbose = FALSE)
    expect_true(all(terra::values(result) == 0))
  }
})

test_that("deterministic kernel smoothing preserves the caller RNG state", {
  f <- refit_spatial_fixture()
  set.seed(773)
  before <- .Random.seed
  for (fft in c(TRUE, FALSE)) {
    kernel_scale.raster(f$rasters[["hab"]], sigma = 19, fft = fft, verbose = FALSE)
    expect_identical(.Random.seed, before)
    kernel_scale.raster(f$rasters[["hab"]], sigma = 19, fft = fft,
      scale_vars = msr_vars(alias = kernel_var("hab")), verbose = FALSE)
    expect_identical(.Random.seed, before)
  }
})

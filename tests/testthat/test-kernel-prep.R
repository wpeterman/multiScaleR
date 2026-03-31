test_that("kernel_prep returns a complete multiScaleR_data object", {
  pts <- terra::vect(cbind(c(5, 7, 9, 11, 13), c(13, 11, 9, 7, 5)))
  r <- terra::rast(matrix(seq_len(20 * 20), nrow = 20))
  names(r) <- "hab"

  out <- kernel_prep(
    pts = pts,
    raster_stack = r,
    max_D = 5,
    kernel = "gaussian",
    verbose = FALSE
  )

  expect_s3_class(out, "multiScaleR_data")
  expect_true(all(c("kernel_dat", "d_list", "raw_cov", "kernel", "shape",
                    "min_D", "max_D", "n_covs", "unit_conv", "sigma", "scl_params") %in% names(out)))
  expect_equal(nrow(out$kernel_dat), length(pts))
  expect_equal(colnames(out$kernel_dat), "hab")
  expect_equal(colnames(out$raw_cov[[1]]), "hab")
  expect_equal(out$kernel, "gaussian")
  expect_equal(out$n_covs, 1)
})

test_that("kernel_prep creates default shape values for expow kernels", {
  pts <- terra::vect(cbind(c(5, 7, 9), c(9, 7, 5)))
  r <- terra::rast(list(a = terra::rast(matrix(runif(20 * 20), 20, 20)),
                        b = terra::rast(matrix(runif(20 * 20), 20, 20))))

  out <- kernel_prep(
    pts = pts,
    raster_stack = r,
    max_D = 5,
    kernel = "expow",
    verbose = FALSE
  )

  expect_equal(out$shape, c(2, 2))
})

test_that("kernel_prep covers sf input plus progress and verbose branches", {
  sim_fix <- make_simulation_fixture()
  pts_sf <- sim_fix$sim$pts

  expect_output(
    out <- kernel_prep(
      pts = pts_sf,
      raster_stack = sim_fix$rs_two,
      max_D = 250,
      kernel = "expow",
      progress = TRUE,
      verbose = TRUE
    ),
    "Calculating weights"
  )

  expect_s3_class(out, "multiScaleR_data")
  expect_equal(length(out$d_list), nrow(pts_sf))
})

test_that("kernel_prep validates key inputs", {
  pts <- terra::vect(cbind(c(5, 7, 9), c(9, 7, 5)))
  r <- terra::rast(matrix(runif(20 * 20), 20, 20))
  names(r) <- "hab"

  expect_error(
    kernel_prep(pts = pts, raster_stack = matrix(1), max_D = 5, verbose = FALSE),
    "SpatRaster"
  )

  expect_error(
    kernel_prep(pts = pts, raster_stack = r, max_D = 5, projected = FALSE, verbose = FALSE),
    "must be projected"
  )

  expect_error(
    kernel_prep(pts = pts, raster_stack = r, max_D = 5, sigma = c(1, 2), verbose = FALSE),
    "Number of sigma values"
  )

  expect_error(
    kernel_prep(pts = pts, raster_stack = r, max_D = 0, verbose = FALSE),
    "max_D"
  )
  expect_error(
    kernel_prep(pts = pts, raster_stack = r, max_D = 5, sigma = -1, verbose = FALSE),
    "sigma"
  )
  expect_error(
    kernel_prep(pts = pts, raster_stack = r, max_D = 5, progress = NA, verbose = FALSE),
    "progress"
  )
})

test_that("print method for multiScaleR_data reports object metadata", {
  fix <- make_core_fixture()

  expect_output(print(fix$kernel_inputs), "spatial covariate")
  expect_output(print(fix$kernel_inputs), "Maximum Distance")
})

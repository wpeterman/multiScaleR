# sim_dat() input handling found in the 0.7.47 help-page review: an sf point
# object, the exponential power kernel without max_D, and rasters whose
# extent does not start at the origin.

sim_dat_raster <- function(xmin = 0, ymin = 0, n = 60, res = 10, crs = "local") {
  set.seed(11)
  r <- terra::rast(nrows = n, ncols = n, xmin = xmin, xmax = xmin + n * res,
                   ymin = ymin, ymax = ymin + n * res, crs = crs)
  terra::values(r) <- stats::rnorm(terra::ncell(r))
  names(r) <- "cov"
  r
}

test_that("an sf point object gives one observation per point", {
  r <- sim_dat_raster(crs = "")
  xy <- cbind(x = seq(200, 400, length.out = 10), y = seq(200, 400, length.out = 10))
  pts <- sf::st_as_sf(data.frame(xy, site = letters[1:10], z = 1:10),
                      coords = c("x", "y"))
  out <- sim_dat(beta = 0.5, sigma = 20, raster_stack = r, n_points = pts,
                 max_D = 80, user_seed = 3)
  expect_length(out$obs, 10L)
  expect_identical(nrow(out$df), 10L)
  expect_identical(nrow(out$pts), 10L)
})

test_that("expow without max_D sets max_D from each layer's sigma and shape", {
  r <- sim_dat_raster()
  out <- sim_dat(beta = 0.5, kernel = "expow", sigma = 20, shape = 3,
                 raster_stack = r, n_points = 12, user_seed = 4)
  expect_length(out$obs, 12L)
  expect_true(all(is.finite(out$df$cov)))
})

test_that("rasters away from the origin place points inside their extent", {
  r <- sim_dat_raster(xmin = 5000, ymin = 9000)
  out <- sim_dat(beta = 0.5, sigma = 20, raster_stack = r, n_points = 12,
                 max_D = 80, user_seed = 5)
  xy <- sf::st_coordinates(out$pts)
  e <- as.vector(terra::ext(r))
  expect_true(all(xy[, 1] > e[1] & xy[, 1] < e[2]))
  expect_true(all(xy[, 2] > e[3] & xy[, 2] < e[4]))
  expect_length(out$obs, 12L)
})

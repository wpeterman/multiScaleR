# kernel_prep() sets min_D, the lower bound of the sigma search, to one cell
# width. It used to take the distance between the first two cells extracted
# for the first point, which spans rows when the top row of that buffer holds
# a single cell.

min_d_raster <- function(resolution) {
  r <- terra::rast(nrows = 100, ncols = 100, xmin = 0, ymin = 0,
                   xmax = 100 * resolution, ymax = 100 * resolution,
                   crs = "EPSG:32617")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- sin(xy[, 1] / (7 * resolution)) + cos(xy[, 2] / (11 * resolution))
  names(r) <- "hab"
  r
}

test_that("min_D is one cell when the first buffer's top row holds one cell", {
  r <- min_d_raster(10)
  # The first buffer reaches 0.05 map units into the row above y = 600 and
  # touches only the cell centered at x = 505 there.
  pts <- sf::st_as_sf(data.frame(x = c(505, 300, 700), y = c(500, 400, 650)),
                      coords = c("x", "y"), crs = 32617)
  inputs <- kernel_prep(pts, r, max_D = 100.05, verbose = FALSE)
  expect_equal(inputs$min_D, 10)
})

test_that("min_D follows the cell size and keeps its earlier value elsewhere", {
  r <- min_d_raster(20)
  pts <- sf::st_as_sf(data.frame(x = c(1010, 600, 1400), y = c(1000, 800, 1300)),
                      coords = c("x", "y"), crs = 32617)
  inputs <- kernel_prep(pts, r, max_D = 300, verbose = FALSE)
  expect_equal(inputs$min_D, 20)
  # Same value as kernel_prep_by_group(), which already used one cell width.
  expect_equal(inputs$min_D, floor(terra::res(r)[1]))
})

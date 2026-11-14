# sim_rast() fields found in the 0.7.47 help-page review: the covariance was
# not circulant-embedded, so half its eigenvalues were negative and clamped,
# and small grids rounded the range to zero cells and returned NaN layers.

lag_correlation <- function(m, h) {
  a <- m[, seq_len(ncol(m) - h)]; b <- m[, h + seq_len(ncol(m) - h)]
  stats::cor(as.vector(a), as.vector(b))
}

test_that("simulated fields follow the stated exponential correlation", {
  range <- 5
  est <- vapply(1:6, function(seed) {
    r <- multiScaleR:::simulate_fft_grf(dim = 120, resolution = 1, range = range,
                                         user_seed = seed)
    m <- terra::as.matrix(r, wide = TRUE)
    c(lag_correlation(m, 1), lag_correlation(m, 5))
  }, numeric(2))
  expect_equal(rowMeans(est), exp(-c(1, 5) / range), tolerance = 0.12)
})

test_that("simulated fields do not repeat across half the grid", {
  # The old construction kept only frequencies whose indices summed to an even
  # number, so shifting the field by half the grid in both directions gave
  # the same field: the top-left and bottom-right quadrants correlated near 1.
  q <- vapply(1:4, function(seed) {
    m <- terra::as.matrix(multiScaleR:::simulate_fft_grf(dim = 100, resolution = 1,
      range = 5, user_seed = seed), wide = TRUE)
    stats::cor(as.vector(m[1:50, 1:50]), as.vector(m[51:100, 51:100]))
  }, numeric(1))
  expect_true(all(abs(q) < 0.3))
})

test_that("small grids give finite layers", {
  r <- sim_rast(dim = 15, resolution = 5, user_seed = 9)
  expect_false(anyNA(terra::values(r)))
  expect_identical(names(r), c("bin1", "bin2", "cont1", "cont2"))
  expect_error(multiScaleR:::simulate_fft_grf(dim = 10, range = 0), "positive")
})

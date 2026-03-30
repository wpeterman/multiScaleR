test_that("kernel distance helpers match expected formulas", {
  expect_equal(
    kernel_dist(sigma = 10, kernel = "gaussian", prob = 0.95),
    round(stats::qnorm(0.975, mean = 0, sd = 10), 2),
    tolerance = 1e-8
  )

  expect_equal(
    kernel_dist(sigma = 10, kernel = "exp", prob = 0.95),
    round(-10 * log(1 - 0.95), 2),
    tolerance = 1e-8
  )

  expect_equal(
    kernel_dist(sigma = 10, kernel = "fixed", prob = 0.5),
    5,
    tolerance = 1e-8
  )

  expect_true(kernel_dist(sigma = 10, kernel = "expow", beta = 2, prob = 0.9) > 0)
})

test_that("kernel distance helpers validate arguments", {
  expect_error(kernel_dist(sigma = 10, kernel = "expow", prob = 0.9), "shape")
  expect_error(kernel_dist(sigma = 10, kernel = "gaussian", prob = 1.2), "between 0 and 1")
  expect_error(plot_kernel(sigma = 10, kernel = "expow"), "beta")
})

test_that("sparse and R kernel weighting implementations agree", {
  d <- c(0, 1, 2)
  dense_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
  sparse_mat <- Matrix::Matrix(dense_mat, sparse = TRUE)

  r_out <- multiScaleR:::scale_type_r(
    d = d,
    kernel = "gaussian",
    sigma = c(1, 2),
    r_stack.df = dense_mat
  )
  sparse_out <- multiScaleR:::scale_type(
    d = d,
    kernel = "gaussian",
    sigma = c(1, 2),
    r_stack.df = sparse_mat
  )

  expect_equal(as.numeric(sparse_out), as.numeric(r_out), tolerance = 1e-8)

  r_wts <- multiScaleR:::scale_type_r(
    d = d,
    kernel = "gaussian",
    sigma = 2,
    r_stack.df = dense_mat[, 1, drop = FALSE],
    output = "wts"
  )
  sparse_wts <- multiScaleR:::scale_type_sparse(
    d = d,
    kernel = "gaussian",
    sigma_ = 2,
    r_stack_df = Matrix::Matrix(dense_mat[, 1, drop = FALSE], sparse = TRUE),
    output = "wts"
  )

  expect_equal(drop(r_wts), drop(sparse_wts), tolerance = 1e-8)
  expect_equal(sum(r_wts), 1, tolerance = 1e-8)
})

test_that("ci function clamps lower bounds and preserves output schema", {
  x <- matrix(c(5, 10, 10, 1), ncol = 2, byrow = TRUE)
  out <- multiScaleR:::ci_func(x, df = 30, min_D = 2, names = c("a", "b"))

  expect_equal(colnames(out), c("Mean", "SE", "2.5%", "97.5%"))
  expect_equal(rownames(out), c("a", "b"))
  expect_true(all(out[, "2.5%"] >= 2, na.rm = TRUE))
  expect_true(all(out[, "97.5%"] >= 2, na.rm = TRUE))
})

test_that("fft convolution preserves constants and missing cells", {
  x <- matrix(1, 5, 5)
  x[3, 3] <- NA_real_
  kernel <- matrix(1, 3, 3)

  out <- multiScaleR:::fft_convolution(x, kernel, fun = "mean", na.rm = TRUE)

  expect_equal(dim(out), dim(x))
  expect_true(is.na(out[3, 3]))
  expect_equal(out[2, 2], 1, tolerance = 1e-8)
  expect_equal(out[4, 4], 1, tolerance = 1e-8)
})

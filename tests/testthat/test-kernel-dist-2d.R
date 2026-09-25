test_that("two-dimensional kernel distances match their closed forms", {
  # Gaussian: radius holding 90% of the planar weight is sqrt(-2 log 0.1) sigma
  expect_equal(kernel_dist(sigma = 300, kernel = "gaussian", dimension = "2d"),
               round(300 * sqrt(-2 * log(0.1)), 2))
  # Negative exponential: (1 + r / sigma) exp(-r / sigma) = 1 - prob at the radius
  r_exp <- multiScaleR:::k_dist(100, prob = 0.9, kernel = "exp", dimension = "2d")
  expect_equal((1 + r_exp / 100) * exp(-r_exp / 100), 0.1, tolerance = 1e-8)
  # Exponential power with beta = 2 is a Gaussian with sigma / sqrt(2);
  # with beta = 1 it is the negative exponential
  expect_equal(
    multiScaleR:::k_dist(100, prob = 0.9, kernel = "expow", beta = 2, dimension = "2d"),
    multiScaleR:::k_dist(100 / sqrt(2), prob = 0.9, kernel = "gaussian", dimension = "2d")
  )
  expect_equal(
    multiScaleR:::k_dist(100, prob = 0.9, kernel = "expow", beta = 1, dimension = "2d"),
    r_exp
  )
  # Fixed buffer: circle holding `prob` of the buffer's area
  expect_equal(kernel_dist(sigma = 100, kernel = "fixed", prob = 0.81, dimension = "2d"), 90)
  # The default stays one-dimensional
  expect_equal(kernel_dist(sigma = 300, kernel = "gaussian"),
               round(qnorm(0.95) * 300, 2))
  expect_error(kernel_dist(sigma = 300, kernel = "gaussian", dimension = "3d"))
})

test_that("the 2D radius holds the stated share of weight on a raster grid", {
  # Weight cells on a fine grid with the package's own kernel weights and check
  # the share of total weight inside the 2D radius
  sigma <- 50
  xy <- expand.grid(x = seq(-400, 400, by = 2), y = seq(-400, 400, by = 2))
  d <- sqrt(xy$x^2 + xy$y^2)
  for (kern in c("gaussian", "exp", "expow")) {
    shape <- if (kern == "expow") 1.5 else NULL
    w <- switch(kern,
                gaussian = exp(-d^2 / (2 * sigma^2)),
                exp      = exp(-d / sigma),
                expow    = exp(-(d / sigma)^shape))
    r2 <- multiScaleR:::k_dist(sigma, prob = 0.9, kernel = kern, beta = shape,
                               dimension = "2d")
    r1 <- multiScaleR:::k_dist(sigma, prob = 0.9, kernel = kern, beta = shape)
    expect_equal(sum(w[d <= r2]) / sum(w), 0.9, tolerance = 0.01)
    expect_lt(sum(w[d <= r1]) / sum(w), 0.85)
  }
})

test_that("summary and print report both distances", {
  fix <- make_core_fixture()
  sum_opt <- summary(fix$opt)

  expect_equal(sum_opt$opt_dist, kernel_dist(fix$opt))
  expect_equal(sum_opt$opt_dist_2d, kernel_dist(fix$opt, dimension = "2d"))
  sum_95 <- summary(fix$opt, prob = 0.95)
  expect_equal(sum_95$opt_dist_2d,
               kernel_dist(fix$opt, prob = 0.95, dimension = "2d"))
  expect_equal(colnames(sum_opt$opt_dist_2d), c("Mean", "2.5%", "97.5%"))
  if (identical(fix$opt$kernel_inputs$kernel, "gaussian")) {
    expect_equal(sum_opt$opt_dist_2d$Mean / sum_opt$opt_dist$Mean,
                 rep(sqrt(-2 * log(0.1)) / qnorm(0.95), nrow(sum_opt$opt_dist)),
                 tolerance = 1e-3)
  }

  printed <- capture.output(print(sum_opt))
  expect_true(any(grepl("1D, along a line", printed, fixed = TRUE)))
  expect_true(any(grepl("2D, radius of the circle", printed, fixed = TRUE)))
  printed_fit <- capture.output(print(fix$opt))
  expect_true(any(grepl("2D, radius of the circle", printed_fit, fixed = TRUE)))

  # Summary objects saved before `opt_dist_2d` existed still print
  old <- sum_opt
  old$opt_dist_2d <- NULL
  printed_old <- capture.output(print(old))
  expect_true(any(grepl("1D: 90% kernel weight along a line", printed_old,
                        fixed = TRUE)))
})

test_that("hard-radius covariates report the radius in both distances", {
  fix <- make_core_fixture()
  opt <- fix$opt
  opt$kernel_inputs$scale_vars$type[] <- "landscape"
  sum_opt <- summary(opt)
  expect_equal(sum_opt$opt_dist_2d[1, "Mean"], round(opt$scale_est$Mean, 2))
  expect_equal(sum_opt$opt_dist_2d, sum_opt$opt_dist)
  expect_true(all(is.na(sum_opt$opt_dist_2d[, c("2.5%", "97.5%")])))
})

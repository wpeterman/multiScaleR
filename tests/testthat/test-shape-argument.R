# `shape` replaces `beta` as the exponential-power shape argument of
# kernel_dist() and plot_kernel(). `beta` stays as a deprecated alias for one
# release.

marker <- function(p) {
  layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomVline"), logical(1)))
  p$layers[[layer]]$data$xintercept
}

test_that("kernel_dist() accepts shape and matches the deprecated beta call", {
  for (dimension in c("2d", "1d")) {
    new <- expect_no_warning(
      kernel_dist(sigma = 100, kernel = "expow", shape = 1.5, prob = 0.95,
                  dimension = dimension)
    )
    expect_warning(
      old <- kernel_dist(sigma = 100, kernel = "expow", beta = 1.5, prob = 0.95,
                         dimension = dimension),
      "`beta` is deprecated in favor of `shape`"
    )
    expect_identical(new, old)
    # Both agree with the internal closed form.
    expect_equal(new, round(multiScaleR:::k_dist(100, 0.95, "expow", beta = 1.5,
                                                  dimension = dimension), 2))
  }
})

test_that("kernel_dist() warns once per call for beta and rejects conflicts", {
  n_warn <- 0
  withCallingHandlers(
    kernel_dist(sigma = 100, kernel = "expow", beta = 2),
    warning = function(w) {
      n_warn <<- n_warn + 1
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(n_warn, 1)

  expect_error(kernel_dist(sigma = 100, kernel = "expow", shape = 2, beta = 3),
               "both supplied with different values")
  # Equal values are not a conflict, but `beta` still warns.
  expect_warning(
    same <- kernel_dist(sigma = 100, kernel = "expow", shape = 2, beta = 2),
    "deprecated"
  )
  expect_identical(same, kernel_dist(sigma = 100, kernel = "expow", shape = 2))
})

test_that("kernel_dist() defaults and messages are unchanged", {
  # Kernels without a shape need neither argument and do not warn.
  expect_no_warning(g <- kernel_dist(sigma = 300, kernel = "gaussian"))
  expect_equal(g, round(300 * sqrt(-2 * log(1 - 0.9)), 2))
  expect_error(kernel_dist(sigma = 10, kernel = "expow"), "shape")
  expect_error(kernel_dist(sigma = 10, kernel = "expow", shape = -1), "`shape` must be > 0")
  # A NULL alias (for example beta = NULL passed programmatically) is not a use.
  expect_no_warning(kernel_dist(sigma = 300, kernel = "gaussian", beta = NULL))
})

test_that("plot_kernel() accepts shape and matches the deprecated beta call", {
  new <- expect_no_warning(plot_kernel(sigma = 300, kernel = "expow", shape = 5))
  expect_warning(old <- plot_kernel(sigma = 300, kernel = "expow", beta = 5),
                 "`beta` is deprecated in favor of `shape`")
  expect_identical(new$data, old$data)
  expect_identical(marker(new), marker(old))
  expect_identical(new$labels$caption, old$labels$caption)
  expect_equal(marker(new), multiScaleR:::k_dist(300, 0.9, "expow", beta = 5))

  # The call used in an existing vignette chunk keeps working.
  expect_warning(
    vig <- plot_kernel(sigma = 300, kernel = "expow", add_label = FALSE, beta = 5),
    "deprecated"
  )
  expect_s3_class(vig, "ggplot")
  expect_identical(vig$data, new$data)
})

test_that("plot_kernel() keeps the positional slot and rejects conflicts", {
  # The third positional argument was `beta` and is now `shape`.
  pos <- expect_no_warning(plot_kernel(0.9, 300, 5, "expow"))
  expect_identical(pos$data, plot_kernel(sigma = 300, kernel = "expow", shape = 5)$data)

  expect_error(plot_kernel(sigma = 300, kernel = "expow", shape = 2, beta = 5),
               "both supplied with different values")
  expect_error(plot_kernel(sigma = 300, kernel = "expow"), "shape")
  expect_error(plot_kernel(sigma = 300, kernel = "expow", shape = 0), "`shape` must be > 0")

  # Kernels without a shape are unaffected.
  expect_no_warning(g <- plot_kernel(sigma = 300, kernel = "gaussian"))
  expect_equal(marker(g), 300 * sqrt(-2 * log(1 - 0.9)))
})

test_that("internal paths that reach kernel_dist() do not raise the beta warning", {
  # Internal callers pass `shape`, so the projection path stays silent.
  r <- terra::rast(matrix(runif(400), 20, 20))
  names(r) <- "x"
  expect_no_warning(
    kernel_scale.raster(raster_stack = r, sigma = 3, shape = 2, kernel = "expow",
                        fft = FALSE, verbose = FALSE)
  )
})

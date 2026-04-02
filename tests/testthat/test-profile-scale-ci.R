test_that("summary uses profile-based sigma intervals when available", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  set.seed(123)

  points <- terra::vect(cbind(c(5, 7, 9, 11, 13),
                              c(13, 11, 9, 7, 5)))
  r1 <- terra::rast(matrix(rnorm(20^2), nrow = 20))
  names(r1) <- "r1"

  kernel_inputs <- kernel_prep(pts = points,
                               raster_stack = r1,
                               max_D = 25,
                               kernel = "gaussian",
                               verbose = FALSE)

  dat <- data.frame(y = rnorm(5), kernel_inputs$kernel_dat)
  mod <- glm(y ~ r1, data = dat)

  opt <- multiScale_optim(fitted_mod = mod,
                          kernel_inputs = kernel_inputs,
                          verbose = FALSE)
  smry <- summary(opt)

  interval_method <- attr(smry$opt_scale, "interval_method")

  expect_identical(unname(interval_method), "profile")
  expect_true(all(is.finite(as.matrix(smry$opt_scale[, c("2.5%", "97.5%")]))))
  expect_true(all(is.finite(as.matrix(smry$opt_dist[, c("2.5%", "97.5%")]))))
})

test_that("summary uses Wald intervals by default and profile intervals on request", {
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

  cache_key <- profile_scale_cache_key(object = opt,
                                       min_D = opt$min_D,
                                       names = row.names(opt$scale_est))
  if (exists(cache_key, envir = .profile_scale_cache, inherits = FALSE)) {
    rm(list = cache_key, envir = .profile_scale_cache)
  }

  smry_wald <- summary(opt)
  interval_method_wald <- attr(smry_wald$opt_scale, "interval_method")

  expect_identical(unname(interval_method_wald), "wald")

  smry_profile <- summary(opt, profile = TRUE)
  interval_method_profile <- attr(smry_profile$opt_scale, "interval_method")

  expect_false(is.null(opt$kernel_inputs$d_list))
  expect_false(is.null(opt$kernel_inputs$raw_cov))
  expect_identical(unname(interval_method_profile), "profile")
  expect_true(all(is.finite(as.matrix(smry_profile$opt_scale[, c("2.5%", "97.5%")]))))
  expect_true(all(is.finite(as.matrix(smry_profile$opt_dist[, c("2.5%", "97.5%")]))))
  expect_true(exists(cache_key, envir = .profile_scale_cache, inherits = FALSE))

  smry_profile_cached <- summary(opt, profile = TRUE)
  expect_equal(smry_profile_cached$opt_scale, smry_profile$opt_scale)
  expect_equal(smry_profile_cached$opt_dist, smry_profile$opt_dist)
})

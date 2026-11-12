# profile_sigma() and summary(profile = TRUE) on the opt-in backends.
#
# profile_glm: the default log grid computed exp(log(min_D)) slightly below
# min_D, so the bounds check rejected the default call. conditional_logit:
# profile_sigma() took the default-refit route, which this backend has no
# context for, and failed after a stray `length.out` warning.

running_profile_glm <- function() {
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700,
                        kernel = "gaussian", verbose = FALSE)
  dat <- data.frame(landscape_counts, inputs$kernel_dat)
  mod <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
  multiScale_optim(mod, inputs, backend = "profile_glm", n_cores = 1, verbose = FALSE)
}

test_that("profile_sigma() defaults work on a profile_glm fit", {
  fit <- running_profile_glm()
  prof <- profile_sigma(fit, verbose = FALSE)
  expect_s3_class(prof, "sigma_profile")
  # The grid spans exactly the fitted bounds.
  expect_identical(range(prof$sigma_grid), c(fit$min_D, fit$max_D))
  expect_true(all(is.finite(prof$profiles$LL)))
  # Each profile peaks next to its estimate.
  for (v in rownames(fit$scale_est)) {
    p <- prof$profiles[prof$profiles$variable == v, ]
    peak <- p$sigma[which.max(p$LL)]
    neighbors <- range(p$sigma[abs(log(p$sigma / fit$scale_est[v, "Mean"])) <
                                 diff(log(prof$sigma_grid[1:2])) * 1.01])
    expect_true(peak >= neighbors[1] && peak <= neighbors[2])
  }

  # Linear spacing also keeps its endpoints.
  lin <- profile_sigma(fit, n_pts = 5, spacing = "linear", verbose = FALSE)
  expect_identical(range(lin$sigma_grid), c(fit$min_D, fit$max_D))

  # A user range outside the bounds keeps an informative error.
  expect_error(profile_sigma(fit, sigma_range = c(10, 600), verbose = FALSE),
               "within the fitted bounds, 20 to 1700")
  expect_error(profile_sigma(fit, sigma_values = c(100, 500, 2000), verbose = FALSE),
               "within the fitted bounds")

  # Profile intervals from summary() are finite and bracket the estimates.
  s <- summary(fit, profile = TRUE)$opt_scale
  expect_true(all(is.finite(as.matrix(s[, c("2.5%", "97.5%")]))))
  expect_true(all(s[, "2.5%"] < s$Mean & s$Mean < s[, "97.5%"]))
})

test_that("the default refit grid also keeps exact endpoints", {
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = surv_pts, raster_stack = land, max_D = 1700,
                        kernel = "gaussian", verbose = FALSE)
  dat <- data.frame(landscape_counts, inputs$kernel_dat)
  mod <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
  fit <- multiScale_optim(mod, inputs, n_cores = 1, verbose = FALSE)
  prof <- profile_sigma(fit, n_pts = 4, verbose = FALSE)
  expect_identical(range(prof$sigma_grid), c(fit$min_D, fit$max_D))
})

test_that("profile_sigma() profiles a conditional_logit fit with native refits", {
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  form <- case_ ~ forest + sl_km + strata(step_id)
  m <- survival::clogit(form, data = fx$df, model = TRUE)
  fit <- multiScale_optim(m, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE)
  sig <- fit$scale_est["forest", "Mean"]

  # No stray warning, and the estimate is the best point on a grid around it.
  expect_no_warning(
    prof <- profile_sigma(fit, sigma_values = sig * c(0.6, 0.8, 1, 1.25, 1.6),
                          verbose = FALSE)
  )
  ll <- prof$profiles$LL
  expect_equal(ll[3], as.numeric(logLik(fit$opt_mod)), tolerance = 1e-8)
  expect_true(all(ll[3] >= ll[-3]))

  # Grid values agree with independent clogit() refits on kernel values
  # prepared at the same sigma (up to distance binning).
  expect_equal(ll[c(1, 5)],
               c(ssf_native_loglik(fx, form, sig * 0.6),
                 ssf_native_loglik(fx, form, sig * 1.6)),
               tolerance = 1e-4)

  # The default grid runs and covers the fitted bounds.
  default <- profile_sigma(fit, n_pts = 4, verbose = FALSE)
  expect_identical(range(default$sigma_grid), c(fit$min_D, fit$max_D))
  expect_true(all(is.finite(default$profiles$LL)))

  # Parallel profiling and out-of-bound grids are refused clearly.
  expect_error(profile_sigma(fit, n_cores = 2, verbose = FALSE),
               "conditional_logit likelihood profiling is serial")
  expect_error(profile_sigma(fit, sigma_values = c(50, 100, 400), verbose = FALSE),
               "within the fitted bounds")
})

test_that("summary(profile = TRUE) on a conditional_logit fit stops with one clear message", {
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  m <- survival::clogit(case_ ~ forest + sl_km + strata(step_id),
                        data = fx$df, model = TRUE)
  fit <- multiScale_optim(m, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE)
  expect_error(summary(fit, profile = TRUE),
               "not available for conditional_logit fits.*profile_sigma\\(\\)")
  # The ordinary summary still works.
  expect_s3_class(summary(fit), "summary_multiScaleR")
})

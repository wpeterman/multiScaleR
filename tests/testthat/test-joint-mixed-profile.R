# Profile-likelihood intervals and sigma profiles for backend = "joint_mixed".
# The fixture fits reuse the frozen annual inputs of test-joint-mixed.R. The
# limits are checked against glmmTMB fits built here, independently of the
# package's own refit function.

skip_joint_mixed <- function() {
  skip_on_cran()
  skip_if_not_installed("RTMB")
  skip_if_not_installed("TMB")
  skip_if_not_installed("glmmTMB")
}

mixed_fixture <- function() {
  f <- readRDS(test_path("fixtures", "joint-mixed-existing.rds"))
  d <- f$data; d$site <- factor(d$site); d$year <- factor(d$year)
  e <- .msr_prepare_kernel_engine(f$inputs$binned, "habitat")
  list(f = f, d = d, e = e,
       form = count ~ temp + year + habitat + (1 | site) + offset(log(effort)))
}

test_that("joint_mixed profile limits sit on the native likelihood cutoff", {
  skip_joint_mixed()
  x <- mixed_fixture()
  d <- x$d
  d$habitat <- .msr_eval_kernel_engine(x$e, log(.23))$value[, 1]
  m <- glmmTMB::glmmTMB(x$form, data = d, family = glmmTMB::nbinom2())
  fit <- multiScale_optim(m, x$f$inputs, backend = "joint_mixed", n_cores = 1,
    profile_control = list(starts = matrix(c(.1, .23, .5), ncol = 1)), verbose = FALSE)
  # A glmmTMB fit with habitat at `sigma`, built outside the package refit.
  native_nll <- function(sigma) {
    d$habitat <- .msr_eval_kernel_engine(fit$mixed_fit$engine,
      log(sigma / x$f$inputs$unit_conv))$value[, 1]
    -as.numeric(stats::logLik(glmmTMB::glmmTMB(x$form, data = d,
      family = glmmTMB::nbinom2())))
  }
  cutoff <- fit$optim_results$value + stats::qchisq(.95, 1) / 2

  s <- summary(fit, profile = TRUE)
  tab <- s$opt_scale
  expect_identical(attr(tab, "interval_method"), c(habitat = "profile"))
  expect_false(any(attr(tab, "bound_limited")))
  limits <- unname(unlist(tab[1, c("2.5%", "97.5%")]))
  expect_lt(limits[1], fit$scale_est$Mean)
  expect_gt(limits[2], fit$scale_est$Mean)
  # Each limit is where the likelihood-ratio statistic reaches the 95% cutoff.
  expect_lt(max(abs(vapply(limits, native_nll, 0) - cutoff)), 1e-3)
  # The SE and the default Wald limits are unchanged.
  expect_equal(tab$SE, fit$scale_est$SE)
  wald <- summary(fit)$opt_scale
  expect_identical(unname(attr(wald, "interval_method")), "log_wald")
  expect_false(isTRUE(all.equal(unlist(wald[1, c("2.5%", "97.5%")]),
                                unlist(tab[1, c("2.5%", "97.5%")]))))
  # Distance limits convert the profile limits.
  expect_equal(unname(unlist(s$opt_dist_2d[1, c("2.5%", "97.5%")])),
    vapply(limits, function(z) k_dist(z, prob = .9, kernel = "gaussian",
      dimension = "2d"), 0), tolerance = 1e-3)
  expect_output(print(s), "Profile-likelihood confidence limits")

  # profile_sigma() refits with the coefficients, random-intercept SD and
  # dispersion re-estimated at each sigma.
  grid <- c(100, fit$scale_est$Mean, 300)
  ps <- profile_sigma(fit, sigma_values = grid, verbose = FALSE)
  expect_s3_class(ps, "sigma_profile")
  expect_equal(ps$sigma_grid, grid)
  expect_equal(ps$profiles$LL, -vapply(grid, native_nll, 0), tolerance = 1e-8)
  expect_equal(ps$profiles$LL[2], -as.numeric(fit$optim_results$value), tolerance = 1e-8)
  expect_true(all(is.finite(ps$profiles$AICc)))
  expect_error(profile_sigma(fit, sigma_values = c(20, 100, 300), verbose = FALSE),
               "within the fitted bounds")
  expect_error(profile_sigma(fit, n_cores = 2, verbose = FALSE), "serial")

  # A stored optimum the profile beats means the search missed a better one.
  shifted <- fit
  shifted$optim_results$value <- shifted$optim_results$value + 1
  expect_error(summary(shifted, profile = TRUE), "better optimum")
  # Fits without the stored refit data cannot be profiled.
  lean <- fit
  lean$mixed_fit$prediction_data <- NULL
  expect_error(summary(lean, profile = TRUE), "lacks the stored data")
  expect_error(profile_sigma(lean, verbose = FALSE), "lacks the stored data")
  # A failed refit is an error, never a likelihood value.
  broken <- fit
  broken$mixed_fit$prediction_data$temp <- NULL
  expect_error(summary(broken, profile = TRUE), "native refit failed")
  expect_error(profile_sigma(broken, sigma_values = grid, verbose = FALSE),
               "every native refit failed")
})

test_that("joint_mixed profile limits at a search bound are labeled", {
  skip_joint_mixed()
  x <- mixed_fixture()
  d <- x$d
  inputs <- x$f$inputs
  inputs$max_D <- 100
  d$habitat <- .msr_eval_kernel_engine(x$e, log(.09))$value[, 1]
  m <- glmmTMB::glmmTMB(x$form, data = d, family = poisson())
  bounded <- multiScale_optim(m, inputs, backend = "joint_mixed", n_cores = 1,
    profile_control = list(starts = matrix(c(.06, .09), ncol = 1)), verbose = FALSE)
  expect_equal(bounded$scale_est$Mean, 100)
  s <- summary(bounded, profile = TRUE)
  tab <- s$opt_scale
  expect_identical(attr(tab, "interval_method"), c(habitat = "profile-bound-limited"))
  expect_identical(unname(attr(tab, "bound_limited")[1, ]), c(FALSE, TRUE))
  expect_equal(tab[1, "97.5%"], 100)
  expect_lt(tab[1, "2.5%"], 100)
  # The SE stays withheld at the bound; the profile still gives a lower limit.
  expect_true(is.na(tab$SE))
  d$habitat <- .msr_eval_kernel_engine(bounded$mixed_fit$engine,
    log(tab[1, "2.5%"] / inputs$unit_conv))$value[, 1]
  nll <- -as.numeric(stats::logLik(glmmTMB::glmmTMB(x$form, data = d, family = poisson())))
  expect_lt(abs(nll - bounded$optim_results$value - stats::qchisq(.95, 1) / 2), 1e-3)
  expect_output(print(s), "not fully resolved")
  expect_output(print(s), "Scale SE withheld for habitat")
})

test_that("joint_mixed profiles reoptimize the other scales", {
  skip_joint_mixed()
  # Two binned radial landscapes, generated as in
  # tools/benchmarks/joint_mixed_wald/design.R, with 200 points in 40 sites.
  set.seed(11)
  n <- 200L
  dist <- seq(.015, 1, length.out = 23L)
  landscape <- function() {
    count <- matrix(sample(1:5, n * length(dist), replace = TRUE), n)
    radial <- outer(rnorm(n), sin(7 * dist)) + outer(rnorm(n), cos(12 * dist)) +
      outer(rnorm(n), dist) + matrix(rnorm(n * length(dist), sd = .25), n)
    list(value = radial * count, count = count)
  }
  l1 <- landscape(); l2 <- landscape()
  inputs <- list(raw_cov = NULL, d_list = NULL, min_D = 50, max_D = 800,
    unit_conv = 1000, kernel = "gaussian",
    binned = list(bin_rep = dist, vsum = list(h1 = l1$value, h2 = l2$value),
      csum = list(h1 = l1$count, h2 = l2$count),
      point_ids = as.character(seq_len(n)), n_points = n))
  engine <- .msr_prepare_kernel_engine(inputs$binned, c("h1", "h2"))
  kernel <- function(log_scales) .msr_eval_kernel_engine(engine, log_scales)$value
  site <- factor(sample(rep(seq_len(40L), length.out = n)))
  truth <- kernel(log(c(.15, .35)))
  dat <- data.frame(z = rnorm(n), site = site)
  eta <- .4 + .8 * truth[, "h1"] + .6 * truth[, "h2"] + .3 * dat$z +
    rnorm(40L, 0, .5)[as.integer(site)]
  dat$y <- rpois(n, exp(eta))
  start <- kernel(log(c(.2, .2)))
  dat$h1 <- start[, "h1"]; dat$h2 <- start[, "h2"]
  form <- y ~ h1 + h2 + z + (1 | site)
  m <- glmmTMB::glmmTMB(form, data = dat, family = poisson())
  fit <- multiScale_optim(m, inputs, backend = "joint_mixed", n_cores = 1, verbose = FALSE)
  tab <- summary(fit, profile = TRUE)$opt_scale
  expect_identical(attr(tab, "interval_method"), c(h1 = "profile", h2 = "profile"))
  expect_true(all(tab[, "2.5%"] < tab$Mean & tab$Mean < tab[, "97.5%"]))

  # At the upper h1 limit, minimizing over h2 with native glmmTMB refits
  # reaches the cutoff, and holding h2 at its estimate does not.
  native_nll <- function(log_scales) {
    k <- kernel(log_scales)
    dat$h1 <- k[, "h1"]; dat$h2 <- k[, "h2"]
    -as.numeric(stats::logLik(glmmTMB::glmmTMB(form, data = dat, family = poisson())))
  }
  z1 <- log(tab["h1", "97.5%"] / 1000)
  k2 <- fit$mixed_fit$log_scale[2]
  best <- stats::optimize(function(z2) native_nll(c(z1, z2)),
    c(k2 - .5, k2 + .5), tol = 1e-5)
  cutoff <- fit$optim_results$value + stats::qchisq(.95, 1) / 2
  expect_lt(abs(best$objective - cutoff), 1e-3)
  expect_gt(native_nll(c(z1, k2)), best$objective)
  # profile_sigma() holds the other scale at its estimate.
  ps <- profile_sigma(fit, sigma_values = c(100, 150, 200), verbose = FALSE)
  expect_equal(ps$profiles$LL[ps$profiles$variable == "h1"],
    -vapply(log(c(.1, .15, .2)), function(z) native_nll(c(z, k2)), 0), tolerance = 1e-8)
})

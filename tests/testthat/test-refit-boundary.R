test_that("refit bound flags use a unit-invariant declared log tolerance", {
  bounds <- function(parameter, lower = rep(1, length(parameter)),
      upper = rep(10, length(parameter)), active = rep(TRUE, length(parameter)),
      multiplier = rep(1, length(parameter)))
    .msr_refit_bounds(parameter, lower, upper, active, multiplier,
      paste0("sigma:x", seq_along(parameter)))
  interior <- bounds(c(2, 5))
  expect_false(interior$active_boundary)
  expect_false(any(interior$parameters$at_lower | interior$parameters$at_upper))
  actual <- bounds(c(1, 10))
  expect_true(actual$active_boundary)
  expect_identical(actual$parameters$at_lower, c(TRUE, FALSE))
  expect_identical(actual$parameters$at_upper, c(FALSE, TRUE))
  expect_identical(bounds(c(exp(.5e-5), 10 * exp(-.5e-5)))$parameters$at_lower,
    c(TRUE, FALSE))
  expect_false(bounds(c(exp(2e-5), 10 * exp(-2e-5)))$active_boundary)
  converted <- bounds(c(1, 10) / 1000, c(1, 1) / 1000, c(10, 10) / 1000,
    multiplier = c(1000, 1000))
  expect_identical(converted$parameters, actual$parameters)
  expect_identical(converted$tolerance_scale, "absolute_log_parameter_gap")
  expect_identical(bounds(0, 0)$parameters$at_lower, TRUE)
  expect_identical(bounds(1e-8, 0)$parameters$at_lower, FALSE)
  expect_false(bounds(c(1, 10), active = c(FALSE, FALSE))$active_boundary)
  expect_true(bounds(c(1, 5), active = c(TRUE, FALSE))$active_boundary)
  expect_false(bounds(c(1, 5), active = c(FALSE, TRUE))$active_boundary)
  expect_error(bounds(NA_real_), "aligned refit")
  expect_error(bounds(1, lower = 10, upper = 1), "aligned refit")
})

test_that("a bound withholds only the bounded parameter's outer Wald SE", {
  hessian <- matrix(c(4, 1, 1, 9), 2)
  raw <- .msr_refit_hessian_se(hessian, c(TRUE, TRUE), c(10, 1))
  bounds <- .msr_refit_bounds(c(1, 5), c(1, 1), c(10, 10), c(TRUE, TRUE),
    c(10, 1), c("sigma:a", "sigma:b"))
  # `a` sits at its lower bound. `b` keeps the SE of its own Hessian block,
  # which holds `a` fixed: sqrt(1 / 9) * 1, not the joint sqrt(4 / 35).
  actual <- .msr_refit_bound_uncertainty(raw, bounds, hessian, c(10, 1))
  expect_true(is.na(actual$SE[1]))
  expect_equal(actual$SE[2], 1 / 3)
  expect_false(isTRUE(all.equal(actual$SE[2], raw$SE[2])))
  expect_identical(actual$diagnostic$conditional_on, "sigma:a")
  expect_identical(actual$diagnostic$status, "parameter_at_optimization_bound")
  expect_identical(actual$diagnostic$curvature_status, "available")
  expect_true(actual$diagnostic$triggered)
  expect_true(all(is.finite(raw$SE)))
  expect_identical(.msr_refit_bound_uncertainty(actual, bounds, hessian, c(10, 1)), actual)
  # Without a Hessian, every smooth SE is withheld.
  bare <- .msr_refit_bound_uncertainty(raw, bounds)
  expect_true(all(is.na(bare$SE)))
  expect_identical(bare$diagnostic$conditional_on, character())
  # Curvature that fails only in the bounded direction leaves `b` usable.
  invalid <- .msr_refit_hessian_se(diag(c(Inf, 9)), c(TRUE, TRUE), c(10, 1))
  rescued <- .msr_refit_bound_uncertainty(invalid, bounds, diag(c(Inf, 9)), c(10, 1))
  expect_identical(rescued$diagnostic$curvature_status, "invalid_hessian")
  expect_equal(rescued$SE, c(NA_real_, 1 / 3))
  # A failed refit anywhere in the Hessian still withholds every smooth SE.
  failed <- .msr_refit_hessian_se(diag(c(1e40, 9)), c(TRUE, TRUE), c(10, 1))
  expect_identical(failed$diagnostic$status, "hessian_includes_failed_refit")
  expect_true(all(is.na(.msr_refit_bound_uncertainty(failed, bounds,
    diag(c(1e40, 9)), c(10, 1))$SE)))
  # An interior parameter outside the smooth mask stays missing.
  hard <- .msr_refit_hessian_se(hessian, c(TRUE, FALSE), c(10, 1))
  hard_bounds <- .msr_refit_bounds(c(1, 5), c(1, 1), c(10, 10), c(TRUE, FALSE),
    c(10, 1), c("sigma:a", "sigma:b"))
  expect_true(all(is.na(.msr_refit_bound_uncertainty(hard, hard_bounds,
    hessian, c(10, 1))$SE)))
  expect_identical(hessian, matrix(c(4, 1, 1, 9), 2))
})

test_that("legacy reporting copies suppress bound Wald uncertainty without changing models", {
  dat <- data.frame(y = c(0, 1, 2, 2, 4, 3), x = seq(-1, 1, length.out = 6))
  model <- glm(y ~ x, data = dat, family = poisson())
  object <- structure(list(scale_est = data.frame(Mean = 10, SE = 2, row.names = "x"),
    shape_est = NULL, min_D = 10, max_D = 100, kernel_inputs = list(kernel = "gaussian",
      unit_conv = 10), optim_results = list(par = 1, hessian = matrix(4, 1)),
    opt_mod = model, diagnostics = list()), class = "multiScaleR")
  original <- object
  view <- .msr_refit_bound_view(object)
  expect_true(is.na(view$scale_est$SE))
  expect_identical(view$scale_est$Mean, object$scale_est$Mean)
  expect_identical(view$opt_mod, model)
  expect_identical(view$optim_results, object$optim_results)
  expect_identical(stats::vcov(view$opt_mod), stats::vcov(model))
  expect_identical(.msr_selection_parameter_count(view), .msr_selection_parameter_count(object))
  expect_identical(object, original)
  expect_true(view$diagnostics$refit_hessian$bounds$active_boundary)
  interior <- object
  interior$optim_results$par <- 4
  interior$scale_est$Mean <- 40
  expect_identical(.msr_refit_bound_view(interior), interior)
  fixed <- object
  fixed$kernel_inputs$kernel <- "fixed"
  expect_identical(.msr_refit_bound_view(fixed), fixed)
  analytic <- object
  analytic$backend <- "profile_glm"
  expect_identical(.msr_refit_bound_view(analytic), analytic)
  incomplete <- object
  incomplete$optim_results <- NULL
  expect_identical(.msr_refit_bound_view(incomplete), incomplete)
  shape <- object
  shape$kernel_inputs$kernel <- "expow"
  shape$optim_results$par <- c(4, .75)
  shape$scale_est$Mean <- 40
  shape$shape_est <- data.frame(Mean = .75, SE = .1, row.names = "x")
  actual <- .msr_refit_bound_view(shape)
  expect_true(all(is.na(c(actual$scale_est$SE, actual$shape_est$SE))))
  expect_identical(actual$diagnostics$refit_hessian$bounds$parameters$at_lower,
    c(FALSE, TRUE))
  # A saved two-scale fit stored joint SEs. The view keeps the interior scale,
  # recomputed from its own block of the stored Hessian, and withholds only the
  # bounded scale and its Wald limits.
  pair <- object
  pair$scale_est <- data.frame(Mean = c(10, 40), SE = c(2, 3), row.names = c("x", "z"))
  pair$optim_results <- list(par = c(1, 4), hessian = matrix(c(4, 1, 1, 9), 2))
  pair$profile_scale_est <- structure(data.frame(Mean = c(10, 40), SE = c(2, 3),
    `2.5%` = c(6, 34), `97.5%` = c(14, 46), row.names = c("x", "z"), check.names = FALSE),
    interval_method = c(x = "wald", z = "wald"))
  viewed <- .msr_refit_bound_view(pair)
  expect_equal(viewed$scale_est$SE, c(NA, 10 / 3))
  expect_identical(viewed$diagnostics$refit_hessian$conditional_on, "sigma:x")
  expect_true(all(is.na(viewed$profile_scale_est["x", c("SE", "2.5%", "97.5%")])))
  expect_equal(unlist(viewed$profile_scale_est["z", c("SE", "2.5%", "97.5%")]),
    c(SE = 3, `2.5%` = 34, `97.5%` = 46))
  expect_identical(.msr_refit_bound_view(viewed), viewed)
})

test_that("bound refits retain native covariance, penalties, profiles and prediction support", {
  # Reuse the existing15-site seed42 fixture. The outer optimizer is mocked;
  # its actual GLM objective/final calls still use the declared lower bound.
  f <- make_core_fixture()
  optim_at_lower <- function(par, fn, ..., lower, upper, method, hessian, control) {
    score <- do.call(fn, c(list(par = lower), list(...)))
    list(par = lower, value = score, convergence = 0L,
      counts = stats::setNames(c(1L, 0L), c("function", "gradient")),
      message = "declared boundary fixture", hessian = diag(4, length(lower)))
  }
  with_mocked_bindings({
    quiet <- capture.output(fit <- multiScale_optim(f$fitted_mod, f$kernel_inputs,
      par = f$opt$optim_results$par, n_cores = 1L, verbose = FALSE))
  }, optim = optim_at_lower)
  expect_false(any(grepl("WARNING!!!|Outer Wald", quiet)))
  expect_equal(fit$scale_est$Mean, rep(fit$min_D, nrow(fit$scale_est)))
  expect_true(all(is.na(fit$scale_est$SE)))
  expect_identical(fit$diagnostics$refit_hessian$status, "parameter_at_optimization_bound")
  expect_identical(fit$diagnostics$refit_hessian$curvature_status, "available")
  expect_identical(fit$optim_results$hessian_unscale, fit$optim_results$hessian)
  expect_true(fit$diagnostics$refit_hessian$bounds$active_boundary)
  original <- fit
  covariance <- stats::vcov(fit$opt_mod)
  count <- attr(stats::logLik(fit$opt_mod), "df") + nrow(fit$scale_est)
  expect_equal(aic_tab(list(fit), AICc = FALSE, verbose = FALSE)$K, count)
  expect_equal(bic_tab(list(fit), verbose = FALSE)$K, count)
  report <- summary(fit)
  expect_true(all(is.na(report$opt_scale[, c("SE", "2.5%", "97.5%")])) )
  expect_identical(stats::vcov(report$fitted_mod), covariance)
  expect_true(any(grepl("optimization bound", capture.output(print(report)), fixed = TRUE)))
  expect_true(any(grepl("optimization bound", capture.output(print(fit)), fixed = TRUE)))
  for (dimension in c("1d", "2d")) {
    expect_true(all(is.na(kernel_dist(fit, dimension = dimension)[, c("2.5%", "97.5%")])) )
  }
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot(fit)
  expect_false(any(vapply(plots, function(p) any(vapply(p$layers,
    function(layer) inherits(layer$geom, "GeomRect"), logical(1))), logical(1))))
  old <- options(multiScaleR.profile_cache_max_bytes = 0)
  on.exit(options(old), add = TRUE)
  endpoints <- cbind(`2.5%` = fit$scale_est$Mean, `97.5%` = 2 * fit$scale_est$Mean)
  rownames(endpoints) <- rownames(fit$scale_est)
  profiled <- with_mocked_bindings(summary(fit, profile = TRUE),
    profile_sigma_intervals = function(...) list(intervals = endpoints,
      method = stats::setNames(rep("profile_bound", nrow(endpoints)), rownames(endpoints))))
  expect_true(all(is.na(profiled$opt_scale$SE)))
  expect_equal(as.matrix(profiled$opt_scale[, c("2.5%", "97.5%")]), endpoints)
  native <- stats::predict(fit$opt_mod, type = "response", se.fit = TRUE)
  actual <- predict(fit, type = "response", se.fit = TRUE, uncertainty = "conditional")
  expect_equal(as.numeric(actual$fit), as.numeric(native$fit), tolerance = 1e-12)
  expect_equal(as.numeric(actual$se.fit), as.numeric(native$se.fit), tolerance = 1e-12)
  contract <- attr(actual, "prediction_contract")
  expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
  expect_true(contract$kernel_support$optimization_bounds$active_boundary)
  expect_identical(contract$kernel_support$extraction_window, fit$diagnostics$max_distance)
  expect_false(contract$kernel_support$establishes_eligibility)
  snapshot <- prepare_prediction(fit, component = "response")
  expect_true(snapshot$prediction_contract$response$kernel_support$optimization_bounds$active_boundary)
  restored <- unserialize(serialize(snapshot, NULL, version = 2))
  expect_identical(restored$prediction_contract, snapshot$prediction_contract)
  expect_equal(predict(restored, se.fit = TRUE, uncertainty = "conditional"),
    predict(snapshot, se.fit = TRUE, uncertainty = "conditional"))
  expect_identical(stats::vcov(fit$opt_mod), covariance)
  expect_identical(fit, original)
})

test_that("one scale at a bound keeps SEs for the other scales in real fits", {
  skip_on_cran()
  # land2 has a true sigma of 500 m, so a 400 m search extent pins it at the
  # upper bound while land1 (true sigma 250 m) stays interior.
  data("landscape_counts", package = "multiScaleR", envir = environment())
  data("surv_pts", package = "multiScaleR", envir = environment())
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = terra::vect(surv_pts), raster_stack = land,
                        max_D = 400, kernel = "gaussian", verbose = FALSE)
  df <- data.frame(landscape_counts, inputs$kernel_dat)
  model <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)

  fit <- suppressWarnings(multiScale_optim(model, inputs, n_cores = 1L, verbose = FALSE))
  h <- fit$diagnostics$refit_hessian
  expect_identical(h$status, "parameter_at_optimization_bound")
  expect_identical(h$conditional_on, "sigma:land2")
  expect_true(is.na(fit$scale_est["land2", "SE"]))
  # The interior SE is the Wald SE of land1's own Hessian entry.
  expected <- sqrt(1 / fit$optim_results$hessian_unscale[1, 1]) * inputs$unit_conv
  expect_equal(fit$scale_est["land1", "SE"], expected)
  report <- summary(fit)
  expect_true(all(is.finite(unlist(report$opt_scale["land1", c("SE", "2.5%", "97.5%")]))))
  expect_true(all(is.na(report$opt_scale["land2", c("SE", "2.5%", "97.5%")])))
  printed <- capture.output(print(report))
  expect_true(any(grepl("At an optimization bound: sigma:land2", printed, fixed = TRUE)))
  expect_true(any(grepl("hold the bounded ones fixed", printed, fixed = TRUE)))
  expect_true(all(is.finite(unlist(kernel_dist(fit)["land1", c("2.5%", "97.5%")]))))

  # profile_glm applies the same rule from joint observed information.
  expect_warning(analytic <- multiScale_optim(model, inputs, backend = "profile_glm",
    n_cores = 1L, verbose = FALSE), "Scale SEs are withheld for: land2")
  expect_false(analytic$uncertainty$available)
  expect_true(is.na(analytic$scale_est["land2", "SE"]))
  expect_equal(analytic$scale_est["land1", "SE"], fit$scale_est["land1", "SE"],
               tolerance = 1e-3)
  expect_identical(analytic$diagnostics$profile_glm$scale_se_conditional_on, "sigma:land2")
  tab <- summary(analytic)$opt_scale
  expect_identical(unname(attr(tab, "interval_method")), c("log_wald", "unavailable"))
  expect_true(any(grepl("hold sigma:land2 fixed", capture.output(print(summary(analytic))),
                        fixed = TRUE)))
  # Joint covariance and joint prediction uncertainty stay withheld.
  expect_error(predict(analytic, se.fit = TRUE), "Joint Wald uncertainty unavailable")
})

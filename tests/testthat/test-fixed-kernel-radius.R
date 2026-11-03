fixed_radius_plot_layer <- function(plot, geom) {
  which(vapply(plot$layers, function(layer) inherits(layer$geom, geom), logical(1)))
}

test_that("fixed kernels classify every estimated scale as a hard radius", {
  specifications <- data.frame(covariate = c("mean", "weighted", "surface", "landscape"),
    type = c("kernel", "surface", "surface", "landscape"),
    weighted = c(FALSE, TRUE, FALSE, FALSE))
  expect_identical(.msr_hard_radius_covariates(specifications), c("surface", "landscape"))
  for (kernel in c("gaussian", "exp", "expow"))
    expect_identical(.msr_hard_radius_covariates(specifications, kernel), c("surface", "landscape"))
  expect_identical(.msr_hard_radius_covariates(specifications, "fixed"), specifications$covariate)
  expect_identical(.msr_hard_radius_covariates(NULL, "fixed", c("old_mean", "old_sq")),
    c("old_mean", "old_sq"))
  expect_identical(.msr_hard_radius_covariates(NULL), character())
})

test_that("saved fixed-kernel fits cannot report legacy Hessian radius limits", {
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  base <- make_core_fixture()$opt
  for (kind in c("kernel", "weighted", "surface", "landscape", "legacy")) {
    fit <- base; fit$kernel_inputs$kernel <- "fixed"; fit$shape_est <- NULL
    fit$scale_est$SE <- fit$scale_est$Mean / 10
    fit$profile_scale_est <- NULL
    if (kind == "legacy") fit$kernel_inputs$scale_vars <- NULL else {
      fit$kernel_inputs$scale_vars$type[] <- if (kind == "weighted") "surface" else kind
      fit$kernel_inputs$scale_vars$weighted[] <- kind == "weighted"
    }
    before <- fit; covariance <- stats::vcov(fit$opt_mod)
    report <- summary(fit)
    radius <- round(fit$scale_est$Mean, 2)
    expect_identical(report$hard_radius, rownames(fit$scale_est))
    expect_true(all(is.na(report$opt_scale$SE)))
    expect_true(all(is.na(report$opt_scale[, c("2.5%", "97.5%")])))
    for (dimension in c("1d", "2d")) {
      direct <- kernel_dist(fit, prob = .81, dimension = dimension)
      expect_equal(direct$Mean, radius)
      expect_true(all(is.na(direct[, c("2.5%", "97.5%")])))
      expect_equal(report[[paste0("opt_dist_", dimension)]]$Mean, radius)
      plots <- plot(fit, prob = .81, dimension = dimension)
      expect_length(plots, nrow(fit$scale_est))
      for (view in plots) {
        expect_length(fixed_radius_plot_layer(view, "GeomRect"), 0L)
        expect_match(view$labels$caption, "Buffer radius")
        expect_equal(view$layers[[fixed_radius_plot_layer(view, "GeomVline")]]$data$xintercept, radius)
      }
    }
    expect_true(any(grepl("Radius", capture.output(print(report)))))
    expect_true(any(grepl("profile = TRUE", capture.output(print(fit)), fixed = TRUE)))
    expect_identical(stats::vcov(fit$opt_mod), covariance)
    expect_identical(fit, before)
  }
})

test_that("fixed-kernel summaries retain profile endpoints and explicit failures", {
  fit <- make_core_fixture()$opt; fit$kernel_inputs$kernel <- "fixed"
  fit$scale_est$SE <- fit$scale_est$Mean / 10
  fit$shape_est <- NULL; fit$profile_scale_est <- NULL
  covariates <- rownames(fit$scale_est); radius <- fit$scale_est$Mean
  old <- options(multiScaleR.profile_cache_max_bytes = 0); on.exit(options(old), add = TRUE)
  for (method in c("profile", "profile_bound")) {
    endpoints <- data.frame(`2.5%` = radius * .8, `97.5%` = radius * 1.3,
      check.names = FALSE, row.names = covariates)
    result <- with_mocked_bindings(summary(fit, profile = TRUE),
      profile_sigma_intervals = function(...) list(intervals = endpoints,
        method = stats::setNames(rep(method, length(covariates)), covariates)))
    expect_true(all(is.na(result$opt_scale$SE)))
    expect_equal(as.matrix(result$opt_scale[, c("2.5%", "97.5%")]), as.matrix(endpoints))
    expect_identical(unname(attr(result$opt_scale, "interval_method")), rep(method, length(covariates)))
    for (dimension in c("1d", "2d")) {
      distance <- result[[paste0("opt_dist_", dimension)]]
      expect_equal(distance$Mean, round(radius, 2))
      expect_equal(as.matrix(distance[, c("2.5%", "97.5%")]), as.matrix(round(endpoints, 2)))
    }
  }
  expect_error(with_mocked_bindings(summary(fit, profile = TRUE),
    profile_sigma_intervals = function(...) stop("retained profile fixture failure")),
    "Profile-likelihood intervals could not be computed: retained profile fixture failure")
})

test_that("fixed-kernel mean and weighted SQ refits omit only radius Hessian SEs", {
  # One actual refit at the declared36-site fixture; no retries or seed search.
  f <- refit_spatial_fixture(kernel = "fixed")
  fit <- multiScale_optim(f$object$opt_mod, f$inputs,
    par = f$object$optim_results$par, n_cores = 1L, backend = "refit", verbose = FALSE)
  expect_identical(fit$optim_results$convergence, 0L)
  expect_true(isTRUE(fit$opt_mod$converged))
  expect_true(all(is.finite(stats::coef(fit$opt_mod))))
  expect_true(all(is.finite(fit$scale_est$Mean)))
  expect_setequal(rownames(fit$scale_est), c("mean", "rough"))
  expect_true(all(is.na(fit$scale_est$SE)))
  expect_identical(fit$diagnostics$refit_hessian$status, "no_active_smooth_parameters")
  expect_false(fit$diagnostics$refit_hessian$triggered)
  expect_identical(fit$diagnostics$refit_hessian$active, c(FALSE, FALSE))
  expect_identical(fit$optim_results$hessian_unscale, fit$optim_results$hessian)
  expect_true(is.matrix(fit$optim_results$hessian))
  expect_null(fit$diagnostics$sigma_precision)
  covariance <- stats::vcov(fit$opt_mod)
  expect_true(all(is.finite(covariance)))
  expect_gt(min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values), 0)
  result <- summary(fit)
  expect_identical(stats::vcov(result$fitted_mod), covariance)
  expected_count <- as.numeric(attr(stats::logLik(fit$opt_mod), "df")) + 2L
  aic <- aic_tab(list(radius = fit), AICc = FALSE, verbose = FALSE)
  bic <- bic_tab(list(radius = fit), verbose = FALSE)
  expect_equal(aic$K, expected_count)
  expect_equal(bic$K, expected_count)
  ll <- as.numeric(stats::logLik(fit$opt_mod))
  expect_equal(aic$AIC, -2 * ll + 2 * expected_count, tolerance = 1e-8)
  expect_equal(bic$BIC, -2 * ll + log(stats::nobs(fit$opt_mod)) * expected_count, tolerance = 1e-8)
})

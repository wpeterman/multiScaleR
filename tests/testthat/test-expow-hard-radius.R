test_that("expow mixtures retain radius transforms without estimating a hard shape", {
  f <- refit_spatial_fixture("expow", weighted = FALSE)
  other <- f$object
  other$shape_est["rough", c("Mean", "SE")] <- NA_real_
  attr(other$shape_est, "not_applicable") <- "rough"
  k <- attr(logLik(other$opt_mod), "df") + 3L
  table <- aic_tab(list(other), AICc = FALSE)
  expect_equal(table$K, k)
  expect_equal(table$AIC, -2 * as.numeric(logLik(other$opt_mod)) + 2 * k)
  altered <- other
  altered$optim_results$par[4] <- 10
  expect_equal(as.numeric(predict(altered, f$data, f$inputs)),
    as.numeric(predict(other, f$data, f$inputs)), tolerance = 1e-12)
  expect_equal(as.numeric(predict(other, f$data, f$inputs, include_fixed = TRUE)),
    as.numeric(predict(other$opt_mod, f$data, type = "response")), tolerance = 1e-10)
  expect_equal(.msr_shape_values(other)[["rough"]], 1)
  recommendation <- .next_run_recommendation(other$scale_est, other$shape_est,
    other$max_D, .empty_diagnostics())
  expect_true(all(is.finite(recommendation$start_par)))
  expect_equal(recommendation$start_shape[["rough"]], 1)
  snap <- prepare_prediction(other, f$data[1:3, ], f$inputs, include_fixed = TRUE)
  expect_equal(as.numeric(predict(snap)),
    as.numeric(predict(other, f$data[1:3, ], f$inputs, include_fixed = TRUE)))
})

test_that("public expow refitting supports hard-radius surface mixtures end to end", {
  f <- refit_spatial_fixture("expow", weighted = FALSE)
  object <- suppressMessages(multiScale_optim(f$object$opt_mod, f$inputs,
    par = f$object$optim_results$par, verbose = FALSE, n_cores = NULL))
  expect_s3_class(object, "multiScaleR")
  expect_true(is.na(object$shape_est["rough", "Mean"]))
  expect_true(is.na(object$shape_est["rough", "SE"]))
  expect_identical(attr(object$shape_est, "not_applicable"), "rough")
  expect_true(all(is.finite(object$optim_results$par)))
  expect_equal(.msr_selection_parameter_count(object),
    attr(logLik(object$opt_mod), "df") + 3L)
  nd <- .model_data(object$opt_mod)[c(7, 19, 3), ]
  expect_equal(as.numeric(predict(object, nd, f$inputs, include_fixed = TRUE)),
    as.numeric(predict(object$opt_mod, nd, type = "response")), tolerance = 1e-10)
  summary <- summary(object)
  expect_equal(summary$opt_dist["rough", "Mean"], round(object$scale_est["rough", "Mean"], 2))
  expect_true(is.na(summary$opt_shape["rough", "Mean"]))
  profiles <- profile_sigma(object, sigma_values = c(18, 28, 38), verbose = FALSE)
  expect_s3_class(profiles, "sigma_profile")
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_true(all(vapply(plot(object), inherits, logical(1), "ggplot")))
  expect_warning(rasters <- kernel_scale.raster(f$rasters, multiScaleR = object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
    offset_values = c(exposure = 1)), "group.*categorical")
  expect_s4_class(rasters, "SpatRaster")
  expect_true(all(c("mean", "rough", "fixed") %in% names(rasters)))
})

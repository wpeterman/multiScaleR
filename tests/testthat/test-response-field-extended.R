test_that("more native response families and spatial covariances retain prediction contracts", {
  skip_if_not_installed("glmmTMB")
  for (family in c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"))
    for (covariance in c("exp", "gau", "mat")) {
      f <- extended_response_field_fixture(family, covariance)
      expect_identical(f$model$fit$convergence, 0L)
      expect_true(f$model$sdr$pdHess)
      expect_true(is.finite(as.numeric(logLik(f$model))))
      expect_false("year" %in% names(f$data))
      expected_K <- 4L + 1L + as.integer(family %in% c("gaussian", "nbinom2", "Gamma"))
      expect_equal(.msr_parameter_count(f$model), expected_K)
      expect_equal(aic_tab(list(native = f$model), AICc = FALSE)$K, expected_K)
      expect_equal(aic_tab(list(field = f$object), AICc = FALSE)$K, expected_K + 2L)
      geometry <- response_field_geometry(f$object, terra::crds(f$points), terra::crs(f$points), batch_size = 2L)
      expect_s3_class(geometry, "multiScaleR_field_geometry")
      rows <- c(7, 42, 7, 105); nd <- f$data[rows, , drop = FALSE]
      keys <- data.frame(row_id = rownames(f$data)); nd$row_id <- keys$row_id[rows]
      for (zero in c(FALSE, TRUE)) {
        args <- if (zero) list(re.form = NA) else list()
        native <- do.call(predict, c(list(object = f$model, newdata = nd, type = "response", se.fit = TRUE), args))
        actual <- predict(f$object, nd, f$inputs, join_by = keys, type = "response", se.fit = TRUE,
          uncertainty = "conditional", native_args = args)
        expect_equal(as.numeric(actual$fit), as.numeric(native$fit), tolerance = 1e-8)
        expect_equal(as.numeric(actual$se.fit), as.numeric(native$se.fit), tolerance = 1e-8)
        if (family == "binomial") expect_true(all(actual$fit >= 0 & actual$fit <= 1))
        contract <- attr(actual, "prediction_contract")
        expect_identical(contract$response_spatial_field$fields[[1L]]$covariance,
          switch(covariance, exp = "exponential", gau = "gaussian", mat = "matern"))
        expect_false(contract$response_spatial_field$kernel_uncertainty_included)
        expect_identical(contract$integration$random_effects, if (zero) "set_to_zero_not_integrated" else "native_method")
      }
      snapshot <- prepare_prediction(f$object, nd, f$inputs, join_by = keys)
      restored <- unserialize(serialize(snapshot, NULL, version = 2L))
      expect_null(restored$opt_mod)
      expect_equal(as.numeric(predict(restored)),
        as.numeric(predict(f$object, nd, f$inputs, join_by = keys)), tolerance = 1e-8)
      # New positions in two-row batches use the same native mean and SE target.
      novel <- nd; novel$field_pos <- glmmTMB::numFactor(novel$x + 7, novel$y_coordinate + 3)
      native <- predict(f$model, novel, type = "response", se.fit = TRUE, allow.new.levels = TRUE)
      thin <- list(object = list(opt_mod = f$model), type = "response", se.fit = TRUE,
        native_args = list(allow.new.levels = TRUE), logical_names = character(), field_geometry = geometry)
      novel$.msr_field_x <- novel$x + 7; novel$.msr_field_y <- novel$y_coordinate + 3
      actual <- .msr_native_raster_eval(thin, novel)
      expect_equal(actual$fit, as.numeric(native$fit), tolerance = 1e-7)
      expect_equal(actual$se.fit, as.numeric(native$se.fit), tolerance = 1e-7)
    }
  cat("RESPONSE_FIELD_EXTENDED families=5 covariances=3 fixed_range=1 native_SE=1\n")
})

test_that("native field effects retain zero-field contracts and guarded joint fitting", {
  skip_if_not_installed("glmmTMB")
  f <- extended_response_field_fixture("nbinom2", "mat")
  path <- tempfile(fileext = ".pdf"); on.exit(unlink(path), add = TRUE)
  grDevices::pdf(path); on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot_marginal_effects(f$object, length.out = 5L)
  expect_true(all(c("mean", "rough", "z") %in% names(curves)))
  contract <- attr(curves, "prediction_contract")
  expect_identical(contract$integration$random_effects, "set_to_zero_not_integrated")
  expect_identical(contract$response_spatial_field$conditioning, "all_response_random_effects_set_to_zero")
  expect_false(contract$response_spatial_field$kernel_uncertainty_included)
  curve <- curves$mean$data
  nd <- f$data[rep(1L, nrow(curve)), , drop = FALSE]
  nd$mean <- (curve$x - f$object$scl_params$mean[["mean"]]) / f$object$scl_params$sd[["mean"]]
  for (name in c("rough", "z", "exposure")) nd[[name]] <- mean(f$data[[name]])
  expected <- predict(f$model, nd, type = "response", re.form = NA, se.fit = TRUE)
  expect_equal(curve$fit, as.numeric(expected$fit), tolerance = 1e-7)
  expect_equal(curve$lwr, as.numeric(expected$fit + qnorm(.025) * expected$se.fit), tolerance = 1e-7)
  expect_match(curves$mean$labels$caption, "set to zero, not integrated", fixed = TRUE)
  expect_false(.msr_glmmTMB_response_fields(f$model)$fit_eligibility_established)
})

test_that("extended spatial families retain native fixed-reference refits and raster bands", {
  skip_if_not_installed("glmmTMB")
  for (pair in list(c("binomial", "mat"), c("nbinom2", "exp"), c("Gamma", "gau"))) {
    f <- extended_response_field_fixture(pair[1L], pair[2L])
    context <- build_opt_context(f$model, f$inputs$raw_cov, scale_vars = f$inputs$scale_vars,
      unit_conv = f$inputs$unit_conv, resolution = f$inputs$resolution,
      n_cols = f$inputs$n_cols, binned = f$inputs$binned, d_list = f$inputs$d_list)
    actual <- kernel_scale_fn(f$par, f$inputs$d_list, f$inputs$raw_cov, f$inputs$kernel,
      f$model, opt_context = context, mod_return = TRUE)
    expect_identical(actual$mod$fit$convergence, 0L)
    expect_true(actual$mod$sdr$pdHess)
    expect_equal(as.numeric(logLik(actual$mod)), as.numeric(logLik(f$model)), tolerance = 1e-5)
    expect_equal(.msr_parameter_count(actual$mod), .msr_parameter_count(f$model))
    expect_identical(actual$mod$modelInfo$reStruc$condReStruc, f$model$modelInfo$reStruc$condReStruc)
    expect_identical(actual$mod$modelInfo$map, f$model$modelInfo$map)
    expect_identical(actual$mod$modelInfo$priors, f$model$modelInfo$priors)
    expect_identical(stats::family(actual$mod)$family, stats::family(f$model)$family)
    expect_identical(stats::family(actual$mod)$link, stats::family(f$model)$link)
    expect_match(attr(actual$mod, "msr_refit_call")$initial_update_error, "provider_family")
    expect_equal(actual$scl_params$mean[c("mean", "rough")], f$object$scl_params$mean[c("mean", "rough")])
    expect_equal(actual$scl_params$sd[c("mean", "rough")], f$object$scl_params$sd[c("mean", "rough")])
  }
  f <- extended_response_field_fixture("nbinom2", "mat")
  geometry <- response_field_geometry(f$object, terra::crds(f$points), terra::crs(f$points), batch_size = 3L)
  whole <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
    site_values = list(z = .2, field_pos = as.character(f$data$field_pos[1L]), field_group = "shared"),
    offset_values = c(exposure = 2))
  raster <- terra::crop(whole, terra::ext(40, 80, 40, 80))
  attr(raster, "multiScaleR_projection") <- attr(whole, "multiScaleR_projection")
  nd <- terra::as.data.frame(raster)
  nd$field_pos <- glmmTMB::numFactor(terra::xyFromCell(raster, seq_len(terra::ncell(raster))))
  nd$field_group <- factor("shared", levels = levels(f$data$field_group))
  expected <- predict(f$model, nd, type = "response", se.fit = TRUE, allow.new.levels = TRUE)
  actual <- predict(f$object, raster, se.fit = TRUE, uncertainty = "conditional",
    native_args = list(allow.new.levels = TRUE), raster_args = list(field_geometry = geometry))
  expect_equal(as.numeric(terra::values(actual)[, 1L]), as.numeric(expected$fit), tolerance = 1e-7)
  expect_equal(as.numeric(terra::values(actual)[, 2L]), as.numeric(expected$se.fit), tolerance = 1e-7)
  expect_identical(attr(actual, "prediction_contract")$response_spatial_field$fields[[1L]]$covariance, "matern")
  expect_false(attr(actual, "prediction_contract")$response_spatial_field$kernel_uncertainty_included)
})

test_that("ordinary lme4 curves attach the same explicit zero-effect contract", {
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  f$data$y <- f$data$y + as.integer(f$data$site)^2
  model <- f$object$opt_mod <- lme4::lmer(y ~ mean + rough + fixed + z + (1 | site),
    data = f$data, REML = FALSE)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot_marginal_effects(f$object, length.out = 5L)
  contract <- attr(curves, "prediction_contract")
  expect_identical(contract$integration$random_effects, "set_to_zero_not_integrated")
  expect_true("observation_noise" %in% contract$uncertainty$omitted)
  expect_true("kernel_parameter_uncertainty" %in% contract$uncertainty$omitted)
  expect_null(contract$response_spatial_field)
  expect_match(curves$rough$labels$caption, "zero, not integrated", fixed = TRUE)
  nd <- f$data[rep(1L, 5L), , drop = FALSE]
  for (name in c("mean", "fixed", "z")) nd[[name]] <- mean(f$data[[name]])
  nd$rough <- (curves$rough$data$x - f$object$scl_params$mean[["rough"]]) / f$object$scl_params$sd[["rough"]]
  expected <- predict(model, nd, type = "response", re.form = NA, se.fit = TRUE)
  expect_equal(curves$rough$data$fit, as.numeric(expected$fit), tolerance = 1e-10)
  expect_equal(curves$rough$data$lwr,
    as.numeric(expected$fit + qnorm(.025) * expected$se.fit), tolerance = 1e-10)
})

test_that("native spatial covariance layouts stay distinct from iid effects and kernel scales", {
  skip_if_not_installed("glmmTMB")
  f <- refit_spatial_fixture()
  xy <- terra::crds(f$points)
  data <- f$data
  data$field_pos <- glmmTMB::numFactor(xy[, 1L], xy[, 2L])
  data$field_group <- factor(rep("shared", nrow(data)))
  for (kind in c("exp", "gau", "mat")) {
    prepared <- glmmTMB::glmmTMB(as.formula(paste0(
      "y ~ mean + rough + ", kind, "(field_pos + 0 | field_group)")),
      data, family = poisson(), doFit = FALSE)
    # Native input-layout metadata wrapper only, never a native fitted model.
    metadata <- structure(list(modelInfo = list(reStruc =
      list(condReStruc = prepared$condReStruc))), class = "glmmTMB")
    info <- .msr_glmmTMB_response_fields(metadata)
    expect_identical(info$status, "native_layout_identified")
    expect_false(info$fit_eligibility_established)
    expect_false(info$range_is_ecological_kernel_scale)
    expect_false(info$kernel_uncertainty_included)
    expect_identical(info$fields[[1L]]$coordinate_variables, "field_pos")
    expect_identical(info$fields[[1L]]$grouping_variables, "field_group")
    expect_equal(info$fields[[1L]]$locations_per_field, 36L)
    expect_equal(info$fields[[1L]]$covariance_parameter_count, if (kind == "mat") 3L else 2L)
    expect_identical(info$fields[[1L]]$covariance,
      switch(kind, exp = "exponential", gau = "gaussian", mat = "matern"))
    metadata$modelInfo$reStruc$condReStruc[[1L]]$blockSize <- 35L
    expect_error(.msr_glmmTMB_response_fields(metadata), "incomplete")
  }
  expect_null(.msr_glmmTMB_response_fields(f$object$opt_mod))
})

test_that("native response fields retain point prediction, uncertainty and parameter penalties", {
  skip_if_not_installed("glmmTMB")
  for (family in c("poisson", "gaussian")) {
    f <- native_response_field_fixture(family)
    model <- f$model
    expect_identical(model$fit$convergence, 0L)
    expect_true(model$sdr$pdHess)
    expect_false("year" %in% names(f$data))
    expect_equal(.msr_parameter_count(model), if (family == "poisson") 6 else 7)
    expect_equal(aic_tab(list(native = model), AICc = FALSE)$K, .msr_parameter_count(model))
    expect_equal(aic_tab(list(spatial = f$object), AICc = FALSE)$K, .msr_parameter_count(model) + 2)
    rows <- c(7, 42, 7, 105)
    nd <- f$data[rows, , drop = FALSE]
    keys <- data.frame(row_id = rownames(f$data))
    nd$row_id <- keys$row_id[rows]
    for (newdata in list(nd, nd[1L, , drop = FALSE])) {
      for (zero in c(FALSE, TRUE)) {
        controls <- if (zero) list(re.form = NA) else list()
        native <- do.call(predict, c(list(object = model, newdata = newdata,
          type = "response", se.fit = TRUE), controls))
        actual <- do.call(predict, c(list(object = f$object, newdata = newdata,
          kernel_inputs = f$inputs, type = "response", se.fit = TRUE,
          uncertainty = "conditional", join_by = keys), controls))
        expect_equal(as.numeric(actual$fit), as.numeric(native$fit), tolerance = 1e-8)
        expect_equal(as.numeric(actual$se.fit), as.numeric(native$se.fit), tolerance = 1e-8)
        contract <- attr(actual, "prediction_contract")
        expect_identical(contract$response_spatial_field$component, "conditional_response")
        expect_false(contract$response_spatial_field$kernel_uncertainty_included)
        expect_identical(contract$response_spatial_field$conditioning, if (zero)
          "all_response_random_effects_set_to_zero" else "native_conditional_spatial_prediction")
        expect_identical(contract$integration$random_effects,
          if (zero) "set_to_zero_not_integrated" else "native_method")
        expect_true("observation_noise" %in% contract$uncertainty$omitted)
      }
    }
    snapshot <- prepare_prediction(f$object, nd, f$inputs, join_by = keys)
    restored <- unserialize(serialize(snapshot, NULL, version = 2))
    expect_null(restored$opt_mod)
    expect_equal(as.numeric(predict(restored)), as.numeric(predict(f$object, nd, f$inputs, join_by = keys)), tolerance = 1e-8)
    expect_identical(restored$prediction_contract$response$response_spatial_field,
      attr(predict(f$object, nd, f$inputs, join_by = keys), "prediction_contract")$response_spatial_field)
    novel <- nd
    novel$field_pos <- glmmTMB::numFactor(novel$x + 7, novel$y_coordinate + 3)
    novel_inputs <- kernel_prep(terra::vect(cbind(novel$x + 7, novel$y_coordinate + 3),
      crs = terra::crs(f$points)), f$rasters, max_D = 40, kernel = "gaussian",
      scale_vars = msr_vars(mean = kernel_var("hab"),
        rough = surface_var("elev", "sq", weighted = TRUE),
        fixed = surface_var("elev", "sa", radius = 25), unused = kernel_var("hab")),
      bin = FALSE, verbose = FALSE)
    rownames(novel) <- rownames(novel_inputs$kernel_dat)
    # The shared transform is already independently tested across kernels.
    # Here it supplies matching covariates at the new physical field positions.
    novel[c("mean", "rough")] <- as.data.frame(
      .msr_refit_spatial_values(f$object, novel_inputs, c("mean", "rough"))$value)
    expected <- predict(model, novel, type = "response", se.fit = TRUE, allow.new.levels = TRUE)
    actual <- predict(f$object, novel, novel_inputs, se.fit = TRUE,
      uncertainty = "conditional", allow.new.levels = TRUE)
    expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-8)
    expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-8)
    expect_identical(attr(actual, "prediction_contract")$response_spatial_field$new_locations,
      "native_new_location_handling_requested")
    missing <- nd; missing$field_pos <- NULL
    expect_error(predict(f$object, missing, f$inputs, join_by = keys), "explicitly contain.*field_pos")
    expect_error(predict(f$object, nd, f$inputs, se.fit = TRUE, join_by = keys), "conditional")
    context <- build_opt_context(model, f$inputs$raw_cov, scale_vars = f$inputs$scale_vars,
      unit_conv = f$inputs$unit_conv, resolution = f$inputs$resolution,
      n_cols = f$inputs$n_cols, binned = f$inputs$binned, d_list = f$inputs$d_list)
    final <- kernel_scale_fn(f$par, f$inputs$d_list, f$inputs$raw_cov, f$inputs$kernel,
      model, opt_context = context, mod_return = TRUE)
    expect_identical(final$mod$modelInfo$reStruc$condReStruc, model$modelInfo$reStruc$condReStruc)
    expect_true(final$mod$sdr$pdHess)
    expect_equal(as.numeric(logLik(final$mod)), as.numeric(logLik(model)), tolerance = 1e-5)
    expect_equal(.msr_parameter_count(final$mod), .msr_parameter_count(model))
    expect_equal(final$scl_params$mean[c("mean", "rough")], f$object$scl_params$mean[c("mean", "rough")])
    expect_equal(final$scl_params$sd[c("mean", "rough")], f$object$scl_params$sd[c("mean", "rough")])
    if (family == "poisson") {
      # Supply the joint backend's required binned mean preparation before
      # asking its unchanged random-structure guard to reject the field.
      joint_inputs <- kernel_prep(f$points, f$rasters, max_D = 40, kernel = "gaussian",
        scale_vars = msr_vars(mean = kernel_var("hab"),
          rough = surface_var("elev", "sq", weighted = TRUE),
          fixed = surface_var("elev", "sa", radius = 25), unused = kernel_var("hab")),
        bin = TRUE, verbose = FALSE)
      expect_error(multiScale_optim(model, joint_inputs,
        backend = "joint_mixed"), "One conditional random intercept")
    }
  }
  cat("NATIVE_RESPONSE_FIELD families=2 spatial_exp=1 no_year=1 fixed_scale_refits=2 joint_guard=retained\n")
})

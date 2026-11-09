test_that("GDR component points retain original learned bases and native uncertainty", {
  f <- gdr_trained_fixture()
  for (name in names(f$models)) {
    model <- f$models[[name]]
    components <- names(.msr_gdr_active_formulas(model))
    for (type in components) {
      raw <- unmarked:::get_covariates(model, type)
      rows <- raw[c(2, 7, 11, 2), , drop = FALSE]
      predictor <- .msr_unmarked_prediction_function(model, component = type)
      for (controls in list(NULL, NA)) for (back in c(TRUE, FALSE)) {
        expected <- gdr_trained_reference(model, type, rows, back, controls)
        actual <- predictor(model, type, rows, backTransform = back, re.form = controls)
        expect_equal(actual, expected, tolerance = 1e-9, info = paste(name, type))
        one <- predictor(model, type, rows[1L, , drop = FALSE], backTransform = back, re.form = controls)
        expect_equal(one, expected[1L, , drop = FALSE], tolerance = 1e-9, ignore_attr = TRUE)
      }
    }
    expect_identical(model@estimates, f$snapshots[[name]]$estimates)
    expect_identical(model@opt, f$snapshots[[name]]$opt)
    expect_identical(model@data, f$snapshots[[name]]$data)
    expect_identical(model@formlist, f$snapshots[[name]]$formlist)
    expect_identical(model@AIC, f$snapshots[[name]]$AIC)
  }
  expect_identical(body(unmarked:::make_mod_matrix), f$provider)
})

test_that("GDR native frames preserve site, primary and removal row layouts", {
  f <- gdr_trained_fixture()
  for (name in names(f$models)) {
    model <- f$models[[name]]
    frame <- model@data[c(9, 2, 15), ]
    predictor <- .msr_unmarked_prediction_function(model)
    for (type in names(.msr_gdr_active_formulas(model))) {
      rows <- unmarked:::get_covariates(frame, type)
      for (controls in list(NULL, NA)) {
        expected <- gdr_trained_reference(model, type, rows, re.form = controls)
        actual <- predictor(model, type, frame, re.form = controls)
        expect_equal(actual, expected, tolerance = 1e-9, info = paste(name, type))
      }
    }
  }
})

test_that("GDR points, snapshots and raster chunks share component recipes", {
  f <- gdr_trained_fixture()
  for (name in names(f$models)) {
    model <- f$models[[name]]
    object <- structure(list(opt_mod = model,
      scl_params = list(mean = c(hab = 0), sd = c(hab = 1))), class = "multiScaleR")
    expect_equal(.msr_parameter_count(model), if (name == "mixed") 10L else 18L)
    for (type in names(.msr_gdr_active_formulas(model))) {
      rows <- unmarked:::get_covariates(model, type)[c(3, 8, 12, 3), , drop = FALSE]
      expected <- gdr_trained_reference(model, type, rows, re.form = NA)
      # These columns are already in native fitted units. Spatial predict()
      # keeps its separate matching-kernel_inputs guard; snapshots expose the
      # public route for already prepared native rows.
      point <- .msr_native_prediction_rows(object, rows, type, TRUE, list(re.form = NA))
      expect_equal(as.matrix(point), as.matrix(expected), tolerance = 1e-9, ignore_attr = TRUE)
      snapshot <- prepare_prediction(object, rows, component = type, native_args = list(re.form = NA))
      restored <- unserialize(serialize(snapshot, NULL))
      expect_equal(predict(restored), expected, tolerance = 1e-9, ignore_attr = TRUE)
      expect_lt(as.numeric(object.size(restored)), 128 * 1024)
      raster_rows <- rows
      logical_names <- .msr_native_raster_logical_names(model, type)
      for (variable in intersect(logical_names, names(raster_rows)))
        raster_rows[[variable]] <- as.numeric(raster_rows[[variable]])
      thin <- list(object = list(opt_mod = model), type = type, se.fit = TRUE,
        native_args = list(re.form = NA), logical_names = logical_names,
        unmarked_predictor = .msr_unmarked_prediction_function(model, re.form = NA, component = type))
      chunks <- lapply(list(1L, 2:4), function(i) .msr_native_raster_eval(thin, raster_rows[i, , drop = FALSE]))
      expect_equal(as.matrix(do.call(rbind, chunks)), as.matrix(expected), tolerance = 1e-9, ignore_attr = TRUE)
      contract <- snapshot$prediction_contract[[type]]
      expect_identical(list(component = contract$component, response = contract$effects$response,
        integration = contract$integration$random_effects, kernel = contract$uncertainty$kernel_parameters,
        noise = contract$uncertainty$observation_noise),
        list(component = type, response = "zero_deviation_reference",
          integration = "set_to_zero_not_integrated", kernel = "fixed_at_estimate", noise = FALSE))
      if (type == "phi") {
        expect_true("wind_good" %in% logical_names)
        bad <- raster_rows; bad$wind_good[1L] <- .2
        expect_error(.msr_native_raster_eval(thin, bad), "logical raster predictor")
      }
    }
  }
})

test_that("GDR effect curves use active component rows and raw modeled values", {
  f <- gdr_trained_fixture()
  old <- options("contrasts"); on.exit(options(old), add = TRUE)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (name in names(f$models)) {
    model <- f$models[[name]]
    object <- structure(list(opt_mod = model,
      scl_params = list(mean = c(hab = 0), sd = c(hab = 1))), class = "multiScaleR")
    for (type in names(.msr_gdr_active_formulas(model))) {
      raw <- unmarked:::get_covariates(model, type)
      variables <- switch(type, lambda = if (name == "mixed") "hab" else c("hab", "land cover"),
        phi = c("wind", "wind_good"), dist = "temperature", rem = if (name == "mixed") character() else "effort")
      defaults <- lapply(raw[variables], function(x) if (is.numeric(x)) mean(x) else x[1L])
      options(contrasts = c("contr.helmert", "contr.poly"))
      plots <- plot_marginal_effects(object, type = type, length.out = 5L)
      expect_setequal(names(plots), variables)
      for (variable in variables) {
        x <- raw[[variable]]
        values <- if (is.numeric(x)) seq(min(x), max(x), length.out = 5L) else
          if (is.factor(x)) x[match(levels(x), as.character(x))] else unique(x)
        rows <- as.data.frame(lapply(defaults, rep, times = length(values)), check.names = FALSE)
        rows[[variable]] <- values
        expected <- gdr_trained_reference(model, type, rows, re.form = NA)
        expect_equal(unname(as.matrix(plots[[variable]]$data[c("fit", "lwr", "upr")])),
          unname(as.matrix(expected[c("Predicted", "lower", "upper")])), tolerance = 1e-8)
        expect_equal(plots[[variable]]$data$x, values, ignore_attr = TRUE)
      }
      contract <- attr(plots, "prediction_contract")
      expect_identical(c(contract$component, contract$integration$random_effects,
        contract$uncertainty$kernel_parameters), c(type, "set_to_zero_not_integrated", "fixed_at_estimate"))
    }
    expect_identical(model@estimates, f$snapshots[[name]]$estimates)
    expect_identical(model@AIC, f$snapshots[[name]]$AIC)
  }
  model <- f$models$mixed
  object <- structure(list(opt_mod = model, scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
  expect_identical(attr(plot_marginal_effects(object, length.out = 3L), "prediction_contract")$component, "lambda")
  expect_error(plot_marginal_effects(object, type = "phi"), "active native GDR component")
})

test_that("GDR scopes and guards preserve unrelated and unavailable metadata", {
  f <- gdr_trained_fixture()
  model <- f$models$TMB
  rows <- unmarked:::get_covariates(model, "lambda")[c(2, 9), c("hab", "land cover")]
  reduced <- model; reduced@data@obsCovs$effort <- NULL
  scoped <- .msr_unmarked_prediction_function(reduced, component = "lambda")
  expect_equal(scoped(reduced, "lambda", rows), gdr_trained_reference(model, "lambda", rows), tolerance = 1e-9)
  expect_error(.msr_unmarked_prediction_function(reduced), "Original unmarked covariates")
  bad <- model; bad@data@siteCovs$hab[1L] <- bad@data@siteCovs$hab[1L] + .1
  expect_error(.msr_unmarked_prediction_function(bad, component = "lambda"), "do not match")
  expect_error(scoped(reduced, "phi", rows), "component-specific")
  expect_error(scoped(reduced, "lambda", reduced@data), "component-specific")
  wrong_level <- rows; wrong_level[["land cover"]] <- "Unknown"
  expect_error(scoped(reduced, "lambda", wrong_level), "new level|factor")
  bad <- model; bad@formlist$phi <- ~wind + offset(wind)
  expect_error(.msr_unmarked_prediction_function(bad, component = "lambda"), "offsets")
  bad <- model; bad@mixture <- "ZIP"
  expect_error(.msr_unmarked_prediction_function(bad), "TMB ZIP")
  mixed <- f$models$mixed
  metadata <- mixed; metadata@formlist$lambda <- ~splines::ns(hab, df = 3) + (1 + scale(hab) | block)
  expect_error(.msr_unmarked_prediction_function(metadata, component = "lambda"), "Learned random-effect")
  zero <- .msr_unmarked_prediction_function(metadata, re.form = NA, component = "lambda")
  points <- unmarked:::get_covariates(mixed, "lambda")[c(2, 9), "hab", drop = FALSE]
  expect_equal(zero(metadata, "lambda", points), gdr_trained_reference(mixed, "lambda", points, re.form = NA), tolerance = 1e-9)
  expect_error(zero(metadata, "lambda", points, re.form = NULL), "zero-effect")
})

test_that("GDR frame missing-input checks and preparation layouts fail explicitly", {
  model <- gdr_trained_fixture()$models$C
  frame <- model@data[c(2, 9), ]
  predictor <- .msr_unmarked_prediction_function(model)
  bad <- frame; bad@yearlySiteCovs$wind[1L] <- NA_real_
  expect_error(predictor(model, "phi", bad), "Missing yearlySiteCovs")
  bad@yDistance[1L, 1L] <- bad@yRemoval[1L, 1L] <- NA_real_
  rows <- unmarked:::get_covariates(bad, "phi")
  expect_equal(predictor(model, "phi", bad), gdr_trained_reference(model, "phi", rows), tolerance = 1e-9)
  original <- .msr_unmarked_local_method
  with_mocked_bindings({
    expect_error(.msr_unmarked_gdr_frame_inputs(model, "lambda", frame, TRUE, NA,
      list(), list(), model@formlist), "layout|evaluation order")
  }, .msr_unmarked_local_method = function(...) function() NULL)
})

test_that("trained GDR chunk contexts reuse checked immutable original covariates", {
  model <- gdr_trained_fixture()$models$TMB
  original <- get("clean_up_covs", asNamespace("unmarked")); count <- 0L
  testthat::with_mocked_bindings({
    predictor <- .msr_unmarked_prediction_function(model, re.form = NA, component = "dist")
    before <- count
    rows <- model@data@yearlySiteCovs[c(2, 9), "temperature", drop = FALSE]
    for (i in seq_len(3L)) expect_true(all(is.finite(predictor(model, "dist", rows)$Predicted)))
    expect_equal(count, before)
    expect_gt(before, 0L)
  }, clean_up_covs = function(...) {count <<- count + 1L; original(...)}, .package = "unmarked")
})

test_that("trained GDR component contexts restore for a native prediction worker", {
  model <- gdr_trained_fixture()$models$mixed
  rows <- unmarked:::get_covariates(model, "dist")[c(2, 9, 12), "temperature", drop = FALSE]
  thin <- list(object = list(opt_mod = model), type = "dist", se.fit = TRUE,
    native_args = list(re.form = NA), logical_names = character(),
    unmarked_predictor = .msr_unmarked_prediction_function(model, re.form = NA, component = "dist"))
  workers <- parallel::makePSOCKcluster(1L)
  on.exit(parallel::stopCluster(workers), add = TRUE)
  cluster_prep(model, workers)
  result <- parallel::clusterCall(workers, function(thin, rows)
    multiScaleR:::.msr_native_raster_eval(thin, rows), thin, rows)[[1L]]
  expect_equal(as.matrix(result), as.matrix(gdr_trained_reference(model, "dist", rows, re.form = NA)), tolerance = 1e-9, ignore_attr = TRUE)
  cat("GDR_TRAINED_PREDICTION components=4 C_TMB=retained rows_and_recipes=verified covariance=native\n")
})

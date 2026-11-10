test_that("native GDR field fits retain normalized targets, counts and every start", {
  f <- gdr_native_field_fixture(); model <- f$model; field <- model@response_field
  expect_s4_class(model, "unmarkedFitGDRField")
  expect_s4_class(model, "unmarkedFitGDR")
  expect_equal(.msr_parameter_count(model), 6L)
  expect_equal(model@AIC, -2 * as.numeric(unmarked::logLik(model)) + 12)
  expect_identical(field$diagnostics$status, "available")
  expect_true(field$diagnostics$joint_precision_positive_definite)
  expect_true(field$diagnostics$conditional_latent_precision_positive_definite)
  expect_lt(field$diagnostics$maximum_absolute_gradient, 1e-3)
  audit <- readRDS(field$audit_file)
  expect_identical(audit$status, "returned")
  expect_length(audit$attempts, 2L)
  expect_true(all(vapply(audit$attempts, function(x) x$eligible, logical(1))))
  expect_identical(audit$annual_field_coordinates, 0L)
  expect_false(audit$range_estimated)
  expect_lt(audit$maximum_start_value_difference, 1e-5)
  expect_equal(gdr_start_values(model)$coordinate_count, 6L)
  index <- which(names(field$point) == "b_lambda")
  physical <- drop(field$lower %*% field$point[index])
  expect_equal(unname(model@estimates@estimates$lambda@estimates[-(1:2)]), unname(physical), tolerance = 1e-12)
  table <- aic_tab(list(iid = f$template$model, spatial = model), AICc = FALSE)
  expect_equal(table$K, c(6, 6))
  expect_match(.msr_selection_family_label(model), "unmarkedFitGDRField", fixed = TRUE)
  # Independent prepared joint target at this actual returned optimizer point.
  gdr_gaussian_require_RTMB()
  reference <- gdr_response_field_fit_fixture()
  context <- .msr_gdr_response_field_context(reference$data, reference$pars,
    reference$coords, 32617, 50000)
  program <- context$prepare(reference$data$Xlambda, reference$data$row_ids$survey)
  p <- model@TMB$env$parList(par = field$point)
  actual <- RTMB::MakeADFun(program$evaluate, p, random = NULL, silent = TRUE)
  native <- TMB::MakeADFun(model@TMB$env$data, p, random = NULL,
    DLL = "unmarked_TMBExports", silent = TRUE)
  key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
  map <- match(key(names(actual$par)), key(names(native$par)))
  expect_false(anyNA(map))
  expect_lt(abs(actual$fn(actual$par) - native$fn(native$par)), 1e-8)
  expect_lt(max(abs(actual$gr(actual$par) - native$gr(native$par)[map])), 1e-7)
  expect_false(context$fit_admission_established)
  cat("GDR_NATIVE_FIELD starts=2 fixed=6 latent=12 native_target=1 RTMB_private_guard=retained\n")
})

test_that("field predictions retain native known-state and zero-effect uncertainty", {
  f <- gdr_native_field_fixture(); rows <- c(3, 45, 3, 159)
  nd <- f$model@data@siteCovs[rows, , drop = FALSE]
  native <- methods::as(f$model, "unmarkedFitGDR")
  for (type in c("lambda", "dist", "rem")) for (zero in c(FALSE, TRUE)) {
    if (type == "rem") {
      nd <- get("get_covariates", asNamespace("unmarked"))(f$model, type)[rows, , drop = FALSE]
    } else nd <- f$model@data@siteCovs[rows, , drop = FALSE]
    controls <- if (zero) list(re.form = NA) else list()
    expected <- do.call(unmarked::predict, c(list(object = native, type = type,
      newdata = nd, backTransform = FALSE), controls))
    actual <- do.call(unmarked::predict, c(list(object = f$model, type = type,
      newdata = nd, backTransform = FALSE), controls))
    for (name in c("Predicted", "SE", "lower", "upper"))
      expect_equal(actual[[name]], expected[[name]], tolerance = 1e-7)
    info <- attr(actual, "response_field_uncertainty")
    expect_equal(info$new_field_residual_predictor_variance, rep(0, length(rows)))
    expect_false(info$SE_includes_new_field_residual)
    expect_true("ecological_scale_uncertainty" %in% info$omitted)
  }
  zero <- unmarked::predict(f$model, type = "lambda", newdata = nd["hab"], re.form = NA)
  expect_true(all(is.finite(zero$SE)))
  expect_identical(attr(zero, "response_field_uncertainty")$point, "zero_effect_reference_not_integrated")
})

test_that("new-location density uncertainty includes distinct field residual variance", {
  f <- gdr_native_field_fixture(); field <- f$model@response_field
  xy <- rbind(f$nodes[2, ], c(490000, 4405000), f$nodes[2, ], c(620000, 4395000))
  nd <- data.frame(hab = c(-.5, 0, .5, 1))
  actual <- unmarked::predict(f$model, type = "lambda", newdata = nd,
    coords = xy, crs = 32617, backTransform = FALSE)
  R <- exp(-as.matrix(dist(f$nodes))/50000)
  K <- exp(-sqrt(outer(xy[, 1], f$nodes[, 1], "-")^2 +
    outer(xy[, 2], f$nodes[, 2], "-")^2)/50000)
  weights <- K %*% solve(R, field$lower)
  G <- matrix(0, 4L, length(field$point))
  G[, names(field$point) == "beta_lambda"] <- cbind(1, nd$hab)
  G[, names(field$point) == "b_lambda"] <- weights
  predicted <- drop(G %*% field$point)
  fitted <- diag(G %*% field$joint_covariance %*% t(G))
  new <- exp(2 * field$point[names(field$point) == "lsigma_lambda"]) *
    pmax(0, 1 - rowSums(K * t(solve(R, t(K)))))
  expect_equal(actual$Predicted, predicted, tolerance = 1e-10)
  expect_equal(actual$SE^2, fitted + new, tolerance = 1e-10)
  info <- attr(actual, "response_field_uncertainty")
  expect_equal(info$fitted_joint_predictor_variance, fitted, tolerance = 1e-10)
  expect_equal(info$new_field_residual_predictor_variance, new, tolerance = 1e-10)
  expect_true(info$SE_includes_new_field_residual)
  expect_gt(min(new[c(2, 4)]), 0)
  permuted <- unmarked::predict(f$model, type = "lambda", newdata = nd[4:1, , drop = FALSE],
    coords = xy[4:1, ], crs = 32617, backTransform = FALSE)
  expect_equal(permuted$Predicted, rev(actual$Predicted), tolerance = 0)
  expect_equal(permuted$SE, rev(actual$SE), tolerance = 0)
  expect_error(unmarked::predict(f$model, type = "lambda", newdata = nd, coords = xy), "declared crs")
  expect_error(unmarked::predict(f$model, type = "lambda", newdata = nd, coords = xy, crs = 4326), "projected CRS")
})

test_that("native GDR fields integrate actual scale refits, snapshots, curves and maps", {
  f <- gdr_native_field_projection_fixture(); object <- f$object
  expect_s4_class(object$opt_mod, "unmarkedFitGDRField")
  expect_equal(object$opt_mod@data@siteCovs$hab, f$field$model@data@siteCovs$hab, tolerance = 1e-10)
  expect_lt(abs(as.numeric(unmarked::logLik(object$opt_mod)) -
    as.numeric(unmarked::logLik(f$field$model))), 1e-5)
  expect_identical(object$opt_mod@response_field$field_id, f$field$model@response_field$field_id)
  expect_equal(aic_tab(list(field = object), AICc = FALSE)$K, 7L)
  expect_identical(f$result$scl_params, list(mean = c(hab = 0), sd = c(hab = 1)))
  nd <- object$opt_mod@data@siteCovs[c(1, 37, 145), , drop = FALSE]
  expected <- unmarked::predict(object$opt_mod, type = "lambda", newdata = nd)
  snapshot <- prepare_prediction(object, nd, component = "lambda")
  restored <- unserialize(serialize(snapshot, NULL, version = 2L))
  expect_null(restored$opt_mod)
  actual <- predict(restored)
  expect_equal(actual$Predicted, expected$Predicted, tolerance = 1e-10)
  expect_equal(actual$SE, expected$SE, tolerance = 1e-10)
  contract <- restored$prediction_contract$lambda$response_spatial_field
  expect_identical(contract$provider, "native_unmarked_TMB")
  expect_false(contract$range_estimated)
  expect_false(contract$kernel_uncertainty_included)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot_marginal_effects(object, type = "lambda", length.out = 5L)
  expect_true("hab" %in% names(curves))
  expect_identical(attr(curves, "prediction_contract")$response_spatial_field$conditioning,
    "all_response_random_effects_set_to_zero")
  whole <- kernel_scale.raster(f$raster, multiScaleR = object, scale_center = TRUE,
    include_fixed = FALSE, match_buffer = TRUE, verbose = FALSE)
  expect_identical(names(whole), "hab")
  raster <- terra::crop(whole, terra::ext(426000, 434000, 4397000, 4403000))
  attr(raster, "multiScaleR_projection") <- attr(whole, "multiScaleR_projection")
  geometry <- response_field_geometry(object, f$field$coords, 32617, batch_size = 3L)
  other <- object; other$backend <- "profile_glm"
  expect_error(response_field_geometry(other, f$field$coords, 32617), "GDR response-field.*refit")
  data <- terra::as.data.frame(raster); xy <- terra::xyFromCell(raster, seq_len(terra::ncell(raster)))
  expected <- unmarked::predict(object$opt_mod, type = "lambda", newdata = data,
    coords = xy, crs = 32617)
  actual <- predict(object, raster, type = "lambda", se.fit = TRUE, uncertainty = "conditional",
    raster_args = list(field_geometry = geometry))
  expect_equal(as.numeric(terra::values(actual)[, 1]), expected$Predicted, tolerance = 1e-8)
  expect_equal(as.numeric(terra::values(actual)[, 2]), expected$SE, tolerance = 1e-8)
  expect_true(any(attr(expected, "response_field_uncertainty")$new_field_residual_predictor_variance > 0))
  expect_error(predict(object, raster, type = "lambda"), "response_field_geometry")
  expect_error(predict(object, raster, type = "lambda", raster_args = list(field_geometry = geometry, cores = 2L)), "cores = 1")
  expect_error(predict(object, raster, type = "lambda", native_args = list(re.form = NA),
    raster_args = list(cores = 2L)), "cores = 1.*zero-effect")
  expect_error(predict(object, raster, type = "dist",
    raster_args = list(field_geometry = geometry)), "lambda, state")
  cat("GDR_FIELD_INTEGRATION kernel_refit=1 snapshot=1 zero_curve=1 coordinate_raster=1 fixed_range_count=0\n")
})

test_that("field target mutations and unqualified inherited routes remain guarded", {
  f <- gdr_native_field_fixture(); model <- f$model
  expect_error(unmarked::update(model), "make_gdr_refitter")
  expect_error(unmarked::simulate(model), "separate field qualification")
  expect_error(unmarked::profile(model), "separate field qualification")
  expect_error(unmarked::confint(model, type = "lambda", method = "profile"), "separate field qualification")
  expect_error(unmarked::ranef(model), "separate spatial GDR qualification")
  expect_error(build_opt_context(model, list(hab = 1)), "make_gdr_refitter")
  bad <- model; bad@response_field$point[1L] <- bad@response_field$point[1L] + .1
  expect_error(unmarked::predict(bad), "context changed")
  bad <- model; bad@K <- bad@K + 1
  expect_error(aic_tab(list(bad)), "context changed")
  bad <- model; bad@response_field$settings$threads <- 2L
  expect_error(unmarked::predict(bad), "context changed: configuration")
  bad <- model; bad@response_field$diagnostics$status <- "unusable"
  expect_error(unmarked::vcov(bad), "context changed: configuration")
  expect_error(response_field_geometry(gdr_native_field_projection_fixture()$object,
    f$coords[160:1, ], 32617), "exact fitted")
  callback <- make_gdr_refitter(model, f$template$starts, gdr_test_audit_path("field-invariants-"),
    tmb_settings = f$template$settings)
  frame <- model@data; frame@yDistance[1, 1] <- frame@yDistance[1, 1] + 1
  expect_error(callback(model, frame, list(covs = "hab")), "histories")
  expect_error(make_gdr_refitter(model, f$template$starts, tempfile(),
    tmb_settings = within(f$template$settings, threads <- 2L)), "retain the original")
})

test_that("retained physical covariance and numeric workers need no saved tape", {
  f <- gdr_native_field_fixture(); model <- f$model; field <- model@response_field
  fixed <- which(startsWith(names(field$point), "beta_"))
  actual <- unmarked::vcov(model)
  expect_equal(unname(actual), unname(field$joint_covariance[fixed, fixed]), tolerance = 1e-12)
  expect_identical(rownames(actual), names(unmarked::coef(model)))
  expect_equal(unmarked::SE(model), sqrt(diag(actual)), tolerance = 1e-12)
  physical <- unmarked::vcov(model, type = "lambda", fixedOnly = FALSE)
  expect_equal(unname(physical), unname(model@estimates@estimates$lambda@covMat), tolerance = 1e-12)
  expect_error(unmarked::vcov(model, method = "nonparboot"), "bootstrap needs separate")
  expect_true(all(vapply(field$bases, function(x) identical(environment(x$terms), baseenv()), logical(1))))
  expect_true(all(vapply(field$forms, function(x) identical(environment(x), baseenv()), logical(1))))
  expect_lt(length(serialize(field, NULL)), 1024^2)
  nd <- model@data@siteCovs[c(2, 50, 151), , drop = FALSE]
  expected <- unmarked::predict(model, type = "lambda", newdata = nd)
  cluster <- parallel::makePSOCKcluster(1L); on.exit(parallel::stopCluster(cluster), add = TRUE)
  installed <- find.package("multiScaleR")
  remote <- parallel::clusterCall(cluster, function(field, data, libraries, installed) {
    .libPaths(libraries); library(multiScaleR)
    stopifnot(identical(normalizePath(find.package("multiScaleR"), winslash = "/"),
      normalizePath(installed, winslash = "/")))
    get(".msr_gdr_field_predict", asNamespace("multiScaleR"))(field, field$forms,
      data, "lambda", TRUE, TRUE, FALSE, .95, NULL, NULL, NULL)
  }, field, nd, .libPaths(), installed)[[1L]]
  expect_equal(remote$Predicted, expected$Predicted, tolerance = 0)
  expect_equal(remote$SE, expected$SE, tolerance = 0)
})

test_that("native fields retain the no-observer target and its five-parameter penalty", {
  f <- gdr_native_field_fixture(); args <- f$template$constructor_args
  map <- gdr_start_values(f$template$model)$map
  keep <- map$block != "lsigma_dist"
  args$distanceformula <- stats::as.formula("~1", env = baseenv())
  args$starts <- args$starts[keep]
  path <- paste0(gdr_test_audit_path("field-no-observer-template-"), ".rds")
  template <- tryCatch(.msr_gdr_TMB_refit(args,
    list(dots = list(control = f$template$settings$control))), error = function(e) {
      .msr_gdr_atomic_audit(list(status = "failed", error = conditionMessage(e),
        optimizer = e$optimizer, covariance = e$covariance), path)
      stop(e)
    })
  .msr_gdr_atomic_audit(list(status = "returned", optimizer = .msr_gdr_optimizer_checkpoint(template@opt),
    diagnostics = attr(template, "msr_gdr_tmb_diagnostics")), path)
  fit <- fit_gdr_response_field(template, f$coords, 32617, 50000,
    list(cold = f$template$starts$cold[keep]), gdr_test_audit_path("field-no-observer-"), f$template$settings)
  expect_equal(.msr_parameter_count(fit), 5L)
  expect_identical(fit@response_field$diagnostics$latent_coordinates, 8L)
  expect_true(fit@response_field$diagnostics$joint_precision_positive_definite)
  expect_identical(fit@response_field$diagnostics$status, "available")
  expect_equal(aic_tab(list(fit), AICc = FALSE)$K, 5L)
  expect_true(all(is.finite(unmarked::predict(fit, type = "dist")$SE)))
})

test_that("field provenance binds code independently of byte compilation", {
  f <- gdr_native_field_fixture()
  before <- .msr_gdr_field_source()
  compiled <- compiler::cmpfun(.msr_gdr_field_fit)
  expect_identical(with_mocked_bindings(.msr_gdr_field_source(),
    .msr_gdr_field_fit = compiled), before)
  changed <- .msr_gdr_field_fit
  body(changed) <- substitute({ ORIGINAL; invisible(NULL) }, list(ORIGINAL = body(changed)))
  expect_false(identical(with_mocked_bindings(.msr_gdr_field_source(),
    .msr_gdr_field_fit = changed), before))
  expect_error(with_mocked_bindings(unmarked::predict(f$model),
    .msr_gdr_field_fit = changed), "context changed: source")
  bad <- f$template$model; bad@data@y[1, 1] <- bad@data@y[1, 1] + 1
  expect_error(fit_gdr_response_field(bad, f$coords, 32617, 50000,
    f$template$starts, gdr_test_audit_path("field-mirror-"), f$template$settings),
    "inherited removal mirror")
  native <- methods::as(f$model, "unmarkedFitGDR")
  expect_equal(unmarked::confint(f$model, type = "lambda", method = "normal"),
    unmarked::confint(native, type = "lambda", method = "normal"), tolerance = 1e-12)
})

test_that("field geometry binds captured inputs independently of compilation", {
  f <- gdr_native_field_fixture(); field <- f$model@response_field
  probe <- field; probe$project <- compiler::cmpfun(field$project)
  expect_identical(.msr_gdr_field_geometry_hash(probe), field$geometry_sha256)
  bad <- f$model; project <- field$project; original <- environment(project)
  scope <- list2env(as.list(original, all.names = TRUE), parent = parent.env(original))
  budget <- scope$budget; environment(budget) <- scope; scope$budget <- budget
  environment(project) <- scope; scope$range <- scope$range + 1
  bad@response_field$project <- project
  expect_error(unmarked::predict(bad), "context changed: geometry")
  expect_identical(.msr_gdr_field_geometry_hash(field), field$geometry_sha256)
})

test_that("field callbacks isolate numeric recipes and guard fit-worker transport", {
  f <- gdr_native_field_fixture(); model <- f$model
  predictor <- .msr_gdr_field_prediction_function(model)
  expect_setequal(ls(environment(predictor), all.names = TRUE), c("field", "forms"))
  expect_identical(parent.env(environment(predictor)), asNamespace("multiScaleR"))
  expect_lt(length(serialize(predictor, NULL)), 1024^2)
  callback <- make_gdr_refitter(model, f$template$starts,
    gdr_test_audit_path("field-isolated-callback-"), tmb_settings = f$template$settings)
  expect_setequal(ls(environment(callback), all.names = TRUE),
    c("field", "cold_starts", "audit_dir", "original_frame", "original_key", "provenance"))
  expect_identical(parent.env(environment(callback)), asNamespace("multiScaleR"))
  expect_lt(length(serialize(callback, NULL)), 2 * 1024^2)
  expect_error(multiScale_optim(model, list(), n_cores = 2), "fit/tape worker transport")
  expect_error(multiScale_optim(model, list(), PSOCK = TRUE), "fit/tape worker transport")
  expect_error(multiScale_optim(model, list(), backend = "profile_glm"), "fit/tape worker transport")
})

test_that("coordinate field batches preserve logical raster column contracts", {
  f <- gdr_native_field_projection_fixture(); object <- f$object
  geometry <- response_field_geometry(object, f$field$coords, 32617, batch_size = 2L)
  xy <- rbind(f$field$nodes[2, ], c(490000, 4405000))
  data <- data.frame(hab = c(0, .5), wet = c(0, 1),
    .msr_field_x = xy[, 1], .msr_field_y = xy[, 2])
  predictor <- .msr_gdr_field_prediction_function(object$opt_mod)
  observer <- function(object, type, newdata, ...) {
    expect_true(is.logical(newdata$wet))
    predictor(object, type, newdata, ...)
  }
  thin <- list(field_geometry = geometry, object = object, type = "lambda", se.fit = TRUE,
    native_args = list(), logical_names = "wet", unmarked_predictor = observer)
  expected <- unmarked::predict(object$opt_mod, type = "lambda", newdata = data,
    coords = xy, crs = 32617)
  actual <- .msr_gdr_field_raster_eval(thin, data)
  expect_equal(actual$fit, expected$Predicted, tolerance = 0)
  expect_equal(actual$se.fit, expected$SE, tolerance = 0)
  data$wet[2L] <- 2
  expect_error(.msr_gdr_field_raster_eval(thin, data), "logical raster predictor.*0/1")
})

test_that("S3 covariance dispatch survives namespace-only loading with S4 exports", {
  cluster <- parallel::makePSOCKcluster(1L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  installed <- find.package("multiScaleR")
  result <- parallel::clusterCall(cluster, function(libraries, installed) {
    .libPaths(libraries)
    stopifnot(requireNamespace("multiScaleR", quietly = TRUE))
    stopifnot(identical(normalizePath(find.package("multiScaleR"), winslash = "/"),
      normalizePath(installed, winslash = "/")))
    fit <- stats::lm(y ~ x, data.frame(x = 1:8, y = c(2, 4, 3, 6, 5, 9, 8, 11)))
    expected <- stats::vcov(fit)
    prototype <- structure(list(backend = "profile_glm",
      uncertainty = list(available = FALSE, conditional_vcov = expected)), class = "multiScaleR")
    list(registered = !is.null(utils::getS3method("vcov", "multiScaleR", optional = TRUE)),
      covariance = stats::vcov(prototype, type = "conditional"), expected = expected,
      attached = "package:multiScaleR" %in% search())
  }, .libPaths(), installed)[[1L]]
  expect_true(result$registered)
  expect_equal(result$covariance, result$expected, tolerance = 1e-12)
  expect_false(result$attached)
})

# Regression tests for model-level refusals that name what failed and what to
# do (audit findings X-10, G-7, G-8, B-4, and C-9).

test_that("response_field_geometry treats an empty CRS like an unknown CRS", {
  skip_if_not_installed("sf")
  # terra::crs() returns "" for a sim_rast() raster.
  for (crs in list("", "  ", NA, "not a crs")) {
    err <- tryCatch(.msr_response_field_crs(crs), error = conditionMessage)
    expect_match(err, "known projected CRS", fixed = TRUE)
    expect_no_match(err, "invalid crs", fixed = TRUE)
  }
  expect_s3_class(.msr_response_field_crs("EPSG:32617"), "crs")
})

test_that("field-map preflight names a prediction raster with an empty CRS", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  f <- native_response_field_fixture("poisson")
  geometry <- response_field_geometry(f$object, terra::crds(f$points),
    terra::crs(f$points))
  # An empty CRS string through the public constructor gets the package message.
  expect_error(response_field_geometry(f$object, terra::crds(f$points), ""),
    "known projected CRS")
  raster <- terra::subset(f$rasters, 1L)
  terra::crs(raster) <- ""
  err <- tryCatch(.msr_response_field_preflight(geometry, f$object, raster, 1L, NULL),
    error = conditionMessage)
  expect_match(err, "prediction raster has no CRS", fixed = TRUE)
  expect_no_match(err, "invalid crs", fixed = TRUE)
})

test_that("joint_mixed names the failing model or kernel requirement", {
  for (pkg in c("RTMB", "TMB", "glmmTMB")) skip_if_not_installed(pkg)
  f <- profile_glm_fixture()
  call_backend <- function(model, kernel) {
    inputs <- f$inputs
    inputs$kernel <- kernel
    .msr_fit_mixed_backend(model, inputs, join_by = NULL, par = NULL,
      start_strategy = "single", n_cores = 1L, refit_fn = NULL,
      control = list(), verbose = FALSE, call = NULL)
  }
  err <- tryCatch(call_backend(f$model, "gaussian"), error = conditionMessage)
  expect_match(err, "requires a glmmTMB model, but the model has class `glm`", fixed = TRUE)
  expect_no_match(err, "exponential")
  mixed <- structure(list(), class = "glmmTMB")
  err <- tryCatch(call_backend(mixed, "expow"), error = conditionMessage)
  expect_match(err, "kernel = \"expow\"", fixed = TRUE)
  expect_no_match(err, "glmmTMB model", fixed = TRUE)
})

test_that("analytic backends name the requested backend for a fixed kernel", {
  f <- profile_glm_fixture()
  inputs <- f$inputs
  inputs$kernel <- "fixed"
  for (backend in c("profile_glm", "joint_nb2")) {
    err <- tryCatch(multiScale_optim(f$model, inputs, backend = backend,
      n_cores = 1, verbose = FALSE), error = conditionMessage)
    expect_match(err, paste0("^", backend, " requires a smooth kernel"))
    expect_match(err, "kernel = \"fixed\"", fixed = TRUE)
  }
})

test_that("scale_inference states its full requirement", {
  fix <- make_core_fixture()
  err <- tryCatch(scale_inference(fix$opt), error = conditionMessage)
  expect_match(err, "NB2 with one spatial covariate", fixed = TRUE)
  expect_match(err, "this fit used the refit backend", fixed = TRUE)
  err <- tryCatch(scale_inference(list()), error = conditionMessage)
  expect_match(err, "`object` is not a multiScaleR fit", fixed = TRUE)
})

test_that("functions needing multiScale_optim() fits point Bayesian fits elsewhere", {
  bayes <- structure(list(), class = "multiScaleR_bayes")
  for (fn in c("profile_sigma", "prepare_prediction", "scale_inference")) {
    err <- tryCatch(get(fn)(bayes), error = conditionMessage)
    expect_match(err, paste0("`", fn, "()` works on fits from `multiScale_optim()`"),
      fixed = TRUE)
    expect_match(err, "posterior_scales()", fixed = TRUE)
    expect_match(err, "predict_scales()", fixed = TRUE)
  }
})

test_that("pscl fits name the prediction argument and its valid values", {
  skip_if_not_installed("pscl")
  z <- make_zeroinfl_fixture()
  expect_identical(.msr_native_prediction_types(z$mod),
                   c("response", "prob", "count", "zero"))
  err <- tryCatch(prepare_prediction(z$obj), error = conditionMessage)
  expect_match(err, "Specify `component =`", fixed = TRUE)
  expect_match(err, "\"response\", \"prob\", \"count\", \"zero\"", fixed = TRUE)
  err <- tryCatch(.msr_native_prediction_rows(z$obj, NULL, "link", FALSE),
    error = conditionMessage)
  expect_match(err, "Prediction type \"link\" (`type` in predict()", fixed = TRUE)
  expect_match(err, "\"response\", \"prob\", \"count\", \"zero\"", fixed = TRUE)
  # Valid and partially matched types reach the native method unchanged.
  expect_equal(unname(.msr_native_prediction_rows(z$obj, NULL, "count", FALSE)),
               unname(predict(z$mod, type = "count")), ignore_attr = TRUE)
  expect_equal(unname(.msr_native_prediction_rows(z$obj, NULL, "zer", FALSE)),
               unname(predict(z$mod, type = "zero")), ignore_attr = TRUE)
  # Classes without readable choices keep a generic description.
  expect_match(.msr_native_type_list(structure(list(), class = "no_such_class")),
               "accepted by the native predict")
})

test_that("conditional_logit accepts any constant-time Surv() spelling", {
  skip_on_cran()
  skip_if_not_installed("survival")
  f <- profile_glm_fixture()
  dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  dat$one <- 1
  # The formula environment supplies `dat` and `n` for the time expressions.
  env <- new.env(parent = asNamespace("survival"))
  env$dat <- dat
  env$n <- nrow(dat)
  fit_with <- function(lhs) {
    form <- stats::as.formula(paste(lhs, "~ hab + z + group + strata(choice)"), env = env)
    model <- survival::coxph(form, data = dat, method = "exact", model = TRUE)
    multiScale_optim(model, f$inputs, backend = "conditional_logit", n_cores = 1,
      verbose = FALSE, profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1)))
  }
  reference <- fit_with("Surv(rep(1, 240L), case_)")
  for (lhs in c("Surv(rep(1, nrow(dat)), case_)", "Surv(rep(1, n), case_)",
                "Surv(one, case_)")) {
    fit <- fit_with(lhs)
    expect_equal(fit$scale_est$Mean, reference$scale_est$Mean, tolerance = 1e-10)
    expect_equal(fit$optim_results$value, reference$optim_results$value, tolerance = 1e-10)
    expect_equal(coef(fit$opt_mod), coef(reference$opt_mod), tolerance = 1e-10)
    expect_equal(predict(fit), predict(reference), tolerance = 1e-10)
    expect_equal(predict(fit, newdata = dat, kernel_inputs = f$inputs),
                 predict(reference, newdata = dat, kernel_inputs = f$inputs),
                 tolerance = 1e-10)
  }
  # Time expressions that are not data columns reuse the retained frame.
  expect_identical(fit_with("Surv(rep(1, n), case_)")$profile_fit$contract$native_data_mode,
                   "retained_model_frame")
  # A transformed case indicator is still refused, with the accepted form.
  err <- tryCatch(fit_with("Surv(rep(1, 240L), case_ == 1)"), error = conditionMessage)
  expect_match(err, "bare 0/1 column in the last argument of Surv()", fixed = TRUE)
  expect_match(err, "Surv(rep(1, nrow(data)), case)", fixed = TRUE)
})

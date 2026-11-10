test_that("GDR numerical policies are explicit plain metadata with unchanged defaults", {
  native <- .msr_gdr_TMB_policy()
  expect_identical(native$random_start, "cached_best")
  expect_null(native$inner_control)
  expect_false(native$record_memory)
  expect_length(native$requested, 0L)
  control <- list(maxit = 1000L, tol = 1e-10, grad.tol = 1e-10, step.tol = 1e-12)
  explicit <- .msr_gdr_TMB_policy(list(random.start = "zero", inner.control = control))
  expect_identical(explicit$inner_control, control)
  expect_identical(explicit$random_start, "zero")
  expect_identical(unserialize(serialize(explicit, NULL, version = 2L)), explicit)
  for (value in list("other", 1, expression(stop("never evaluate")), structure("zero", payload = new.env())))
    expect_error(.msr_gdr_TMB_policy(list(random.start = value)), "random.start")
  for (value in list(list(maxit = 10L), modifyList(control, list(maxit = 1.5)),
      modifyList(control, list(tol = Inf)), modifyList(control, list(step.tol = 0))))
    expect_error(.msr_gdr_TMB_policy(list(inner.control = value)), "inner.control")
  bad <- control; attr(bad$tol, "payload") <- new.env()
  expect_error(.msr_gdr_TMB_policy(list(inner.control = bad)), "inner.control")
  expect_error(.msr_gdr_TMB_policy(list(record_memory = NA)), "record_memory")
})

test_that("completed GDR optimizer points survive a declared covariance interruption", {
  f <- gdr_tmb_refit_fixture("P")
  directory <- gdr_test_audit_path("gdr-pre-covariance-")
  callback <- make_gdr_refitter(f$model, f$starts, directory, tmb_settings = f$settings)
  before <- gdr_guard_state(f$model)
  failure <- with_mocked_bindings(tryCatch(callback(f$model, f$data,
    list(return_model = TRUE, covs = "hab")), error = identity),
    sdreport = function(...) stop("declared covariance interruption"), .package = "TMB")
  expect_s3_class(failure, "msr_gdr_refit_error")
  expect_length(failure$audit$attempts, 2L)
  expect_identical(gdr_guard_state(f$model), before)
  contract <- readRDS(file.path(directory, "contract.rds"))
  packets <- list()
  for (i in seq_along(failure$audit$attempts)) {
    attempt <- failure$audit$attempts[[i]]
    expect_match(attempt$error, "declared covariance interruption")
    expect_false(attempt$eligible)
    expect_true(file.exists(attempt$optimizer_point_file))
    expect_true(file.exists(attempt$workflow$covariance_request_file))
    progress <- readRDS(attempt$workflow$covariance_request_file)
    expect_false(progress$inference_eligible)
    expect_identical(progress$point_sha256, digest::digest(file = attempt$optimizer_point_file, algo = "sha256"))
    expect_identical(progress$phase$name, "covariance_started")
    expect_gte(attempt$workflow$phases$covariance_failed$elapsed_seconds, 0)
    point <- readRDS(attempt$optimizer_point_file); packets[[i]] <- point
    expect_s3_class(point, "msr_gdr_optimizer_point")
    expect_false(methods::is(point, "unmarkedFit"))
    expect_false(point$inference_eligible)
    expect_identical(point$covariance_status, "not_started")
    expect_false(point$models_or_tapes_saved)
    expect_false(point$latent_states_saved)
    expect_false(point$covariance_saved)
    expect_equal(point$optimizer, attempt$optimizer, tolerance = 0)
    expect_identical(point$coordinates$position, seq_along(attempt$optimizer$par))
    expect_identical(point$coordinates$native_display_label, names(attempt$optimizer$par))
    expect_gt(anyDuplicated(point$coordinates$native_display_label), 0L)
    expect_identical(point$finite_cutoff$storage, typeof(f$model@TMB$env$data$K))
    expect_identical(point$finite_cutoff$value, f$model@TMB$env$data$K)
    expect_identical(point$binding$contract_sha256, contract$contract_sha256)
    expect_identical(point$binding$input_context_sha256, failure$audit$input_context_sha256)
    expect_identical(point$policy$random_start, "cached_best")
    expect_identical(names(point$phases), c("construction_started", "construction_complete", "optimization_complete"))
    expect_gte(point$native_preparation$elapsed_seconds, 0)
    expect_identical(point$native_preparation$started$name, "native_preparation_started")
    expect_true(all(vapply(point$phases, function(x) is.finite(x$elapsed_seconds) && x$elapsed_seconds >= 0, logical(1))))
    expect_false(any(vapply(point$phases, function(x) x$memory$requested, logical(1))))
    expect_lt(file.info(attempt$optimizer_point_file)$size, 20000)
    plain <- function(x) if (is.list(x)) all(vapply(x, plain, logical(1))) else
      !is.function(x) && !is.environment(x) && typeof(x) != "externalptr"
    expect_true(plain(point))
  }
  expect_identical(packets[[1L]]$native_data_sha256, packets[[2L]]$native_data_sha256)
})

test_that("native zero-inner policy retains joint inference, counts and independent references", {
  f <- gdr_tmb_refit_fixture("mixed")
  settings <- f$settings
  settings$random.start <- "zero"
  settings$inner.control <- list(maxit = 2500L, tol = 1e-10, grad.tol = 1e-10, step.tol = 1e-12)
  settings$record_memory <- requireNamespace("ps", quietly = TRUE)
  directory <- gdr_test_audit_path("gdr-explicit-inner-")
  callback <- make_gdr_refitter(f$model, f$starts, directory, tmb_settings = settings)
  model <- callback(f$model, f$data, list(return_model = TRUE, covs = "hab"))
  audit <- attr(model, "msr_gdr_refit_audit"); contract <- readRDS(audit$contract_file)
  expect_s4_class(model, "unmarkedFitGDR")
  expect_length(audit$attempts, 3L)
  expect_equal(.msr_parameter_count(model), 6L)
  expect_identical(contract$controls$policy$requested, settings[c("random.start", "inner.control", "record_memory")])
  for (attempt in audit$attempts) {
    d <- attempt$covariance
    expect_true(attempt$eligible)
    expect_true(d$conditional_latent_stationary)
    expect_true(d$covariance_at_current_latent_mode)
    expect_true(d$current_fixed_coordinates_at_optimizer)
    expect_lt(d$returned_value_gap, 1e-7)
    expect_lt(d$conditional_latent_maximum_absolute_gradient, 1e-3)
    expect_lt(d$conditional_latent_Newton_decrement, 1e-3)
    expect_false(d$coefficient_hessians_omitted)
    expect_identical(d$workflow$policy$random_start, "zero")
    expect_identical(d$workflow$policy$inner_control, settings$inner.control)
    expect_true(all(vapply(d$workflow$phases, function(p) identical(p$memory$requested, settings$record_memory), logical(1))))
    if (settings$record_memory) {
      expect_true(all(vapply(d$workflow$phases, function(p) p$memory$available, logical(1))))
      expect_true(all(vapply(d$workflow$phases, function(p) p$memory$rss_bytes > 0, logical(1))))
    }
    expect_true(all(vapply(d$workflow$phases, function(p) grepl("not_phase_peak", p$memory$scope, fixed = TRUE), logical(1))))
    expect_true(file.exists(attempt$optimizer_point_file))
  }
  # Fresh native objective at the same fixed point, using the literal policy.
  # This is a software reference, not a claim of unique mode/integral accuracy.
  random <- c("b_lambda", "b_dist")
  model@TMB$gr(model@opt$par)
  current <- model@TMB$env$last.par
  parameters <- model@TMB$env$parList(par = current)
  parameters[random] <- lapply(parameters[random], function(x) rep(0, length(x)))
  native <- TMB::MakeADFun(.msr_gdr_TMB_preparation(model), parameters,
    random = random, random.start = expression(rep(0, length(random))),
    inner.control = settings$inner.control, silent = TRUE, DLL = "unmarked_TMBExports")
  expect_identical(names(native$par), names(model@opt$par))
  expect_lt(abs(native$fn(model@opt$par) - model@opt$value), 1e-7)
  expect_equal(as.numeric(native$gr(model@opt$par)),
    as.numeric(model@TMB$gr(model@opt$par)), tolerance = 1e-7)
  reference <- TMB::sdreport(native, par.fixed = model@opt$par, getJointPrecision = TRUE)
  actual <- TMB::sdreport(model@TMB, par.fixed = model@opt$par, getJointPrecision = TRUE)
  expect_equal(as.matrix(reference$cov.fixed), as.matrix(actual$cov.fixed), tolerance = 1e-6)
  expect_equal(as.matrix(reference$jointPrecision), as.matrix(actual$jointPrecision), tolerance = 1e-6)
  expect_true(all(is.finite(unmarked::SE(model))))
  d <- attr(model, "msr_gdr_tmb_diagnostics")
  score <- model@TMB$env$f(model@TMB$env$last.par, order = 1L)[model@TMB$env$random]
  precision <- model@TMB$env$spHess(model@TMB$env$last.par, random = TRUE)
  expect_equal(d$conditional_latent_maximum_absolute_gradient, max(abs(score)), tolerance = 1e-7)
  expect_equal(d$conditional_latent_Newton_decrement,
    sqrt(max(0, sum(score * as.numeric(solve(precision, score))))), tolerance = 1e-7)
  expect_equal(aic_tab(list(model), AICc = FALSE, verbose = FALSE)$K, 6L)
  cat("GDR_TMB_WORKFLOW inner_policy=explicit native_joint_stationarity=1 pre_covariance_points=typed full_SE=retained priors=unchanged\n")
})

test_that("positive curvature alone cannot admit a nonstationary or mismatched mode", {
  env <- new.env(parent = emptyenv()); env$random <- 2L; env$last.par <- c(0, .3)
  env$spHess <- function(par, random) Matrix::Diagonal(1L, 2)
  env$f <- function(par, order) c(0, .1)
  result <- list(opt = list(par = 0, value = 0),
    TMB = list(env = env, fn = function(par) 0, gr = function(par) 0),
    sdr = list(pdHess = TRUE, cov.fixed = matrix(1), par.fixed = 0, par.random = .3,
      jointPrecision = Matrix::Diagonal(2L)))
  d <- .msr_gdr_TMB_diagnostics(result)
  expect_true(d$conditional_latent_precision_positive_definite)
  expect_false(d$conditional_latent_stationary)
  expect_identical(d$status, "unusable")
  env$f <- function(par, order) c(0, 0)
  result$sdr$par.random <- .2
  d <- .msr_gdr_TMB_diagnostics(result)
  expect_true(d$conditional_latent_stationary)
  expect_false(d$covariance_at_current_latent_mode)
  expect_identical(d$status, "unusable")
  result$sdr$par.random <- .3; result$TMB$fn <- function(par) .01
  d <- .msr_gdr_TMB_diagnostics(result)
  expect_gt(d$returned_value_gap, 1e-7)
  expect_identical(d$status, "unusable")
  result$TMB$fn <- function(par) 0
  env$f <- function(par, order) stop("declared joint gradient failure")
  d <- .msr_gdr_TMB_diagnostics(result)
  expect_match(d$conditional_stationarity_error, "declared joint gradient failure")
  expect_identical(d$status, "unusable")
  env$f <- function(par, order) c(0, 0)
  env$last.par <- c(.1, .3)
  d <- .msr_gdr_TMB_diagnostics(result)
  expect_false(d$current_fixed_coordinates_at_optimizer)
  expect_identical(d$status, "unusable")
})

test_that("zero-random-effect parents retain zero latent diagnostics and native counts", {
  for (kind in c("P", "NB")) {
    f <- gdr_tmb_refit_fixture(kind)
    d <- attr(f$model, "msr_gdr_tmb_diagnostics")
    expect_identical(d$latent_coordinates, 0L)
    expect_true(d$conditional_latent_stationary)
    expect_equal(d$conditional_latent_maximum_absolute_gradient, 0)
    expect_equal(d$conditional_latent_Newton_decrement, 0)
    expect_true(d$covariance_at_current_latent_mode)
    expect_equal(.msr_parameter_count(f$model), if (kind == "P") 4L else 5L)
    data <- .msr_gdr_TMB_preparation(f$model)
    parameters <- f$model@TMB$env$parList(par = f$model@TMB$env$last.par)
    zero <- .msr_gdr_TMB_adfun(data$model, data[names(data) != "model"], parameters, NULL,
      .msr_gdr_TMB_policy(list(random.start = "zero")))
    expect_lt(abs(zero$fn(f$model@opt$par) - f$model@opt$value), 1e-7)
    expect_equal(as.numeric(zero$gr(f$model@opt$par)),
      as.numeric(f$model@TMB$gr(f$model@opt$par)), tolerance = 1e-7)
  }
})

test_that("field fitting and refitting retain typed points and their fixed-range target", {
  f <- gdr_native_field_fixture(); model <- f$model
  field <- .msr_gdr_field_guard(model)
  files <- list.files(dirname(field$audit_file), pattern = "-attempt-[12]-point[.]rds$", full.names = TRUE)
  expect_length(files, 2L)
  for (file in files) {
    point <- readRDS(file)
    expect_s3_class(point, "msr_gdr_optimizer_point")
    expect_false(point$inference_eligible)
    expect_identical(point$binding$field_id, field$field_id)
    expect_identical(point$binding$source, field$source)
    expect_identical(point$policy$random_start, "cached_best")
    expect_false(point$covariance_saved)
  }
  expect_equal(.msr_parameter_count(model), 6L)
  expect_identical(model@response_field$definition$range, 50000)
  expect_false(model@response_field$definition$range_estimated)
  changed <- f$template$settings; changed$random.start <- "zero"
  expect_error(make_gdr_refitter(model, f$template$starts,
    gdr_test_audit_path("field-policy-change-"), tmb_settings = changed), "retain")
  settings <- f$template$settings
  settings$random.start <- "zero"
  settings$inner.control <- list(maxit = 2500L, tol = 1e-10, grad.tol = 1e-10, step.tol = 1e-12)
  zero <- fit_gdr_response_field(f$template$model, f$coords, 32617, 50000,
    f$template$starts, gdr_test_audit_path("native-field-zero-policy-"), settings)
  expect_lt(abs(as.numeric(unmarked::logLik(zero)) - as.numeric(unmarked::logLik(model))), 1e-6)
  expect_equal(unmarked::coef(zero), unmarked::coef(model), tolerance = 1e-3)
  expect_equal(.msr_parameter_count(zero), 6L)
  expect_false(zero@response_field$definition$range_estimated)
  expect_identical(zero@response_field$settings$policy$requested,
    settings[c("random.start", "inner.control")])
  callback <- make_gdr_refitter(zero, f$template$starts,
    gdr_test_audit_path("field-zero-policy-refitter-"), tmb_settings = settings)
  expect_true(is.function(callback))
})

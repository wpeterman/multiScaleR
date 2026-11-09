test_that("native TMB refitters retain actual fits, start order, counts and joint covariance", {
  for (kind in c("P", "NB", "mixed")) {
    f <- gdr_tmb_refit_fixture(kind)
    expect_identical(as.integer(f$model@opt$convergence), 0L, info = kind)
    if (kind != "mixed") {
      native_C <- gdr_refit_fixture(kind)$model
      expect_lt(abs(as.numeric(unmarked::logLik(f$model)) - as.numeric(unmarked::logLik(native_C))), 1e-4)
      expect_equal(unmarked::coef(f$model), unmarked::coef(native_C), tolerance = 1e-3)
    }
    before <- gdr_guard_state(f$model)
    refit <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path(paste0("tmb-refit-", kind, "-")),
      tmb_settings = f$settings)
    point <- refit(f$model, f$data, list(return_model = FALSE, covs = "hab"))
    expect_s4_class(point, "unmarkedFitGDR")
    expect_false(is.null(point@TMB))
    audit <- attr(point, "msr_gdr_refit_audit")
    expect_length(audit$attempts, 2L)
    expect_true(audit$objective_only)
    expect_true(all(vapply(audit$attempts, function(a) a$covariance$sdreport_computed, logical(1))))
    expect_false(any(vapply(audit$attempts, function(a) a$covariance$coefficient_hessians_omitted, logical(1))))
    expect_true(all(vapply(audit$attempts, function(a) a$covariance$pdHess, logical(1))))
    expect_true(all(vapply(audit$attempts, function(a) all(is.finite(a$gradient)), logical(1))))
    contract <- readRDS(audit$contract_file)
    expect_identical(contract$controls$engine, "TMB")
    expect_identical(contract$controls$dots$control, f$settings$control)
    expect_match(contract$provider$TMB$native_TMB_DLL_sha256, "^[[:xdigit:]]{64}$")
    expect_identical(contract$original_coordinates, gdr_start_values(f$model))
    final <- refit(f$model, f$data, list(return_model = TRUE, covs = "hab"))
    expect_identical(final@call[[1L]], quote(unmarked::gdistremoval))
    expect_identical(final@call$control, f$settings$control)
    last <- tail(attr(final, "msr_gdr_refit_audit")$attempts, 1L)[[1L]]
    expect_identical(last$name, "final")
    expect_true(last$SE_requested)
    expect_identical(last$covariance$status, "available")
    expect_true(last$covariance$fixed_covariance_positive_definite)
    expect_true(last$covariance$joint_precision_positive_definite)
    expect_true(last$covariance$conditional_latent_precision_positive_definite)
    expect_true(last$covariance$covariance_at_optimizer)
    expect_equal(.msr_parameter_count(final), if (kind == "mixed") 6L else 4L + (kind == "NB"))
    expect_identical(last$constructor_coordinates, .msr_gdr_refit_start(final, contract$controls))
    expect_identical(final@data, f$data)
    expect_identical(final@sitesRemoved, f$model@sitesRemoved)
    expect_gte(as.numeric(unmarked::logLik(final)), as.numeric(unmarked::logLik(point)) - 1e-5)
    expect_identical(gdr_guard_state(f$model), before)
    for (a in attr(final, "msr_gdr_refit_audit")$attempts) {
      expect_true(file.exists(a$request_file))
      expect_identical(readRDS(a$audit_file), a)
    }
    table <- aic_tab(list(final), AICc = FALSE)
    expect_equal(table$K, .msr_parameter_count(final))
    expect_equal(table$AIC, -2 * as.numeric(unmarked::logLik(final)) + 2 * table$K)
    if (kind == "mixed") {
      expect_equal(last$covariance$latent_coordinates, 12L)
      expect_equal(contract$original_coordinates$random_modes_omitted, c(lambda = 8L, dist = 4L, rem = 0L))
    }
  }
  cat("GDR_TMB_REFIT native_P_NB_mixed=3 fixed_SD_counts=preserved joint_covariance=native\n")
})

test_that("TMB conditional curvature uses the current mode rather than its initial default", {
  # An API witness, not a fit: spHess's default point has negative curvature,
  # while the explicitly supplied current conditional mode is positive.
  env <- new.env(parent = emptyenv())
  env$random <- 2L
  env$last.par <- c(0, .3)
  env$spHess <- function(par = c(0, 0), random = TRUE) {
    Matrix::Diagonal(1L, if (identical(par, c(0, .3))) 2 else -2)
  }
  result <- list(opt = list(par = 0), TMB = list(env = env, gr = function(par) 0),
    sdr = list(pdHess = TRUE, cov.fixed = matrix(1), par.fixed = 0,
      jointPrecision = Matrix::Diagonal(2L)))
  diagnostic <- .msr_gdr_TMB_diagnostics(result)
  expect_identical(diagnostic$status, "available")
  expect_true(diagnostic$conditional_latent_precision_positive_definite)
  expect_identical(diagnostic$conditional_curvature_point, "native_current_mode")
})

test_that("explicit TMB settings, target bindings and existing native guards fail closed", {
  f <- gdr_tmb_refit_fixture("P")
  expect_error(make_gdr_refitter(f$model, f$starts, tempfile()), "explicitly contain")
  for (setting in list(list(), within(f$settings, threads <- 2L),
      within(f$settings, control <- list(maxit = 0L)), within(f$settings, method <- "SANN")))
    expect_error(make_gdr_refitter(f$model, f$starts, tempfile(), tmb_settings = setting))
  expect_error(make_gdr_refitter(gdr_refit_fixture("P")$model, f$starts, tempfile(),
    tmb_settings = f$settings), "only to native TMB")
  wrong <- f$model; wrong@K <- wrong@K + 1
  expect_error(make_gdr_refitter(wrong, f$starts, tempfile(), tmb_settings = f$settings), "target differs")
  wrong <- f$model; wrong@data@yDistance[1, 1] <- wrong@data@yDistance[1, 1] + 1
  expect_error(make_gdr_refitter(wrong, f$starts, tempfile(), tmb_settings = f$settings))
  wrong <- f$model; wrong@mixture <- "ZIP"
  expect_error(make_gdr_refitter(wrong, f$starts, tempfile(), tmb_settings = f$settings), "TMB ZIP")
  wrong <- f$model; wrong@formlist$lambda <- ~hab + offset(hab)
  expect_error(make_gdr_refitter(wrong, f$starts, tempfile(), tmb_settings = f$settings), "offsets")
  refit <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path("tmb-data-"), tmb_settings = f$settings)
  wrong <- f$data; wrong@yRemoval[1, 1] <- wrong@yRemoval[1, 1] + 1
  error <- tryCatch(refit(f$model, wrong, list(covs = "hab")), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_length(error$audit$attempts, 0L)
  expect_match(conditionMessage(error), "changed histories")
  expect_true(file.exists(error$audit_file))
})

test_that("native TMB uses explicit positive iteration budgets and retains failures", {
  f <- gdr_tmb_refit_fixture("P")
  settings <- f$settings; settings$control$maxit <- 1L
  refit <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path("tmb-iterations-"), tmb_settings = settings)
  error <- tryCatch(refit(f$model, f$data, list(covs = "hab")), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_length(error$audit$attempts, 2L)
  expect_true(any(vapply(error$audit$attempts, function(a) !a$converged, logical(1))))
  expect_true(all(file.exists(vapply(error$audit$attempts, `[[`, character(1), "audit_file"))))
})

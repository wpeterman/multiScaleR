test_that("native TMB S4 results work with all four kernel refit paths", {
  f <- gdr_tmb_refit_fixture("P")
  frame <- f$model@data; n <- nrow(frame@siteCovs)
  distances <- rep(list(c(.1, .4, .9)), n)
  cells <- lapply(frame@siteCovs$hab, function(x) matrix(x + c(-.2, .1, .4), ncol = 1,
    dimnames = list(NULL, "hab")))
  names(cells) <- names(distances) <- rownames(frame@siteCovs)
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    callback <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path(paste0("tmb-", kernel, "-")),
      tmb_settings = f$settings)
    context <- build_opt_context(f$model, cells, refit_fn = callback)
    context$covariate_transform <- "raw"
    par <- if (kernel == "expow") c(.3, 1.5) else .3
    expected <- .msr_kernel_cov_w(par, distances, cells, kernel, context)[, 1L]
    point <- kernel_scale_fn(par, distances, cells, kernel, f$model, opt_context = context)
    final <- kernel_scale_fn(par, distances, cells, kernel, f$model, opt_context = context, mod_return = TRUE)
    expect_true(is.finite(point))
    expect_equal(point, -as.numeric(unmarked::logLik(final$mod)), tolerance = 1e-5)
    expect_s4_class(final$mod, "unmarkedFitGDR")
    expect_equal(unname(final$mod@data@siteCovs$hab), unname(expected), tolerance = 1e-12)
    expect_identical(final$scl_params, list(mean = c(hab = 0), sd = c(hab = 1)))
    expect_identical(final$mod@data@yDistance, frame@yDistance)
    expect_identical(final$mod@data@yRemoval, frame@yRemoval)
    expect_equal(.msr_parameter_count(final$mod), 4L)
    expect_identical(.msr_gdr_variance_status(final$mod, TRUE)$status, "available")
    prediction <- unmarked::predict(final$mod, type = "lambda",
      newdata = final$mod@data@siteCovs[c(1, 20, 50), , drop = FALSE])
    expect_true(all(is.finite(prediction$Predicted)))
    expect_true(all(is.finite(prediction$SE) & prediction$SE > 0))
  }
  # One real native callback in a fresh worker verifies provider/DLL restoration.
  callback <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path("tmb-worker-"), tmb_settings = f$settings)
  worker <- gdr_refit_worker_probe(callback, f$model, f$data)
  expect_equal(worker$K, 4L)
  audit <- readRDS(worker$audit_file)
  expect_true(audit$provider_check$files_rechecked)
  expect_match(audit$provider_check$TMB$native_TMB_DLL_sha256, "^[[:xdigit:]]{64}$")
  cat("GDR_TMB_KERNELS kernels=4 native_S4=retained worker_provider=verified\n")
})

test_that("unusable highest TMB attempts retain optimizer evidence and cannot lose to a lower fit", {
  f <- gdr_tmb_refit_fixture("P")
  counter <- 0L
  # A declared constructor stub tests the audit selection contract only.
  # Its returned original model is not a fitted update or scientific reference.
  constructor <- function(args, settings) {
    counter <<- counter + 1L
    if (counter == 1L) stop(structure(list(message = "declared post-optimizer covariance failure",
      call = NULL, optimizer = list(par = f$model@opt$par, value = -1e9,
        convergence = 0L), covariance = NULL), class = c("msr_gdr_tmb_error", "error", "condition")))
    f$model
  }
  testthat::with_mocked_bindings({
    callback <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path("tmb-high-failure-"), tmb_settings = f$settings)
    error <- tryCatch(callback(f$model, f$model@data), error = identity)
    expect_s3_class(error, "msr_gdr_refit_error")
    expect_match(conditionMessage(error), "highest declared native likelihood")
    expect_length(error$audit$attempts, 2L)
    expect_identical(error$audit$selected, 1L)
    expect_equal(error$audit$attempts[[1L]]$optimizer$value, -1e9)
    expect_equal(error$audit$attempts[[1L]]$logLik, 1e9)
    expect_false(error$audit$attempts[[1L]]$eligible)
    expect_match(error$audit$attempts[[1L]]$error, "post-optimizer covariance failure")
  }, .msr_gdr_TMB_refit = constructor, .package = "multiScaleR")
})

test_that("restored TMB file mismatches remain failures on every request", {
  f <- gdr_tmb_refit_fixture("P")
  callback <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path("tmb-files-"), tmb_settings = f$settings)
  environment(callback)$counter$provider_pid <- 0L
  environment(callback)$source$TMB$native_TMB_DLL_sha256 <- "declared-changed-provider"
  for (i in 1:2) {
    error <- tryCatch(callback(f$model, f$data), error = identity)
    expect_s3_class(error, "msr_gdr_refit_error")
    expect_match(conditionMessage(error), "native TMB files changed")
    expect_true(error$audit$provider_check$files_rechecked)
    expect_length(error$audit$attempts, 0L)
    expect_identical(environment(callback)$counter$provider_pid, 0L)
  }
})

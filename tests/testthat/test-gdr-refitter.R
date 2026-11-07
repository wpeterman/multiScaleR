test_that("C refitters preserve evaluated settings, starts, penalties and final uncertainty", {
  for (mixture in c("P", "NB", "ZIP")) {
    fixture <- gdr_refit_fixture(mixture)
    model <- fixture$model
    expect_identical(as.integer(model@opt$convergence), 0L)
    before <- gdr_guard_state(model)
    path <- gdr_test_audit_path(paste0("gdr-refit-", mixture, "-"))
    refit <- make_gdr_refitter(model, fixture$starts, path,
      provenance = list(fixture = mixture, seed = fixture$seed))
    contract <- readRDS(file.path(path, "contract.rds"))
    expect_identical(contract$controls$engine, "C")
    expect_equal(contract$controls$K, 400)
    expect_true(contract$controls$se)
    expect_identical(contract$controls$dots$control, list(maxit = 2500, reltol = 1e-10))
    expect_match(contract$provider$native_DLL_sha256, "^[[:xdigit:]]{64}$")
    objective <- refit(model, fixture$data, list(return_model = FALSE, covs = "hab"))
    audit <- attr(objective, "msr_gdr_refit_audit")
    expect_identical(audit$status, "returned")
    expect_length(audit$attempts, 3L)
    expect_true(all(vapply(audit$attempts, function(x) !x$SE_requested, logical(1))))
    expect_null(objective@opt$hessian)
    expect_equal(as.numeric(unmarked::logLik(objective)),
      max(vapply(audit$attempts, `[[`, numeric(1), "logLik")))
    expect_true(all(file.exists(vapply(audit$attempts, `[[`, character(1), "audit_file"))))
    expect_identical(readRDS(audit$path)$status, "returned")
    expect_match(audit$input_sha256, "^[[:xdigit:]]{64}$")
    expect_equal(.msr_parameter_count(objective), 4L + as.integer(mixture != "P"))
    final <- refit(model, fixture$data, list(return_model = TRUE, covs = "hab"))
    final_audit <- attr(final, "msr_gdr_refit_audit")
    expect_length(final_audit$attempts, 4L)
    expect_true(final_audit$attempts[[4L]]$SE_requested)
    expect_identical(final_audit$attempts[[4L]]$covariance$status, "available")
    expect_true(is.matrix(final@opt$hessian))
    expect_gte(as.numeric(unmarked::logLik(final)), as.numeric(unmarked::logLik(objective)) - 1e-5)
    expect_identical(names(unmarked::coef(final)), names(unmarked::coef(model)))
    expect_equal(.msr_parameter_count(final), .msr_parameter_count(model))
    expect_identical(final@data, fixture$data)
    expect_identical(gdr_guard_state(model), before)
    # The callback and its compact contract retain working PSOCK/RDS semantics.
    expect_type(unserialize(serialize(refit, NULL)), "closure")
    if (mixture == "P") {
      worker <- gdr_refit_worker_probe(refit, model, fixture$data)
      expect_equal(worker$logLik, as.numeric(unmarked::logLik(objective)), tolerance = 1e-6)
      expect_equal(worker$K, 4L)
      expect_true(file.exists(worker$audit_file))
      worker_audit <- readRDS(worker$audit_file)
      expect_true(worker_audit$provider_check$files_rechecked)
      expect_identical(worker_audit$provider_check$native_DLL_sha256, contract$provider$native_DLL_sha256)
    }
  }
})


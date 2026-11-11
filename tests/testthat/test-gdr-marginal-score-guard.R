test_that("generic native GDR eligibility rejects excessive finite marginal scores", {
  f <- gdr_tmb_refit_fixture("P")
  original <- attr(f$model, "msr_gdr_tmb_diagnostics")
  # Selection-only stubs preserve a real fitted template and its covariance.
  # They do not refit data or supply another scientific likelihood reference.
  for (score in c(1e-3, .01)) {
    constructor <- function(args, settings, checkpoint = NULL) {
      model <- f$model
      diagnostic <- original
      diagnostic$gradient <- rep(score, length(model@opt$par))
      diagnostic$maximum_absolute_gradient <- score
      attr(model, "msr_gdr_tmb_diagnostics") <- diagnostic
      model
    }
    with_mocked_bindings({
      callback <- make_gdr_refitter(f$model, f$starts,
        gdr_test_audit_path("generic-score-rejected-"), tmb_settings = f$settings)
      error <- tryCatch(callback(f$model, f$model@data), error = identity)
      expect_s3_class(error, "msr_gdr_refit_error")
      expect_match(conditionMessage(error), "highest declared native likelihood")
      expect_length(error$audit$attempts, 2L)
      for (attempt in error$audit$attempts) {
        expect_true(attempt$converged)
        expect_identical(attempt$covariance$status, "available")
        expect_false(attempt$eligible)
        expect_false(attempt$marginal_score$eligible)
        expect_equal(attempt$marginal_score$maximum_absolute_gradient, score)
        expect_identical(attempt$marginal_score$threshold, 1e-3)
      }
    }, .msr_gdr_TMB_refit = constructor, .package = "multiScaleR")
  }
  expect_identical(attr(f$model, "msr_gdr_tmb_diagnostics"), original)
})

test_that("qualified zero-score TMB stubs preserve counts and final uncertainty", {
  f <- gdr_tmb_refit_fixture("P")
  original <- attr(f$model, "msr_gdr_tmb_diagnostics")
  constructor <- function(args, settings, checkpoint = NULL) {
    model <- f$model
    diagnostic <- original
    diagnostic$gradient <- rep(0, length(model@opt$par))
    diagnostic$maximum_absolute_gradient <- 0
    attr(model, "msr_gdr_tmb_diagnostics") <- diagnostic
    model
  }
  with_mocked_bindings({
    callback <- make_gdr_refitter(f$model, f$starts,
      gdr_test_audit_path("generic-score-qualified-"), tmb_settings = f$settings)
    model <- callback(f$model, f$model@data)
    audit <- attr(model, "msr_gdr_refit_audit")
    expect_s4_class(model, "unmarkedFitGDR")
    expect_length(audit$attempts, 3L)
    expect_true(all(vapply(audit$attempts, function(z) z$eligible && z$marginal_score$eligible, logical(1))))
    expect_equal(.msr_parameter_count(model), 4L)
    expect_identical(unmarked::SE(model), unmarked::SE(f$model))
    expect_identical(model@data, f$model@data)
  }, .msr_gdr_TMB_refit = constructor, .package = "multiScaleR")
  expect_identical(attr(f$model, "msr_gdr_tmb_diagnostics"), original)
})

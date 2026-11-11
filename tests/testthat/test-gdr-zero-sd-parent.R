test_that("exact zero-SD parents preserve the conditional finite-K target", {
  gdr_gaussian_require_RTMB()
  f <- gdr_gaussian_fixture(TRUE, TRUE, 2L)
  baseline <- .msr_gdr_hash(list(f$data, f$pars))
  lambda_zero <- .msr_gdr_gaussian_zero_sd_parent(f$data, f$pars, "lambda")
  expected_data <- f$data
  expected_data$block <- integer()
  fit <- RTMB::MakeADFun(lambda_zero$program$evaluate, lambda_zero$parameters, random = NULL, silent = TRUE)
  expect_lt(abs(fit$fn(fit$par) - gdr_gaussian_numeric(expected_data, f$pars)), 1e-9)
  expect_identical(lambda_zero$parameters$b_dist, f$pars$b_dist)
  expect_identical(lambda_zero$parameters$lsigma_dist, f$pars$lsigma_dist)

  zero <- .msr_gdr_gaussian_zero_sd_parent(f$data, f$pars, c("lambda", "dist"))
  native <- gdr_gaussian_fixture(FALSE, FALSE, 2L)$captured
  ref <- TMB::MakeADFun(c(list(model = native$model), native$data), zero$parameters,
    random = NULL, DLL = "unmarked_TMBExports", silent = TRUE)
  parent <- RTMB::MakeADFun(zero$program$evaluate, zero$parameters, random = NULL, silent = TRUE)
  point <- parent$fn(parent$par)
  expect_lt(abs(point - ref$fn(ref$par)), 1e-9)
  expect_identical(zero$program$latent_coordinate_count, 0L)
  expect_identical(zero$program$fixed_coordinate_count, 13L)
  expect_false(zero$fit_admission_established)
  expect_false(zero$uncertainty_evaluated)

  # The positive-SD Laplace path should approach the exact parent. The
  # boundary itself is evaluated through the native no-field reference above.
  program <- .msr_gdr_gaussian_program(f$data, f$pars)
  near <- f$pars
  near$b_lambda[] <- near$b_dist[] <- 0
  error <- vapply(c(1e-3, 1e-4, 1e-5), function(sd) {
    near$lsigma_lambda <- near$lsigma_dist <- log(sd)
    almost <- RTMB::MakeADFun(program$evaluate, near, random = f$captured$random, silent = TRUE)
    abs(almost$fn(almost$par) - point)
  }, numeric(1))
  expect_true(all(is.finite(error)))
  expect_true(all(error[-1L] < 0.05 * head(error, -1L)))
  expect_lt(error[3L], 1e-6)
  expect_identical(.msr_gdr_hash(list(f$data, f$pars)), baseline)
  expect_error(.msr_gdr_gaussian_zero_sd_parent(f$data, f$pars, c("lambda", "lambda")), "existing")
  expect_error(.msr_gdr_gaussian_zero_sd_parent(f$data, f$pars, "rem"), "existing")
})

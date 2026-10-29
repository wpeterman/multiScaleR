test_that("metadata separates posterior computation from hierarchical likelihood", {
  describe <- multiScaleR:::.msr_scale_inference_metadata
  for (method in c("prepare", "nuts", "joint_laplace", "marginal_laplace", "importance_corrected")) {
    info <- describe("posterior", method)
    expect_identical(info$target, "posterior")
    expect_identical(info$population_priors, "proper")
    expect_identical(info$objective, "negative_log_posterior")
    expect_identical(info$uncertainty$intervals,
      if (method == "prepare") "not_provided" else "marginal_credible")
    expect_identical(info$uncertainty$calibration, "not_established")
    expect_false(info$metadata_establishes_eligibility)
  }
  for (method in c("laplace", "factorized")) {
    info <- describe("marginal_likelihood", method)
    expect_identical(info$population_priors, "none")
    expect_identical(info$objective, "negative_marginal_log_likelihood")
    expect_identical(info$integration$population_parameters, "optimized")
    expect_identical(info$uncertainty$population_parameters, "not_propagated")
    expect_identical(info$uncertainty$intervals, "not_provided")
  }
  expect_error(describe("posterior", "factorized"), "incompatible")
  expect_error(describe("marginal_likelihood", "nuts"), "incompatible")
  expect_error(describe("unknown", "prepare"), "Unknown")
})

test_that("prepared and rejected objects retain legacy fields and diagnostics", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  prepared <- prepare_bayesian_fixture(f)
  expect_identical(prepared$inference, "prepared")
  expect_identical(prepared$inference_metadata$target, "posterior")
  expect_identical(prepared$inference_metadata$solver, "none")
  expect_null(prepared$draws)
  expect_null(prepared$diagnostics)
  expect_identical(tail(names(prepared), 1), "inference_metadata")
  # All-zero responses have a known unbounded count-intercept likelihood.
  f$data$count[] <- 0
  spec <- multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    list(habitat = ~1 + (1 | year)), "poisson", "key", 200)
  for (fit in list(multiScaleR:::.msr_marginal_fit(spec), multiScaleR:::.msr_factorized_fit(spec))) {
    expect_false(fit$diagnostics$passed)
    expect_false(fit$uncertainty$available)
    expect_length(fit$attempts, 0)
    expect_identical(fit$inference_metadata$target, "marginal_likelihood")
    expect_identical(fit$inference_metadata$diagnostics_field, "diagnostics")
    expect_false(fit$inference_metadata$metadata_establishes_eligibility)
    expect_identical(tail(names(fit), 1), "inference_metadata")
  }
})

test_that("older stored fits remain usable without metadata", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f))
  current <- posterior_scales(fit)
  current_prediction <- predict(fit, summary = FALSE)
  fit$inference_metadata <- NULL
  expect_identical(posterior_scales(fit), current)
  expect_identical(predict(fit, summary = FALSE), current_prediction)
  expect_equal(dim(current_prediction), c(8, nrow(f$data)))
  expect_identical(summary(fit)$inference, "test fixture")
})

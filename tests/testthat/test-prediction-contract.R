test_that("posterior contracts separate means, replicates and zeroed effects", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"))
  fit$diagnostics <- list(passed = FALSE, reasons = "retained diagnostic failure")
  before <- fit
  set.seed(31); state <- .Random.seed
  mean <- predict(fit)
  contract <- attr(mean, "prediction_contract")
  expect_identical(contract$inference_target, "posterior")
  expect_identical(contract$quantity, "expected_response")
  expect_false(contract$uncertainty$observation_noise)
  expect_identical(contract$uncertainty$intervals, "conditional_posterior_quantiles")
  expect_identical(contract$diagnostics_status, "failed")
  expect_false(contract$metadata_establishes_eligibility)
  expect_true(contract$joint_parameter_map_complete)
  expect_identical(.Random.seed, state)
  expect_null(attr(predict(fit, summary = FALSE), "prediction_contract"))
  expect_identical(attr(predict(fit, type = "link"), "prediction_contract")$quantity,
                   "linear_predictor")
  simulated <- simulate(fit, nsim = 2, seed = 41)
  predictive <- attr(simulated, "prediction_contract")
  expect_identical(predictive$quantity, "replicate_response")
  expect_true(predictive$uncertainty$observation_noise)
  expect_identical(predictive$uncertainty$intervals, "posterior_predictive_replicates")
  expect_identical(.Random.seed, state)
  reference <- attr(predict(fit, re.form = NA, scale_re.form = NA), "prediction_contract")
  expect_identical(reference$effects$scale$known, "zero_deviation_reference")
  expect_identical(reference$effects$response$known, "zero_deviation_reference")
  expect_identical(reference$uncertainty$intervals, "posterior_reference_quantiles")
  expect_identical(fit, before)
})

test_that("new-group contracts preserve joint effects, row order and unused responses", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"))
  nd <- f$data[c(2, 1, 2), ]
  nd$year <- c("future", "y1", "future")
  nd$site <- c("new_site", "s1", "new_site")
  p <- predict(fit, nd, kernel_inputs = f$inputs, new_scale_levels = "sample",
    new_response_levels = "zero", seed = 91)
  contract <- attr(p, "prediction_contract")
  expect_identical(contract$effects$scale$unseen, "sample_once_per_group_and_joint_vector")
  expect_identical(contract$effects$response$unseen, "zero_deviation_reference")
  expect_true(contract$effects$scale$unseen_present)
  expect_true(contract$effects$response$unseen_present)
  expect_identical(contract$uncertainty$intervals, "new_group_predictive_quantiles")
  expect_identical(p$point_id, as.character(nd$key))
  changed <- nd; changed$count <- c(NA, -1, 1e6)
  expect_identical(p, predict(fit, changed, kernel_inputs = f$inputs,
    new_scale_levels = "sample", new_response_levels = "zero", seed = 91))
  draws <- predict(fit, nd, kernel_inputs = f$inputs, summary = FALSE,
    new_scale_levels = "sample", new_response_levels = "zero", seed = 91)
  expect_identical(draws[, 1], draws[, 3])
  expect_null(attr(draws, "prediction_contract"))
  scale <- predict_scales(fit, nd, new_levels = "sample", seed = 91)
  expect_identical(attr(scale, "prediction_contract")$quantity, "kernel_scale")
  expect_identical(attr(scale, "prediction_contract")$effects$response$known, "not_applicable")
  expect_identical(attr(predict_scales(fit, nd, re.form = NA),
    "prediction_contract")$effects$scale$known, "zero_deviation_reference")
  permutation <- c(3, 1, 2)
  reordered <- predict(fit, nd[permutation, ], kernel_inputs = f$inputs,
    new_scale_levels = "sample", new_response_levels = "zero", seed = 91)
  expect_identical(attr(reordered, "prediction_contract"), contract)
  expect_identical(reordered$median[order(permutation)], p$median)
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(predict(restored, nd, kernel_inputs = f$inputs,
    new_scale_levels = "sample", new_response_levels = "zero", seed = 91), p)
  expect_error(predict(fit, nd, kernel_inputs = f$inputs), "New scale group levels")
})

test_that("older fits describe missing provenance without manufacturing eligibility", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f, "gaussian"))
  expected <- predict(fit, summary = FALSE)
  fit$inference_metadata <- NULL
  restored <- unserialize(serialize(fit, NULL))
  expect_identical(predict(restored, summary = FALSE), expected)
  info <- attr(predict(restored), "prediction_contract")
  expect_identical(info$integration, list(recorded = FALSE, prediction_vectors = "stored_joint_draws"))
  expect_identical(info$diagnostics_status, "unavailable")
  expect_true(info$joint_parameter_map_complete)
  context <- multiScaleR:::.msr_bayes_prediction_prepare(restored, NULL, NULL,
    NULL, "error", "error", TRUE, TRUE)
  restored$draws$response_z <- NULL
  incomplete <- multiScaleR:::.msr_bayes_prediction_contract(restored, context,
    "expected_response", "response")
  expect_false(incomplete$joint_parameter_map_complete)
  expect_false(incomplete$metadata_establishes_eligibility)
})

test_that("analytic contracts preserve point values and covariance restrictions", {
  f <- profile_glm_fixture(n = 120)
  fit <- profile_glm_fit(f, profile_control = list(n_starts = 3))
  expected <- predict(fit)
  snapshot <- prepare_prediction(fit)
  expect_identical(predict(snapshot), expected)
  info <- snapshot$prediction_contract$response
  expect_identical(info$inference_target, "likelihood")
  expect_identical(info$quantity, "expected_response")
  expect_identical(info$uncertainty$intervals, "mean_confidence_when_eligible")
  expect_false(info$uncertainty$observation_noise)
  expect_true(info$fixed_inputs)
  old <- snapshot; old$prediction_contract <- NULL
  expect_identical(predict(unserialize(serialize(old, NULL))), expected)
  bad <- fit; bad$uncertainty$available <- FALSE
  bad$uncertainty$reasons <- "retained weak identification"
  weak <- prepare_prediction(bad)
  expect_false(weak$prediction_contract$response$uncertainty$joint_available)
  expect_identical(predict(weak), expected)
  expect_error(predict(weak, se.fit = TRUE), "Joint Wald uncertainty unavailable")
  expect_identical(predict(weak, se.fit = TRUE, uncertainty = "conditional"),
    predict(bad, se.fit = TRUE, uncertainty = "conditional"))
  expect_error(prepare_prediction(fit, max_bytes = as.numeric(object.size(snapshot)) - 1),
    "exceed max_bytes")
  conditional <- fit; conditional$backend <- "conditional_logit"
  choice <- multiScaleR:::.msr_analytic_prediction_contract(conditional, "response")
  expect_identical(choice$quantity, "choice_probability")
  expect_identical(choice$inference_target, "conditional_likelihood")
  expect_identical(choice$uncertainty$intervals, "not_provided")
  expect_false(choice$uncertainty$joint_available)
})

test_that("native contracts retain components and never equate zeroing with integration", {
  describe <- multiScaleR:::.msr_native_prediction_contract
  for (zero in list(NA, local(~0))) {
    info <- describe(NULL, "glmmTMB", "response", list(re.form = zero))
    expect_identical(info$effects$response, "zero_deviation_reference")
    expect_identical(info$integration$random_effects, "set_to_zero_not_integrated")
    expect_identical(info$uncertainty$kernel_parameters, "fixed_at_estimate")
  }
  expect_identical(describe(NULL, "glmmTMB", "response", list())$effects$response,
                   "native_conditional_effects")
  expect_identical(describe(NULL, "glmmTMB", "zprob", list())$quantity, "native_component")
  for (type in c("psi", "col", "ext", "det", "state")) {
    info <- describe(NULL, "unmarked", type, list(backTransform = FALSE))
    expect_identical(info$component, type)
    expect_identical(info$quantity, "native_component")
    expect_identical(info$output_scale, "native_link")
    expect_identical(info$uncertainty$intervals, "native_conditional_confidence")
    expect_false(info$uncertainty$observation_noise)
    expect_false(info$metadata_establishes_eligibility)
  }
})

test_that("native snapshots retain exact native output and survive metadata-free reload", {
  skip_if_not_installed("unmarked")
  f <- make_unmarked_fixture()
  object <- f$obj
  snapshot <- prepare_prediction(object, component = "state")
  expected <- unmarked::predict(object$opt_mod, type = "state")
  expect_identical(snapshot$prediction_contract$state$uncertainty$standard_errors, "provided")
  point <- predict(object, type = "state")
  expect_equal(as.numeric(point), expected$Predicted)
  expect_identical(attr(point, "prediction_contract")$uncertainty$standard_errors, "not_requested")
  table <- predict(object, type = "state", se.fit = TRUE, uncertainty = "conditional")
  expect_identical(attr(table, "prediction_contract")$uncertainty$standard_errors, "provided")
  expect_identical(predict(snapshot), expected)
  expect_identical(snapshot$prediction_contract$state$component, "state")
  expect_null(snapshot$opt_mod)
  expect_null(snapshot$kernel_inputs)
  expect_lt(as.numeric(object.size(snapshot)), 64 * 1024^2)
  snapshot$prediction_contract <- NULL
  expect_identical(predict(unserialize(serialize(snapshot, NULL))), expected)
  expect_error(predict(snapshot, uncertainty = "joint"), "native prediction tables")
})

test_that("factorized quantiles remain conditional on the supplied population vector", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  spec <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp,
    f$inputs, f$data, list(habitat = ~1), "gaussian", "key", 200)
  p <- spec$parameters; p$scale_z <- numeric()
  result <- multiScaleR:::.msr_factorized_reconstruct(spec, p)
  info <- result$prediction_contract
  expect_identical(info$inference_target, "likelihood")
  expect_identical(info$uncertainty$population_parameters, "conditioned_on_supplied_vector")
  expect_identical(info$uncertainty$intervals, "deterministic_at_supplied_population")
  expect_contains(info$uncertainty$omitted, "population_parameter_uncertainty")
  expect_false(info$uncertainty$observation_noise)
  expect_true(result$conditional_only)
  expect_identical(info$diagnostics_status, "unavailable")
})

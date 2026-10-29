test_that("neutral preparation has no population priors or Bayesian dependency", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  formula <- count ~ habitat + temp + offset(log(effort))
  scales <- list(habitat = ~1 + (1 | year))
  neutral <- multiScaleR:::.msr_scale_prepare(formula, f$inputs, f$data,
    scales, "nbinom2", "key", 200)
  expect_null(neutral$priors)
  expect_null(neutral$scales$habitat$prior_mean)
  expect_null(neutral$scales$habitat$prior_sd)
  expect_equal(neutral$parameters$log_scale_sd, log(.2))
  expect_identical(neutral$observation_ids, f$data$key)
  expect_identical(environment(neutral$scales$habitat$formula), baseenv())
  local_mocked_bindings(
    .msr_bayes_prepare = function(...) stop("Bayesian preparation was called"),
    .msr_bayes_priors = function(...) stop("Bayesian priors were called"),
    .package = "multiScaleR")
  likelihood <- multiScaleR:::.msr_marginal_prepare(formula, f$inputs, f$data,
    scales, "nbinom2", "key", 200)
  expect_identical(likelihood$V, neutral$V)
  expect_identical(likelihood$center, neutral$center)
  expect_identical(likelihood$spread, neutral$spread)
  expect_null(likelihood$priors)
})

test_that("Bayesian attachment preserves neutral designs and explicit prior starts", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3, annual = FALSE)
  formula <- count ~ habitat + temp + (1 | site)
  scales <- list(habitat = ~0 + year)
  neutral <- multiScaleR:::.msr_scale_prepare(formula, f$inputs, f$data,
    scales, "gaussian", "key", 300)
  priors <- multiScaleR:::.msr_bayes_priors(list(scale_median = 300,
    scale_coef_sd = .7, response_group_sd = 1.6, dispersion_log_mean = -2), "habitat", f$inputs)
  attached <- multiScaleR:::.msr_bayes_attach_priors(neutral, priors)
  expect_null(neutral$priors)
  expect_null(neutral$scales$habitat$prior_mean)
  for (field in c("rows", "V", "C", "delta", "center", "spread", "X_fixed", "response"))
    expect_identical(attached[[field]], neutral[[field]])
  expect_equal(attached$scales$habitat$prior_mean, rep(log(.3), 3))
  expect_equal(attached$parameters$scale_beta, rep(log(.3), 3))
  expect_equal(attached$parameters$log_response_sd, log(.8))
  expect_identical(attached$parameters$log_dispersion, -2)
  expect_identical(tail(names(attached), 2), c("priors", "parameters"))
  expect_identical(multiScaleR:::.msr_bayes_attach_priors(attached, priors), attached)
  named <- prepare_bayesian_fixture(f, priors = list(scale_group_sd = c(scale = .8),
    response_group_sd = c(response = 1.6)))
  expect_identical(named$parameters$log_scale_sd, c(scale = log(.4)))
  expect_identical(named$parameters$log_response_sd, c(response = log(.8)))
})

test_that("neutral reference scales align by name and retain pooled row identity", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  f$inputs$binned$vsum$other <- f$inputs$binned$vsum$habitat[, 15:1]
  f$inputs$binned$csum$other <- f$inputs$binned$csum$habitat[, 15:1]
  formula <- count ~ habitat + other + temp
  scales <- list(habitat = ~1 + (1 | year), other = ~1 + (1 | site))
  original <- multiScaleR:::.msr_scale_prepare(formula, f$inputs, f$data,
    scales, "poisson", "key", c(other = 350, habitat = 200))
  index <- rev(seq_len(nrow(f$data)))
  shuffled <- multiScaleR:::.msr_scale_prepare(formula, f$inputs, f$data[index, ],
    scales, "poisson", "key", c(habitat = 200, other = 350))
  expect_equal(shuffled$center, original$center, tolerance = 1e-14)
  expect_equal(shuffled$spread, original$spread, tolerance = 1e-14)
  expect_identical(shuffled$V$habitat, original$V$habitat[index, ])
  expect_equal(original$parameters$scale_beta, log(c(.2, .35)))
  expect_error(multiScaleR:::.msr_scale_reference(c(100, 200), c("habitat", "other")), "name every")
  expect_error(multiScaleR:::.msr_scale_reference(NA_real_, "habitat"), "positive distances")
})

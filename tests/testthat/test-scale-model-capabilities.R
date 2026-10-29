test_that("capabilities distinguish structural support from inference targets", {
  table <- multiScaleR:::.msr_scale_capabilities()
  expect_equal(nrow(table), 7)
  expect_false(anyDuplicated(table[c("target", "method")]) > 0L)
  expect_false(any(table$public_interface[table$target == "marginal_likelihood"]))
  expect_true(all(is.infinite(table$max_kernel_predictors[table$target == "posterior"])))
  expect_equal(table$required_scale_sd_count[table$method %in%
    c("marginal_laplace", "importance_corrected")], c(1L, 1L))
  expect_identical(table$response_effects[table$method == "factorized"],
    "gaussian_within_scale_group")
  for (i in seq_len(nrow(table))) {
    info <- multiScaleR:::.msr_scale_inference_metadata(table$target[i], table$method[i])
    expect_identical(info$target, table$target[i])
    expect_identical(info$method, table$method[i])
    expect_false(info$metadata_establishes_eligibility)
    expect_identical(info$uncertainty$calibration, "not_established")
  }
  select <- multiScaleR:::.msr_scale_capability
  expect_error(select("posterior", "factorized"), "incompatible")
  expect_error(select("marginal_likelihood", "nuts"), "incompatible")
  expect_error(select("unknown", "nuts"), "Unknown")
  expect_error(select("posterior", c("nuts", "prepare")), "incompatible")
  expect_error(select("posterior", NA_character_), "incompatible")
})

test_that("supported families, kernels and zero-SD faces pass without changing models", {
  check <- multiScaleR:::.msr_scale_check_capability
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 4, n_year = 3, kernel = kernel)
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    for (scale_formula in list(~1, ~1 + (1 | year))) {
      spec <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp,
        f$inputs, f$data, list(habitat = scale_formula), family, "key", 200)
      expect_identical(check(spec, "marginal_likelihood", "laplace"), spec)
      expect_identical(check(spec, "marginal_likelihood", "factorized"), spec)
      if (length(spec$parameters$log_scale_sd)) {
        reduced <- multiScaleR:::.msr_marginal_reduce(spec, "scale")
        expect_identical(check(reduced, "marginal_likelihood", "factorized"), reduced)
      }
    }
    b <- multiScale_bayes(count ~ habitat + temp + (1 | site), f$inputs, f$data,
      list(habitat = ~1 + (1 | year)), family, "key", method = "prepare")$specification
    for (method in c("prepare", "nuts", "joint_laplace", "marginal_laplace", "importance_corrected"))
      expect_identical(check(b, "posterior", method), b)
    expect_identical(check(b, "marginal_likelihood", "laplace"), b)
    if (family != "gaussian")
      expect_error(check(b, "marginal_likelihood", "factorized"), "Gaussian response")
    else expect_error(check(b, "marginal_likelihood", "factorized"), "coupling")
  }
})

test_that("marginal posterior requires one total scale SD, rather than one kernel", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  # Reversing the bins supplies a distinct landscape, rather than a collinear
  # rescaling of habitat that the existing design-rank check correctly rejects.
  columns <- rev(seq_len(ncol(f$inputs$binned$vsum$habitat)))
  f$inputs$binned$vsum$cover <- f$inputs$binned$vsum$habitat[, columns]
  f$inputs$binned$csum$cover <- f$inputs$binned$csum$habitat[, columns]
  prepare <- function(scales) multiScale_bayes(count ~ habitat + cover + (1 | site),
    f$inputs, f$data, scales, row_id = "key", method = "prepare")
  one <- prepare(list(habitat = ~1 + (1 | year), cover = ~1))
  two <- prepare(list(habitat = ~1 + (1 | year), cover = ~1 + (1 | year)))
  none <- prepare(list(habitat = ~1, cover = ~1))
  check <- multiScaleR:::.msr_scale_check_capability
  for (fit in list(one, two, none)) for (method in c("prepare", "nuts", "joint_laplace"))
    expect_identical(check(fit$specification, "posterior", method), fit$specification)
  for (method in c("marginal_laplace", "importance_corrected")) {
    expect_identical(check(one$specification, "posterior", method), one$specification)
    expect_error(check(two$specification, "posterior", method), "exactly one scale-group SD")
    expect_error(check(none$specification, "posterior", method), "exactly one scale-group SD")
  }
  for (method in c("laplace", "factorized"))
    expect_error(check(one$specification, "marginal_likelihood", method), "one kernel")
  # Unsupported public requests must fail before dispatch or any optional solver.
  local_mocked_bindings(.msr_bayes_marginal = function(...) stop("solver was called"),
    .package = "multiScaleR")
  expect_error(multiScale_bayes(count ~ habitat + cover, f$inputs, f$data,
    scales = list(habitat = ~1 + (1 | year), cover = ~1 + (1 | year)),
    row_id = "key", method = "marginal"), "exactly one scale-group SD")
})

test_that("factorization checks nesting by observed membership rather than group names", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  f$data$replicate <- interaction(f$data$year, rep(c("a", "a", "b", "b"), 3), drop = TRUE)
  f$data$alias <- f$data$year
  prepare <- function(response, scales = list(habitat = ~1 + (1 | year)))
    multiScaleR:::.msr_marginal_prepare(response, f$inputs, f$data, scales, "gaussian", "key", 200)
  nested <- prepare(count ~ habitat + (1 | replicate))
  alias <- prepare(count ~ habitat + (1 | alias))
  crossed <- prepare(count ~ habitat + (1 | site))
  shared <- prepare(count ~ habitat + (1 | site), list(habitat = ~1))
  check <- multiScaleR:::.msr_scale_check_capability
  for (spec in list(nested, alias, shared)) {
    expect_identical(check(spec, "marginal_likelihood", "factorized"), spec)
    expect_identical(multiScaleR:::.msr_factorized_validate(spec), spec)
  }
  expect_identical(check(crossed, "marginal_likelihood", "laplace"), crossed)
  expect_error(check(crossed, "marginal_likelihood", "factorized"), "coupling")
  zero <- multiScaleR:::.msr_marginal_reduce(crossed, "response")
  expect_identical(check(zero, "marginal_likelihood", "factorized"), zero)
  orphan <- zero; orphan$parameters$response_z <- 0
  expect_error(check(orphan, "marginal_likelihood", "factorized"), "no corresponding grouping term")
  unsupported <- nested; unsupported$family <- "custom"
  expect_error(check(unsupported, "marginal_likelihood", "factorized"), "Gaussian response")
  # Generic family/kernel validation remains with the neutral preparer. Private
  # optimizer refinements can inspect structure without a complete data model.
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    family = "custom", row_id = "key", method = "prepare"), "arg")
  f$inputs$kernel <- "expow"
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", method = "prepare"), "binned Gaussian or exponential")
})

test_that("private likelihood solvers reject crossed terms before optimization", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  fit <- prepare_bayesian_fixture(f,
    scales = list(habitat = ~1 + (1 | year) + (1 | site)))
  for (method in c("prepare", "nuts", "joint_laplace"))
    expect_identical(multiScaleR:::.msr_scale_check_capability(fit$specification,
      "posterior", method), fit$specification)
  local_mocked_bindings(.msr_marginal_attempt = function(...) stop("optimizer was called"),
    .package = "multiScaleR")
  expect_error(multiScaleR:::.msr_marginal_fit(fit$specification), "one kernel and at most one scale")
  expect_error(multiScaleR:::.msr_factorized_fit(fit$specification), "one kernel and at most one scale")
  expect_error(multiScaleR:::.msr_marginal_prepare(count ~ habitat, f$inputs, f$data,
    scales = list(habitat = ~1 + (1 | year) + (1 | site)), row_id = "key"),
    "private marginal prototype supports")
})

test_that("importance proposal density includes the log-SD Jacobian", {
  node <- list(mode = c(scale_z = 0, log_scale_sd = 0), root = matrix(1),
    covariance = matrix(1), hessian = matrix(1))
  grid <- multiScaleR:::.msr_bayes_marginal_weights(c(0, .5, 1), c(0, 0, 0))
  grid$nodes <- rep(list(node), 3); grid$prior_sd <- .5
  samples <- cbind(c(-1, 0, 2), log(c(.1, .5, 1)))
  expected <- dt(samples[, 1] / 1.2, df = 7, log = TRUE) - log(1.2) +
    dnorm(exp(samples[, 2]), 0, .5, log = TRUE) + log(2) + samples[, 2]
  expect_equal(multiScaleR:::.msr_bayes_importance_logq(samples, grid, 2), expected,
    tolerance = 1e-12)
  expect_equal(multiScaleR:::.msr_bayes_logadd(c(-Inf, -1000), c(-Inf, -1000)),
    c(-Inf, -1000 + log(2)))
  density <- function(x) exp(multiScaleR:::.msr_bayes_importance_density(
    matrix(x, ncol = 1), list(mode = 0, hessian = matrix(1))))
  expect_equal(integrate(density, -Inf, Inf)$value, 1, tolerance = 1e-6)
})

test_that("importance correction recovers a skewed posterior and preserves RNG", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("posterior", minimum_version = "1.6.0")
  parameters <- list(beta = numeric(), scale_beta = numeric(), scale_z = log(2),
    log_scale_sd = log(.5), response_z = numeric(), log_response_sd = numeric(),
    log_dispersion = numeric())
  # exp(scale_z) is Gamma(shape=2, rate=1), independently of a half-normal SD.
  objective <- function(p) exp(p$scale_z) - 2 * p$scale_z + lgamma(2) -
    RTMB::dnorm(exp(p$log_scale_sd), 0, .5, log = TRUE) - log(2) - p$log_scale_sd
  obj <- RTMB::MakeADFun(objective, parameters, silent = TRUE)
  ctl <- list(starts = 1L, iter.max = 1000L, max_grid = 65L,
    cdf_tol = .005, quantile_tol = .02, logz_tol = .01,
    importance_max = 16384L, importance_refinements = 0L, nested_response = FALSE)
  evaluate <- function(x, warm) {
    tau <- if (x == 0) .5e-8 else if (x == 1) 4 else .5 * qnorm((1 + x) / 2)
    multiScaleR:::.msr_bayes_marginal_node(obj, 2, log(tau), .5, warm, ctl)
  }
  grid <- multiScaleR:::.msr_bayes_marginal_grid(evaluate, .5, ctl)
  grid$hyper_index <- 2L
  out <- list(parameters = parameters, specification = list(parameters = parameters,
    scales = list(), response = list(random = list())), marginal = grid,
    diagnostics = list(passed = TRUE))
  set.seed(193); before <- .Random.seed
  fit <- multiScaleR:::.msr_bayes_importance(out, obj, 2000L, 914L, ctl)
  expect_identical(fit$inference_metadata$method, "importance_corrected")
  expect_identical(tail(names(fit), 1), "inference_metadata")
  expect_identical(fit$inference_metadata$integration$population_parameters, "joint_importance_sampling")
  expect_identical(.Random.seed, before)
  expect_true(fit$diagnostics$passed)
  expect_equal(sum(fit$importance$weights), 1, tolerance = 1e-12)
  candidate <- fit$importance$candidates
  weight <- fit$importance$weights
  expect_equal(sum(weight * candidate[, 1]), digamma(2), tolerance = .025)
  expect_equal(sum(weight * exp(candidate[, 1])), 2, tolerance = .035)
  expect_equal(multiScaleR:::.msr_bayes_weighted_quantile(exp(candidate[, 1]),
    weight, c(.025, .5, .975)), qgamma(c(.025, .5, .975), 2), tolerance = .04)
  expect_equal(sum(weight * exp(candidate[, 2])), .5 * sqrt(2 / pi), tolerance = .015)
  expect_true(all(fit$importance$resample_ids %in% seq_len(nrow(candidate))))
  expect_equal(fit$draws$scale_z[, 1], candidate[fit$importance$resample_ids, 1])
  # A maximum is a cap, whereas an explicit minimum must retain the larger
  # population even when the original small batch already passes.
  larger <- ctl; larger$importance_min <- larger$importance_max <- 32768L
  more <- multiScaleR:::.msr_bayes_importance(out, obj, 2000L, 914L, larger)
  expect_true(more$diagnostics$passed)
  expect_equal(nrow(more$importance$candidates), 32768L)
  expect_equal(more$diagnostics$min_candidates, 32768L)
  expect_false(more$importance$history$budget_met[1])
  expect_true(tail(more$importance$history$budget_met, 1))
  explicit <- ctl; explicit$importance_min <- 16384L; explicit$importance_df <- 7
  same <- multiScaleR:::.msr_bayes_importance(out, obj, 2000L, 914L, explicit)
  expect_identical(same$importance$candidates, fit$importance$candidates)
  expect_identical(same$importance$weights, fit$importance$weights)
  expect_identical(same$draws, fit$draws)
  heavy <- ctl; heavy$importance_df <- 5
  alternative <- multiScaleR:::.msr_bayes_importance(out, obj, 2000L, 914L, heavy)
  expect_true(alternative$diagnostics$passed)
  expect_equal(alternative$diagnostics$proposal_df, 5)
  expect_equal(sum(alternative$importance$weights * exp(alternative$importance$candidates[, 1])),
    2, tolerance = .035)
  expect_true(all(vapply(alternative$importance$adapted_proposal$components,
    function(x) identical(x$df, 5), logical(1))))
})

test_that("batched predictors and conditional response curvature match the posterior", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 6, n_year = 3)
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    prepared <- prepare_bayesian_fixture(f, family = family)
    spec <- prepared$specification
    p <- prepared$parameters
    p$beta <- seq_along(p$beta) / 10
    p$log_response_sd[] <- log(.4)
    p$scale_z <- c(-.3, .1, .4)
    samples <- matrix(rep(unlist(p), each = 3), nrow = 3)
    samples[, 1] <- samples[, 1] + c(-.1, 0, .1)
    blocks <- multiScaleR:::.msr_bayes_parameter_blocks(samples, p)
    eta <- multiScaleR:::.msr_bayes_batch_eta(spec, blocks)
    for (i in 1:3) {
      q <- lapply(blocks, function(x) unname(x[i, ]))
      expect_equal(eta[i, ], multiScaleR:::.msr_bayes_components(spec, q,
        include_response = FALSE)$eta, tolerance = 1e-12)
    }
    geometry <- multiScaleR:::.msr_bayes_response_geometry(spec, blocks)
    obj <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(spec, ad = TRUE), p, silent = TRUE)
    random <- which(names(obj$par) == "response_z")
    for (i in 1:3) {
      q <- samples[i, ]; q[random] <- geometry$mode[i, ]
      expect_lt(max(abs(obj$gr(q)[random])), 1e-5)
      expect_equal(diag(obj$he(q)[random, random]), unname(geometry$hessian[i, ]),
        tolerance = 1e-9)
    }
  }
})

test_that("importance controls fail explicitly for unsupported settings", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture()
  call_fit <- function(control) multiScale_bayes(count ~ habitat, f$inputs, f$data,
    scales = list(habitat = ~1 + (1 | year)), row_id = "key",
    method = "marginal", control = control)
  expect_error(call_fit(list(correction = "other")), "correction")
  expect_error(call_fit(list(importance_max = 100)), "importance_max")
  expect_error(call_fit(list(importance_max = "65536")), "numeric")
  expect_error(call_fit(list(importance_min = 100)), "importance_min")
  expect_error(call_fit(list(importance_min = "65536")), "numeric")
  expect_error(call_fit(list(importance_min = 131072, importance_max = 65536)), "must not exceed")
  expect_error(call_fit(list(importance_refinements = 3)), "importance_refinements")
  expect_error(call_fit(list(importance_df = 6)), "importance_df")
  expect_error(call_fit(list(importance_df = "5")), "numeric")
  expect_error(call_fit(list(importance_df = NA)), "importance_df")
  expect_error(call_fit(list(nested_response = NA)), "nested_response")
  expect_error(call_fit(list(nested_response = TRUE)), "requires correction")
})

test_that("conditional proposal fallbacks stay normalized and retain every row", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  prepared <- prepare_bayesian_fixture(f, family = "gaussian")
  spec <- prepared$specification
  samples <- matrix(rep(unlist(prepared$parameters), each = 2), nrow = 2)
  blocks <- multiScaleR:::.msr_bayes_parameter_blocks(samples, prepared$parameters)
  # An overflowing conditional curvature falls back to a normalized prior
  # proposal for that row. The importance target still decides its contribution.
  blocks$log_response_sd[1, ] <- 400
  geometry <- suppressWarnings(multiScaleR:::.msr_bayes_response_geometry(spec, blocks))
  expect_equal(nrow(geometry$mode), 2)
  expect_equal(geometry$fallbacks, 1)
  expect_equal(unname(geometry$mode[1, ]), rep(0, 3))
  expect_equal(unname(geometry$hessian[1, ]), rep(1, 3))
  expect_true(all(is.finite(geometry$mode)))
  expect_true(all(is.finite(geometry$hessian)))
})

test_that("corrected repeated-site nodes avoid dense site covariance storage", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("posterior", minimum_version = "1.6.0")
  f <- bayesian_scale_fixture(n_site = 6, n_year = 3)
  prepared <- prepare_bayesian_fixture(f, family = "gaussian")
  p <- prepared$parameters
  p$beta <- c(.5, .3, 1.5)
  p$log_response_sd[] <- log(.3)
  f$data$count <- multiScaleR:::.msr_bayes_components(prepared$specification, p)$eta +
    .3 * sin(seq_len(nrow(f$data)))
  set.seed(196); before <- .Random.seed
  fit <- suppressWarnings(multiScale_bayes(
    count ~ habitat + temp + (1 | site) + offset(log(effort)), f$inputs, f$data,
    family = "gaussian", scales = list(habitat = ~1 + (1 | year)), row_id = "key",
    method = "marginal", seed = 831,
    control = list(correction = "importance", importance_max = 16384,
                   importance_refinements = 0)))
  expect_identical(.Random.seed, before)
  total <- sum(lengths(fit$parameters))
  outer <- total - length(fit$parameters$response_z)
  expect_equal(length(fit$marginal$nodes[[1]]$mode), outer)
  expect_equal(nrow(fit$marginal$nodes[[1]]$covariance), outer - 1)
  expect_equal(length(fit$marginal$parameter_indices), outer)
  expect_false(any(names(fit$marginal$nodes[[1]]$mode) == "response_z"))
  expect_equal(dim(fit$draws$response_z), c(2000L, 6L))
  expect_true(fit$importance$conditional_response)
  expect_equal(nrow(posterior_scales(fit)), 3)
  expect_equal(posterior_scales(unserialize(serialize(fit, NULL))), posterior_scales(fit))
})

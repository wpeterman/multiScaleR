test_that("conditional response coefficient curvature matches the full posterior", {
  skip_if_not_installed("RTMB")
  fixture <- bayesian_scale_fixture(n_site = 5, n_year = 3)
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- fixture
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    prepared <- prepare_bayesian_fixture(f, family = family)
    spec <- prepared$specification; parameters <- prepared$parameters
    parameters$scale_z[] <- c(-.3, .2, .5)
    parameters$response_z[] <- seq_along(parameters$response_z) / 10
    component <- multiScaleR:::.msr_bayes_components(spec, parameters)
    X <- cbind(spec$X_fixed, do.call(cbind, component$kernel))
    offset <- spec$offset
    for (term in spec$response$random)
      offset <- offset + parameters$response_z[term$z_index][term$index] *
        exp(parameters$log_response_sd[term$sd_index])
    dispersion <- if (length(parameters$log_dispersion)) parameters$log_dispersion[1] else 0
    geometry <- multiScaleR:::.msr_bayes_beta_geometry(X, spec$y, spec$trials, offset,
      family, dispersion, spec$priors$beta_sd, parameters$beta)
    expect_false(geometry$fallback)
    parameters$beta <- geometry$mode
    obj <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(spec, ad = TRUE), parameters, silent = TRUE)
    beta <- which(names(obj$par) == "beta")
    expect_equal(unname(obj$he(obj$par)[beta, beta]), unname(geometry$hessian), tolerance = 1e-8)
    expect_lt(max(abs(obj$gr(obj$par)[beta])), 1e-5)
  }
})

test_that("Gaussian coefficient proposals have the exact conditional density", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  prepared <- prepare_bayesian_fixture(f, family = "gaussian")
  spec <- prepared$specification; parameters <- prepared$parameters
  obj <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(spec, ad = TRUE), parameters, silent = TRUE)
  beta <- which(names(obj$par) == "beta"); keep <- setdiff(seq_along(obj$par), beta)
  outer <- matrix(rep(obj$par[keep], each = 8), nrow = 8)
  set.seed(3791)
  samples <- multiScaleR:::.msr_bayes_beta_candidates(outer, keep, spec, parameters$beta)
  replay <- multiScaleR:::.msr_bayes_beta_candidates(outer, keep, spec, parameters$beta,
    given = samples$samples)
  expect_equal(samples$logq, replay$logq, tolerance = 1e-10)
  difference <- vapply(seq_len(nrow(outer)), function(i) -obj$fn(samples$samples[i, ]), 0) - samples$logq
  expect_lt(diff(range(difference)), 1e-8)
  expect_equal(samples$fallbacks, 0)
  # Numerical failure uses a normalized prior proposal and keeps every row.
  dispersion <- match(which(names(obj$par) == "log_dispersion"), keep)
  outer[, dispersion] <- -400
  bad <- suppressWarnings(multiScaleR:::.msr_bayes_beta_candidates(outer, keep, spec, parameters$beta))
  expect_equal(bad$fallbacks, nrow(outer))
  expect_equal(nrow(bad$samples), nrow(outer))
  expect_true(all(is.finite(bad$samples)))
  expect_equal(bad$logq, rowSums(dnorm(bad$samples[, beta, drop = FALSE],
    sd = spec$priors$beta_sd, log = TRUE)), tolerance = 1e-12)
  expect_error(multiScaleR:::.msr_bayes_beta_candidates(outer, keep, spec, parameters$beta,
    given = matrix(0, 1, 1)), "layout")
})

test_that("latent-prior defense has a normalized density including SD Jacobians", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  prepared <- prepare_bayesian_fixture(f)
  parameters <- prepared$parameters
  flat <- unlist(parameters, use.names = FALSE)
  names <- rep(names(parameters), lengths(parameters))
  keep <- which(names != "beta"); hyper <- which(names == "log_scale_sd")
  prepared$marginal <- list(hyper_index = hyper, weights = 1,
    nodes = list(list(mode = flat, covariance = diag(length(flat) - 1))))
  grid <- multiScaleR:::.msr_bayes_latent_prior_proposal(prepared, keep)
  set.seed(5812)
  draws <- multiScaleR:::.msr_bayes_importance_candidates(grid, match(hyper, keep), 40000)
  full <- matrix(0, nrow(draws), length(flat)); full[, keep] <- draws
  b <- multiScaleR:::.msr_bayes_parameter_blocks(full, parameters)
  p <- prepared$specification$priors
  # A different normalized, lighter-tailed target tests both the proposal's
  # sampling law and its density against known independent moments.
  logp <- rowSums(dnorm(b$scale_z, sd = .8, log = TRUE)) +
    rowSums(dnorm(b$response_z, sd = .8, log = TRUE))
  for (sm in grid$scales) for (j in seq_along(sm$beta_index))
    logp <- logp + dnorm(b$scale_beta[, sm$beta_index[j]], sm$prior_mean[j],
      .8 * sm$prior_sd[j], log = TRUE)
  logp <- logp + rowSums(dnorm(exp(b$log_scale_sd), 0, .8 * p$scale_group_sd, log = TRUE) +
    log(2) + b$log_scale_sd)
  logp <- logp + rowSums(dnorm(exp(b$log_response_sd), 0, .8 * p$response_group_sd, log = TRUE) +
    log(2) + b$log_response_sd) + dnorm(b$log_dispersion[, 1], log = TRUE)
  logq <- multiScaleR:::.msr_bayes_importance_logq(draws, grid, match(hyper, keep))
  ratio <- exp(logp - logq); weights <- ratio / sum(ratio)
  expect_equal(mean(ratio), 1, tolerance = .035)
  expect_equal(sum(weights * exp(b$log_scale_sd[, 1])),
    .8 * p$scale_group_sd * sqrt(2 / pi), tolerance = .012)
  expect_equal(sum(weights * exp(b$log_response_sd[, 1])),
    .8 * p$response_group_sd * sqrt(2 / pi), tolerance = .02)
  expect_equal(sum(weights * b$log_dispersion[, 1]), 0, tolerance = .03)
  expect_equal(sum(weights * b$scale_z[, 1]^2), .8^2, tolerance = .03)
})

test_that("conditional coefficient controls reject contradictory settings", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture()
  call_fit <- function(control) multiScale_bayes(count ~ habitat, f$inputs, f$data,
    scales = list(habitat = ~1 + (1 | year)), row_id = "key", method = "marginal", control = control)
  expect_error(call_fit(list(conditional_beta = NA)), "conditional_beta")
  expect_error(call_fit(list(conditional_beta = TRUE)), "requires correction")
  expect_error(call_fit(list(correction = "importance", conditional_beta = TRUE,
    nested_response = TRUE)), "nested_response = FALSE")
})

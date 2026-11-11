test_that("slice sampling fits nonannual scale formulas and reuses posterior methods", {
  f <- bayesian_scale_fixture(n_site = 8L, n_year = 3L, kernel = "exp")
  d <- f$data
  d$region <- factor(rep(rep(c("valley", "ridge", "upland"), length.out = 8), 3))
  d$cover_type <- factor(rep(c("open", "forest"), length.out = nrow(d)))
  d$success <- as.integer(d$count >= 3L)
  d$failure <- 2L - d$success
  b <- f$inputs$binned
  f$inputs$binned$vsum$forest <- b$vsum$habitat[, rev(seq_along(b$bin_rep))]
  f$inputs$binned$csum$forest <- b$csum$habitat[, rev(seq_along(b$bin_rep))]
  response <- cbind(success, failure) ~ habitat * cover_type + forest + (1 | site)
  scales <- list(habitat = ~ temp + (1 | region), forest = ~1)
  run <- function(...) multiScale_bayes(response, f$inputs, d, scales = scales,
    row_id = "key", family = "binomial", method = "slice", ...)
  before <- .Random.seed
  fit <- suppressWarnings(run(chains = 2L, iter = 45L, warmup = 15L, seed = 18L))
  expect_identical(.Random.seed, before)
  same <- suppressWarnings(run(chains = 2L, iter = 45L, warmup = 15L, seed = 18L))
  expect_identical(fit$chains, same$chains)
  expect_equal(dim(fit$chains), c(30L, 2L, sum(lengths(fit$parameters))))
  expect_false(identical(fit$chains[, 1L, ], fit$chains[, 2L, ]))
  expect_true(all(is.finite(fit$chains)))
  expect_null(fit$sampler)
  expect_identical(fit$inference_metadata$solver, "slice")
  expect_identical(fit$inference_metadata$integration$response_effects, "joint_sampling")
  expect_equal(nrow(posterior_scales(fit)), nrow(d) + 1L)
  prediction <- predict(fit, newdata = d, kernel_inputs = f$inputs, draws = 4)
  expect_equal(prediction$mean, predict(fit, draws = 4)$mean)
  expect_false(any(prediction$new_scale_level | prediction$new_response_level))
  expect_equal(dim(simulate(fit, nsim = 3, seed = 12)), c(nrow(d), 3L))
  expect_error(run(chains = 0), "chains")
  expect_error(run(control = list(adapt_delta = 0.95)), "empty control")
  if (requireNamespace("posterior", quietly = TRUE)) {
    expect_true(diagnose_bayes_scales(fit)$available)
    expect_true(diagnose_bayes_probabilities(fit, list(forest = 200))$available)
    expect_false(fit$diagnostics$passed)
  }
})

test_that("slice sampling runs with no optional Bayesian packages", {
  local_mocked_bindings(.msr_bayes_package_status = function(...) list(passed = FALSE),
    .package = "multiScaleR")
  f <- bayesian_scale_fixture(n_site = 6L, n_year = 1L)
  expect_true(check_bayes_setup("slice")$ready)
  for (family in c("nbinom2", "poisson")) {
    fit <- multiScale_bayes(count ~ habitat + offset(log(effort)), f$inputs, f$data,
      family = family, row_id = "key", method = "slice", chains = 1L, iter = 12L, seed = 55)
    expect_true(all(is.finite(fit$chains)))
    expect_true(is.na(fit$diagnostics$passed))
    expect_match(fit$diagnostics$reason, "posterior")
  }
})

test_that("shared-scale bin caching preserves the posterior kernel calculation", {
  for (kernel in c("gaussian", "exp")) {
    f <- bayesian_scale_fixture(n_site = 6L, n_year = 3L, annual = FALSE, kernel = kernel)
    fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key",
      method = "prepare")
    s <- fit$specification
    p <- fit$parameters
    bins <- .msr_bayes_slice_shared_bins(s)
    expect_gt(length(bins[[1L]]), 1L)
    for (log_s in c(-500, -8, -2, 0, 7, 500)) {
      p$scale_beta <- log_s
      expect_equal(.msr_bayes_slice_kernel(s, p, 1L, bins),
        .msr_bayes_components(s, p)$kernel$habitat, tolerance = 1e-13)
    }
  }
})

test_that("the slice scale posterior agrees with independent Gaussian quadrature", {
  f <- bayesian_scale_fixture(n_site = 24L, n_year = 1L)
  f$data$count <- rnorm(nrow(f$data), sd = .45)
  formula <- count ~ habitat + offset(log(effort))
  priors <- list(scale_intercept_sd = .65, beta_sd = 1.5,
    dispersion_log_mean = log(.45), dispersion_log_sd = .25)
  prepare <- multiScale_bayes(formula, f$inputs, f$data, row_id = "key", priors = priors,
    family = "gaussian", method = "prepare")
  s <- prepare$specification
  p <- prepare$parameters
  # Integrate the independent normal regression priors out analytically.
  # The resulting Gaussian marginal has covariance s^2 I + beta_sd^2 X X'.
  residual <- s$y - s$offset
  scale_grid <- seq(p$scale_beta - 3, p$scale_beta + 3, length.out = 151)
  sd_grid <- seq(p$log_dispersion - 1.5, p$log_dispersion + 1.5, length.out = 101)
  density <- vapply(scale_grid, function(a) {
    p$scale_beta <- a
    raw <- .msr_bayes_components(s, p)$kernel$habitat
    X <- cbind(1, raw)
    covariance <- priors$beta_sd^2 * tcrossprod(X)
    log_joint <- vapply(sd_grid, function(b) {
      C <- chol(covariance + diag(exp(2 * b), length(residual)))
      -sum(log(diag(C))) - .5 * sum(forwardsolve(t(C), residual)^2) +
        dnorm(b, priors$dispersion_log_mean, priors$dispersion_log_sd, log = TRUE)
    }, numeric(1))
    sum(exp(log_joint)) * dnorm(a, s$parameters$scale_beta, priors$scale_intercept_sd)
  }, numeric(1))
  cumulative <- cumsum(density) / sum(density)
  expected <- stats::approx(cumulative, scale_grid, xout = c(.1, .5, .9))$y
  fit <- suppressWarnings(multiScale_bayes(formula, f$inputs, f$data, row_id = "key",
    priors = priors, family = "gaussian", method = "slice", chains = 2L,
    iter = 2200L, warmup = 200L, seed = 361L))
  actual <- unname(stats::quantile(fit$draws$scale_beta, c(.1, .5, .9)))
  expect_lt(max(abs(expected - actual)), .10)
})


test_that("the Gaussian reference leaves the binomial coefficient target unchanged", {
  set.seed(3141)
  y <- c(1, 4, 6)
  trials <- c(8, 9, 10)
  state <- list(X = matrix(c(-.5, 1, 2), ncol = 1), offset = c(-.2, .1, .3))
  prior_sd <- 1.5
  loglik <- function(eta) sum(dbinom(y, trials, plogis(eta), log = TRUE))
  # Deliberately shift and narrow the Gaussian reference to test its correction.
  geometry <- list(mode = 0.4, R = matrix(4, 1, 1), fallback = FALSE)
  grid <- seq(-5, 5, length.out = 2001)
  density <- vapply(grid, function(b)
    exp(loglik(state$offset + drop(state$X * b))) * dnorm(b, sd = prior_sd), 0)
  expected <- sum(grid * density) / sum(density)
  beta <- 0
  state$ll <- loglik(state$offset)
  draws <- numeric(5000)
  for (i in seq_along(draws)) {
    update <- .msr_bayes_slice_beta(beta, state, geometry, prior_sd, loglik)
    draws[i] <- beta <- update$beta
    state$ll <- update$likelihood
  }
  expect_lt(abs(mean(draws[-seq_len(500)]) - expected), .07)
})

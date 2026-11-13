test_that("event assessment respects overlap and simulation uncertainty", {
  n <- 10000L; w <- rep(1 / n, n)
  first <- rep(c(TRUE, FALSE), n / 2L)
  second <- rep(c(TRUE, FALSE, FALSE, FALSE), n / 4L)
  assess <- .msr_laplace_audit_event
  expect_identical(assess(first, first, w, .05, TRUE)$status, "within_tolerance")
  expect_identical(assess(first, second, w, .05, TRUE)$status, "discrepancy")
  expect_identical(assess(second, first, w, .05, TRUE)$status, "discrepancy")
  expect_identical(assess(first, first, w, .05, FALSE)$status, "unresolved")
  expect_identical(assess(first, rep(TRUE, n), w, .05, TRUE)$status, "unresolved")
  expect_identical(assess(first, c(rep(TRUE, n / 2L), rep(FALSE, n / 2L)),
    w, .05, TRUE)$status, "unresolved")
  expect_identical(assess(first, first, w, .001, TRUE)$status, "unresolved")
})

test_that("the diagnostic preserves the original fit and caller RNG", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("posterior", minimum_version = "1.6.0")
  f <- bayesian_scale_fixture(n_site = 9L, n_year = 3L)
  d <- f$data
  fit <- multiScale_bayes(count ~ habitat, f$inputs, d, family = "gaussian",
    row_id = "key", method = "laplace", laplace_draws = 2500L, seed = 817L)
  original <- fit
  before <- .Random.seed
  result <- diagnose_bayes_laplace(fit, list(habitat = c(150, 250, 350)), candidates = 2048L, seed = 84L)
  expect_identical(before, .Random.seed)
  expect_identical(fit, original)
  expect_identical(result$raw_diagnostics, fit$diagnostics)
  expect_false(result$metadata_establishes_eligibility)
  expect_true(result$available, info = result$reason)
  expect_equal(result$target_evaluation[["evaluations"]], 2048)
  expect_equal(result$probabilities$approximation,
    vapply(c(150, 250, 350), function(t)
      mean(fit$draws$scale_beta[, 1] <= log(t / fit$specification$unit_conv)), 0))
  expect_identical(result, diagnose_bayes_laplace(fit, list(habitat = c(150, 250, 350)),
    candidates = 2048L, seed = 84L))
  expect_error(diagnose_bayes_laplace(fit, list(habitat = c(150, 150))), "distinct finite")

  # Failed density evaluations remain explicit, with the original approximation.
  local_mocked_bindings(.msr_bayes_importance_target = function(samples, obj, spec)
    list(log_density = rep(NA_real_, nrow(samples)), diagnostics = c(unresolved = nrow(samples))),
    .package = "multiScaleR")
  invalid <- diagnose_bayes_laplace(fit, list(habitat = 150), candidates = 2000L, seed = 84L)
  expect_false(invalid$available)
  expect_match(invalid$reason, "No candidate was dropped")
  expect_identical(invalid$probabilities$status, "unresolved")
  expect_true(is.finite(invalid$probabilities$approximation))
  expect_true(is.na(invalid$probabilities$reference))
})

test_that("fresh weighted scale reference agrees with independent quadrature", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("posterior", minimum_version = "1.6.0")
  f <- bayesian_scale_fixture(n_site = 18L, n_year = 1L)
  f$data$count <- rnorm(nrow(f$data), sd = .45)
  priors <- list(scale_median = 250, scale_intercept_sd = .65, beta_sd = 1.5,
    dispersion_log_mean = log(.45), dispersion_log_sd = .25)
  fit <- multiScale_bayes(count ~ habitat + offset(log(effort)), f$inputs, f$data,
    row_id = "key", priors = priors, family = "gaussian", method = "laplace",
    laplace_draws = 5000L, seed = 361L)
  spec <- fit$specification
  b <- f$inputs$binned
  scale <- seq(log(.25) - 3.5, log(.25) + 3.5, length.out = 281L)
  log_sd <- seq(log(.45) - 1.5, log(.45) + 1.5, length.out = 121L)
  y <- f$data$count - log(f$data$effort)
  normalized <- function(a) {
    w <- exp(-.5 * (b$bin_rep / exp(a))^2)
    drop(b$vsum$habitat %*% w / (b$csum$habitat %*% w))
  }
  ref <- normalized(log(.25))
  density <- vapply(scale, function(a) {
    x <- (normalized(a) - mean(ref)) / sd(ref)
    V <- priors$beta_sd^2 * tcrossprod(cbind(1, x))
    v <- vapply(log_sd, function(s) {
      R <- chol(V + diag(exp(2 * s), length(y)))
      -.5 * sum(forwardsolve(t(R), y)^2) - sum(log(diag(R))) +
        dnorm(s, log(.45), .25, log = TRUE)
    }, numeric(1))
    sum(exp(v)) * dnorm(a, log(.25), .65)
  }, numeric(1))
  reference <- stats::approx(scale, cumsum(density) / sum(density),
    xout = log(c(.15, .25, .35)))$y
  result <- diagnose_bayes_laplace(fit, list(habitat = c(150, 250, 350)), seed = 187L)
  expect_true(result$available, info = result$reason)
  expect_lt(max(abs(result$probabilities$reference - reference)), .045)
})

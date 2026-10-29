test_that("correlated proposal sampling and density agree in whitened coordinates", {
  H <- matrix(c(1, .9, .9, 1), 2)
  geometry <- multiScaleR:::.msr_bayes_laplace_geometry(H)
  node <- c(list(mode = c(0, 0, 0), hessian = H), geometry[c("root", "covariance")])
  grid <- multiScaleR:::.msr_bayes_marginal_weights(c(0, .5, 1), c(0, 0, 0))
  grid$nodes <- rep(list(node), 3); grid$prior_sd <- .5
  set.seed(338)
  samples <- multiScaleR:::.msr_bayes_importance_candidates(grid, 3L, 50000L)
  u <- samples[, 1:2]
  white <- u %*% t(chol(H)) / 1.2
  # Both shared- and coordinatewise-scale t components have these exact
  # univariate marginals. Dividing by coordinate scales after correlating
  # normals breaks this property and the joint proposal density.
  expect_equal(colMeans(white^2), rep(7/5, 2), tolerance = .08)
  expect_equal(colMeans(white <= qt(.975, 7)), rep(.975, 2), tolerance = .004)
  logq <- multiScaleR:::.msr_bayes_importance_logq(samples, grid, 3L)
  logp <- .5 * as.numeric(determinant(H, logarithm = TRUE)$modulus) - log(2*pi) -
    rowSums(u * (u %*% H)) / 2 +
    dnorm(exp(samples[, 3]), 0, .5, log = TRUE) + log(2) + samples[, 3]
  ratio <- exp(logp - logq)
  expect_equal(mean(ratio), 1, tolerance = .02)
  weights <- ratio / sum(ratio)
  expect_equal(colSums(u * weights), c(0, 0), tolerance = .025)
  expect_equal(crossprod(u * sqrt(weights)), solve(H), tolerance = .03)
})

test_that("heavier-tail proposals and mixed degrees of freedom have exact densities", {
  H <- matrix(c(1, .7, .7, 1), 2)
  geometry <- multiScaleR:::.msr_bayes_laplace_geometry(H)
  node <- c(list(mode = c(0, 0, 0), hessian = H), geometry[c("root", "covariance")])
  grid <- multiScaleR:::.msr_bayes_marginal_weights(c(0, .5, 1), c(0, 0, 0))
  grid$nodes <- rep(list(node), 3); grid$prior_sd <- .5; grid$df <- 5
  original <- grid; original$df <- 7
  mixture <- list(components = list(original, grid), weights = c(.3, .7))
  set.seed(9147)
  samples <- multiScaleR:::.msr_bayes_importance_candidates(grid, 3L, 50000L, mixture)
  u <- samples[, 1:2]
  logq <- multiScaleR:::.msr_bayes_importance_logq(samples, grid, 3L, mixture)
  logp <- .5 * as.numeric(determinant(H, logarithm = TRUE)$modulus) - log(2*pi) -
    rowSums(u * (u %*% H)) / 2 +
    dnorm(exp(samples[, 3]), 0, .5, log = TRUE) + log(2) + samples[, 3]
  ratio <- exp(logp - logq); weights <- ratio / sum(ratio)
  expect_equal(mean(ratio), 1, tolerance = .02)
  expect_equal(crossprod(u * sqrt(weights)), solve(H), tolerance = .035)
  adapted <- multiScaleR:::.msr_bayes_importance_adapt(samples, logp - logq, grid, 3L)
  expect_identical(adapted$df, 5)
  point <- cbind(0, 0, log(5))
  grid$tail_fraction <- .01
  expect_true(is.finite(multiScaleR:::.msr_bayes_importance_logq(point, grid, 3L)))
})

test_that("the log-SD proposal tail extends beyond numerical prior-CDF inversion", {
  node <- list(mode = c(0, 0), root = matrix(1), covariance = matrix(1), hessian = matrix(1))
  grid <- multiScaleR:::.msr_bayes_marginal_weights(c(0, .5, 1), c(0, 0, 0))
  grid$nodes <- rep(list(node), 3); grid$prior_sd <- .5; grid$tail_fraction <- .01
  sample <- cbind(c(-1, 0, 2), log(c(.1, .5, 5)))
  h <- sample[, 2]
  prior_h <- 2 * dnorm(exp(h), 0, .5) * exp(h)
  expected <- dt(sample[, 1] / 1.2, 7, log = TRUE) - log(1.2) +
    log(.99 * prior_h + .01 * dnorm(h, log(1), .6))
  expect_equal(multiScaleR:::.msr_bayes_importance_logq(sample, grid, 2L), expected,
    tolerance = 1e-12)
  expect_true(all(is.finite(expected)))
})

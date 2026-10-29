warm_benchmark_helpers <- function() {
  root <- normalizePath(file.path('..', '..'), winslash = '/', mustWork = FALSE)
  path <- file.path(root, 'tools/benchmarks/joint_scale_assessment/phase13/warm.R')
  skip_if_not(file.exists(path), 'Development benchmark helper is not packaged')
  e <- new.env(parent = globalenv()); sys.source(path, e); e
}

test_that('warm conditional coefficients match cold NB2 and Poisson solves', {
  e <- warm_benchmark_helpers()
  set.seed(523)
  X <- cbind(1, rnorm(100), rnorm(100)); colnames(X) <- c('intercept', 'a', 'b')
  offset <- log(runif(100, .5, 2)); weights <- runif(100, .3, 2)
  y <- rnbinom(100, mu = exp(drop(X %*% c(.5, .3, -.2)) + offset), size = .4)
  norm <- .msr_nb2_normalizer(y, weights)
  w <- e$warm_coefficients()
  for (alpha in c(0, 1e-12, 1e-8, .001, .5, 3, 10, 0)) {
    z <- w$solve(X, y, offset, weights, alpha, norm)
    cold <- .msr_nb2_coefficients(X, y, offset, weights, alpha, norm)
    expect_equal(z$beta, cold$beta, tolerance = 1e-6)
    expect_equal(z$terms$value, cold$terms$value, tolerance = 1e-10)
    expect_lte(max(abs(crossprod(X, z$terms$score))), 1e-6)
  }
  expect_gt(unname(w$state$counts['warm']), 0)
  w$state$beta <- rep(Inf, ncol(X))
  z <- w$solve(X, y, offset, weights, .5, norm)
  expect_equal(unname(w$state$counts['fallback']), 1)
  expect_equal(unname(w$state$failures['invalid warm coefficients']), 1)
  w$reset(); expect_null(w$state$beta)
  expect_equal(unname(w$state$counts['resets']), 1)
  Xbad <- cbind(X, X[, 2])
  w$state$beta <- rep(0, ncol(Xbad))
  expect_error(w$solve(Xbad, y, offset, weights, .5, norm), 'rank-deficient')
  expect_null(w$state$beta)
  expect_equal(unname(w$state$counts['fallback']), 2)
  # A deliberately exhausted warm budget must invoke the unchanged cold solve.
  short <- e$warm_coefficients(maxit = 0L); short$state$beta <- rep(50, ncol(X))
  z <- short$solve(X, y, offset, weights, .5, norm)
  expect_equal(unname(short$state$counts['fallback']), 1)
  expect_equal(z$beta, .msr_nb2_coefficients(X, y, offset, weights, .5, norm)$beta)
})

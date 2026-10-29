policy18 <- function() {
  path <- file.path('..', '..', 'tools/benchmarks/joint_scale_assessment/phase18/policy.R')
  skip_if_not(file.exists(path), 'Development-only benchmark policy')
  e <- new.env(parent = globalenv()); sys.source(path, e); e
}
test_that('stable information reporting preserves raw evidence and safeguards', {
  e <- policy18()
  a <- list(available = FALSE, boundary = TRUE, weak = TRUE, multiple_modes = TRUE,
    poisson_limit = FALSE, poisson_compatible = FALSE,
    reasons = c('parameter at an optimization bound',
      'nonpositive or ill-conditioned observed information', 'unknown diagnostic'))
  b <- a; b$reasons[2] <- 'log-kernel standard error exceeds one'
  x <- e$stable_information18(a); y <- e$stable_information18(b)
  expect_identical(x$raw, a); expect_identical(y$raw, b)
  expect_identical(x[-1], y[-1]); expect_false(x$available)
  expect_true('unknown diagnostic' %in% x$categories)
  b$poisson_compatible <- TRUE
  expect_false(identical(x[-1], e$stable_information18(b)[-1]))
})
test_that('cold fallback rejects rank deficiency and preserves successful reference', {
  e <- policy18(); X <- cbind(1, seq(-1, 1, length.out = 20))
  y <- rep(c(0, 1, 2, 3), 5); w <- rep(1, 20); offset <- rep(0, 20)
  norm <- .msr_nb2_normalizer(y, w)
  a <- .msr_nb2_coefficients(X, y, offset, w, .5, norm)
  expect_identical(e$cold_coefficients18(X, y, offset, w, .5, norm), a)
  expect_error(e$cold_coefficients18(cbind(X, X[, 2]), y, offset, w, .5, norm),
    'rank-deficient')
})

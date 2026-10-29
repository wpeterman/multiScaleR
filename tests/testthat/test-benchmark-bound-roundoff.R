roundoff_benchmark_helpers <- function() {
  root <- normalizePath(file.path('..', '..'), winslash = '/', mustWork = FALSE)
  path <- file.path(root, 'tools/benchmarks/joint_scale_assessment/phase15/bounds.R')
  skip_if_not(file.exists(path), 'Development benchmark helper is not packaged')
  e <- new.env(parent = globalenv()); sys.source(path, e); e
}
test_that('only generated starts receive auditable roundoff projection', {
  e <- roundoff_benchmark_helpers(); eps <- .Machine$double.eps
  lo <- c(-1, 0); hi <- c(2, 3)
  x <- c(-1 - eps, 3 + 4 * eps)
  a <- e$check_generated_start(x, lo, hi, TRUE)
  expect_true(a$valid); expect_identical(a$original, x)
  expect_identical(a$adjusted, c(-1, 3)); expect_identical(a$changed, 1:2)
  expect_equal(a$magnitude, 4 * eps)
  b <- e$check_generated_start(x, lo, hi, FALSE)
  expect_false(b$valid); expect_identical(b$adjusted, x)
  expect_length(b$changed, 0)
  a <- e$check_generated_start(c(0, 3), lo, hi, TRUE)
  expect_true(a$valid); expect_length(a$changed, 0)
  expect_identical(a$adjusted, c(0, 3))
})
test_that('material and nonfinite violations remain failures', {
  e <- roundoff_benchmark_helpers()
  for (x in list(c(-1.001, 0), c(0, 3.001), c(NA_real_, 1), c(Inf, 1))) {
    a <- e$check_generated_start(x, c(-1, 0), c(2, 3), TRUE)
    expect_false(a$valid); expect_identical(a$adjusted, x)
    expect_length(a$changed, 0)
  }
  # A distant upper bound cannot inflate tolerance near a lower bound.
  a <- e$check_generated_start(-1 - 1e-8, -1, 1e8, TRUE)
  expect_false(a$valid)
  expect_match(a$reason, 'roundoff tolerance')
  # Mixed tiny/material violations reject the whole start without partial repair.
  x <- c(-1 - .Machine$double.eps, 3.001)
  a <- e$check_generated_start(x, c(-1, 0), c(2, 3), TRUE)
  expect_false(a$valid); expect_identical(a$adjusted, x)
})

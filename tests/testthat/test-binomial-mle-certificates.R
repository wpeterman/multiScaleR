test_that("finite geometry certificates are exact and bound to the represented target", {
  skip_if_not_installed("gmp")
  skip_if_not_installed("lpSolve")
  X <- cbind(intercept = 1, x = c(-2, -1, 1, 2)); y <- c(0, 1, 0, 1)
  result <- .msr_binomial_mle_status(X, y)
  expect_identical(result$status, "finite")
  expect_identical(result$qualification, "exact_original_design_certificate")
  expect_true(.msr_binomial_verify_certificate(X, y, rep(1, 4), result$certificate)$valid)
  changed <- result$certificate; changed$corrected_weights[[1L]] <- "0"
  expect_false(.msr_binomial_verify_certificate(X, y, rep(1, 4), changed)$valid)
  changed <- result$certificate
  row <- setdiff(seq_along(changed$weights), changed$corrected_rows)[[1L]]
  changed$weights[[row]] <- changed$weights[[row]] + .01
  expect_false(.msr_binomial_verify_certificate(X, y, rep(1, 4), changed)$valid)
  expect_false(.msr_binomial_verify_certificate(X, 1 - y, rep(1, 4), result$certificate)$valid)
  changedX <- X; changedX[1, 2] <- changedX[1, 2] + 1e-12
  expect_false(.msr_binomial_verify_certificate(changedX, y, rep(1, 4), result$certificate)$valid)
  expect_false(.msr_binomial_verify_certificate(X, y, rep(1, 4), list())$valid)
})

test_that("exact inverse witnesses handle noninteger designs and zero pivots", {
  skip_if_not_installed("gmp")
  for (X in list(matrix(c(0, .1, .3, .7), 2),
      rbind(c(0, 0, .3), c(.2, 0, 1), c(.7, .1, 1)))) {
    inverse <- .msr_binomial_qinverse(X)
    expect_false(is.null(inverse))
    expect_true(all(gmp::`%*%`(gmp::as.bigq(X), inverse) == diag(nrow(X))))
  }
  expect_null(.msr_binomial_qinverse(matrix(c(.1, .1, .2, .2), 2)))
})

test_that("paired outcome bases establish exact rank without numerical LP admission", {
  skip_if_not_installed("gmp")
  X <- cbind(1, rep(c(-2, -1, 1, 2), each = 2)); y <- rep(c(0, 1), 4)
  failed <- function(...) stop("The paired shortcut must not call the LP provider")
  for (scale in list(c(1, 1), c(2^-400, 2^400), c(2^400, 2^-400))) {
    design <- sweep(X, 2, scale, "*")
    result <- .msr_binomial_mle_status(design, y, .solver = failed)
    expect_identical(result$status, "finite")
    expect_identical(result$certificate$kind, "paired_basis")
    expect_true(.msr_binomial_verify_certificate(design, y, rep(1, 8), result$certificate)$valid)
    expect_false(result$selection_eligible)
    changed <- result$certificate; changed$failure_rows <- changed$success_rows
    expect_false(.msr_binomial_verify_certificate(design, y, rep(1, 8), changed)$valid)
  }
  # No intercept is required by the geometry contract.
  design <- matrix(rep(c(-2, -1, 1, 2), each = 2), ncol = 1)
  expect_identical(.msr_binomial_mle_status(design, y, .solver = failed)$status, "finite")
  expect_identical(.msr_binomial_mle_status(X[8:1, ], y[8:1], .solver = failed)$status, "finite")
  expect_identical(.msr_binomial_mle_status(cbind(X, -X[, 2]), y)$status, "unknown")
})

test_that("non-strict sum maximizers cannot mislabel complete separation", {
  skip_if_not_installed("gmp")
  skip_if_not_installed("lpSolve")
  fixtures <- list(
    list(X = cbind(1, c(-1, 0, 1)), y = c(0, 1, 1), subtype = "complete"),
    list(X = cbind(1, c(-1, 0, 0, 1)), y = c(0, 0, 1, 1), subtype = "quasi_complete"),
    list(X = matrix(c(-2, -1, 1, 2), ncol = 1), y = c(0, 0, 1, 1), subtype = "complete"))
  for (f in fixtures) {
    result <- .msr_binomial_mle_status(f$X, f$y)
    expect_identical(result$status, "nonexistent")
    expect_identical(result$separation, f$subtype)
    expect_true(.msr_binomial_verify_certificate(f$X, f$y, rep(1, length(f$y)), result$certificate)$valid)
    changed <- result$certificate; changed$direction <- rep("0", ncol(f$X))
    expect_false(.msr_binomial_verify_certificate(f$X, f$y, rep(1, length(f$y)), changed)$valid)
    for (scale in list(rep(2^400, ncol(f$X)), seq_len(ncol(f$X)) / 7)) {
      transformed <- sweep(f$X, 2, scale, "*")
      other <- .msr_binomial_mle_status(transformed, f$y)
      expect_identical(other$status, "nonexistent")
      expect_true(.msr_binomial_verify_certificate(transformed, f$y, rep(1, length(f$y)), other$certificate)$valid)
    }
  }
})

test_that("false successful and malformed LP returns cannot certify finite existence", {
  skip_if_not_installed("gmp")
  X <- cbind(1, c(-2, -1, 1, 2)); y <- c(0, 0, 1, 1)
  false_success <- function(objective.in, ...) list(status = 0L,
    solution = rep(1 / length(objective.in), length(objective.in)))
  result <- .msr_binomial_mle_status(X, y, .solver = false_success)
  expect_identical(result$status, "unknown")
  for (bad in list(list(status = "bad"), list(status = NA_integer_), list(status = 0L, solution = NA_real_),
      list(status = 0L, solution = numeric()))) {
    solver <- local({ value <- bad; function(...) value })
    expect_identical(.msr_binomial_mle_status(X, y, .solver = solver)$status, "unknown")
  }
  cat("\nBINOMIAL_CERTIFICATES original_design=1 exact_arithmetic=1 public_admission=0\n")
})

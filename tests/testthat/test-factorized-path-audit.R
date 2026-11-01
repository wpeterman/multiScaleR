test_that("unrestricted annual-SD paths retain native densities and selected zero faces", {
  cases <- 0L; maximum <- 0
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (response_zero in if (response) c(FALSE, TRUE) else FALSE) {
          cases <- cases + 1L
          s <- mixture_fixture(kernel, family, FALSE, response)$spec
          if (family == "binomial") s$trials[] <- 4
          if (response_zero) s <- .msr_marginal_reduce(s, "response")
          p <- private_sq_parameters(s); p$beta[3L] <- .12
          p$scale_z <- p$response_z <- numeric()
          before <- list(s, p, .Random.seed)
          full <- .msr_factorized_path_limits(s, p)
          zero <- p; zero$log_scale_sd <- numeric()
          boundary <- .msr_factorized_path_limits(s, zero)
          expect_true(full$applicable)
          expect_false(full$candidate_zero_sd)
          expect_true(boundary$candidate_zero_sd)
          expect_identical(full$references, boundary$references)
          expected <- numeric()
          for (poisson in if (family == "nbinom2") c(FALSE, TRUE) else FALSE) {
            g <- s; q <- p
            if (poisson) {g$family <- "poisson"; q$log_dispersion <- numeric()}
            broad <- mixture_native_nll(g, q, "broad")
            nearest <- mixture_native_nll(g, q, "nearest")
            expected <- c(expected, vapply(pnorm(c(-.75, 0, .75)), function(w)
              mixture_native_total(broad, nearest, w), 0))
          }
          actual <- unname(vapply(full$references, `[[`, 0, "objective"))
          maximum <- max(maximum, abs(actual - expected))
          expect_equal(actual, unname(expected), tolerance = 1e-10, info = paste(kernel, family, response, response_zero))
          expect_identical(list(s, p, .Random.seed), before)
          fixed <- .msr_marginal_reduce(s, "scale")
          absent <- .msr_factorized_path_limits(fixed, zero)
          expect_false(absent$applicable)
          expect_length(absent$references, 0L)
        }
  expect_equal(cases, 12L)
  cat(sprintf("\nANNUAL_SD_PATH cases=12 native_error=%.12g\n", maximum))
})

test_that("path scope rejects coupling and priors without inventing annual directions", {
  s <- mixture_fixture("gaussian", "gaussian", FALSE, TRUE)$spec
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  bad <- s; bad$priors <- list(beta = 1)
  expect_error(.msr_factorized_path_limits(bad, p), "neutral")
  bad <- s; bad$response$random[[1L]]$index[] <- 1L
  expect_error(.msr_factorized_path_limits(bad, p), "coupling")
  bad <- s; bad$scales[[1L]]$X <- cbind(1, 1)
  expect_false(.msr_factorized_path_limits(bad, p)$applicable)
  bad <- p; bad$scale_z <- 0
  expect_true(all(vapply(.msr_factorized_path_limits(s, bad)$references,
    function(x) !x$completed && grepl("population layout", x$error), TRUE)))
  raw <- mixture_fixture("gaussian", "gaussian", TRUE)$spec
  expect_error(.msr_factorized_path_limits(raw, private_sq_parameters(raw)), "Weighted SQ")
  expect_identical(.msr_factorized_path_limits(s, p)$direction, "original_unrestricted_model")
})

test_that("path witnesses preserve old starts and close the audit on partial failure", {
  s <- mixture_fixture("gaussian", "gaussian", FALSE)$spec
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  native <- .msr_factorized_path_limits(s, p)
  candidate <- min(vapply(native$references, `[[`, 0, "objective")) + .01
  expect_lt(candidate, 100)
  # Only optimizer callbacks are mocked. Real native path leaves remain active.
  # Deliberately singular old references must still retain every attempted start.
  testthat::local_mocked_bindings(
    .msr_factorized_adfun = function(spec, parameters, ...) {
      par <- unlist(parameters, use.names = FALSE)
      env <- new.env(); env$last.par <- par
      env$parList <- function(x, par) parameters
      list(par = par, env = env, fn = function(x) 100,
        gr = function(x) rep(0, length(x)), he = function(x) matrix(0, length(x), length(x)))
    },
    .msr_factorized_objective = function(...) function(parameters) 100,
    .msr_marginal_polish = function(...) list(passed = FALSE), .package = "multiScaleR")
  audit <- .msr_factorized_limit_audit(s, p, candidate)
  expect_true(audit$completed)
  expect_true(any(grepl("^annual-sd-infinite-", audit$witnesses)))
  expect_length(audit$references$`broad-gaussian`$attempts, 3L)
  expect_length(audit$references$`nearest-gaussian`$attempts, 3L)
  expect_true(all(vapply(audit$references$`broad-gaussian`$attempts, function(a) !a$curvature$passed, TRUE)))
  expect_identical(audit$references[names(native$references)], native$references)
  # Partial failure is retained alongside successful paths and all old attempts.
  group_limit <- .msr_factorized_group_mixture_limit
  testthat::local_mocked_bindings(.msr_factorized_group_mixture_limit = function(group, parameters, ratio) {
    if (ratio == 0) stop("mock path failure")
    group_limit(group, parameters, ratio)
  }, .package = "multiScaleR")
  failed <- .msr_factorized_limit_audit(s, p, candidate)
  expect_false(failed$completed)
  expect_match(failed$references$`annual-sd-infinite-gaussian-zero`$error, "mock path failure")
  expect_true(failed$references$`annual-sd-infinite-gaussian-negative`$completed)
  expect_length(failed$references$`broad-gaussian`$attempts, 3L)
  zero <- p; zero$log_scale_sd <- numeric()
  boundary <- .msr_factorized_limit_audit(s, zero, candidate)
  expect_true(boundary$infinite_sd_path$candidate_zero_sd)
  expect_identical(boundary$references$`annual-sd-infinite-gaussian-negative`,
    failed$references$`annual-sd-infinite-gaussian-negative`)
})

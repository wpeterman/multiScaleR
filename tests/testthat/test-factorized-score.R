test_that("conditional score audits agree with full AD derivatives", {
  skip_if_not_installed("RTMB")
  for (family in c("gaussian", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
    s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
      list(habitat = ~1 + (1 | year)), family, "key", 200)
    p <- s$parameters; p$beta <- c(.2, -.1, .4)
    a <- multiScaleR:::.msr_factorized_adfun(s, p, rel.tol = 1e-10, step = .5)
    b <- multiScaleR:::.msr_factorized_score_adfun(s, p, rel.tol = 1e-10, step = .5)
    expect_equal(b$fn(b$par), as.numeric(a$fn(a$par)), tolerance = 1e-6)
    expect_equal(b$gr(b$par), as.numeric(a$gr(a$par)), tolerance = 1e-4)
    H <- a$he(a$par)
    expect_lt(max(abs(b$he(b$par) - H)) / max(1, max(abs(H))), 1e-5)
  }
})

test_that("score auditing is explicit and preserves default and boundary behavior", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 6, n_year = 2)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1), "gaussian", "key", 200)
  p <- s$parameters; p$beta <- c(.5, .2, 2)
  s$y <- multiScaleR:::.msr_bayes_components(s, p)$eta + .15 * sin(seq_along(s$y))
  before <- .Random.seed
  a <- multiScaleR:::.msr_factorized_fit(s, list(starts = 1L))
  b <- multiScaleR:::.msr_factorized_fit(s, list(starts = 1L, cold_audit = "score"))
  expect_identical(.Random.seed, before)
  expect_identical(a$control$cold_audit, "ad")
  expect_identical(b$control$cold_audit, "score")
  expect_true(a$diagnostics$passed && b$diagnostics$passed)
  expect_equal(a$log_likelihood, b$log_likelihood, tolerance = 1e-10)
  expect_equal(a$annual, b$annual, tolerance = 1e-10)
  expect_identical(b$attempts[[1L]]$cold_audit, "score")
  expect_false(b$uncertainty$available)
  for (bad in list("unknown", NA_character_, character(), c("ad", "score"), 1))
    expect_error(multiScaleR:::.msr_factorized_fit(s, list(cold_audit = bad)), "cold_audit")
})

test_that("a failed score audit remains a retained failed attempt", {
  skip_if_not_installed("RTMB")
  f <- bayesian_scale_fixture(n_site = 6, n_year = 2)
  s <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp, f$inputs, f$data,
    list(habitat = ~1), "gaussian", "key", 200)
  p <- s$parameters; p$beta <- c(.5, .2, 2)
  s$y <- multiScaleR:::.msr_bayes_components(s, p)$eta + .15 * sin(seq_along(s$y))
  local_mocked_bindings(.msr_factorized_score_adfun = function(...) stop("injected score failure"),
    .package = "multiScaleR")
  fit <- multiScaleR:::.msr_factorized_fit(s, list(starts = 1L, cold_audit = "score"))
  expect_false(fit$diagnostics$passed)
  expect_length(fit$attempts, 1)
  expect_identical(fit$attempts[[1L]]$error_stage, "cold integration")
  expect_match(fit$attempts[[1L]]$error, "injected score failure")
})

test_that("one-group joint-score tape retains native response and normal density constants", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        s <- private_sq_fixture(kernel, family, response = response)$spec
        if (family == "binomial") s$trials[] <- 4
        g <- .msr_factorized_group(s, 2L); p <- private_sq_parameters(g)
        p$beta[3] <- .12; p$scale_z <- p$response_z <- numeric()
        before <- g; tape <- .msr_factorized_group_score_tape(g, p)
        expect_length(tape$z, 1L)
        expect_identical(names(tape$obj$par)[tape$outer], rep(names(p), lengths(p)))
        for (z in c(-1.5, 0, 1.2)) {
          v <- tape$obj$par; v[tape$z] <- z
          q <- p; q$scale_z <- z
          expected <- sq_group_native_response(g, q) - dnorm(z, log = TRUE)
          expect_equal(as.numeric(tape$obj$fn(v)), expected, tolerance = 1e-9)
          expect_true(all(is.finite(tape$obj$gr(v))))
          expect_true(all(is.finite(tape$obj$he(v))))
        }
        expect_identical(g, before)
        expect_error(.msr_factorized_score_adfun(s), "Weighted SQ")
      }
})

test_that("conditional signed-score moments match integrated SQ nuisance derivatives", {
  skip_if_not_installed("RTMB")
  # Twenty fixed-parameter identities, not a model-fitting or sampler panel.
  cases <- expand.grid(kernel = c("gaussian", "exp"),
    kind = c("gaussian", "poisson", "binomial", "nbinom2", "response"), sd = c(.3, 2))
  for (i in seq_len(nrow(cases))) {
    kind <- as.character(cases$kind[i]); family <- if (kind == "response") "gaussian" else kind
    s <- private_sq_fixture(as.character(cases$kernel[i]), family, response = kind == "response")$spec
    if (family == "binomial") s$trials[] <- 4
    g <- .msr_factorized_group(s, 2L); p <- private_sq_parameters(g)
    p$beta[3] <- .12; p$log_scale_sd[] <- log(cases$sd[i])
    p$scale_z <- p$response_z <- numeric()
    before <- g; tape <- .msr_factorized_group_score_tape(g, p)
    normalizer <- sq_group_native_integral(g, p)
    bound <- if (family == "gaussian") -length(g$y) * (p$log_dispersion + .5 * log(2 * pi)) else 0
    cache <- new.env(parent = emptyenv())
    node <- function(z, order = 0L) {
      key <- sprintf("%a", z)
      record <- if (exists(key, cache, inherits = FALSE)) get(key, cache, inherits = FALSE) else list()
      v <- tape$obj$par; v[tape$z] <- z
      if (is.null(record$value)) record$value <- as.numeric(tape$obj$fn(v))
      if (order >= 1L && is.null(record$score)) record$score <- as.numeric(tape$obj$gr(v))[tape$outer]
      if (order >= 2L && is.null(record$hessian)) record$hessian <- tape$obj$he(v)[tape$outer, tape$outer, drop = FALSE]
      assign(key, record, cache); record
    }
    weight <- function(z) {
      # Only a floating-point zero of the valid envelope skips a tail node.
      if (exp(bound + dnorm(z, log = TRUE) + normalizer) == 0) return(0)
      exp(-node(z)$value + normalizer)
    }
    integrate_moment <- function(fn, absolute) integrate(function(z) vapply(z, fn, 0),
      -Inf, Inf, rel.tol = 1e-9, abs.tol = absolute,
      subdivisions = 500L, stop.on.error = TRUE)$value
    expect_equal(integrate_moment(weight, 0), 1, tolerance = 1e-8)
    center <- node(0, 1L)$score
    difference <- vapply(seq_along(center), function(j) integrate_moment(function(z) {
      w <- weight(z); if (w == 0) return(0)
      w * (node(z, 1L)$score[j] - center[j])
    }, 1e-10), 0)
    score <- center + difference
    H <- matrix(0, length(center), length(center))
    for (j in seq_along(center)) for (k in seq_len(j)) {
      H[j, k] <- H[k, j] <- integrate_moment(function(z) {
        w <- weight(z); if (w == 0) return(0)
        a <- node(z, 2L); delta <- a$score - center
        w * (a$hessian[j, k] - delta[j] * delta[k])
      }, 1e-8) + difference[j] * difference[k]
    }
    fn <- .msr_factorized_group_integral(g, TRUE)
    obj <- RTMB::MakeADFun(function(q) {
      part <- fn(q); -log(part$value) - part$shift
    }, p, silent = TRUE)
    expect_identical(names(obj$par), names(tape$obj$par)[tape$outer])
    expect_equal(as.numeric(obj$fn(obj$par)), normalizer, tolerance = 1e-6)
    expect_equal(score, as.numeric(obj$gr(obj$par)), tolerance = 2e-5)
    expect_equal(unname(H), unname(obj$he(obj$par)), tolerance = 2e-5)
    native <- function(x) sq_group_native_integral(g, private_sq_unpack(x, p))
    expect_equal(score, private_sq_fd_gradient(native, obj$par), tolerance = 2e-5)
    expect_equal(unname(H), private_sq_fd_hessian(native, obj$par), tolerance = 2e-5)
    expect_identical(g, before)
    expect_error(.msr_factorized_score_adfun(s), "Weighted SQ")
  }
})

test_that("joint-score group tapes enforce layout and scientific RAM before native tape entry", {
  skip_if_not_installed("RTMB")
  s <- private_sq_fixture(response = TRUE)$spec; g <- .msr_factorized_group(s, 1L)
  p <- private_sq_parameters(g); p$scale_z <- p$response_z <- numeric()
  expect_error(.msr_factorized_group_score_tape(s, p), "one recoded")
  bad <- g; bad$observation_ids <- s$observation_ids
  expect_error(.msr_factorized_group_score_tape(bad, p), "observation metadata")
  unintegrated <- p; unintegrated$response_z <- 0
  expect_error(.msr_factorized_group_score_tape(g, unintegrated), "integrated response")
  launched <- 0L
  for (ram in c(NA_real_, 1024^3)) with_mocked_bindings({
    with_mocked_bindings({
      expect_error(.msr_factorized_group_score_tape(g, p),
        if (is.na(ram)) "Cannot detect RAM" else "RAM budget")
    }, MakeADFun = function(...) { launched <<- launched + 1L; stop("Tape reached") }, .package = "RTMB")
  }, .msr_total_ram_bytes = function() ram, .package = "multiScaleR")
  expect_equal(launched, 0L)
  expect_error(.msr_factorized_score_adfun(s), "Weighted SQ")
})

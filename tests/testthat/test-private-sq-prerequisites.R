test_that("exact SQ limits include all nearest ties and keep constant support exact", {
  buffers <- list(shifted = list(v = 1e12 + c(0, 2, 4), d = c(1, 1, 2)),
    constant = list(v = rep(7, 4), d = c(0, 0, 1, 2)))
  expect_equal(.msr_scale_surface_limit(buffers, "nearest"), c(shifted = 1, constant = 0))
  expect_equal(.msr_scale_surface_limit(buffers, "broad"),
    c(shifted = sqrt(8 / 3), constant = 0), tolerance = 1e-12)
  expect_error(.msr_scale_surface_limit(list(), "broad"), "nonempty")
  expect_error(.msr_scale_surface_limit(list(list(v = c(1, NA), d = c(0, 1))), "nearest"), "finite")
  expect_error(.msr_scale_surface_limit(list(list(v = 1:2, d = -1:0)), "nearest"), "nonnegative")
  expect_error(.msr_scale_surface_limit(list(list(v = 1:2, d = 0)), "broad"), "matching")
  expect_error(.msr_scale_surface_limit(buffers, "invalid"), "arg")
  for (kernel in c("gaussian", "exp")) {
    f <- private_sq_fixture(kernel)
    for (limit in c("broad", "nearest")) {
      b <- f$spec$surface$buffers$rough
      direct <- vapply(b, function(x) {
        v <- if (limit == "nearest") x$v[x$d == min(x$d)] else x$v
        sqrt(mean((v - mean(v))^2))
      }, 0)
      exact <- .msr_scale_surface_limit(b, limit)
      expect_equal(exact, direct, tolerance = 1e-12)
      extreme <- if (limit == "broad") 1e20 else 1e-20
      expect_equal(.msr_scale_surface_raw(b, rep(log(extreme), length(b)), kernel),
        exact, tolerance = 1e-10)
    }
  }
})

test_that("pure SQ fixed-kernel reductions preserve neutral normalized density and reference maps", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- private_sq_fixture(kernel, family)
    s <- f$spec
    if (family == "binomial") s$trials[] <- 4
    before <- s
    for (limit in c("broad", "nearest")) for (poisson in
        if (family == "nbinom2") c(FALSE, TRUE) else FALSE) {
      # This internal reduction constructs a fixed density witness, not a fit.
      reduced <- .msr_factorized_kernel_limit(s, limit, poisson)
      raw <- .msr_scale_surface_limit(s$surface$buffers$rough, limit)
      expect_identical(names(reduced$V), "rough")
      expect_identical(names(reduced$C), "rough")
      expect_identical(names(reduced$delta), "rough")
      expect_equal(as.numeric(reduced$V$rough), raw)
      expect_equal(as.numeric(reduced$C$rough), rep(1, length(s$y)))
      expect_length(reduced$surface$buffers, 0)
      expect_identical(reduced$center, s$center)
      expect_identical(reduced$spread, s$spread)
      expect_identical(reduced$rows, s$rows)
      expect_identical(reduced$observation_ids, s$observation_ids)
      expect_identical(reduced$offset, s$offset)
      expect_identical(reduced$trials, s$trials)
      expect_null(reduced$priors)
      expect_false(reduced$fixed_kernel_limit$finite_scale_estimate)
      expect_length(reduced$parameters$scale_beta, 0)
      expect_length(reduced$parameters$scale_z, 0)
      expect_length(reduced$parameters$log_scale_sd, 0)
      p <- reduced$parameters; p$beta <- c(.4, -.2, .7)
      if (length(p$log_dispersion)) p$log_dispersion[] <- if (family == "gaussian") log(.8) else log(3)
      eta <- drop(s$X_fixed %*% p$beta[seq_len(ncol(s$X_fixed))]) + s$offset +
        p$beta[ncol(s$X_fixed) + 1L] * (raw - s$center[[1]]) / s$spread[[1]]
      expected <- -sum(switch(reduced$family,
        gaussian = dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE),
        poisson = dpois(s$y, exp(eta), log = TRUE),
        binomial = dbinom(s$y, s$trials, plogis(eta), log = TRUE),
        nbinom2 = dnbinom(s$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE)))
      expect_equal(.msr_factorized_objective(reduced, FALSE)(p), expected, tolerance = 1e-10)
      obj <- .msr_factorized_adfun(reduced, p)
      expect_equal(as.numeric(obj$fn(obj$par)), expected, tolerance = 1e-9)
      expect_true(all(is.finite(obj$gr(obj$par))))
      expect_true(all(is.finite(obj$he(obj$par))))
      # Independent response-density derivatives in the surviving coordinates.
      X <- cbind(s$X_fixed, (raw - s$center[[1]]) / s$spread[[1]])
      mu <- if (reduced$family == "binomial") plogis(eta) else exp(eta)
      score_eta <- switch(reduced$family,
        gaussian = (eta - s$y) / exp(2 * p$log_dispersion),
        poisson = mu - s$y, binomial = s$trials * mu - s$y,
        nbinom2 = exp(p$log_dispersion) * (mu - s$y) / (exp(p$log_dispersion) + mu))
      weight <- switch(reduced$family,
        gaussian = rep(exp(-2 * p$log_dispersion), length(eta)),
        poisson = mu, binomial = s$trials * mu * (1 - mu),
        nbinom2 = exp(p$log_dispersion) * mu * (exp(p$log_dispersion) + s$y) /
          (exp(p$log_dispersion) + mu)^2)
      score <- drop(crossprod(X, score_eta))
      H <- crossprod(X, X * weight)
      if (reduced$family == "gaussian") {
        ss <- sum((s$y - eta)^2) * exp(-2 * p$log_dispersion)
        H <- rbind(cbind(H, -2 * score), c(-2 * score, 2 * ss))
        score <- c(score, length(eta) - ss)
      } else if (reduced$family == "nbinom2") {
        size <- exp(p$log_dispersion)
        first <- digamma(s$y + size) - digamma(size) + log(size) + 1 -
          log(size + mu) - (s$y + size) / (size + mu)
        second <- trigamma(s$y + size) - trigamma(size) + 1 / size -
          1 / (size + mu) - (mu - s$y) / (size + mu)^2
        mixed <- drop(crossprod(X, size * mu * (mu - s$y) / (size + mu)^2))
        H <- rbind(cbind(H, mixed), c(mixed, -size * sum(first) - size^2 * sum(second)))
        score <- c(score, -size * sum(first))
      }
      expect_equal(as.numeric(obj$gr(obj$par)), unname(score), tolerance = 1e-8)
      expect_equal(unname(obj$he(obj$par)), unname(H), tolerance = 1e-8)
    }
    expect_identical(s, before)
    expect_error(.msr_factorized_limit_spec(s, "broad"), "Weighted SQ")
    expect_error(.msr_scale_check_capability(s, "marginal_likelihood", "laplace"), "Weighted SQ")
  }
})

test_that("kernel independence removes annual SD only and preserves Gaussian response covariance", {
  skip_if_not_installed("RTMB")
  for (annual in c(FALSE, TRUE)) {
    f <- private_sq_fixture(response = TRUE, annual = annual)
    s <- f$spec
    for (limit in c("broad", "nearest")) {
      reduced <- .msr_factorized_kernel_limit(s, limit)
      expect_identical(reduced$parameters$log_response_sd, s$parameters$log_response_sd)
      expect_identical(reduced$response, s$response)
      p <- reduced$parameters; p$beta <- c(.4, -.2, .7)
      p$log_dispersion[] <- log(.8); p$log_response_sd[] <- log(.3)
      eta <- .msr_bayes_components(reduced, p, include_response = FALSE)$eta
      group <- reduced$response$random[[1]]$index
      covariance <- diag(exp(2 * p$log_dispersion), length(eta)) +
        outer(group, group, "==") * exp(2 * p$log_response_sd)
      ch <- chol(covariance); residual <- s$y - eta
      expected <- .5 * (length(eta) * log(2 * pi) + 2 * sum(log(diag(ch))) +
        sum(backsolve(ch, residual, transpose = TRUE)^2))
      expect_equal(.msr_factorized_objective(reduced)(p), expected, tolerance = 1e-10)
      expect_length(reduced$parameters$log_scale_sd, 0)
    }
  }
})

test_that("SQ group subsetting recomputes memory from the retained rows", {
  f <- private_sq_fixture(); s <- f$spec; before <- s
  for (i in seq_along(s$scales[[1]]$random[[1]]$levels)) {
    rows <- which(s$scales[[1]]$random[[1]]$index == i)
    group <- .msr_factorized_group(s, i)
    expected <- sum(vapply(s$surface$buffers$rough[rows], function(b) length(b$v), 0L))
    expect_equal(group$surface$memory$observation_cells, expected)
    expect_equal(group$surface$memory$estimated_AD_bytes, expected * 512)
    expect_identical(group$surface$buffers$rough, s$surface$buffers$rough[rows])
    expect_equal(length(group$y), length(rows))
  }
  expect_identical(s, before)
})

test_that("every SQ differentiation entry rejects unresolved or insufficient RAM before tapes", {
  skip_if_not_installed("RTMB")
  s <- private_sq_fixture()$spec
  launched <- 0L
  constructors <- list(
    function() .msr_marginal_objective(s, TRUE),
    function() .msr_marginal_adfun(s),
    function() .msr_marginal_inner_audit(s, s$parameters),
    function() .msr_factorized_objective(s, TRUE),
    function() .msr_factorized_scale_objective(s, TRUE),
    function() .msr_factorized_adfun(s),
    function() .msr_factorized_score_adfun(s),
    function() .msr_factorized_conditional_cache(s))
  for (ram in c(NA_real_, 1024^3)) with_mocked_bindings({
    with_mocked_bindings({
      for (constructor in constructors) expect_error(constructor(),
        if (is.na(ram)) "Cannot detect RAM" else "RAM budget")
    }, MakeADFun = function(...) { launched <<- launched + 1L; stop("Tape reached") },
       MakeTape = function(...) { launched <<- launched + 1L; stop("Tape reached") },
       .package = "RTMB")
    expect_true(is.function(.msr_marginal_objective(s, FALSE)))
  }, .msr_total_ram_bytes = function() ram,
     # Isolate admission only; no raw-SQ fit or tape is constructed by this test.
     .msr_factorized_validate = function(spec) invisible(spec), .package = "multiScaleR")
  expect_equal(launched, 0)
  expect_error(.msr_factorized_validate(s), "Weighted SQ")
})

test_that("a singular nearest-SQ design is not promoted to finite curvature", {
  skip_if_not_installed("RTMB")
  s <- private_sq_fixture()$spec
  s$surface$buffers$rough <- lapply(s$surface$buffers$rough, function(b) { b$d <- seq_along(b$d); b })
  reduced <- .msr_factorized_kernel_limit(s, "nearest")
  raw <- as.numeric(reduced$V$rough)
  design <- cbind(reduced$X_fixed, (raw - s$center[[1]]) / s$spread[[1]])
  expect_true(qr(design)$rank < ncol(design))
  obj <- .msr_factorized_adfun(reduced)
  expect_false(.msr_marginal_curvature(obj$he(obj$par))$passed)
})

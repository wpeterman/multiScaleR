test_that("neutral raw SQ joint density and every derivative block match an independent target", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
    for (annual in c(FALSE, TRUE)) {
      f <- private_sq_fixture(kernel, family, annual = annual)
      s <- f$spec; if (family == "binomial") s$trials[] <- 4
      before <- s; p <- private_sq_parameters(s)
      expect_null(s$priors)
      ref <- private_sq_reference(s, p)
      got <- .msr_bayes_components(s, p)
      expect_equal(got$eta, ref$eta, tolerance = 1e-12)
      expect_equal(as.numeric(got$kernel$rough), ref$kernel, tolerance = 1e-12)
      expect_equal(got$scales$rough, ref$physical, tolerance = 1e-12)
      # Evaluate admission before entering RTMB; no optimizer or latent integration.
      objective <- .msr_marginal_objective(s, TRUE)
      obj <- RTMB::MakeADFun(objective, p, silent = TRUE)
      native <- function(x) private_sq_joint_reference(s, private_sq_unpack(x, p))
      expect_identical(names(obj$par), rep(names(p), lengths(p)))
      expect_equal(as.numeric(obj$fn(obj$par)), native(obj$par), tolerance = 1e-9)
      expect_equal(.msr_marginal_objective(s)(p), native(obj$par), tolerance = 1e-10)
      expect_equal(as.numeric(obj$gr(obj$par)), private_sq_fd_gradient(native, obj$par), tolerance = 2e-5)
      expect_equal(unname(obj$he(obj$par)), private_sq_fd_hessian(native, obj$par), tolerance = 2e-5)
      expect_identical(s, before)
      expect_error(.msr_scale_check_capability(s, "marginal_likelihood", "laplace"), "Weighted SQ")
      expect_error(.msr_scale_check_capability(s, "marginal_likelihood", "factorized"), "Weighted SQ")
    }
})

test_that("raw SQ likelihood and AD coordinates are invariant to physical distance re-expression", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- private_sq_fixture(kernel, family)
    s <- f$spec; if (family == "binomial") s$trials[] <- 4
    p <- private_sq_parameters(s)
    moved <- f$inputs
    moved$d_list <- lapply(moved$d_list, function(d) d * moved$unit_conv)
    moved$unit_conv <- 1
    other <- .msr_scale_prepare(f$formula, moved, f$data, f$scales, family, "key", reference_scale = 250)
    if (family == "binomial") other$trials[] <- 4
    q <- p; q$scale_beta <- p$scale_beta + log(s$unit_conv)
    expect_equal(other$center, s$center, tolerance = 1e-12)
    expect_equal(other$spread, s$spread, tolerance = 1e-12)
    expect_identical(other$observation_ids, s$observation_ids)
    expect_equal(.msr_bayes_components(other, q)$scales, .msr_bayes_components(s, p)$scales, tolerance = 1e-12)
    expect_equal(.msr_marginal_objective(other)(q), .msr_marginal_objective(s)(p), tolerance = 1e-10)
    a <- .msr_marginal_objective(s, TRUE); b <- .msr_marginal_objective(other, TRUE)
    oa <- RTMB::MakeADFun(a, p, silent = TRUE); ob <- RTMB::MakeADFun(b, q, silent = TRUE)
    expect_equal(as.numeric(ob$gr(ob$par)), as.numeric(oa$gr(oa$par)), tolerance = 1e-9)
    expect_equal(unname(ob$he(ob$par)), unname(oa$he(oa$par)), tolerance = 1e-8)
  }
})

test_that("private group metadata and Gaussian response density follow reordered observations", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) {
    f <- private_sq_fixture(kernel, response = TRUE)
    order <- c(10, 3, 7, 1, 12, 5, 8, 2, 11, 4, 9, 6)
    s <- .msr_scale_prepare(f$formula, f$inputs, f$data[order, , drop = FALSE],
      f$scales, "gaussian", "key", reference_scale = 250)
    before <- s; p <- private_sq_parameters(s); values <- numeric(3)
    for (i in seq_len(3)) {
      rows <- which(s$scales[[1L]]$random[[1L]]$index == i)
      g <- .msr_factorized_group(s, i)
      expect_identical(g$rows, s$rows[rows])
      expect_identical(g$observation_ids, s$observation_ids[rows])
      expect_identical(g$data, s$data[rows, , drop = FALSE])
      expect_identical(g$y, s$y[rows])
      expect_identical(g$offset, s$offset[rows])
      expect_identical(g$trials, s$trials[rows])
      expect_identical(g$surface$buffers$rough, s$surface$buffers$rough[rows])
      expect_equal(nrow(g$data), length(g$y))
      expect_equal(g$surface$memory$observation_cells, sum(vapply(g$surface$buffers$rough, function(b) length(b$v), 0)))
      q <- p; q$scale_z <- p$scale_z[i]
      observed <- sort(unique(s$response$random[[1L]]$index[rows]))
      q$response_z <- p$response_z[observed]
      expect_equal(.msr_marginal_objective(g)(q), private_sq_joint_reference(g, q), tolerance = 1e-10)
      eta <- private_sq_reference(g, q, response = FALSE)$eta
      index <- g$response$random[[1L]]$index
      covariance <- diag(exp(2 * q$log_dispersion), length(eta)) +
        outer(index, index, "==") * exp(2 * q$log_response_sd)
      ch <- chol(covariance)
      expected <- .5 * (length(eta) * log(2 * pi) + 2 * sum(log(diag(ch))) +
        sum(backsolve(ch, g$y - eta, transpose = TRUE)^2))
      values[i] <- .msr_factorized_response_nll(g, q)
      expect_equal(values[i], expected, tolerance = 1e-10)
      q$response_z <- numeric()
      .msr_scale_surface_ram(g)
      obj <- RTMB::MakeADFun(function(z) sum(.msr_factorized_response_nll(g, z, TRUE)), q, silent = TRUE)
      expect_equal(as.numeric(obj$fn(obj$par)), expected, tolerance = 1e-9)
      expect_true(all(is.finite(obj$gr(obj$par))))
      expect_true(all(is.finite(obj$he(obj$par))))
    }
    expect_equal(sum(values), sum(.msr_factorized_response_nll(s, p)), tolerance = 1e-10)
    expect_identical(s, before)
    expect_error(.msr_factorized_validate(s), "Weighted SQ")
  }
})

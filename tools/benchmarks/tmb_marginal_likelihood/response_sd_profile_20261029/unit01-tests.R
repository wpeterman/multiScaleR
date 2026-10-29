response_sd_profile_fixture <- function() {
  f <- bayesian_scale_fixture(n_site = 6L, n_year = 5L, seed = 713403)
  prepare <- function(data) multiScaleR:::.msr_marginal_prepare(
    count ~ habitat + temp + (1 | year) + offset(log(effort)), f$inputs, data,
    list(habitat = ~1), "gaussian", "key", 200)
  spec <- prepare(f$data)
  p <- spec$parameters; p$beta <- c(.5, .25, 2); p$scale_beta[] <- log(.2)
  p$response_z[] <- c(-1, -.3, 0, .4, .9); p$log_response_sd[] <- log(.35)
  p$log_dispersion[] <- log(.2)
  eta <- multiScaleR:::.msr_bayes_components(spec, p)$eta
  n <- nrow(f$data)
  residual <- qnorm((((seq_len(n) * 7L - 1L) %% n) + .5) / n)
  f$data$count <- eta + .2 * residual
  list(spec = prepare(f$data), parameters = p, inputs = f$inputs, data = f$data)
}

response_profile_covariance_nll <- function(spec, p, sd) {
  eta <- multiScaleR:::.msr_bayes_components(spec, p, include_response = FALSE)$eta
  Z <- model.matrix(~0 + factor(spec$response$random[[1L]]$index))
  covariance <- diag(exp(2 * p$log_dispersion), length(spec$y)) + sd^2 * tcrossprod(Z)
  residual <- spec$y - eta
  length(residual) / 2 * log(2 * pi) + sum(log(diag(chol(covariance)))) +
    drop(crossprod(residual, solve(covariance, residual))) / 2
}

test_that("response SD traces retain exact covariance targets, starts and fixed-SD witnesses", {
  skip_if_not_installed("RTMB")
  f <- response_sd_profile_fixture()
  before <- serialize(f$spec, NULL); rng <- .Random.seed
  trace <- multiScaleR:::.msr_factorized_response_sd_profile(f$spec, c(0, .2, .5),
    initial_parameters = list(f$parameters))
  expect_s3_class(trace, "msr_response_sd_profile")
  expect_identical(serialize(f$spec, NULL), before)
  expect_identical(.Random.seed, rng)
  expect_identical(trace$curve$sd, c(0, .2, .5))
  expect_length(trace$points, 3L)
  expect_true(any(trace$curve$available & trace$curve$sd > 0))
  expect_false(trace$global_fit_eligible)
  expect_false(trace$uncertainty$available)
  expect_false(trace$diagnostics$interval_available)
  expect_false(any(c("confidence_level", "cutoff", "interval", "p_value") %in% names(trace)))
  expect_identical(trace$response_group, "year")
  for (point in trace$points) {
    expect_length(point$attempts, 4L)
    expect_false(point$global_fit_eligible)
    for (a in point$attempts) {
      expect_false(a$global_fit_eligible)
      expect_equal(a$initial$log_response_sd, if (point$sd == 0) numeric() else log(point$sd))
      if (isTRUE(a$audited_value)) {
        expect_equal(a$objective, response_profile_covariance_nll(f$spec, a$parameters, point$sd), tolerance = 1e-8)
        expect_lt(a$cold_error, 1e-6)
        expect_lt(a$reference_error, 1e-6)
        expect_true(all(is.finite(a$full_gradient)))
        expect_true(all(is.finite(a$full_hessian)))
        expect_equal(a$parameters$log_response_sd, if (point$sd == 0) numeric() else log(point$sd))
      }
    }
    if (!is.null(point$limits)) {
      expect_true(point$limits$completed)
      expect_equal(point$limits$fixed_sd, point$sd)
      for (limit in point$limits$references) {
        expect_length(limit$attempts, 3L)
        for (a in limit$attempts) {
          expect_false(a$global_fit_eligible)
          expect_equal(a$initial$log_response_sd, if (point$sd == 0) numeric() else log(point$sd))
          if (!is.null(a$parameters)) expect_equal(a$parameters$log_response_sd,
            if (point$sd == 0) numeric() else log(point$sd))
        }
      }
    }
  }
  # Check all nuisance directions against an independent normalized covariance
  # density, and its analytic beta block, at a retained positive-SD optimum.
  point <- trace$points[[which(trace$curve$available & trace$curve$sd > 0)[1L]]]
  a <- point$attempts[[point$selected]]
  template <- a$parameters
  unpack <- function(x) {
    p <- template; at <- 0L
    for (nm in names(p)) {
      n <- length(p[[nm]])
      if (nm != "log_response_sd") {
        p[[nm]] <- if (n) unname(x[at + seq_len(n)]) else numeric(); at <- at + n
      }
    }
    p
  }
  h <- 1e-5
  reference_score <- vapply(seq_along(a$par), function(j) {
    hi <- lo <- a$par; hi[j] <- hi[j] + h; lo[j] <- lo[j] - h
    (response_profile_covariance_nll(f$spec, unpack(hi), point$sd) -
      response_profile_covariance_nll(f$spec, unpack(lo), point$sd)) / (2 * h)
  }, 0)
  expect_lt(max(abs(reference_score - a$full_gradient[a$free_indices])), 3e-5)
  kernel <- multiScaleR:::.msr_bayes_components(f$spec, a$parameters, include_response = FALSE)$kernel
  X <- cbind(f$spec$X_fixed, do.call(cbind, kernel))
  Z <- model.matrix(~0 + factor(f$spec$response$random[[1L]]$index))
  covariance <- diag(exp(2 * a$parameters$log_dispersion), length(f$spec$y)) + point$sd^2 * tcrossprod(Z)
  beta <- which(names(a$par) == "beta")
  expect_equal(unname(a$nuisance_hessian[beta, beta]), unname(crossprod(X, solve(covariance, X))), tolerance = 1e-8)
  expect_gt(abs(a$full_gradient[a$fixed_indices]), 1e-4)
  # Exact zero is block removal, not a small positive optimization floor.
  zero <- multiScaleR:::.msr_response_profile_problem(f$spec, 0, f$parameters)
  p0 <- zero$parameters(zero$par)
  expect_length(p0$log_response_sd, 0L)
  expect_length(p0$response_z, 0L)
  expect_length(zero$fixed, 0L)
  expect_equal(zero$fn(zero$par), response_profile_covariance_nll(f$spec, p0, 0), tolerance = 1e-9)
  expect_equal(response_profile_covariance_nll(f$spec, p0, 1e-8),
    response_profile_covariance_nll(f$spec, p0, 0), tolerance = 1e-8)
})

test_that("profile scope and identification failures stop before fitting", {
  f <- response_sd_profile_fixture()
  profile <- multiScaleR:::.msr_factorized_response_sd_profile
  for (grid in list(numeric(), .2, c(0, 0), c(0, NA), c(0, Inf), c(0, -.1), matrix(0), "0"))
    expect_error(profile(f$spec, grid), "sd")
  expect_error(profile(f$spec, c(0, .2), list(cutoff = 3.84)), "control")
  expect_error(profile(f$spec, c(0, .2), list(starts = 0)), "starts")
  expect_error(profile(f$spec, c(0, .2), initial_parameters = list(list(beta = 1))), "blocks")
  posterior <- f$spec; posterior$priors <- list(scale_median = 200)
  expect_error(profile(posterior, c(0, .2)), "unpenalized")
  counts <- f$spec; counts$family <- "nbinom2"
  expect_error(profile(counts, c(0, .2)), "Gaussian")
  random_scale <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + (1 | year),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  expect_error(profile(random_scale, c(0, .2)), "shared")
  constant <- f$spec; constant$y <- constant$offset + as.numeric(factor(constant$response$random[[1]]$index))
  result <- profile(constant, c(0, .2))
  expect_false(result$diagnostics$trace_complete)
  expect_false(any(result$curve$available))
  expect_true(all(vapply(result$points, function(p) length(p$attempts) == 0L, TRUE)))
  expect_match(result$curve$reason[1L], "unbounded")
})

test_that("conditional selection preserves lower ineligible and competing starts", {
  f <- response_sd_profile_fixture()
  calls <- 0L
  local_mocked_bindings(.msr_response_profile_attempt = function(spec, sd, initial, iter.max) {
    calls <<- calls + 1L; i <- (calls - 1L) %% 3L + 1L
    list(eligible_nuisance = i != 2L, stationary = TRUE, audited_value = TRUE,
      objective = c(0, -1, .5)[i], par = c(0, 0), nuisance_hessian = diag(2), parameters = initial,
      fixed_sd = sd, initial = initial, global_fit_eligible = FALSE)
  }, .msr_response_profile_limits = function(...) list(completed = TRUE, witnesses = character()), .package = "multiScaleR")
  lower <- multiScaleR:::.msr_factorized_response_sd_profile(f$spec, c(0, .2))
  expect_false(any(lower$curve$available))
  expect_true(all(vapply(lower$points, function(p) p$better_ineligible_candidate, TRUE)))
  expect_true(all(vapply(lower$points, function(p) length(p$attempts) == 3L, TRUE)))
  calls <- 0L
  local_mocked_bindings(.msr_response_profile_attempt = function(spec, sd, initial, iter.max) {
    calls <<- calls + 1L; i <- (calls - 1L) %% 3L + 1L
    list(eligible_nuisance = TRUE, stationary = TRUE, audited_value = TRUE,
      objective = c(0, 1, .5)[i], par = c(if (i == 2L) 2 else 0, 0),
      nuisance_hessian = diag(2), parameters = initial, fixed_sd = sd, global_fit_eligible = FALSE)
  }, .package = "multiScaleR")
  competing <- multiScaleR:::.msr_factorized_response_sd_profile(f$spec, c(0, .2))
  expect_false(any(competing$curve$available))
  expect_true(all(vapply(competing$points, function(p) p$competing_nuisance_optima, TRUE)))
})

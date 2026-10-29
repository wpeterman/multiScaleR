annual_sd_profile_fixture <- function() {
  f <- bayesian_scale_fixture(n_site = 6L, n_year = 5L, seed = 713403)
  prepare <- function(data) multiScaleR:::.msr_marginal_prepare(
    count ~ habitat + temp + offset(log(effort)), f$inputs, data,
    list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  spec <- prepare(f$data)
  p <- spec$parameters; p$beta <- c(.5, .25, 2); p$scale_beta[] <- log(.2)
  p$scale_z[] <- c(-1, -.3, 0, .4, .9); p$log_scale_sd[] <- log(.35)
  p$log_dispersion[] <- log(.2)
  eta <- multiScaleR:::.msr_bayes_components(spec, p)$eta
  n <- nrow(f$data)
  residual <- qnorm((((seq_len(n) * 7L - 1L) %% n) + .5) / n)
  f$data$count <- eta + .2 * residual
  list(spec = prepare(f$data), parameters = p, inputs = f$inputs, data = f$data)
}

# Independent base-R integration of Gaussian observation densities and the
# standard-normal group density. Centering changes neither limits nor measure.
annual_profile_reference <- function(spec, parameters, sd) {
  index <- spec$scales[[1L]]$random[[1L]]$index
  if (sd == 0) {
    p <- parameters; p$scale_z <- rep(0, length(unique(index))); p$log_scale_sd <- 0
    eta <- multiScaleR:::.msr_bayes_components(spec, p)$eta
    return(-sum(dnorm(spec$y, eta, exp(p$log_dispersion), log = TRUE)))
  }
  sum(vapply(seq_along(unique(index)), function(group) {
    g <- multiScaleR:::.msr_factorized_group(spec, group)
    log_density <- function(z) {
      p <- parameters; p$log_scale_sd <- log(sd); p$scale_z <- z
      eta <- multiScaleR:::.msr_bayes_components(g, p)$eta
      sum(dnorm(g$y, eta, exp(p$log_dispersion), log = TRUE)) + dnorm(z, log = TRUE)
    }
    shift <- -optimize(function(z) -log_density(z), c(-12, 12))$objective
    bound <- -length(g$y) * (parameters$log_dispersion + .5 * log(2 * pi))
    # Skip only values whose Gaussian density upper bound already rounds to
    # zero. This avoids 0*Inf kernel arithmetic without imposing a tail cutoff.
    value <- integrate(function(z) vapply(z, function(x) {
      if (exp(bound + dnorm(x, log = TRUE) - shift) == 0) return(0)
      exp(log_density(x) - shift)
    }, 0),
      -Inf, Inf, rel.tol = 1e-9, subdivisions = 1000L)$value
    -shift - log(value)
  }, 0))
}

test_that("annual SD traces preserve normalized fixed targets and exact limits", {
  skip_if_not_installed("RTMB")
  f <- annual_sd_profile_fixture()
  before <- serialize(f$spec, NULL); rng <- .Random.seed
  trace <- multiScaleR:::.msr_factorized_scale_sd_profile(f$spec, c(0, .2, .5),
    initial_parameters = list(f$parameters))
  observer <- getOption("multiScaleR.annual_profile_test_observer")
  if (is.function(observer)) observer(trace)
  expect_s3_class(trace, "msr_scale_sd_profile")
  expect_identical(serialize(f$spec, NULL), before)
  expect_identical(.Random.seed, rng)
  expect_identical(trace$curve$sd, c(0, .2, .5))
  expect_length(trace$points, 3L)
  expect_false(trace$global_fit_eligible)
  expect_false(trace$uncertainty$available)
  expect_false(trace$diagnostics$interval_available)
  expect_identical(trace$scale_group, "year")
  expect_identical(trace$units, "log kernel scale deviations")
  audited <- list()
  for (point in trace$points) {
    expect_length(point$attempts, 4L)
    expect_false(point$global_fit_eligible)
    for (a in point$attempts) {
      expect_false(a$global_fit_eligible)
      expect_equal(a$initial$log_scale_sd, if (point$sd == 0) numeric() else log(point$sd))
      if (isTRUE(a$audited_value)) {
        expect_equal(a$parameters$log_scale_sd, if (point$sd == 0) numeric() else log(point$sd))
        expect_lt(a$cold_error, 1e-6)
        expect_lt(a$reference_error, 1e-6)
        expect_true(all(is.finite(a$full_gradient)))
        expect_true(all(is.finite(a$full_hessian)))
        audited[[length(audited) + 1L]] <- list(sd = point$sd, attempt = a)
      }
    }
    if (!is.null(point$limits)) {
      expect_true(point$limits$completed)
      expect_equal(point$limits$fixed_sd, point$sd)
      for (limit in point$limits$references) {
        expect_length(limit$attempts, 3L)
        expect_identical(limit$sd_coordinate, "structurally_omitted_kernel_limit")
        for (a in limit$attempts) {
          expect_false(a$global_fit_eligible)
          expect_length(a$initial$log_scale_sd, 0L)
          expect_identical(a$sd_coordinate, "structurally_omitted_kernel_limit")
          if (!is.null(a$parameters)) expect_length(a$parameters$log_scale_sd, 0L)
        }
      }
    }
  }
  expect_gt(length(audited), 0)
  for (point in trace$points) {
    accepted <- Filter(function(a) isTRUE(a$audited_value), point$attempts)
    if (length(accepted)) expect_equal(accepted[[1L]]$objective,
      annual_profile_reference(f$spec, accepted[[1L]]$parameters, point$sd), tolerance = 1e-7)
  }
  # Finite-difference all nuisance directions at a retained positive-SD vector.
  positive <- Filter(function(x) x$sd > 0, audited)
  expect_gt(length(positive), 0)
  if (length(positive)) {
    a <- positive[[1L]]$attempt; sd <- positive[[1L]]$sd
    problem <- multiScaleR:::.msr_fixed_sd_profile_problem(f$spec, sd, a$parameters, "scale")
    reference <- function(x) annual_profile_reference(f$spec, problem$parameters(x), sd)
    h <- 1e-5
    score <- vapply(seq_along(a$par), function(j) {
      hi <- lo <- a$par; hi[j] <- hi[j] + h; lo[j] <- lo[j] - h
      (reference(hi) - reference(lo)) / (2 * h)
    }, 0)
    expect_lt(max(abs(score - a$full_gradient[a$free_indices])), 1e-4)
    expect_gt(abs(a$full_gradient[a$fixed_indices]), 1e-4)
  }
  zero <- multiScaleR:::.msr_fixed_sd_profile_problem(f$spec, 0, f$parameters, "scale")
  p0 <- zero$parameters(zero$par)
  expect_length(p0$scale_z, 0L); expect_length(p0$log_scale_sd, 0L)
  expect_equal(zero$fn(zero$par), annual_profile_reference(f$spec, p0, 0), tolerance = 1e-9)
  expect_equal(annual_profile_reference(f$spec, p0, 1e-8),
    annual_profile_reference(f$spec, p0, 0), tolerance = 1e-7)
  for (kernel in c("broad", "nearest")) {
    limit <- multiScaleR:::.msr_factorized_limit_spec(f$spec, kernel)
    values <- vapply(c(0, .2, .5), function(sd) {
      p <- f$parameters; p$scale_beta <- numeric()
      problem <- multiScaleR:::.msr_fixed_sd_profile_problem(limit, sd, p, "scale", allow_irrelevant = TRUE)
      expect_length(problem$fixed, 0L)
      expect_length(problem$parameters(problem$par)$log_scale_sd, 0L)
      problem$fn(problem$par)
    }, 0)
    expect_equal(values, rep(values[1], 3), tolerance = 0)
  }
})

test_that("annual profiles reject unsupported scopes before fitting", {
  f <- annual_sd_profile_fixture(); profile <- multiScaleR:::.msr_factorized_scale_sd_profile
  for (grid in list(numeric(), .2, c(0, 0), c(0, -1), c(0, NA), matrix(0)))
    expect_error(profile(f$spec, grid), "sd")
  expect_error(profile(f$spec, c(0, .2), list(cutoff = 3.84)), "control")
  expect_error(profile(f$spec, c(0, .2), initial_parameters = list(list(beta = 1))), "blocks")
  posterior <- f$spec; posterior$priors <- list(beta_sd = 2)
  expect_error(profile(posterior, c(0, .2)), "unpenalized")
  counts <- f$spec; counts$family <- "nbinom2"
  expect_error(profile(counts, c(0, .2)), "Gaussian")
  response <- multiScaleR:::.msr_marginal_prepare(count ~ habitat + temp + (1 | year),
    f$inputs, f$data, list(habitat = ~1 + (1 | year)), "gaussian", "key", 200)
  expect_error(profile(response, c(0, .2)), "no response group")
  expect_error(multiScaleR:::.msr_fixed_sd_profile_problem(f$spec, .2, f$parameters, "scale",
    allow_irrelevant = TRUE), "exact kernel limit")
})

test_that("annual conditional selection retains lower and competing candidates", {
  f <- annual_sd_profile_fixture(); calls <- 0L
  local_mocked_bindings(.msr_fixed_sd_profile_attempt = function(spec, sd, initial, iter.max, ...) {
    calls <<- calls + 1L; i <- (calls - 1L) %% 3L + 1L
    list(eligible_nuisance = i != 2L, audited_value = TRUE, objective = c(0, -1, .5)[i],
      par = c(0, 0), nuisance_hessian = diag(2), parameters = initial, global_fit_eligible = FALSE)
  }, .msr_scale_profile_limits = function(...) list(completed = TRUE, witnesses = character()), .package = "multiScaleR")
  result <- multiScaleR:::.msr_factorized_scale_sd_profile(f$spec, c(0, .2))
  expect_false(any(result$curve$available))
  expect_true(all(vapply(result$points, function(p) p$better_ineligible_candidate, TRUE)))
  calls <- 0L
  local_mocked_bindings(.msr_fixed_sd_profile_attempt = function(spec, sd, initial, iter.max, ...) {
    calls <<- calls + 1L; i <- (calls - 1L) %% 3L + 1L
    list(eligible_nuisance = TRUE, audited_value = TRUE, objective = c(0, 1, .5)[i],
      par = c(if (i == 2L) 2 else 0, 0), nuisance_hessian = diag(2), parameters = initial, global_fit_eligible = FALSE)
  }, .package = "multiScaleR")
  result <- multiScaleR:::.msr_factorized_scale_sd_profile(f$spec, c(0, .2))
  expect_false(any(result$curve$available))
  expect_true(all(vapply(result$points, function(p) p$competing_nuisance_optima, TRUE)))
  local_mocked_bindings(.msr_fixed_sd_profile_attempt = function(spec, sd, initial, iter.max, ...)
    list(eligible_nuisance = TRUE, objective = 0, par = c(0, 0), nuisance_hessian = diag(2), parameters = initial),
    .msr_scale_profile_limits = function(...) list(completed = FALSE, witnesses = character()), .package = "multiScaleR")
  result <- multiScaleR:::.msr_factorized_scale_sd_profile(f$spec, c(0, .2))
  expect_false(any(result$curve$available)); expect_match(result$curve$reason[1], "audit failed")
  local_mocked_bindings(.msr_scale_profile_limits = function(...) list(completed = TRUE, witnesses = "broad"), .package = "multiScaleR")
  result <- multiScaleR:::.msr_factorized_scale_sd_profile(f$spec, c(0, .2))
  expect_false(any(result$curve$available)); expect_match(result$curve$reason[1], "witnesses")
})

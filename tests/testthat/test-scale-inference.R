testthat::skip_if_not_installed("MASS")

test_that("full profiles retain distinct regions and outside brackets", {
  x <- c(1, 2, 3, 4, 5, 6, 7)
  z <- .msr_scale_components(x, c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE))
  expect_equal(z$lower, c(1, 5))
  expect_equal(z$upper, c(2, 6))
  expect_equal(z$upper_outside, c(3, 7))
  expect_identical(z$lower_bound_limited, c(TRUE, FALSE))
  expect_null(.msr_scale_components(x, c(TRUE, NA, rep(FALSE, 5))))
  expect_equal(nrow(.msr_scale_components(x, rep(FALSE, 7))), 0L)
})

test_that("NB2 profile inference refines likelihood crossings without mutating fits", {
  f <- profile_glm_fixture(n = 180)
  x <- suppressWarnings(profile_glm_fit(f))
  before <- suppressWarnings(serialize(x, NULL))
  a <- scale_inference(x, grid_points = 13, n_starts = 4)
  expect_identical(suppressWarnings(serialize(x, NULL)), before)
  expect_s3_class(a, "msr_scale_inference")
  expect_lte(a$reference$value, x$optim_results$value + 1e-7)
  expect_equal(range(a$grid$sigma), a$bounds)
  expect_true(any(a$grid$profile_accept))
  setup <- .msr_inference_setup(x, 4)
  for (j in seq_len(nrow(a$profile_set))) for (side in c("lower", "upper")) {
    if (!a$profile_set[[paste0(side, "_bound_limited")]][j]) {
      p <- .msr_inference_refit(setup, fixed_scale =
        log(a$profile_set[[side]][j] / setup$unit))
      expect_equal(2 * (p$value - a$reference$value), qchisq(.95, 1), tolerance = 1e-5)
    }
  }
  # Independent fixed-design NB2 fit verifies the zero-effect likelihood.
  ref <- MASS::glm.nb(y ~ z + group + offset(log(exposure)), data = f$data)
  expect_equal(a$diagnostics$zero_effect_deviance,
    2 * (-as.numeric(logLik(ref)) - a$reference$value), tolerance = 1e-6)
  expect_error(scale_inference(x, sigma_values = 900), "within the fitting bounds")
  expect_error(scale_inference(x, bootstrap = 1, seed = 1), "too small")
  expect_error(scale_inference(x, bootstrap_values = 250), "requires")
})

test_that("weak-effect inference finds the boundary and reports an unresolved scale", {
  f <- profile_glm_fixture(n = 300, seed = 43008, effect = .05)
  x <- suppressWarnings(profile_glm_fit(f))
  a <- scale_inference(x, grid_points = 17, n_starts = 4)
  expect_gt(a$reference$improvement, .03)
  expect_equal(a$estimate, 800, tolerance = 1e-5)
  expect_identical(a$diagnostics$profile_status, "unresolved")
  expect_true(a$diagnostics$zero_effect_compatible)
  expect_true(all(a$grid$profile_accept))
  expect_equal(a$profile_set$lower, 50, tolerance = 1e-8)
  expect_equal(a$profile_set$upper, 800, tolerance = 1e-8)
  expect_true(a$profile_set$lower_bound_limited && a$profile_set$upper_bound_limited)
})

test_that("bootstrap refits are reproducible, audit failures, and preserve RNG", {
  f <- profile_glm_fixture(n = 80)
  x <- suppressWarnings(profile_glm_fit(f))
  set.seed(919)
  state <- .Random.seed
  a <- scale_inference(x, sigma_values = 250, level = .8, bootstrap = 9,
                       bootstrap_values = 250, seed = 71, n_starts = 2)
  expect_identical(.Random.seed, state)
  b <- scale_inference(x, sigma_values = 250, level = .8, bootstrap = 9,
                       bootstrap_values = 250, seed = 71, n_starts = 2)
  expect_equal(a$bootstrap_records, b$bootstrap_records)
  expect_equal(nrow(a$bootstrap_records), 9)
  expect_true(all(a$bootstrap_records$ok))
  expect_true(all(a$bootstrap_records$constrained_nll >=
                    a$bootstrap_records$full_nll - 1e-5))
  p <- a$grid$bootstrap_p[a$grid$sigma == 250]
  observed <- a$grid$lr[a$grid$sigma == 250]
  expect_equal(p, (1 + sum(a$bootstrap_records$lr >= observed - 1e-8)) / 10)
  expect_null(a$bootstrap_set)
  expect_identical(a$diagnostics$bootstrap_status, "pointwise_only")
  full <- scale_inference(x, sigma_values = c(50, 800), level = .5,
    bootstrap = 3, seed = 16, n_starts = 2)
  expect_true(is.data.frame(full$bootstrap_set))
  expect_false(anyNA(full$grid$bootstrap_p))
  expect_true(all(is.na(.msr_bootstrap_tail(c(1, NA), 2, .95))))
  expect_equal(.msr_bootstrap_tail(rep(0, 19), 3, .95)["p"], c(p = .05))
  weighted <- x
  weighted$profile_fit$contract$objective$weights[1] <- 2
  expect_error(scale_inference(weighted, bootstrap = 19, seed = 1), "unit observation")
  bad <- x
  bad$profile_fit$contract$objective$family <- "poisson"
  expect_error(scale_inference(bad), "NB2 with one spatial")
})

test_that("failed bootstrap pairs remain visible and dispersion cache is exact", {
  f <- profile_glm_fixture(n = 60)
  x <- suppressWarnings(profile_glm_fit(f))
  setup <- .msr_inference_setup(x, 2)
  obj <- setup$obj
  for (theta in c(.5, 20, 2)) {
    par <- c(log(.25), log(theta))
    fresh <- .msr_profile_glm_objective(obj$engine, obj$y, obj$family,
      obj$X_fixed, obj$offset, obj$weights, obj$trials)
    expect_equal(obj$fn(par), fresh$fn(par), tolerance = 1e-12)
    expect_equal(obj$gr(par), fresh$gr(par), tolerance = 1e-12)
  }
  reference <- list(par = x$profile_fit$par, value = x$optim_results$value)
  local_mocked_bindings(.msr_inference_refit = function(...) stop("forced refit failure"))
  result <- .msr_scale_bootstrap(setup, reference, reference, 2)
  expect_equal(nrow(result), 2L)
  expect_false(any(result$ok))
  expect_true(all(grepl("forced refit failure", result$error)))
  expect_true(all(is.na(.msr_bootstrap_tail(result$lr, 1, .95))))
})

test_that("zero-effect diagnostics recover oscillating reduced NB2 fits", {
  for (seed in c(733072L, 733090L)) {
    # Initial conditional GLM iterations can warn for these difficult data.
    f <- suppressWarnings(profile_glm_fixture(n = 100, seed = seed, theta = .5))
    x <- suppressWarnings(profile_glm_fit(f))
    setup <- .msr_inference_setup(x, 2)
    obj <- setup$obj
    value <- suppressWarnings(.msr_nb2_zero_effect(setup))
    # Independent joint optimization uses numerical derivatives, without
    # conditional NB2 GLM iterations. glm.nb also fails on the second case.
    beta <- glm.fit(obj$X_fixed, obj$y, offset = obj$offset,
                    family = poisson())$coefficients
    likelihood <- function(par) -sum(dnbinom(obj$y,
      mu = exp(drop(obj$X_fixed %*% head(par, -1)) + obj$offset),
      size = exp(tail(par, 1)), log = TRUE))
    initial <- nlminb(c(beta, log(.3)), likelihood,
      lower = c(rep(-Inf, length(beta)), tail(setup$lower, 1)),
      upper = c(rep(Inf, length(beta)), tail(setup$upper, 1)),
      control = list(iter.max = 1000L, eval.max = 2000L, rel.tol = 1e-12))
    independent <- optim(initial$par, likelihood, method = "BFGS",
      control = list(maxit = 2000L, reltol = 1e-11))
    expect_equal(independent$convergence, 0L)
    expect_equal(value, independent$value, tolerance = 1e-7)
    score <- vapply(seq_along(independent$par), function(j) {
      h <- rep(0, length(independent$par)); h[j] <- 1e-5
      (likelihood(independent$par + h) - likelihood(independent$par - h)) / 2e-5
    }, numeric(1))
    expect_lt(max(abs(score)), 1e-3)
    a <- suppressWarnings(scale_inference(x, sigma_values = c(50, 250, 800), n_starts = 2))
    expect_true(is.finite(a$diagnostics$zero_effect_deviance))
    expect_lte(max(a$grid$lr), a$diagnostics$zero_effect_deviance + 1e-5)
  }
})

test_that("no-intercept NB2 profiles recover oscillating conditional GLM fits", {
  f <- profile_glm_fixture(n = 100, seed = 9822)
  f$model <- glm(y ~ hab - 1 + offset(log(exposure)), data = f$data,
    family = MASS::negative.binomial(2))
  x <- suppressWarnings(profile_glm_fit(f))
  setup <- .msr_inference_setup(x, 2)
  for (theta in c(.3, 2, 20)) {
    p <- c(log(.05), log(theta))
    value <- suppressWarnings(setup$obj$evaluate(p))
    X <- value$design
    # A one-dimensional derivative-free optimization is an independent check.
    objective <- function(beta) -sum(dnbinom(setup$obj$y,
      mu = exp(drop(X * beta) + setup$obj$offset), size = theta, log = TRUE))
    independent <- optimize(objective, c(-20, 20), tol = 1e-10)
    expect_equal(value$value, independent$objective, tolerance = 1e-8)
    expect_equal(unname(value$beta), independent$minimum, tolerance = 1e-6)
    numeric_score <- (suppressWarnings(setup$obj$fn(p + c(1e-5, 0))) -
                       suppressWarnings(setup$obj$fn(p - c(1e-5, 0)))) / 2e-5
    expect_equal(value$gradient[1], numeric_score, tolerance = 1e-5)
  }
  a <- suppressWarnings(scale_inference(x, grid_points = 5, n_starts = 2))
  expect_true(is.finite(a$diagnostics$zero_effect_deviance))
  expect_equal(range(a$grid$sigma), c(50, 800), tolerance = 1e-8)
})

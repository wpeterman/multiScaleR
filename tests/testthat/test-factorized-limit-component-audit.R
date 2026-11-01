test_that("optimized raw exact-limit witnesses agree with normalized native Gaussian ML", {
  skip_if_not_installed("RTMB")
  s <- path_component_fixture("gaussian", "gaussian", FALSE, "tied", "canonical")
  p <- private_sq_parameters(s); p$beta[3L] <- .12
  p$scale_z <- p$response_z <- numeric()
  before <- list(s, p, .Random.seed)
  candidate <- .msr_factorized_density_component(s)(p)
  audit <- .msr_factorized_limit_component_audit(s, p, candidate)
  expect_true(audit$completed)
  expect_identical(audit[c("component_only", "target", "finite_scale_estimate", "uncertainty_available")],
    list(component_only = TRUE, target = "marginal_likelihood_limit",
      finite_scale_estimate = FALSE, uncertainty_available = FALSE))
  maximum <- 0
  for (kernel in c("broad", "nearest")) {
    a <- audit$references[[paste(kernel, "gaussian", sep = "-")]]
    expect_length(a$attempts, 3L)
    expect_true(all(vapply(a$attempts, function(x) isTRUE(x$audited), TRUE)))
    # Independent cell-based endpoint, then ordinary Gaussian GLM ML solution.
    raw <- vapply(s$surface$buffers[[1L]], function(b) {
      v <- if (kernel == "nearest") b$v[b$d == min(b$d)] else b$v
      sqrt(mean((v - mean(v))^2))
    }, 0)
    k <- (raw - s$center[[1L]]) / s$spread[[1L]]
    native <- stats::glm.fit(cbind(s$X_fixed, k), s$y,
      offset = s$offset, family = stats::gaussian())
    sigma <- sqrt(mean((s$y - native$fitted.values)^2))
    optimum <- -sum(stats::dnorm(s$y, native$fitted.values, sigma, log = TRUE))
    expect_lt(abs(a$objective - optimum), 1e-6)
    for (attempt in a$attempts) {
      reference <- sum(mixture_native_nll(s, attempt$parameters, kernel))
      maximum <- max(maximum, abs(attempt$objective - reference))
      expect_lt(max(abs(c(attempt$objective, attempt$ad_objective) - reference)), 1e-6)
      expect_identical(lengths(attempt$parameters)[c("scale_beta", "scale_z", "log_scale_sd", "response_z")],
        c(scale_beta = 0L, scale_z = 0L, log_scale_sd = 0L, response_z = 0L))
    }
  }
  expect_identical(audit$infinite_sd_path$direction, "original_unrestricted_model")
  expect_identical(list(s, p, .Random.seed), before)
  cat(sprintf("\nOPTIMIZED_LIMIT_COMPONENT preparations=1 endpoint_attempts=6 native_error=%.12g\n", maximum))
})

test_that("component orchestration retains starts, zero faces, Poisson coordinates and failures", {
  # These are software orchestration mocks, not accepted optimizer references.
  # Real ordinary point values and exact path leaves remain active. Deliberate
  # nonconvergence and singular curvature must survive as feasible witnesses.
  seen <- list()
  testthat::local_mocked_bindings(.msr_factorized_adfun = function(spec, parameters, ...) {
    seen[[length(seen) + 1L]] <<- c(raw = length(spec$surface$buffers),
      scale_beta = length(spec$parameters$scale_beta))
    par <- unlist(parameters, use.names = FALSE)
    env <- new.env(); env$last.par <- par
    env$parList <- function(x, par) private_sq_unpack(x, parameters)
    ordinary <- .msr_factorized_density_component(spec)
    list(par = par, env = env,
      fn = function(x) ordinary(private_sq_unpack(x, parameters)),
      gr = function(x) rep(0, length(x)),
      he = function(x) matrix(0, length(x), length(x)))
  }, .msr_marginal_polish = function(...) list(passed = FALSE), .package = "multiScaleR")
  testthat::local_mocked_bindings(nlminb = function(start, objective, gradient, control) {
    stopifnot(identical(control, list(iter.max = 250L, eval.max = 750L, rel.tol = 1e-12)))
    list(par = start, objective = objective(start), convergence = 1L,
      message = "deliberate software nonconvergence")
  }, .package = "stats")
  base_density_factory <- .msr_factorized_audit_density
  density_calls <- character()
  testthat::local_mocked_bindings(.msr_factorized_audit_density = function(spec, ad = FALSE,
      rel.tol = NULL, step = NULL, coordinate = c("standard", "physical")) {
    density_calls <<- c(density_calls, match.arg(coordinate))
    base_density_factory(spec, ad, rel.tol, step, coordinate)
  }, .package = "multiScaleR")
  extra <- c("component_only", "target", "finite_scale_estimate", "uncertainty_available")
  mean <- mixture_fixture("gaussian", "gaussian", FALSE)$spec
  p <- private_sq_parameters(mean); p$scale_z <- p$response_z <- numeric()
  strict <- .msr_factorized_limit_audit(mean, p, 1000)
  component <- .msr_factorized_limit_component_audit(mean, p, 1000)
  expect_identical(strict, component[setdiff(names(component), extra)])
  cases <- 0L
  for (zero in c(FALSE, TRUE)) {
    s <- path_component_fixture(if (zero) "exp" else "gaussian", "nbinom2", FALSE,
      if (zero) "unique" else "tied", if (zero) "contrast" else "canonical")
    p <- private_sq_parameters(s); p$beta[3L] <- .12
    p$scale_z <- p$response_z <- numeric()
    p$log_scale_sd <- if (zero) numeric() else log(9)
    before <- list(s, p, .Random.seed)
    previous_calls <- length(density_calls)
    audit <- .msr_factorized_limit_component_audit(s, p, 1000)
    called <- density_calls[seq.int(previous_calls + 1L, length(density_calls))]
    expect_identical(as.integer(table(factor(called, levels = c("standard", "physical")))),
      c(13L, if (zero) 0L else 1L))
    if (!zero) expect_true(exp(p$log_scale_sd) >= 8)
    expect_true(audit$completed)
    expect_identical(audit$infinite_sd_path$candidate_zero_sd, zero)
    expect_length(audit$references, 11L)
    expect_length(audit$witnesses, 11L)
    expected_path <- .msr_factorized_path_component(s, p)
    expect_identical(audit$references[names(expected_path$references)], expected_path$references)
    for (label in c("broad-nbinom2", "broad-poisson", "nearest-nbinom2", "nearest-poisson")) {
      a <- audit$references[[label]]; cases <- cases + 1L
      expect_length(a$attempts, 3L)
      expect_true(all(vapply(a$attempts, function(x) isTRUE(x$audited) &&
        identical(x$optimizer$convergence, 1L) && !x$curvature$passed && !x$refinement$passed, TRUE)))
      initial <- a$attempts[[1L]]$initial
      expect_equal(a$attempts[[2L]]$initial, initial + .05 * sin(seq_along(initial) + 2L), tolerance = 0)
      expect_equal(a$attempts[[3L]]$initial, initial + .05 * sin(seq_along(initial) + 3L), tolerance = 0)
    }
    poisson <- if (zero) .msr_marginal_reduce(s, "scale") else s
    poisson$family <- "poisson"
    # The source template still has NB2 dispersion here; only the explicit
    # family-reduction callback removes it. Generic density remains strict.
    q <- p; q$log_dispersion <- numeric()
    expect_error(.msr_factorized_density_component(poisson)(q), "population blocks")
    standard <- .msr_factorized_audit_density(poisson, FALSE, 1e-10, .5)
    physical <- .msr_factorized_audit_density(poisson, coordinate = "physical")
    expect_equal(audit$references$`fitted-scales-poisson`$objective, standard(q), tolerance = 1e-10)
    expect_equal(audit$references$`fitted-scales-poisson`$coordinate_error,
      if (zero) 0 else abs(standard(q) - physical(q)), tolerance = 1e-10)
    expect_identical(attr(standard, "component")[c("rel.tol", "step")], list(rel.tol = 1e-10, step = .5))
    expect_identical(attr(physical, "component")[c("rel.tol", "step")], list(rel.tol = 1e-10, step = 2))
    expect_identical(list(s, p, .Random.seed), before)
  }
  expect_identical(cases, 8L)
  # Retain failures of two endpoint factories alongside completed starts and
  # infinite-SD paths. Also retain a failed independent coordinate comparison.
  s <- path_component_fixture("gaussian", "nbinom2", FALSE, "tied", "canonical")
  p <- private_sq_parameters(s); p$beta[3L] <- .12
  p$scale_z <- p$response_z <- numeric(); p$log_scale_sd <- log(9)
  kernel_factory <- .msr_factorized_kernel_component
  density_factory <- .msr_factorized_audit_density
  testthat::local_mocked_bindings(.msr_factorized_kernel_component = function(spec, kernel, ...) {
    if (kernel == "nearest") stop("mock exact endpoint failure")
    kernel_factory(spec, kernel, ...)
  }, .msr_factorized_audit_density = function(spec, ad = FALSE, rel.tol = NULL, step = NULL,
      coordinate = c("standard", "physical")) {
    if (match.arg(coordinate) == "physical") stop("mock physical coordinate failure")
    density_factory(spec, ad, rel.tol, step, coordinate)
  }, .package = "multiScaleR")
  failed <- .msr_factorized_limit_component_audit(s, p, 1000)
  expect_false(failed$completed)
  expect_length(failed$references, 11L)
  expect_identical(vapply(failed$references[c("nearest-nbinom2", "nearest-poisson")], `[[`, "", "error"),
    c(`nearest-nbinom2` = "mock exact endpoint failure", `nearest-poisson` = "mock exact endpoint failure"))
  expect_match(failed$references$`fitted-scales-poisson`$error, "mock physical coordinate failure")
  expect_length(failed$references$`broad-nbinom2`$attempts, 3L)
  expect_length(failed$references$`broad-poisson`$attempts, 3L)
  expect_true(all(vapply(failed$references[grep("^annual-sd-infinite-", names(failed$references))],
    function(x) isTRUE(x$completed), TRUE)))
  expect_false(any(c("nearest-nbinom2", "nearest-poisson", "fitted-scales-poisson") %in% failed$witnesses))
  expect_true(length(seen) > 0L && all(vapply(seen, function(x) all(x == 0L), TRUE)))
})

test_that("explicit component auditing preserves strict admission and neutral source restrictions", {
  s <- path_component_fixture("gaussian", "gaussian", TRUE, "tied", "canonical")
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  expect_error(.msr_factorized_limit_audit(s, p, 100), "Weighted SQ")
  expect_error(.msr_factorized_adfun(s, p), "Weighted SQ")
  expect_error(.msr_factorized_fit(s), "Weighted SQ")
  expect_error(.msr_factorized_audit_density(s, ad = TRUE), "ordinary numeric")
  bad <- s; bad$priors <- list(beta = 1)
  expect_error(.msr_factorized_limit_component_audit(bad, p, 100), "neutral")
  bad <- s; bad$response$random[[1L]]$index[] <- 1L
  expect_error(.msr_factorized_limit_component_audit(bad, p, 100), "coupling")
  bad <- s; bad$rows <- bad$rows[-1L]
  expect_error(.msr_factorized_limit_component_audit(bad, p, 100), "metadata")
  expect_error(.msr_factorized_limit_component_audit(s, p, Inf), "finite numeric candidate")
  expect_identical(names(formals(.msr_factorized_limit_audit)), c("spec", "parameters", "candidate_nll"))
})

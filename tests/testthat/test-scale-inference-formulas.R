# Bounded software tests using existing fixtures, not a calibration study.
inference_product_object <- function(boundary = FALSE, warm = FALSE, weighted = FALSE) {
  f <- if (weighted) formula_weighted_fixture(surface_only = TRUE) else
    formula_integration_fixture("nbinom2")
  if (weighted) f$model <- glm(f$formula, data = f$data,
    family = MASS::negative.binomial(2))
  contract <- .msr_profile_contract(f$model, f$inputs)
  control <- .msr_profile_control(list(n_starts = 3L,
    coefficient_start = if (warm) "warm" else "cold"), joint_nb2 = boundary)
  if (boundary) {
    obj <- contract$objective
    contract$objective <- .msr_joint_nb2_objective(obj$engine, obj$y, obj$family,
      obj$X_fixed, obj$offset, obj$weights, obj$trials,
      coefficient_start = control$coefficient_start,
      response_design = obj$response_design)
  } else if (warm) {
    contract$objective <- .msr_profile_nb2_continuation(contract$objective)
  }
  builder <- if (boundary) .msr_joint_nb2_bounds else .msr_profile_starts
  bounds <- builder(contract, f$inputs, NULL, control, f$model)
  object <- structure(list(backend = if (boundary) "joint_nb2" else "profile_glm",
    opt_mod = f$model, kernel_inputs = f$inputs,
    profile_fit = c(list(contract = contract, control = control,
      par = bounds$starts[nrow(bounds$starts), ]), bounds)), class = "multiScaleR")
  list(object = object, data = f$data, formula = formula(f$model))
}

test_that("NB2 scale refits retain products in cold and warm coordinates", {
  for (weighted in c(FALSE, TRUE)) for (boundary in c(FALSE, TRUE)) for (warm in c(FALSE, TRUE)) {
    f <- inference_product_object(boundary, warm, weighted)
    original <- f$object$profile_fit$contract$objective
    setup <- .msr_inference_setup(f$object, 2L)
    expect_identical(setup$obj$response_design, original$response_design)
    expect_false(identical(environment(setup$obj$fn), environment(original$fn)))
    for (theta in c(.5, 2, 6, if (boundary) Inf)) {
      p <- c(log(.3), if (boundary) log1p(1 / theta) else log(theta))
      expect_equal(setup$obj$fn(p), original$fn(p), tolerance = 1e-12)
      expect_equal(setup$obj$gr(p), original$gr(p), tolerance = 1e-12)
    }
    at <- log(.3)
    result <- .msr_inference_refit(setup, fixed_scale = at)
    kd <- .msr_eval_kernel_engine(setup$obj$engine, at)
    nd <- f$data
    nm <- setup$obj$engine$covariates[1L]
    nd[[nm]] <- kd$value[, nm]
    native <- MASS::glm.nb(f$formula, data = nd,
      control = glm.control(epsilon = 1e-11, maxit = 100L))
    expect_lt(abs(result$value + as.numeric(logLik(native))), 1e-5)
    expect_true(any(result$records$eligible))
    native_zero <- MASS::glm.nb(y ~ group + z + offset(log(exposure)), data = nd,
      control = glm.control(epsilon = 1e-11, maxit = 100L))
    expect_lt(abs(.msr_nb2_zero_effect(setup) + as.numeric(logLik(native_zero))), 1e-5)
  }
})

test_that("bootstrap target reconstruction keeps the full native coefficient basis", {
  audit <- function(boundary, warm, weighted) {
    f <- inference_product_object(boundary, warm, weighted)
    setup <- .msr_inference_setup(f$object, 2L)
    p <- c(log(.3), if (boundary) log1p(.5) else log(2))
    reference <- list(par = p, value = setup$obj$fn(p))
    calls <- 0L
    # Inspect every reconstructed target. Optimization itself is covered by
    # real constrained fits above and the existing end-to-end bootstrap tests.
    local_mocked_bindings(.msr_inference_refit = function(setup, obj, starts,
                                                       fixed_scale = NULL) {
      calls <<- calls + 1L
      stopifnot(identical(obj$response_design, setup$obj$response_design))
      kd <- .msr_eval_kernel_engine(obj$engine, p[1L])
      X <- .msr_response_design_matrix(obj$response_design, kd$value)
      native <- glm.fit(X, obj$y, offset = obj$offset, weights = obj$weights,
        family = MASS::negative.binomial(2),
        control = glm.control(epsilon = 1e-11, maxit = 100L))
      nll <- -sum(dnbinom(obj$y, mu = native$fitted.values, size = 2, log = TRUE))
      stopifnot(isTRUE(native$converged), abs(obj$fn(p) - nll) < 1e-5)
      list(par = p, value = nll, records = data.frame(eligible = TRUE))
    }, .package = "multiScaleR")
    # Two pairs only audit software plumbing, not bootstrap accuracy.
    result <- .msr_scale_bootstrap(setup, reference, reference, 2L)
    expect_equal(calls, 4L)
    expect_equal(nrow(result), 2L)
    expect_true(all(result$ok), info = paste(result$error, collapse = "; "))
    expect_equal(result$lr, c(0, 0))
    expect_equal(result$starts_attempted, c(2L, 2L))
  }
  set.seed(173L)
  for (weighted in c(FALSE, TRUE)) for (boundary in c(FALSE, TRUE)) for (warm in c(FALSE, TRUE))
    audit(boundary, warm, weighted)
})

test_that("the public product assessment retains all starts and its scale null", {
  f <- formula_integration_fixture("nbinom2")
  fit <- multiScale_optim(f$model, f$inputs, backend = "profile_glm",
    n_cores = 1, verbose = FALSE,
    profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1L)))
  plan <- fit$profile_fit$contract$objective$response_design
  assessment <- scale_inference(fit, sigma_values = c(50, 250, 800), n_starts = 2L)
  expect_identical(fit$profile_fit$contract$objective$response_design, plan)
  expect_s3_class(assessment, "msr_scale_inference")
  expect_lte(assessment$reference$value, fit$optim_results$value + 1e-7)
  expect_equal(range(assessment$grid$sigma), assessment$bounds)
  expect_identical(assessment$diagnostics$bootstrap_status, "not_requested")
  expect_equal(ncol(fit$profile_fit$contract$objective$X_fixed), 3L)
  expect_false(any(grepl("hab", colnames(fit$profile_fit$contract$objective$X_fixed), fixed = TRUE)))
  expect_null(assessment$bootstrap_set)
})

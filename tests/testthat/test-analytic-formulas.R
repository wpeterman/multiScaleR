test_that("profile products preserve native family likelihoods and full joint scores", {
  for (family in c("poisson", "gaussian", "binomial", "nbinom2")) {
    f <- formula_integration_fixture(family, two = TRUE)
    contract <- .msr_profile_contract(f$model, f$inputs)
    obj <- contract$objective
    expect_s3_class(obj$response_design, "msr_response_design")
    par <- c(log(c(.3, .4)), if (family == "nbinom2") log(2))
    a <- obj$evaluate(par)
    nd <- f$data; kd <- .msr_eval_kernel_engine(obj$engine, par[seq_len(2)])
    for (nm in obj$engine$covariates) nd[[nm]] <- kd$value[, nm]
    reference <- glm(formula(f$model), data = nd, family = f$model$family,
      control = glm.control(epsilon = 1e-11, maxit = 100))
    expect_equal(a$value, -as.numeric(logLik(reference)), tolerance = 1e-9)
    expect_equal(a$beta, coef(reference), tolerance = 1e-9)
    numerical <- vapply(seq_along(par), function(j) {
      x <- y <- par; x[j] <- x[j] + 1e-5; y[j] <- y[j] - 1e-5
      (obj$fn(x) - obj$fn(y)) / 2e-5
    }, 0)
    expect_equal(obj$gr(par), numerical, tolerance = 1e-5, ignore_attr = TRUE)
    joint <- .msr_profile_joint(obj, par)
    numerical <- vapply(seq_along(joint$par), function(j) {
      x <- y <- joint$par; x[j] <- x[j] + 1e-5; y[j] <- y[j] - 1e-5
      (joint$fn(x) - joint$fn(y)) / 2e-5
    }, 0)
    expect_equal(joint$gr(joint$par), numerical, tolerance = 1e-5, ignore_attr = TRUE)
  }
})

test_that("expanded profile fitting retains predictions, covariance and snapshots", {
  f <- formula_integration_fixture()
  starts <- matrix(c(.08, .25, .6), ncol = 1)
  fit <- multiScale_optim(f$model, f$inputs, backend = "profile_glm", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = starts))
  expect_equal(exp(fit$profile_fit$starts), starts, ignore_attr = TRUE)
  expect_equal(unname(predict(fit)), unname(predict(fit$opt_mod, type = "response")), tolerance = 1e-8)
  conditional <- predict(fit, type = "link", se.fit = TRUE, uncertainty = "conditional")
  reference <- predict(fit$opt_mod, type = "link", se.fit = TRUE)
  # Native IRLS and numerically differentiated joint information agree at the
  # established numerical SE accuracy; this is not bitwise covariance identity.
  cat(sprintf("FORMULA_CONDITIONAL_SE max_absolute=%.12g max_relative=%.12g\n",
    max(abs(conditional$se.fit - reference$se.fit)),
    max(abs(conditional$se.fit / reference$se.fit - 1))))
  expect_equal(unname(conditional$se.fit), unname(reference$se.fit), tolerance = 1e-6)
  snapshot <- prepare_prediction(fit)
  expect_identical(predict(snapshot), predict(fit))
  expect_identical(predict(snapshot, se.fit = TRUE, uncertainty = "conditional"),
    predict(fit, se.fit = TRUE, uncertainty = "conditional"))
  if (fit$uncertainty$available) expect_identical(predict(snapshot, se.fit = TRUE), predict(fit, se.fit = TRUE))
  else expect_error(predict(fit, se.fit = TRUE))
  order <- rev(seq_len(nrow(f$data)))
  expect_equal(unname(predict(fit, f$data[order, ], f$inputs)),
    unname(predict(fit)[order]), tolerance = 1e-10)
  plan <- fit$profile_fit$contract$objective$response_design
  d <- .msr_profile_prediction_design(fit, NULL, NULL)
  beta <- fit$uncertainty$par[fit$uncertainty$beta_index]
  objective <- fit$profile_fit$contract$objective
  k <- fit$profile_fit$par[seq_len(length(objective$engine$parameter_names))]
  eta <- function(k) drop(.msr_response_design_matrix(plan,
    .msr_eval_kernel_engine(objective$engine, k)$value) %*% beta)
  expect_equal(drop(d$ds * d$kernel_sensitivity), (eta(k + 1e-5) - eta(k - 1e-5)) / 2e-5,
    tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(.msr_selection_parameter_count(fit), length(coef(fit$opt_mod)) + 1L)
})

test_that("expanded conditional models retain native cold audits and choice predictions", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture(); dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- as.formula("case_ ~ hab * group + z + offset(log(exposure)) + strata(choice)", env = asNamespace("survival"))
  m <- survival::clogit(form, dat, model = TRUE, x = TRUE, y = TRUE)
  contract <- .msr_clogit_contract(m, f$inputs); obj <- .msr_clogit_objective(contract)
  par <- c(log(.3), rep(.15, obj$nb))
  numerical <- vapply(seq_along(par), function(j) {
    a <- b <- par; a[j] <- a[j] + 1e-5; b[j] <- b[j] - 1e-5
    (obj$fn(a) - obj$fn(b)) / 2e-5
  }, 0)
  expect_equal(obj$gr(par), numerical, tolerance = 1e-6, ignore_attr = TRUE)
  native <- .msr_clogit_native(contract, log(.3))
  expect_equal(obj$fn(c(log(.3), unname(coef(native)[contract$original_columns]))),
    -as.numeric(logLik(native)), tolerance = 1e-10)
  fit <- withCallingHandlers(multiScale_optim(m, f$inputs, backend = "conditional_logit",
    n_cores = 1, verbose = FALSE, profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1))),
    msr_clogit_fit_error = function(e) print(e$start_results))
  expect_true(any(fit$profile_fit$start_results$eligible))
  expect_true(all(is.na(fit$scale_est$SE)))
  probability <- predict(fit)
  expect_equal(as.numeric(rowsum(probability, dat$choice)), rep(1, 80), tolerance = 1e-12)
  expect_equal(predict(fit, dat, f$inputs), probability, tolerance = 1e-12)
  expect_identical(predict(prepare_prediction(fit)), probability)
  expect_identical(predict(prepare_prediction(fit), se.fit = TRUE, uncertainty = "conditional"),
    predict(fit, se.fit = TRUE, uncertainty = "conditional"))
  expect_error(predict(fit, se.fit = TRUE), "conditional SEs")
})

test_that("active terms exclude removed kernels and GAM targets are not linearized", {
  f <- profile_glm_fixture(); f$inputs$binned$vsum$cover <- f$inputs$binned$vsum$hab
  f$inputs$binned$csum$cover <- f$inputs$binned$csum$hab
  f$data$cover <- f$data$hab
  m <- glm(y ~ . - cover - exposure, f$data, family = poisson())
  contract <- .msr_profile_contract(m, f$inputs)
  expect_identical(contract$covariates, "hab")
  expect_null(contract$objective$response_design)
  skip_if_not_installed("mgcv")
  m <- mgcv::gam(as.formula("y ~ hab + s(z,k=5) + group + offset(log(exposure))",env = asNamespace("mgcv")),
    data = f$data, family = poisson(), method = "ML")
  expect_error(.msr_profile_contract(m, f$inputs), "GAM smooth and penalty targets")
  expect_s3_class(m, "gam")
})


test_that("conditional coefficient steps reuse kernels without bypassing cold refits", {
  skip_if_not_installed("survival")
  f <- profile_glm_fixture(); dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- as.formula("Surv(rep(1, 240), case_) ~ hab * group + z + strata(choice)",env=asNamespace("survival"))
  m <- survival::coxph(form, dat, model=TRUE,x=TRUE,y=TRUE,ties="exact")
  contract <- .msr_clogit_contract(m, f$inputs)
  original <- .msr_eval_kernel_engine; calls <- 0L
  local_mocked_bindings(.msr_eval_kernel_engine = function(...) {
    calls <<- calls + 1L; original(...)
  }, .package = "multiScaleR")
  obj <- .msr_clogit_objective(contract)
  p <- c(log(.3), rep(.1, obj$nb))
  obj$fn(p); obj$gr(p)
  expect_identical(calls, 1L)
  p[2] <- p[2] + .01; obj$fn(p); obj$gr(p)
  expect_identical(calls, 1L)
  p[1] <- p[1] + .01; obj$fn(p)
  expect_identical(calls, 2L)
  .msr_clogit_native(contract, p[1])
  expect_identical(calls, 3L)
})

test_that("expanded scale assessments retain the complete response plan", {
  f <- formula_integration_fixture("nbinom2")
  contract <- .msr_profile_contract(f$model, f$inputs)
  control <- .msr_profile_control(list(n_starts = 3L))
  bounds <- .msr_profile_starts(contract, f$inputs, NULL, control, f$model)
  object <- structure(list(backend = "profile_glm", opt_mod = f$model,
    kernel_inputs = f$inputs, profile_fit = c(list(contract = contract,
    control = control, par = bounds$starts[1L, ]), bounds)), class = "multiScaleR")
  setup <- .msr_inference_setup(object, 2L)
  expect_identical(setup$obj$response_design, contract$objective$response_design)
  expect_equal(setup$obj$fn(object$profile_fit$par),
    contract$objective$fn(object$profile_fit$par), tolerance = 1e-12)
  expect_equal(setup$obj$gr(object$profile_fit$par),
    contract$objective$gr(object$profile_fit$par), tolerance = 1e-12)
  expect_true(any(setup$starts[, 1L] == setup$lower[1L]))
  expect_true(any(setup$starts[, 1L] == setup$upper[1L]))
})

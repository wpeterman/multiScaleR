# Bounded native/AD software checks on existing inputs. No sampling study.
test_that("kernel expression leaves retain native algebra and zero derivatives", {
  dat <- data.frame(y = c(1, 2, 1, 3, 2, 4), h = c(0, -.4, .2, .8, -.7, .5),
    k = c(.3, -.2, .5, .9, .7, -.6), z = seq_len(6) / 6,
    group = factor(rep(c("a", "b"), 3)))
  forms <- list(y ~ I(h^2) * group + I(h*k) + z,
    y ~ h:I(h^2) + I(h + 1) + k + z,
    y ~ 0 + I(h^2):group + I(h/(h^2 + 1)) + z,
    y ~ exp(h) + log(h^2 + 1) + sqrt(k^2 + 1) + z,
    y ~ sin(h) + cos(k) + tan(h) + z,
    y ~ base::I(h^2) + base::exp(k) + z)
  for (form in forms) {
    tt <- terms(form); frame <- model.frame(tt, dat)
    design <- .msr_compile_response_design(tt, frame, c("h", "k"))
    kernel <- as.matrix(dat[c("h", "k")])
    X <- model.matrix(tt, frame)
    expect_equal(.msr_response_design_matrix(design, kernel), X,
      tolerance = 1e-12, ignore_attr = TRUE)
    beta <- seq(-.2, .2, length.out = ncol(X))
    expect_equal(.msr_response_design_eta(design, kernel, beta), drop(X %*% beta),
      tolerance = 1e-12, ignore_attr = TRUE)
    sensitivity <- .msr_response_kernel_sensitivity(design, kernel, beta)
    for (j in seq_len(ncol(kernel))) {
      hi <- lo <- kernel; hi[, j] <- hi[, j] + 1e-6; lo[, j] <- lo[, j] - 1e-6
      numerical <- drop((.msr_response_design_matrix(design, hi) -
        .msr_response_design_matrix(design, lo)) %*% beta) / 2e-6
      expect_lt(max(abs(sensitivity[, j] - numerical)), 1e-6)
    }
    updated <- .msr_response_update_frame(design, frame, kernel)
    expect_equal(model.matrix(tt, updated), X, tolerance = 1e-12, ignore_attr = TRUE)
    if (requireNamespace("RTMB", quietly = TRUE)) {
      obj <- RTMB::MakeADFun(function(p) sum(.msr_response_design_eta(design,
        p$kernel, p$beta)^2), list(kernel = kernel, beta = beta), silent = TRUE)
      score <- vapply(seq_along(obj$par), function(j) {
        hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
        (obj$fn(hi) - obj$fn(lo)) / 2e-5
      }, 0)
      expect_lt(max(abs(as.numeric(obj$gr(obj$par)) - score)), 2e-5)
    }
  }
})

test_that("expression eta supports complete AD scores and curvature", {
  skip_if_not_installed("RTMB")
  dat <- data.frame(y = 1:4, h = c(0, -.4, .2, .8), k = c(.3, -.2, .5, .9),
                    group = factor(c("a", "b", "a", "b")))
  tt <- terms(y ~ I(h^2) * group + I(h*k) + sin(h))
  design <- .msr_compile_response_design(tt, model.frame(tt, dat), c("h", "k"))
  parameters <- list(kernel = as.matrix(dat[c("h", "k")]),
    beta = seq(-.2, .2, length.out = design$p))
  obj <- RTMB::MakeADFun(function(p) sum(.msr_response_design_eta(design,
    p$kernel, p$beta)^2), parameters, silent = TRUE)
  numerical <- vapply(seq_along(obj$par), function(j) {
    hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
    (obj$fn(hi) - obj$fn(lo)) / 2e-5
  }, 0)
  expect_lt(max(abs(as.numeric(obj$gr(obj$par)) - numerical)), 2e-5)
  H <- vapply(seq_along(obj$par), function(j) {
    hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
    (as.numeric(obj$gr(hi)) - as.numeric(obj$gr(lo))) / 2e-5
  }, numeric(length(obj$par)))
  expect_lt(max(abs(obj$he(obj$par) - H)), 1e-4)
  local_mocked_bindings(.msr_total_ram_bytes = function() NA_real_, .package = "multiScaleR")
  expect_error(.msr_response_design_eta(design, parameters$kernel, parameters$beta), "Cannot detect RAM")
})

test_that("expression profiles preserve all native family likelihoods and joint scores", {
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- formula_integration_fixture(family, two = TRUE)
    response <- formula(f$model)[[2L]]
    form <- as.formula(call("~", response, quote(I(hab^2) * group + I(hab * cover) +
      sin(hab) + z + offset(log(exposure)))))
    f$model <- glm(form, f$data, family = f$model$family)
    contract <- .msr_profile_contract(f$model, f$inputs); obj <- contract$objective
    p <- c(log(c(.3, .4)), if (family == "nbinom2") log(2))
    a <- obj$evaluate(p); nd <- f$data
    kernel <- .msr_eval_kernel_engine(obj$engine, p[1:2])
    nd$hab <- kernel$value[, "hab"]; nd$cover <- kernel$value[, "cover"]
    native <- glm(form, nd, family = f$model$family,
      control = glm.control(epsilon = 1e-11, maxit = 100))
    expect_lt(abs(a$value + as.numeric(logLik(native))), 1e-7)
    rebuilt <- .msr_profile_conditional_model(contract, a, f$model, p[1:2])
    expect_equal(model.matrix(rebuilt), model.matrix(native), tolerance = 1e-10, ignore_attr = TRUE)
    joint <- .msr_profile_joint(obj, p)
    numerical <- vapply(seq_along(joint$par), function(j) {
      hi <- lo <- joint$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
      (joint$fn(hi) - joint$fn(lo)) / 2e-5
    }, 0)
    expect_lt(max(abs(joint$gr(joint$par) - numerical)), 2e-5)
  }
})

test_that("prepared expressions preserve reference units, complete priors and AD targets", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (surface in c(FALSE, TRUE))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
      f <- prepared_formula_fixture(kernel, surface)
      if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
      form <- count ~ I(habitat^2) * year + I(habitat * other) + poly(temp, 2) + offset(log(effort))
      fit <- multiScale_bayes(form, f$inputs, f$data, f$scales, family,
        row_id = "key", method = "prepare", priors = list(scale_median = 250))
      spec <- fit$specification; p <- spec$parameters
      p$beta[] <- seq(-.2, .2, length.out = length(p$beta))
      p$scale_z[] <- seq(-.1, .2, length.out = length(p$scale_z))
      native <- prepared_formula_native_eta(spec, p, f$data)
      eta <- .msr_bayes_components(spec, p)$eta
      expect_equal(eta, native$eta, tolerance = 1e-12, ignore_attr = TRUE)
      expect_identical(spec$beta_names, colnames(native$X))
      expect_identical(spec$priors, .msr_bayes_priors(list(scale_median = 250),
        spec$covariates, f$inputs))
      likelihood <- switch(family,
        gaussian = dnorm(spec$y, eta, exp(p$log_dispersion), log = TRUE),
        poisson = dpois(spec$y, exp(eta), log = TRUE),
        binomial = dbinom(spec$y, spec$trials, plogis(eta), log = TRUE),
        nbinom2 = dnbinom(spec$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE))
      expect_equal(.msr_bayes_objective(spec)(p),
        -sum(likelihood) - prepared_formula_prior(spec, p), tolerance = 1e-10)
      obj <- RTMB::MakeADFun(.msr_bayes_objective(spec, TRUE), p, silent = TRUE)
      expect_equal(as.numeric(obj$fn(obj$par)), .msr_bayes_objective(spec)(p), tolerance = 1e-10)
      score <- vapply(seq_along(obj$par), function(j) {
        hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
        (obj$fn(hi) - obj$fn(lo)) / 2e-5
      }, 0)
      expect_lt(max(abs(as.numeric(obj$gr(obj$par)) - score)), 2e-5)
      if (surface) expect_error(.msr_scale_check_capability(spec,
        "marginal_likelihood", "factorized"), "Weighted SQ")
    }
})

test_that("expression admission rejects nonlocal functions, fixed context and masked bases", {
  expect_error(.msr_validate_kernel_expressions(terms(y ~ poly(h, 2)), "h"), "Unsupported")
  expect_error(.msr_validate_kernel_expressions(terms(y ~ I(h*z)), "h"), "kernel inputs")
  env <- new.env(parent = baseenv()); env$exp <- function(x) x
  form <- as.formula("y ~ exp(h)", env = env)
  expect_error(.msr_validate_kernel_expressions(terms(form), "h"), "original base")
})

test_that("expression scale and expow-shape scores preserve native kernel semantics", {
  for (kernel in c("exp", "expow")) {
    f <- formula_integration_fixture()
    f$inputs$kernel <- kernel
    model <- glm(y ~ I(hab^2) + hab * group + z + offset(log(exposure)),
      data = f$data, family = poisson())
    contract <- .msr_profile_contract(model, f$inputs); obj <- contract$objective
    p <- c(log(.3), if (kernel == "expow") log(1.8))
    kd <- .msr_eval_kernel_engine(obj$engine, p)
    nd <- f$data; nd$hab <- kd$value[, "hab"]
    native <- glm(formula(model), nd, family = poisson(),
      control = glm.control(epsilon = 1e-11, maxit = 100L))
    expect_lt(abs(obj$fn(p) + as.numeric(logLik(native))), 1e-7)
    score <- vapply(seq_along(p), function(j) {
      hi <- lo <- p; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
      (obj$fn(hi) - obj$fn(lo)) / 2e-5
    }, 0)
    expect_lt(max(abs(obj$gr(p) - score)), 2e-5)
  }
})


test_that("numeric names do not mask functions in dynamic expressions", {
  context <- new.env(parent = baseenv())
  context$exp <- c(-.3, 0, .2, .5)
  dat <- data.frame(y = 1:4, exp = context$exp)
  form <- y ~ exp(exp)
  environment(form) <- context
  tt <- terms(form)
  design <- .msr_compile_response_design(tt, model.frame(tt, dat), "exp")
  kernel <- as.matrix(dat["exp"])
  expect_equal(.msr_response_design_matrix(design, kernel),
    model.matrix(tt, model.frame(tt, dat)), ignore_attr = TRUE)
  expect_equal(as.numeric(.msr_response_kernel_sensitivity(design, kernel, c(.2, .3))),
    .3 * base::exp(dat$exp))
  context$exp <- function(x) x + 10
  expect_error(.msr_compile_response_design(tt, model.frame(tt, dat), "exp"),
    "original base function binding")
})

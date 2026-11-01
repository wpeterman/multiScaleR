test_that("native mixed expression likelihoods retain complete AD scores", {
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  f <- readRDS(test_path("fixtures", "joint-mixed-existing.rds"))
  dat <- f$data; dat$site <- factor(dat$site); dat$year <- factor(dat$year)
  engine <- .msr_prepare_kernel_engine(f$inputs$binned, "habitat")
  kp <- log(.23)
  dat$habitat <- .msr_eval_kernel_engine(engine, kp)$value[, 1L]
  for (family in c("poisson", "nbinom2")) {
    model <- glmmTMB::glmmTMB(count ~ habitat + I(habitat^2) + year + poly(temp, 2) +
      (1 | site) + offset(log(effort)), dat,
      family = if (family == "poisson") poisson() else glmmTMB::nbinom2())
    args <- .msr_native_mixed_design(model, engine, "habitat", kp, rownames(dat))
    expect_true(!is.null(args$response_design$expression_programs))
    obj <- do.call(.msr_mixed_rtmb, c(args, list(fixed_scale = TRUE)))
    expect_lt(abs(obj$fn(obj$par) + as.numeric(logLik(model))), 1e-5)
    at <- obj$par + .02
    numerical <- vapply(seq_along(at), function(j) {
      hi <- lo <- at; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
      (obj$fn(hi) - obj$fn(lo)) / 2e-5
    }, 0)
    expect_lt(max(abs(as.numeric(obj$gr(at)) - numerical)), 2e-5)
  }
})

test_that("public expression fits retain native prediction, covariance, curves and K", {
  f <- formula_integration_fixture()
  model <- glm(y ~ hab * group + I(hab^2) + z + offset(log(exposure)),
    f$data, family = poisson())
  fit <- multiScale_optim(model, f$inputs, backend = "profile_glm", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1L)))
  expect_equal(nrow(fit$profile_fit$start_results), 3L)
  expect_true(any(fit$profile_fit$start_results$eligible))
  expect_equal(unname(predict(fit)), unname(predict(fit$opt_mod, type = "response")), tolerance = 1e-8)
  actual <- predict(fit, type = "link", se.fit = TRUE, uncertainty = "conditional")
  expected <- predict(fit$opt_mod, type = "link", se.fit = TRUE)
  expect_equal(actual$fit, expected$fit, tolerance = 1e-8)
  expect_equal(actual$se.fit, expected$se.fit, tolerance = 1e-6)
  snapshot <- prepare_prediction(fit)
  expect_identical(predict(snapshot), predict(fit))
  expect_identical(predict(snapshot, type = "link", se.fit = TRUE, uncertainty = "conditional"), actual)
  order <- rev(seq_len(nrow(f$data)))
  expect_equal(unname(predict(fit, f$data[order, ], f$inputs)), rev(unname(predict(fit))), tolerance = 1e-10)
  expect_equal(.msr_selection_parameter_count(fit), as.numeric(attr(logLik(fit$opt_mod), "df")) + 1L)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  plots <- tryCatch(suppressMessages(plot_marginal_effects(fit, length.out = 5)),
    finally = grDevices::dev.off())
  expect_setequal(names(plots), c("hab", "group", "z"))
  expect_true(all(vapply(plots, function(x) all(is.finite(x$data$fit)), logical(1))))
})


test_that("conditional GLM metadata retains every native formula offset", {
  for (family in c("poisson", "gaussian", "binomial", "nbinom2")) {
    for (expression in c(FALSE, TRUE)) {
      f <- formula_integration_fixture(family)
      form <- if (expression) update.formula(formula(f$model), . ~ . + I(hab^2)) else formula(f$model)
      model <- glm(form, f$data, family = f$model$family)
      contract <- .msr_profile_contract(model, f$inputs)
      kp <- log(.3)
      fixed <- contract$objective$evaluate(c(kp, if (family == "nbinom2") log(2)))
      native <- .msr_profile_conditional_model(contract, fixed, model, kp)
      expect_identical(native$offset, contract$objective$offset)
      expect_equal(model.matrix(native), native$x, tolerance = 1e-12, ignore_attr = TRUE)
      expected <- drop(native$x %*% coef(native)) + native$offset
      expect_equal(unname(predict(native, type = "link", se.fit = TRUE)$fit),
        unname(expected), tolerance = 1e-10)
      expect_equal(unname(predict(native, type = "response", se.fit = TRUE)$fit),
        unname(native$family$linkinv(expected)), tolerance = 1e-10)
    }
  }
})

test_that("saved conditional logistic fits retain response and strata for refits", {
  skip_if_not_installed("survival")
  f <- make_core_fixture()
  n <- nrow(f$kernel_inputs$kernel_dat)
  dat <- data.frame(case_ = rep(c(TRUE, FALSE, FALSE), length.out = n),
                    cont1 = f$kernel_inputs$kernel_dat$cont1,
                    stratum = factor(rep(seq_len(n / 3), each = 3)))
  coxph <- survival::coxph
  strata <- survival::strata
  form <- stats::as.formula("case_ ~ cont1 + strata(stratum)",
                            env = asNamespace("survival"))
  fit <- survival::clogit(form, data = dat,
                         model = TRUE, x = TRUE, y = TRUE)
  # Use a standalone formula environment, avoiding pkgload's transient test
  # environment in the serialized object.
  environment(fit$terms) <- asNamespace("survival")
  environment(fit$formula) <- asNamespace("survival")
  environment(attr(fit$model, "terms")) <- asNamespace("survival")
  fit <- unserialize(serialize(fit, NULL))
  # Exercise reconstruction when the generic extractor returns the response
  # under a synthetic name, as happens for saved survival model frames.
  incomplete <- dat
  incomplete$case_ <- NULL
  local_mocked_bindings(get_data = function(...) incomplete)
  recovered <- .model_data(fit)
  expect_equal(as.numeric(recovered$case_), as.numeric(dat$case_))
  expect_equal(recovered$stratum, dat$stratum)
  ctx <- build_opt_context(fit, f$kernel_inputs$raw_cov)
  expect_true(all(names(dat) %in% names(ctx$data_template)))
  result <- kernel_scale_fn(40 / f$kernel_inputs$unit_conv,
                            f$kernel_inputs$d_list, f$kernel_inputs$raw_cov,
                            "gaussian", fit, mod_return = TRUE, opt_context = ctx)
  expect_s3_class(result$mod, "coxph")
  expect_true(is.finite(as.numeric(logLik(result$mod))))
  expect_identical(result$mod$method, "exact")
  expect_equal(result$mod$y, fit$y)
  expect_equal(model.frame(result$mod)[["strata(stratum)"]],
               model.frame(fit)[["strata(stratum)"]])
})

test_that("the final native Cox refit saves a supported nested basis", {
  skip_if_not_installed("survival")
  f <- make_core_fixture()
  n <- nrow(f$kernel_inputs$kernel_dat)
  dat <- data.frame(case_ = rep(c(TRUE, FALSE, FALSE), length.out = n),
    cont1 = f$kernel_inputs$kernel_dat$cont1,
    z = sin(seq_len(n) * .73),
    stratum = factor(rep(seq_len(n / 3), each = 3)))
  coxph <- survival::coxph
  strata <- survival::strata
  form <- case_ ~ cont1 + I(base::scale(z)) + strata(stratum)
  fit <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE,
    na.action = stats::na.omit)
  local_mocked_bindings(get_data = function(...) dat)
  ctx <- build_opt_context(fit, f$kernel_inputs$raw_cov)
  result <- kernel_scale_fn(40 / f$kernel_inputs$unit_conv,
    f$kernel_inputs$d_list, f$kernel_inputs$raw_cov, "gaussian", fit,
    mod_return = TRUE, opt_context = ctx)
  actual <- result$mod
  expect_false(.msr_needs_nested_metadata(terms(actual)))
  nd <- dat[c(2L, 5L), , drop = FALSE]
  nd$z <- c(1, 2)
  expected <- coef(actual)[1L] * nd$cont1 +
    coef(actual)[2L] * (nd$z - mean(dat$z)) / sd(dat$z)
  expect_equal(as.numeric(predict(actual, nd, reference = "zero")),
    expected, tolerance = 1e-10)
})

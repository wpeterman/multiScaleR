test_that("conditional weighted SQ preserves native likelihood and prediction contracts", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture(n = 120)
  b <- f$inputs$binned
  f$inputs$scale_vars <- msr_vars(hab = kernel_var("hab"),
                                rough = surface_var("hab", "sq", weighted = TRUE))
  f$inputs$d_list <- rep(list(b$bin_rep), nrow(f$data))
  f$inputs$raw_cov <- lapply(seq_len(nrow(f$data)), function(i)
    data.frame(hab = b$vsum$hab[i, ] / b$csum$hab[i, ]))
  e <- .msr_prepare_surface_engine(f$inputs, c("hab", "rough"), seq_len(nrow(f$data)))
  kd <- .msr_eval_kernel_engine(e, log(c(.25, .4)))
  dat <- f$data
  dat$rough <- kd$value[, "rough"]
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3), each = 3)
  form <- stats::as.formula("case_ ~ hab + rough + z + strata(choice)",
                            env = asNamespace("survival"))
  model <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE)
  c <- .msr_clogit_contract(model, f$inputs)
  o <- .msr_clogit_objective(c)
  p <- c(log(c(.25, .4)), rep(.1, o$nb))
  numerical <- vapply(seq_along(p), function(i) {
    a <- b <- p; a[i] <- a[i] + 1e-5; b[i] <- b[i] - 1e-5
    (o$fn(a) - o$fn(b)) / 2e-5
  }, 0)
  expect_equal(o$gr(p), numerical, tolerance = 1e-6, ignore_attr = TRUE)
  native <- .msr_clogit_native(c, log(c(.25, .4)))
  beta <- coef(native)[c(c$fixed_columns, c$covariates)]
  expect_equal(o$fn(c(log(c(.25, .4)), beta)), -as.numeric(logLik(native)), tolerance = 1e-10)
  fit <- multiScale_optim(model, f$inputs, backend = "conditional_logit", n_cores = 1,
                         profile_control = list(n_starts = 3), verbose = FALSE)
  expect_equal(nrow(fit$profile_fit$start_results), 3)
  stored <- fit$profile_fit$contract$objective$engine
  expect_null(stored$buffers)
  expect_null(stored$means)
  full <- fit
  full$profile_fit$contract$objective$engine <- .msr_prepare_surface_engine(
    fit$kernel_inputs, stored$covariates, stored$rows,
    stored$standardize, stored$implementation)
  expect_lt(as.numeric(object.size(fit)), as.numeric(object.size(full)))
  expect_equal(predict(fit, se.fit = TRUE, uncertainty = "conditional"),
               predict(full, se.fit = TRUE, uncertainty = "conditional"))
  expect_true(all(is.na(fit$scale_est$SE)))
  a <- predict(fit, se.fit = TRUE, uncertainty = "conditional")
  restored <- unserialize(serialize(fit, NULL))
  restored$profile_fit$contract$objective$engine$implementation <- "stream"
  expect_identical(restored$profile_fit$start_results, fit$profile_fit$start_results)
  expect_equal(predict(restored, se.fit = TRUE, uncertainty = "conditional"),
               a, tolerance = 1e-12)
  expect_equal(predict(restored, se.fit = TRUE, uncertainty = "conditional"),
               a, tolerance = 1e-12)
  b <- predict(fit, newdata = dat, kernel_inputs = f$inputs,
               se.fit = TRUE, uncertainty = "conditional")
  expect_equal(a, b, tolerance = 1e-10)
  prepared <- prepare_prediction(fit)
  expect_equal(predict(prepared, se.fit = TRUE, uncertainty = "conditional"), a)
  expect_equal(predict(prepare_prediction(fit, dat, f$inputs)), a$fit)
  expect_equal(predict(unserialize(serialize(prepared, NULL)), type = "link"),
               predict(fit, type = "link"))
  expect_null(prepared$kernel_inputs$raw_cov)
  expect_null(prepared$opt_mod)
  expect_lt(as.numeric(object.size(prepared)), as.numeric(object.size(fit)))
  expect_error(predict(prepared, se.fit = TRUE), "conditional SEs")
  expect_error(predict(prepared, interval = "confidence"), "conditional SEs")
  expect_error(predict(prepared, newdata = dat), "fixed inputs")
  expect_error(prepare_prediction(fit, max_bytes = 1), "max_bytes")
  changed_fit <- fit
  changed_fit$profile_fit$par[1] <- changed_fit$profile_fit$par[1] + .2
  expect_equal(predict(prepared), a$fit)
  expect_false(isTRUE(all.equal(predict(changed_fit), a$fit)))
  reversed <- predict(fit, newdata = dat[nrow(dat):1, ], kernel_inputs = f$inputs)
  expect_equal(unname(reversed), rev(unname(a$fit)), tolerance = 1e-10)
  expect_equal(predict(prepare_prediction(fit, dat[nrow(dat):1, ], f$inputs)), reversed)
  units <- f$inputs
  units$unit_conv <- units$unit_conv / 2
  units$d_list <- lapply(units$d_list, function(d) d * 2)
  units$binned$bin_rep <- units$binned$bin_rep * 2
  expect_equal(predict(fit, newdata = dat, kernel_inputs = units),
               predict(fit), tolerance = 1e-10)
  expect_equal(predict(prepare_prediction(fit, dat, units)), a$fit, tolerance = 1e-10)
  missing_raw <- f$inputs; missing_raw$raw_cov <- NULL
  expect_error(predict(fit, newdata = dat, kernel_inputs = missing_raw), "raw buffers")
  changed <- f$inputs; changed$scale_vars$source[2] <- "different"
  expect_error(predict(fit, newdata = dat, kernel_inputs = changed), "specifications")
  expect_error(prepare_prediction(fit, dat, changed), "specifications")
  expect_error(predict(fit, se.fit = TRUE), "conditional SEs")
  hard <- f$inputs; hard$scale_vars$weighted[2] <- FALSE
  expect_error(.msr_clogit_contract(model, hard), "hard-radius")
})

# print() and summary() regressions from the integration audit of 2026-11-24
# (findings C-1, X-6, C-7, G-5, U-6, X-7). Each case starts from a real fit.

audit_hab_inputs <- function() {
  data(pts, package = "multiScaleR", envir = environment())
  data(count_data, package = "multiScaleR", envir = environment())
  hab <- terra::rast(system.file("extdata", "hab.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                        kernel = "gaussian", verbose = FALSE)
  list(inputs = inputs, data = count_data)
}

test_that("a fit made through do.call() prints a compact call quickly", {
  a <- audit_hab_inputs()
  m <- glm(y ~ hab, family = poisson, data = a$data)
  fit <- do.call(multiScale_optim, list(fitted_mod = m, kernel_inputs = a$inputs,
                                        verbose = FALSE))
  # The stored call holds the evaluated model, inputs, and function.
  expect_true(is.function(fit$call[[1L]]))
  elapsed <- system.time(printed <- capture.output(print(fit)))[["elapsed"]]
  expect_lt(elapsed, 1)
  call_text <- paste(printed[seq_len(grep("Kernel used", printed)[1L])], collapse = " ")
  expect_match(call_text, "multiScale_optim(fitted_mod = <glm>", fixed = TRUE)
  expect_match(call_text, "kernel_inputs = <multiScaleR_data>", fixed = TRUE)
  expect_match(call_text, "verbose = FALSE", fixed = TRUE)
  printed <- capture.output(print(summary(fit)))
  expect_true(any(grepl("fitted_mod = <glm>", printed, fixed = TRUE)))

  # An ordinary call prints exactly as before.
  plain <- multiScale_optim(m, a$inputs, verbose = FALSE)
  expect_identical(capture.output(.msr_print_call(plain$call)),
                   capture.output(print(plain$call)))
})

test_that("a Hessian touched by a failed refit is explained in print and summary", {
  a <- audit_hab_inputs()
  fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = a$data),
                          a$inputs, verbose = FALSE)
  # Record the failed-refit outcome that multiScale_optim() stores; see
  # .msr_refit_hessian_se(). A real failure needs a refit that breaks near
  # the optimum, which no small public fixture produces reliably.
  fit$scale_est$SE <- NA_real_
  fit$diagnostics$refit_hessian$status <- "hessian_includes_failed_refit"
  fit$diagnostics$refit_hessian$curvature_status <- "hessian_includes_failed_refit"
  for (printed in list(capture.output(print(fit)),
                       capture.output(print(summary(fit))))) {
    at <- grep("refit failed while the Hessian", printed, fixed = TRUE)
    expect_length(at, 1L)
    # The note sits with the scale table, before the fitted model.
    expect_lt(at, grep("Fitted Model", printed, fixed = TRUE)[1L])
    expect_true(any(grepl("summary(x, profile = TRUE)", printed, fixed = TRUE)))
  }
})

test_that("nlme::lme fits are summarized, plotted, and converted to distances", {
  skip_if_not_installed("nlme")
  a <- audit_hab_inputs()
  dat <- a$data
  dat$grp <- factor(rep(letters[1:8], length.out = nrow(dat)))
  m <- nlme::lme(y ~ hab, random = ~ 1 | grp, data = dat)
  fit <- multiScale_optim(m, a$inputs, verbose = FALSE)
  expect_s3_class(fit$opt_mod, "lme")
  # insight returns one residual df per fixed term for lme; the scale limits
  # use N minus the number of fixed coefficients, as for gls.
  expect_identical(.msr_residual_df(fit$opt_mod),
                   nrow(dat) - length(nlme::fixef(fit$opt_mod)))
  s <- summary(fit)
  se <- fit$scale_est["hab", "SE"]
  q <- stats::qt(.975, nrow(dat) - 2L)
  expect_equal(unname(unlist(s$opt_scale["hab", c("2.5%", "97.5%")])),
               fit$scale_est["hab", "Mean"] + c(-1, 1) * q * se, tolerance = 1e-10)
  expect_true(all(is.finite(unlist(kernel_dist(fit)))))
  expect_no_error(capture.output(print(s)))
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit))
})

test_that("summary() and plot() unwrap amt fits under the default refit backend", {
  skip_on_cran()
  skip_if_not_installed("amt")
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; strata <- survival::strata
  wrapped <- amt::fit_clogit(fx$df, case_ ~ forest + sl_km + strata(step_id),
                             model = TRUE)
  fit <- multiScale_optim(wrapped, fx$inputs, n_cores = 1, verbose = FALSE)
  expect_s3_class(fit$opt_mod, "fit_clogit")
  s <- summary(fit)
  expect_s3_class(s$fitted_mod, "coxph")
  expect_true(is.finite(s$opt_scale["forest", "SE"]))
  expect_no_error(capture.output(print(s)))
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fit))
  # Generics read the native model inside the wrapper.
  expect_equal(coef(fit), coef(fit$opt_mod$model))
  expect_equal(AIC(fit), aic_tab(list(fit), AICc = FALSE)$AIC, tolerance = 1e-10)
})

test_that("conditional_logit print and summary point to profile_sigma()", {
  skip_on_cran()
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; strata <- survival::strata
  m <- survival::clogit(case_ ~ forest + sl_km + strata(step_id),
                        data = fx$df, model = TRUE)
  fit <- multiScale_optim(m, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE,
                          profile_control = list(starts = matrix(0.3, ncol = 1)))
  for (printed in list(capture.output(print(fit)),
                       capture.output(print(summary(fit))))) {
    at <- grep("profile_sigma(x)", printed, fixed = TRUE)
    expect_length(at, 1L)
    # Directly under the scale table, before the distances.
    expect_lt(at, grep("Scale of Effect -- Distance", printed, fixed = TRUE)[1L])
  }
})

test_that("profile_glm summaries name withheld scale SEs under the table", {
  f <- profile_glm_fixture()
  # The true scale is 250; a 120 search limit puts sigma on the upper bound.
  f$inputs$max_D <- 120
  expect_warning(fit <- profile_glm_fit(f), "optimization bound")
  expect_true(is.na(fit$scale_est["hab", "SE"]))
  for (printed in list(capture.output(print(summary(fit))),
                       capture.output(print(fit)))) {
    at <- grep("Scale SEs and limits are withheld for sigma:hab", printed, fixed = TRUE)
    expect_length(at, 1L)
    expect_match(printed[at], "parameter at an optimization bound", fixed = TRUE)
    expect_lt(at, grep("Scale of Effect -- Distance", printed, fixed = TRUE)[1L])
  }
})

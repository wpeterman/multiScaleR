test_that("joint mixed models retain native predictions and independent starts", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("TMB")
  skip_if_not_installed("glmmTMB")
  # Frozen existing annual example inputs, shared-scale frequentist model.
  # No simulation or Bayesian sampling occurs in this test.
  f <- readRDS(test_path("fixtures", "joint-mixed-existing.rds"))
  d <- f$data; d$site <- factor(d$site); d$year <- factor(d$year)
  e <- multiScaleR:::.msr_prepare_kernel_engine(f$inputs$binned, "habitat")
  d$habitat <- multiScaleR:::.msr_eval_kernel_engine(e, log(.23))$value[, 1]
  for (zi in c(FALSE, TRUE)) {
    m <- glmmTMB::glmmTMB(count ~ temp + year + habitat + (1 | site) + offset(log(effort)),
      data = d, family = if (zi) poisson() else glmmTMB::nbinom2(),
      ziformula = if (zi) ~1 else ~0,
      contrasts = if (zi) list(year = "contr.sum") else NULL)
    fit <- multiScale_optim(m, f$inputs, backend = "joint_mixed", n_cores = 1,
      profile_control = list(starts = matrix(c(.1, .23, .5), ncol = 1)), verbose = FALSE)
    expect_s3_class(fit, "multiScaleR")
    expect_s3_class(fit$opt_mod, "glmmTMB")
    expect_s3_class(summary(fit), "summary_multiScaleR")
    expect_equal(fit$diagnostics$joint_mixed$eligible_starts, 3L)
    expect_true(all(is.na(fit$scale_est$SE)))
    expect_true(is.finite(kernel_dist(fit)[1, 1]))
    expect_error(summary(fit, profile = TRUE), "not validated")
    expect_error(profile_sigma(fit), "not validated")
    expect_error(vcov(fit), "unavailable")
    expect_equal(vcov(fit, type = "conditional"), vcov(fit$opt_mod))
    expect_equal(predict(fit), predict(fit$opt_mod, type = "response"))
    expect_equal(predict(fit, se.fit = TRUE, uncertainty = "conditional"),
                 predict(fit$opt_mod, type = "response", se.fit = TRUE))
    expect_error(predict(fit, se.fit = TRUE), "conditional")
    expect_error(predict(fit, interval = "confidence"), "conditional")
    revdata <- d[nrow(d):1, ]; revdata$habitat <- 999
    expect_equal(unname(predict(fit, newdata = revdata, kernel_inputs = f$inputs)),
                 rev(unname(predict(fit))), tolerance = 1e-8)
    expect_equal(predict(fit, re.form = NA), predict(fit$opt_mod, type = "response", re.form = NA))
    prepared <- prepare_prediction(fit)
    expect_equal(predict(prepared), predict(fit))
    expect_equal(predict(prepared, type = "link", se.fit = TRUE, uncertainty = "conditional"),
                 predict(fit, type = "link", se.fit = TRUE, uncertainty = "conditional"))
    expect_equal(predict(unserialize(serialize(prepared, NULL))), predict(fit))
    expect_null(prepared$opt_mod)
    expect_null(prepared$kernel_inputs)
    expect_error(predict(prepared, se.fit = TRUE), "conditional")
    expect_error(predict(prepared, interval = "confidence"), "conditional")
    expect_error(predict(prepared, re.form = NA), "fixed inputs")
    expect_error(prepare_prediction(fit, max_bytes = 1), "max_bytes")
    zero_re <- prepare_prediction(fit, native_args = list(re.form = NA))
    expect_equal(predict(zero_re), predict(fit, re.form = NA))
    novel <- revdata
    novel$site <- factor(rep("new_site", nrow(novel)))
    new_groups <- prepare_prediction(fit, novel, f$inputs,
      native_args = list(allow.new.levels = TRUE))
    expect_equal(predict(new_groups), predict(fit, novel, f$inputs, allow.new.levels = TRUE))
    capture_native <- function(expr) {
      messages <- character()
      value <- withCallingHandlers(expr, warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
      list(value = value, warnings = messages)
    }
    ordinary <- capture_native(predict(fit, novel, f$inputs))
    warning_snapshot <- capture_native(prepare_prediction(fit, novel, f$inputs))$value
    expect_equal(capture_native(predict(warning_snapshot)), ordinary)
    reordered <- prepare_prediction(fit, revdata, f$inputs)
    expect_equal(predict(reordered), predict(fit, revdata, f$inputs))
    changed_units <- f$inputs
    changed_units$unit_conv <- changed_units$unit_conv / 2
    changed_units$binned$bin_rep <- changed_units$binned$bin_rep * 2
    expect_equal(predict(prepare_prediction(fit, revdata, changed_units)), predict(reordered))
    bad_new <- revdata; bad_new$effort <- NULL
    expect_error(prepare_prediction(fit, bad_new, f$inputs), "explicitly")
    expect_error(prepare_prediction(fit, native_args = list(newparams = list())), "supported")
    if (zi) {
      probability <- prepare_prediction(fit, component = "zprob")
      expect_equal(predict(probability), predict(fit$opt_mod, type = "zprob"))
      expect_equal(predict(probability, type = "zprob", se.fit = TRUE, uncertainty = "conditional"),
                   predict(fit$opt_mod, type = "zprob", se.fit = TRUE))
      expect_error(predict(probability, type = "response"), "not prepared")
    }
    expect_error(predict(fit, newdata = d), "together")
    bad <- d; bad$effort <- NULL
    expect_error(predict(fit, newdata = bad, kernel_inputs = f$inputs), "explicitly")
    expect_error(multiScale_optim(m, f$inputs, backend = "joint_mixed",
      profile_control = list(coefficient_start = "warm")), "only starts")
    expect_error(multiScale_optim(m, f$inputs, backend = "joint_mixed", n_cores = 2), "serial")
    bad_model <- m; bad_model$modelInfo$map <- list(theta = factor(NA))
    expect_error(multiScale_optim(bad_model, f$inputs, backend = "joint_mixed"), "maps")
    bad_model <- m; bad_model$modelInfo$reTrms$cond$cnms[[1]] <- c("(Intercept)", "temp")
    expect_error(multiScale_optim(bad_model, f$inputs, backend = "joint_mixed"), "random intercept")
    bad_model <- m; bad_model$sdr$pdHess <- FALSE
    expect_error(multiScale_optim(bad_model, f$inputs, backend = "joint_mixed"), "not eligible")
    bad_inputs <- f$inputs; bad_inputs$binned$point_ids <- rev(bad_inputs$binned$point_ids)
    # Valid identifiers alone do not certify the intended spatial mapping.
    # Duplicate identifiers must fail before fitting rather than silently join.
    bad_inputs$binned$point_ids[1] <- bad_inputs$binned$point_ids[2]
    expect_error(multiScale_optim(m, bad_inputs, backend = "joint_mixed"), "unique")
  }
})

test_that("joint mixed failures preserve every attempted start", {
  expect_error(profile_sigma(1), "fitted.*multiScaleR")
  z <- multiScaleR:::.msr_native_joint(list(p = 1), "habitat", character(),
    matrix(log(c(.1, .23, .5)), ncol = 1), log(.04), log(1),
    function(kp) stop("intentional native failure"))
  expect_length(z$outcomes, 3L)
  expect_true(is.na(z$selected))
  expect_null(z$model)
  expect_true(all(vapply(z$outcomes, function(x) identical(x$error, "intentional native failure"), TRUE)))
})

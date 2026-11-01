test_that("actual amt iSSF objects retain movement corrections and conditional predictions", {
  skip_if_not_installed("amt")
  # Re-express an existing numerical fixture without adding a research panel.
  f <- profile_glm_fixture()
  dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$step_id_ <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  dat$log_sl_ <- dat$z
  dat$sl_ <- exp(dat$log_sl_)
  dat$cos_ta_ <- .8 * (dat$group == "b") - .4
  dat$ta_ <- acos(dat$cos_ta_)
  attr(dat, "sl_") <- amt::make_gamma_distr(shape = 2, scale = 1)
  attr(dat, "ta_") <- amt::make_vonmises_distr(kappa = 2)
  form <- as.formula(paste("case_ ~ hab + log_sl_ + cos_ta_ +",
    "offset(log(exposure)) + strata(step_id_)"), env = asNamespace("survival"))
  wrapped <- amt::fit_issf(dat, form, model = TRUE, x = TRUE, y = TRUE,
    more = list(marker = "original proposal"))
  original <- serialize(wrapped, NULL)
  starts <- matrix(c(.08, .25, .6), ncol = 1)
  fit <- withCallingHandlers(multiScale_optim(wrapped, f$inputs, backend = "conditional_logit",
    n_cores = 1, verbose = FALSE, profile_control = list(starts = starts)),
    msr_clogit_fit_error = function(e) print(e$start_results))
  expect_identical(serialize(wrapped, NULL), original)
  expect_identical(fit$opt_mod$sl_, wrapped$sl_)
  expect_identical(fit$opt_mod$ta_, wrapped$ta_)
  expect_identical(fit$opt_mod$more, wrapped$more)
  expect_identical(fit$profile_fit$contract$native_data_mode, "retained_model_frame")
  expect_equal(exp(fit$profile_fit$starts), starts, ignore_attr = TRUE)
  expect_true(any(fit$profile_fit$start_results$eligible))
  expect_equal(coef(fit$opt_mod), coef(fit$opt_mod$model))
  expect_equal(AIC(fit$opt_mod), AIC(fit$opt_mod$model))
  expect_true(all(is.na(fit$scale_est$SE)))
  expect_false(fit$uncertainty$available)

  # amt corrects the original proposal once, using the selected coefficients.
  reference <- wrapped; reference$model <- fit$opt_mod$model
  step_correction <- function(object) {
    warnings <- character()
    value <- withCallingHandlers(amt::update_sl_distr(object), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    list(value = value, warnings = warnings)
  }
  sl <- step_correction(fit$opt_mod)
  expect_equal(sl, step_correction(reference))
  # This valid log-step-only model leaves the proposal's scale unchanged.
  expect_length(sl$warnings, 1L)
  expect_match(sl$warnings, "sl_.*did not appear")
  expect_equal(amt::update_ta_distr(fit$opt_mod), amt::update_ta_distr(reference))
  expect_identical(fit$opt_mod$sl_, wrapped$sl_)

  probability <- predict(fit)
  expect_equal(as.numeric(rowsum(probability, dat$step_id_)), rep(1, 80), tolerance = 1e-12)
  expect_equal(predict(unserialize(serialize(fit, NULL))), probability)
  nd <- dat
  nd$hab <- .msr_eval_kernel_engine(fit$profile_fit$contract$objective$engine,
    fit$profile_fit$par)$value[, "hab"]
  x1 <- nd[c(1, 2), , drop = FALSE]; x2 <- nd[3, , drop = FALSE]
  rss <- amt::log_rss(fit$opt_mod, x1 = x1, x2 = x2, ci = "se")
  eta1 <- predict(fit$opt_mod$model, x1, type = "lp", reference = "zero")
  eta2 <- predict(fit$opt_mod$model, x2, type = "lp", reference = "zero")
  expect_equal(unname(rss$df$log_rss), unname(eta1 - eta2[1]), tolerance = 1e-12)
  # Native intervals cover coefficients conditionally on the selected kernel.
  X1 <- model.matrix(fit$opt_mod$model, data = x1)
  X2 <- model.matrix(fit$opt_mod$model, data = x2)
  delta <- sweep(X1, 2, X2[1, ])
  se <- sqrt(rowSums((delta %*% vcov(fit$opt_mod$model)) * delta))
  expect_equal(rss$df$upr - rss$df$log_rss, unname(qnorm(.975) * se), tolerance = 1e-12)
  expect_equal(rss$df$log_rss - rss$df$lwr, unname(qnorm(.975) * se), tolerance = 1e-12)
  expect_equal(.msr_selection_parameter_count(fit), length(coef(fit$opt_mod$model)) + 1L)
  expect_equal(bic_tab(list(fit), mod_names = "iSSF")$BIC,
    -2 * as.numeric(logLik(fit$opt_mod$model)) +
      log(attr(logLik(fit$opt_mod$model), "nobs")) * .msr_selection_parameter_count(fit),
    tolerance = 1e-12)
})

test_that("saved Cox frames preserve fixed transformations and compound offsets", {
  skip_if_not_installed("amt")
  f <- profile_glm_fixture()
  dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$step_id_ <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  dat$duration <- rep(c(1, 1.2, 1.6), length.out = nrow(dat))
  form <- as.formula(paste("case_ ~ hab + poly(z, 2) + group +",
    "offset(log(exposure * duration)) + strata(step_id_)"), env = asNamespace("survival"))
  m <- amt::fit_clogit(dat, form, model = TRUE, x = TRUE, y = TRUE)$model
  contract <- .msr_clogit_contract(m, f$inputs)
  expect_identical(contract$native_data_mode, "retained_model_frame")
  native <- .msr_clogit_native(contract, log(.3))
  dat$hab <- .msr_eval_kernel_engine(contract$engine, log(.3))$value[, "hab"]
  expected <- survival::coxph(formula(m), dat, method = m$method,
    model = TRUE, x = TRUE, y = TRUE,
    control = survival::coxph.control(iter.max = 100L, eps = 1e-10))
  expect_equal(coef(native), coef(expected), tolerance = 1e-12)
  expect_equal(logLik(native), logLik(expected), tolerance = 1e-12)
  expect_equal(vcov(native), vcov(expected), tolerance = 1e-12)
  expect_equal(model.matrix(native), model.matrix(expected), tolerance = 1e-12)
  # Predict outside the training polynomial values using its original basis.
  nd <- dat; nd$z <- nd$z + .35
  expect_equal(predict(native, nd, type = "lp", reference = "zero"),
    predict(expected, nd, type = "lp", reference = "zero"), tolerance = 1e-12)
  design <- model.matrix(contract$terms, nd, contrasts.arg = contract$contrasts)
  expect_equal(design[, contract$original_columns, drop = FALSE],
    model.matrix(expected, data = nd), tolerance = 1e-12, ignore_attr = TRUE)
  broken <- m; rownames(broken$model)[1] <- "unknown row"
  expect_error(.msr_clogit_contract(broken, f$inputs), "row names")
  broken <- m; attr(broken$terms, "predvars") <- NULL
  expect_error(.msr_clogit_contract(broken, f$inputs), "saved prediction transformations")
})

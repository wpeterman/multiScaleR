# Deterministic software comparisons, not an ecological recovery study.
test_that("compiled choice normalization preserves stable probabilities and RNG", {
  reference <- function(eta, groups) {
    p <- lp <- numeric(length(eta))
    for (i in groups) {
      shift <- max(eta[i]); w <- exp(eta[i] - shift)
      p[i] <- w / sum(w); lp[i] <- eta[i] - shift - log(sum(w))
    }
    list(probability = p, log_probability = lp)
  }
  groups <- list(c(5L, 2L, 9L), c(7L, 1L, 3L, 4L), c(8L, 6L))
  for (eta in list(rep(0, 9), c(-1000, 1001, 0, -1001, 1000, 999, 2, -999, 1),
                   cos(seq_len(9)), c(-1e308, 1e308, rep(0, 7)))) {
    got <- .msr_clogit_softmax(eta, groups)
    expect_equal(got, reference(eta, groups), tolerance = 1e-14)
    expect_equal(vapply(groups, function(i) sum(got$probability[i]), 0), rep(1, 3),
      tolerance = 1e-14)
  }
  expect_identical(.msr_clogit_softmax(numeric(), list()), reference(numeric(), list()))
  expect_equal(.msr_clogit_softmax(1:3, list(c(3, 1), 2))$probability[2], 1)
  set.seed(991L); before <- .Random.seed
  .msr_clogit_softmax(c(0, 1), list(1:2))
  expect_identical(.Random.seed, before)
  expect_error(.msr_clogit_softmax(c(0, Inf), list(1:2)), "finite")
  expect_error(.msr_clogit_softmax(c(0, 1), list(c(1, 1), 2)), "overlap")
  expect_error(.msr_clogit_softmax(c(0, 1), list(1)), "every")
  expect_error(.msr_clogit_softmax(c(0, 1), list(c(1, 3))), "invalid")
  expect_error(.msr_clogit_softmax(c(0, 1), list(c(1, 1.5))), "invalid")
  expect_error(.msr_clogit_softmax(c(0, 1), list(integer(), 1:2)), "empty")
})

ssf_projection_fixture <- function(weighted) {
  f <- if (weighted) formula_weighted_fixture(surface_only = TRUE) else profile_glm_fixture()
  dat <- f$data
  if (weighted) {
    nm <- "rough"
  } else {
    # A colon inside a quoted name is a variable name, not a cross term.
    nm <- "hab:itat"
    names(f$inputs$binned$vsum) <- names(f$inputs$binned$csum) <- nm
    names(dat)[names(dat) == "hab"] <- nm
  }
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- as.formula(paste0("case_ ~ `", nm,
    "` * group + poly(z,2) + offset(log(exposure)) + strata(choice)"),
    env = asNamespace("survival"))
  model <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE)
  list(model = model, data = dat, inputs = f$inputs, covariate = nm)
}

test_that("step-selection link and relative-score projections do not require sets", {
  skip_if_not_installed("survival")
  for (weighted in c(FALSE, TRUE)) {
    f <- ssf_projection_fixture(weighted)
    fit <- multiScale_optim(f$model, f$inputs, backend = "conditional_logit",
      n_cores = 1, verbose = FALSE,
      profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1L)))
    expect_equal(nrow(fit$profile_fit$start_results), 3L)
    expect_true(any(fit$profile_fit$start_results$eligible))
    native_data <- f$data
    engine <- if (weighted) .msr_prepare_surface_engine(f$inputs, f$covariate,
      seq_len(nrow(f$data))) else .msr_prepare_kernel_engine(f$inputs$binned,
        f$covariate, kernel = f$inputs$kernel)
    native_data[[f$covariate]] <- .msr_eval_kernel_engine(engine,
      fit$profile_fit$par[seq_along(engine$parameter_names)])$value[, 1L]
    native <- .analysis_model(fit$opt_mod)
    lp <- predict(native, native_data, type = "lp", reference = "zero", se.fit = TRUE)
    risk <- predict(native, native_data, type = "risk", reference = "zero")
    # Survival centers formula offsets even with reference="zero". The
    # conditional backend preserves its unshifted formula-offset contract.
    offset_mean <- mean(model.offset(model.frame(native)))
    expected_link <- lp$fit + offset_mean
    expected_risk <- exp(expected_link)
    nd <- f$data[, setdiff(names(f$data), c("choice", "case_", "y")), drop = FALSE]
    link <- predict(fit, nd, f$inputs, type = "link", se.fit = TRUE, uncertainty = "conditional")
    score <- predict(fit, nd, f$inputs, type = "risk", se.fit = TRUE, uncertainty = "conditional")
    expect_equal(unname(link$fit), unname(expected_link), tolerance = 1e-10)
    expect_equal(unname(link$se.fit), unname(lp$se.fit), tolerance = 1e-10)
    expect_equal(unname(score$fit), unname(risk * exp(offset_mean)), tolerance = 1e-10)
    expect_equal(unname(score$se.fit), unname(expected_risk * lp$se.fit), tolerance = 1e-10)
    expect_error(predict(fit, nd, f$inputs), "choice-set IDs")
    expect_error(predict(fit, nd, f$inputs, type = "risk", se.fit = TRUE), "conditional")
    snap <- prepare_prediction(fit, nd, f$inputs)
    expect_identical(predict(snap, type = "link", se.fit = TRUE, uncertainty = "conditional"), link)
    expect_identical(predict(snap, type = "risk", se.fit = TRUE, uncertainty = "conditional"), score)
    expect_error(predict(snap), "choice-set IDs")
    expect_false(snap$prediction_contract$response$available)
    expect_identical(snap$prediction_contract$risk$component, "relative_selection_score")
    expect_false(snap$prediction_contract$risk$uncertainty$joint_available)
    expect_identical(snap$prediction_contract$risk$uncertainty$risk_standard_errors,
      "delta_method_risk_times_link_SE")
    expect_false("supplied_choice_sets" %in% snap$prediction_contract$link$conditioning)
    bad <- nd; bad$choice <- seq_len(nrow(bad))
    expect_equal(predict(fit, bad, f$inputs, type = "risk"), score$fit)
    bad_snap <- prepare_prediction(fit, bad, f$inputs)
    expect_identical(predict(bad_snap, type = "link"), link$fit)
    expect_error(predict(bad_snap), "two alternatives")
    expect_error(predict(fit, bad, f$inputs), "two alternatives")
    missing <- nd; missing$z <- NULL
    expect_error(prepare_prediction(fit, missing, f$inputs), "missing formula")
    complete <- prepare_prediction(fit, f$data, f$inputs)
    expect_equal(predict(complete), predict(fit), tolerance = 1e-12)
    expect_equal(predict(fit, type = "r"), predict(fit), tolerance = 1e-12)
    expect_equal(.msr_selection_parameter_count(fit), length(coef(native)) + 1L)
  }
})

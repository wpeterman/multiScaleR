# Observe allocation lifetime without measuring or claiming peak RAM savings.
curve_reference_observe_basis_lifetime <- function(object, model, reference) {
  predictor <- .msr_predict_gam; observed <- list()
  watch <- function(object, ...) {
    args <- list(...)
    if (identical(args$type, "link")) {
      caller <- parent.frame()
      observed[[length(observed) + 1L]] <<- vapply(c("fitted_design", "reference_design"),
        exists, logical(1), envir = caller, inherits = FALSE)
    }
    do.call(predictor, c(list(object = object), args))
  }
  result <- with_mocked_bindings(.msr_curve_reference_context(object, model, reference),
    .msr_predict_gam = watch)
  expect_length(observed, 1L)
  expect_identical(unname(observed[[1L]]), c(FALSE, FALSE))
  expect_identical(.msr_predict_gam, predictor)
  result
}

# Supplied raw inputs establish fitted feature equivalence, not an inverse of
# evaluated transformations. All targets below use existing software fixtures.
test_that("expired GLM curves verify raw references with quotes, bases and offsets", {
  f <- refit_spatial_fixture(); dat <- f$data
  names(dat)[names(dat) == "z"] <- "z value|raw"
  dat$a <- seq_len(nrow(dat)) %% 2L == 0L
  dat$b <- seq_len(nrow(dat)) %% 5L == 0L
  dat$unused <- seq_len(nrow(dat))
  bare <- glm(y ~ mean + rough + poly(`z value|raw`, 2) * group +
    I(fixed^2) + I(a | b) + offset(log(exposure) + .1 * a),
    data = dat, family = poisson())
  model <- glm(y ~ mean + rough + stats::poly(`z value|raw`, 2) * group +
    I(fixed^2) + I(a | b) + offset(log(exposure) + .1 * a),
    data = dat, family = poisson())
  model$call$data <- quote(expired_raw_curve_data)
  environment(model$formula) <- asNamespace("stats")
  attr(model$terms, ".Environment") <- asNamespace("stats")
  attr(attr(model$model, "terms"), ".Environment") <- asNamespace("stats")
  object <- f$object; object$opt_mod <- model; original <- model
  reference <- dat; row.names(reference) <- row.names(model$model)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  expect_error(suppressWarnings(plot_marginal_effects(object, length.out = 7)),
    "Raw inputs|recover|raw")
  plots <- plot_marginal_effects(object, length.out = 7,
    reference_data = reference[rev(seq_len(nrow(reference))), , drop = FALSE])
  expect_identical(object$opt_mod, original)
  expect_true(all(c("z value|raw", "fixed", "group", "a", "b") %in% names(plots)))
  expect_false(any(c("exposure", "unused", "site") %in% names(plots)))
  vars <- c("mean", "rough", "z value|raw", "fixed", "group", "a", "b", "exposure")
  defaults <- lapply(dat[vars], function(v) if (is.numeric(v)) mean(v) else v[1L])
  nd <- as.data.frame(lapply(defaults, rep, times = 7), check.names = FALSE)
  nd[["z value|raw"]] <- seq(min(dat[["z value|raw"]]), max(dat[["z value|raw"]]), length.out = 7)
  expected <- predict(bare, nd, type = "link", se.fit = TRUE)
  expect_equal(plots[["z value|raw"]]$data$x, nd[["z value|raw"]], tolerance = 1e-10)
  expect_equal(plots[["z value|raw"]]$data$fit, unname(exp(expected$fit)), tolerance = 1e-8)
  expect_equal(plots[["z value|raw"]]$data$lwr,
    unname(exp(expected$fit + qnorm(.025) * expected$se.fit)), tolerance = 1e-8)
  expect_equal(plots[["z value|raw"]]$data$upr,
    unname(exp(expected$fit + qnorm(.975) * expected$se.fit)), tolerance = 1e-8)
  info <- attr(plots, "reference_data_contract")
  expect_true(info$fitted_feature_equivalence)
  expect_true(info$native_design_and_offsets_verified)
  expect_false(info$original_raw_values_uniquely_established)
  expect_false(info$metadata_establishes_eligibility)
  expect_identical(info$alignment, "explicit_fitted_row_IDs")
  one <- plot_marginal_effects(object, length.out = 1, reference_data = reference)
  point <- nd[1L, , drop = FALSE]
  point[["z value|raw"]] <- min(dat[["z value|raw"]])
  expect_equal(one[["z value|raw"]]$data$fit,
    unname(predict(bare, point, type = "response")), tolerance = 1e-8)
})

test_that("raw curve references reject row, feature, contrast and offset changes", {
  f <- refit_spatial_fixture(); dat <- f$data
  model <- lm(y ~ mean + rough + poly(z, 2) * group + I(fixed^2) +
    offset(log(exposure)), data = dat)
  object <- f$object; object$opt_mod <- model
  reference <- dat; row.names(reference) <- row.names(model$model)
  context <- .msr_curve_reference_context(object, model, reference)
  expect_identical(row.names(context$data), row.names(model$model))
  expect_identical(context$model$coefficients, model$coefficients)
  automatic <- as.data.frame(lapply(reference, identity))
  expect_error(.msr_curve_reference_context(object, model, automatic), "explicit fitted row IDs")
  wrong <- reference[-1L, , drop = FALSE]
  expect_error(.msr_curve_reference_context(object, model, wrong), "row IDs")
  wrong <- reference; row.names(wrong)[1L] <- "unknown"
  expect_error(.msr_curve_reference_context(object, model, wrong), "row IDs")
  wrong <- reference; attr(wrong, "row.names") <- rep("duplicate", nrow(wrong))
  expect_error(.msr_curve_reference_context(object, model, wrong), "unique explicit")
  wrong <- reference; wrong$z <- NULL
  expect_error(.msr_curve_reference_context(object, model, wrong), "missing raw inputs")
  for (variable in c("z", "fixed", "exposure")) {
    wrong <- reference; wrong[[variable]][1L] <- wrong[[variable]][1L] + .1
    expect_error(.msr_curve_reference_context(object, model, wrong), "fitted feature|offset")
  }
  for (bad in c(NA_real_, Inf, NaN)) {
    wrong <- reference; wrong$z[1L] <- bad
    expect_error(.msr_curve_reference_context(object, model, wrong), "finite scalar")
  }
  wrong <- reference; wrong$group <- factor("new", levels = c(levels(dat$group), "new"))
  expect_error(.msr_curve_reference_context(object, model, wrong), "unknown factor levels")
  wrong <- reference; contrasts(wrong$group) <- contr.sum(nlevels(wrong$group))
  expect_error(.msr_curve_reference_context(object, model, wrong), "contrasts differ")
  external <- lm(y ~ mean + rough + poly(z, 2), data = dat, offset = log(exposure))
  expect_error(.msr_curve_reference_context(object, external, reference), "outside the model formula")
  missing_frame <- model; missing_frame$model <- NULL
  expect_error(.msr_curve_reference_context(object, missing_frame, reference), "retained evaluated")
  nested <- lm(y ~ mean + I(scale(z)^2), data = dat)
  expect_error(.msr_curve_reference_context(object, nested, reference), "nested learned")
  explicit <- lm(y ~ mean + I(scale(z, center = 0, scale = 2)^2), data = dat)
  expect_true(.msr_curve_reference_context(object, explicit, reference)$contract$fitted_feature_equivalence)
  expect_error(plot_marginal_effects(object, reference_data = reference, component = 1), "only to multi-component")
  bayes <- structure(list(), class = "multiScaleR_bayes")
  expect_error(plot_marginal_effects(bayes, reference_data = reference), "do not redefine Bayesian")
  unsupported <- object; unsupported$opt_mod <- structure(list(), class = "glmmTMB")
  expect_error(plot_marginal_effects(unsupported, reference_data = reference), "not supported")
})

test_that("raw references preserve dot expansion, removed terms and nonunique inverses", {
  dat <- data.frame(y = c(1, 3, 2, 7, 4, 6, 8, 5, 9, 4, 6, 2),
    z = c(-2, -1.7, -1.2, -.8, -.4, 0, .3, .7, 1, 1.5, 2.2, 3),
    group = factor(rep(c("a", "b"), 6)), discard = seq_len(12))
  model <- lm(y ~ . - z - discard + poly(z, 2), data = dat)
  object <- structure(list(opt_mod = model, scl_params = list(mean = numeric(), sd = numeric())),
    class = "multiScaleR")
  reference <- dat; row.names(reference) <- row.names(model$model)
  reference$extra <- 20
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(object, reference_data = reference, length.out = 7)
  expect_setequal(names(plots), c("group", "z"))
  expected <- predict(model, data.frame(z = plots$z$data$x, group = dat$group[1L],
    discard = mean(dat$discard)), se.fit = TRUE)
  expect_equal(plots$z$data$fit, unname(expected$fit), tolerance = 1e-10)
  expect_equal(plots$z$data$lwr, unname(expected$fit + qnorm(.025) * expected$se.fit), tolerance = 1e-10)
  square <- lm(y ~ I(z^2), data = dat)
  reflected <- reference; reflected$z <- -reflected$z
  context <- .msr_curve_reference_context(object, square, reflected)
  expect_true(context$contract$fitted_feature_equivalence)
  expect_false(context$contract$original_raw_values_uniquely_established)
})

test_that("raw references verify only retained fitted rows after native NA exclusion", {
  f <- refit_spatial_fixture(); dat <- f$data
  dat$z[3L] <- NA_real_
  models <- list(
    lm(y ~ mean + base::scale(z), data = dat, na.action = na.exclude),
    glm(y ~ mean + base::scale(z), data = dat, family = poisson(), na.action = na.exclude))
  if (requireNamespace("mgcv", quietly = TRUE))
    models[[3L]] <- mgcv::gam(y ~ base::scale(z) + s(mean, k = 4), data = dat,
      method = "REML", na.action = na.exclude)
  for (model in models) {
    reference <- dat[row.names(model$model), , drop = FALSE]
    row.names(reference) <- row.names(model$model)
    context <- if (inherits(model, "gam"))
      curve_reference_observe_basis_lifetime(f$object, model, reference) else
        .msr_curve_reference_context(f$object, model, reference)
    expect_identical(row.names(context$data), row.names(model$model))
    expect_equal(context$contract$fitted_rows, nrow(dat) - 1L)
    expect_true(context$contract$native_design_and_offsets_verified)
    whole <- dat; row.names(whole) <- row.names(dat)
    expect_error(.msr_curve_reference_context(f$object, model, whole), "row IDs")
  }
})

test_that("analytic raw curves rebuild selected kernels without changing the fitted target", {
  f <- profile_glm_fixture()
  model <- glm(y ~ hab * group + poly(z, 2) + offset(log(exposure)),
    data = f$data, family = poisson())
  contract <- .msr_profile_contract(model, f$inputs)
  solved <- contract$objective$evaluate(log(.3))
  native <- .msr_profile_conditional_model(contract, solved, model, log(.3))
  object <- structure(list(backend = "profile_glm", opt_mod = native,
    kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = log(.3)),
    scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
  reference <- f$data; row.names(reference) <- row.names(native$model)
  context <- .msr_curve_reference_context(object, native, reference)
  expect_equal(as.numeric(context$data$hab), as.numeric(native$model$hab), tolerance = 1e-10)
  expect_identical(context$contract$kernel_source, "verified_selected_preparation")
  before <- object
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(object, reference_data = reference, length.out = 7)
  nd <- data.frame(hab = mean(context$data$hab), group = reference$group[1L],
    z = seq(min(reference$z), max(reference$z), length.out = 7),
    exposure = mean(reference$exposure))
  expected <- predict(native, nd, type = "link", se.fit = TRUE)
  expect_equal(plots$z$data$fit, unname(exp(expected$fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$lwr, unname(exp(expected$fit + qnorm(.025) * expected$se.fit)), tolerance = 1e-8)
  expect_identical(object, before)
})

test_that("retained conditional Cox curves verify raw bases, strata and compound offsets", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture(); dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  contrasts(dat$group) <- contr.sum(nlevels(dat$group))
  form <- as.formula("case_ ~ hab * group + poly(z,2) + offset(log(exposure) + .1*z) + strata(choice)",
    env = asNamespace("survival"))
  model <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE)
  contract <- .msr_clogit_contract(model, f$inputs)
  native <- .msr_clogit_native(contract, log(.3))
  contract$objective <- list(engine = contract$engine, X_fixed = contract$X_fixed,
    offset = contract$offset, response_design = contract$response_design)
  native$call$data <- quote(expired_cox_curve_data)
  object <- structure(list(backend = "conditional_logit", opt_mod = native,
    kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = log(.3)),
    scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
  reference <- dat; row.names(reference) <- row.names(native$model)
  context <- .msr_curve_reference_context(object, native, reference)
  expect_identical(ncol(model.matrix(context$model)), length(coef(native)))
  wrong <- reference; wrong$choice[1L] <- wrong$choice[1L] + 1L
  expect_error(.msr_curve_reference_context(object, native, wrong), "fitted feature|level")
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (response in c(FALSE, TRUE)) {
    plots <- plot_marginal_effects(object, reference_data = reference,
      length.out = 7, link = response)
    expect_false(any(c("choice", "exposure", "case_") %in% names(plots)))
    nd <- data.frame(hab = mean(context$data$hab), group = dat$group[1L],
      z = seq(min(dat$z), max(dat$z), length.out = 7),
      exposure = mean(dat$exposure), choice = dat$choice[1L])
    expected <- predict(native, nd, type = "lp", reference = "zero", se.fit = TRUE)
    expected_fit <- if (response) exp(expected$fit) else expected$fit
    expected_lower <- expected$fit + qnorm(.025) * expected$se.fit
    if (response) expected_lower <- exp(expected_lower)
    expect_equal(plots$z$data$fit, unname(expected_fit), tolerance = 1e-8)
    expect_equal(plots$z$data$lwr, unname(expected_lower), tolerance = 1e-8)
  }
})

test_that("multi-component GAM curves select native link and response columns explicitly", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture(); dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  dat$category <- (seq_len(nrow(dat)) %/% 3L) %% 3L
  make <- function(kind, qualified) {
    sc <- if (qualified) "base::scale(z)" else "scale(z)"
    pol <- if (qualified) "stats::poly(exposure,2)" else "poly(exposure,2)"
    first <- as.formula(paste(if (kind == "gaulss") "y ~" else "category ~",
      "rough * group +", sc, "+ s(mean,k=4) + offset(log(exposure))"), env = asNamespace("mgcv"))
    second <- as.formula(paste("~ fixed +", sc, "+", pol), env = asNamespace("mgcv"))
    mgcv::gam(list(first, second), data = dat,
      family = if (kind == "gaulss") mgcv::gaulss() else mgcv::multinom(K = 2), method = "REML")
  }
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (kind in c("gaulss", "multinom")) {
    bare <- make(kind, FALSE); model <- make(kind, TRUE)
    object <- f$object; object$opt_mod <- model; original <- model
    reference <- dat; row.names(reference) <- row.names(model$model)
    checked <- curve_reference_observe_basis_lifetime(object, model, reference)
    expect_true(checked$contract$native_design_and_offsets_verified)
    predictor <- .msr_predict_gam
    for (bad in c(.01, Inf, NaN)) {
      calls <- 0L
      wrong_basis <- function(object, ...) {
        args <- list(...); value <- do.call(predictor, c(list(object = object), args))
        if (identical(args$type, "lpmatrix")) {
          calls <<- calls + 1L
          if (calls == 2L) value[nrow(value), ncol(value)] <- value[nrow(value), ncol(value)] + bad
        }
        value
      }
      expect_error(with_mocked_bindings(.msr_curve_reference_context(object, model, reference),
        .msr_predict_gam = wrong_basis), "does not reproduce the fitted native design or basis")
      expect_identical(calls, 2L)
    }
    expect_error(plot_marginal_effects(object, reference_data = reference), "explicit native output selection")
    expect_error(plot_marginal_effects(object, component = 0, reference_data = reference), "component")
    expect_error(plot_marginal_effects(object, component = 3, reference_data = reference), "linear predictor")
    for (response in c(FALSE, TRUE)) {
      components <- if (response && kind == "multinom") 1:3 else 1:2
      for (column in components) {
        plots <- plot_marginal_effects(object, reference_data = reference,
          component = column, link = response, length.out = 7)
        vars <- c("rough", "group", "z", "mean", "fixed", "exposure")
        defaults <- lapply(dat[vars], function(v) if (is.numeric(v)) mean(v) else v[1L])
        nd <- as.data.frame(lapply(defaults, rep, times = 7))
        nd$z <- seq(min(dat$z), max(dat$z), length.out = 7)
        expected <- predict(bare, nd, type = if (response) "response" else "link", se.fit = TRUE)
        expect_equal(plots$z$data$fit, unname(expected$fit[, column]), tolerance = 1e-8)
        expect_equal(plots$z$data$lwr,
          unname(expected$fit[, column] + qnorm(.025) * expected$se.fit[, column]), tolerance = 1e-8)
        expect_equal(plots$z$data$upr,
          unname(expected$fit[, column] + qnorm(.975) * expected$se.fit[, column]), tolerance = 1e-8)
        expect_identical(attr(plots, "prediction_contract")$native_output$column, column)
        expect_identical(attr(plots, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
        if (response) expect_identical(attr(plots, "prediction_contract")$quantity, "native_component")
        else if (column == 2L) expect_false(any(c("rough", "mean", "group") %in% names(plots)))
      }
    }
    expect_error(plot_marginal_effects(object, component = 4, link = TRUE,
      reference_data = reference), "available native GAM output")
    expect_identical(object$opt_mod, original)
  }
})

test_that("ordinary retained clogit curves exponentiate native Cox scores", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture(); dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- as.formula("case_ ~ hab + z + strata(choice)", env = asNamespace("survival"))
  model <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE)
  contract <- .msr_clogit_contract(model, f$inputs)
  native <- .msr_clogit_native(contract, log(.3))
  contract$objective <- list(engine = contract$engine, X_fixed = contract$X_fixed,
    offset = contract$offset, response_design = contract$response_design)
  object <- structure(list(backend = "conditional_logit", opt_mod = native,
    kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = log(.3)),
    scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(object, link = TRUE, length.out = 7)
  nd <- data.frame(hab = mean(native$model$hab),
    z = seq(min(native$model$z), max(native$model$z), length.out = 7),
    choice = dat$choice[1L])
  expected <- predict(native, nd, type = "lp", reference = "zero", se.fit = TRUE)
  expect_false("choice" %in% names(plots))
  expect_equal(plots$z$data$fit, unname(exp(expected$fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$lwr,
    unname(exp(expected$fit + qnorm(.025) * expected$se.fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$upr,
    unname(exp(expected$fit + qnorm(.975) * expected$se.fit)), tolerance = 1e-8)
  expect_identical(plots$z$labels$y, "Relative risk")
})

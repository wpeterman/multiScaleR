# Bounded native-provider conformance, using the existing spatial fixture.
# These tests qualify prediction geometry, not recovery or interval calibration.
test_that("native qualified fixed transforms preserve fitted targets and all row routes", {
  f <- refit_spatial_fixture()
  dat <- f$data; dat$y <- dat$y + as.integer(dat$site)^2
  make <- function(kind, qualified) {
    expr <- if (qualified) "base::scale(z)" else "scale(z)"
    rhs <- paste("mean + rough * group + fixed +", expr)
    if (kind == "gam") rhs <- paste("rough * group + fixed +", expr, "+ s(mean,k=4)")
    if (kind %in% c("lmer", "glmmTMB")) rhs <- paste(rhs, "+ (1 | site)")
    form <- as.formula(paste("y ~", rhs), env = if (kind == "gam") asNamespace("mgcv") else environment())
    switch(kind, lm = lm(form, data = dat), glm = glm(form, data = dat, family = poisson()),
      gam = mgcv::gam(form, data = dat, method = "ML"),
      lmer = lme4::lmer(form, data = dat, REML = FALSE),
      glmmTMB = glmmTMB::glmmTMB(form, data = dat, family = gaussian()))
  }
  for (kind in c("lm", "glm", "gam", "lmer", "glmmTMB")) {
    pkg <- switch(kind, gam = "mgcv", lmer = "lme4", glmmTMB = "glmmTMB", "stats")
    if (!requireNamespace(pkg, quietly = TRUE)) next
    bare <- make(kind, FALSE); qualified <- make(kind, TRUE)
    original <- qualified
    repaired <- .msr_native_fixed_metadata(qualified)
    expect_identical(qualified, original)
    expect_identical(coef(repaired), coef(qualified))
    expect_identical(fitted(repaired), fitted(qualified))
    expect_identical(logLik(repaired), logLik(qualified))
    expect_equal(vcov(repaired), vcov(qualified), tolerance = 0)
    expect_identical(.msr_native_fixed_metadata(repaired), repaired)
    if (kind == "glmmTMB") {
      expect_true(qualified$sdr$pdHess)
      expect_identical(repaired$sdr$pdHess, qualified$sdr$pdHess)
    }
    object <- f$object; object$opt_mod <- qualified
    mixed <- kind %in% c("lmer", "glmmTMB")
    controls <- if (mixed) list(re.form = NA) else list()
    nd <- dat[c(5, 23, 17), , drop = FALSE]
    native <- function(data, se = FALSE, type = "response")
      do.call(stats::predict, c(list(object = bare, newdata = data, type = type),
        if (se) list(se.fit = TRUE) else list(), controls))
    live <- function(data, se = FALSE, type = "response")
      do.call(stats::predict, c(list(object = object, newdata = data,
        kernel_inputs = f$inputs, include_fixed = TRUE, type = type,
        se.fit = se, uncertainty = "conditional"), controls))
    for (rows in list(nd, nd[1, , drop = FALSE], transform(nd, z = .2))) {
      expect_equal(as.numeric(live(rows)), as.numeric(native(rows)), tolerance = 1e-8)
      if (kind != "lmer") {
        actual <- live(rows, TRUE); expected <- native(rows, TRUE)
        expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-8)
        expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-8)
      }
    }
    expect_equal(as.numeric(live(nd, type = "link")),
      as.numeric(native(nd, type = if (kind == "lm") "response" else "link")), tolerance = 1e-8)
    snap <- prepare_prediction(object, nd, f$inputs, include_fixed = TRUE, native_args = controls)
    restored <- unserialize(serialize(snap, NULL))
    expect_null(restored$opt_mod)
    expect_equal(as.numeric(predict(restored)), as.numeric(native(nd)), tolerance = 1e-8)
    if (kind != "lmer")
      expect_equal(as.numeric(predict(restored, se.fit = TRUE, uncertainty = "conditional")$se.fit),
        as.numeric(native(nd, TRUE)$se.fit), tolerance = 1e-8)
    expect_error(predict(object, nd, f$inputs, se.fit = TRUE), "conditional")
    scenario <- list(group = "b", z = .2)
    if (mixed) scenario$site <- levels(dat$site)[1L]
    grid <- kernel_scale.raster(f$rasters, multiScaleR = object,
      scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE,
      site_values = scenario, verbose = FALSE)
    map <- do.call(stats::predict, c(list(object = object, newdata = grid), controls))
    expected_rows <- dat; expected_rows$group <- factor("b", levels = levels(dat$group)); expected_rows$z <- .2
    expect_equal(terra::extract(map, f$points, ID = FALSE)[[1]],
      as.numeric(native(expected_rows)), tolerance = 1e-8)
  }
})

test_that("native effect curves retain fixed normalization when raw inputs are retained", {
  f <- refit_spatial_fixture(); dat <- f$data
  bare <- glm(y ~ mean + rough * group + fixed + scale(z) + z:fixed,
    data = dat, family = poisson())
  qualified <- glm(y ~ mean + rough * group + fixed + base::scale(z) + z:fixed,
    data = dat, family = poisson())
  expect_null(extract_namespace(qualified))
  explicit <- qualified; explicit$call[[1L]] <- quote(stats::glm)
  expect_identical(extract_namespace(explicit), "stats")
  explicit$call[[1L]] <- quote(stats:::glm)
  expect_identical(extract_namespace(explicit), "stats")
  expect_identical(namespace(qualified), "stats")
  incomplete <- qualified; incomplete$model <- NULL
  expect_error(.msr_native_fixed_metadata(incomplete), "retained evaluated training frame")
  complete <- .msr_native_fixed_metadata(qualified); complete$model <- NULL
  expect_identical(.msr_native_fixed_metadata(complete), complete)
  object <- f$object; object$opt_mod <- qualified
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(object, length.out = 7)
  vars <- c("mean", "rough", "fixed", "group", "z")
  defaults <- lapply(dat[vars], function(v) if (is.numeric(v)) mean(v) else v[1L])
  nd <- as.data.frame(lapply(defaults, rep, times = 7))
  nd$z <- seq(min(dat$z), max(dat$z), length.out = 7)
  expected <- predict(bare, nd, type = "link", se.fit = TRUE)
  expect_equal(plots$z$data$x, nd$z, tolerance = 1e-10)
  expect_equal(plots$z$data$fit, unname(exp(expected$fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$lwr, unname(exp(expected$fit + qnorm(.025) * expected$se.fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$upr, unname(exp(expected$fit + qnorm(.975) * expected$se.fit)), tolerance = 1e-8)
})

test_that("legacy analytic full fits recover metadata without changing the fitted design", {
  f <- profile_glm_fixture()
  model <- glm(y ~ hab + base::scale(z) + offset(log(exposure)), data = f$data, family = poisson())
  contract <- .msr_profile_contract(model, f$inputs)
  solved <- contract$objective$evaluate(log(.3))
  native <- .msr_profile_conditional_model(contract, solved, model, log(.3))
  object <- structure(list(backend = "profile_glm", opt_mod = native,
    kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = log(.3))), class = "multiScaleR")
  nd <- f$data[c(8, 2, 19), , drop = FALSE]
  correct <- .msr_profile_prediction_design(object, nd, f$inputs)
  legacy <- object; legacy$profile_fit$contract$terms <- terms(model)
  old_terms <- legacy$profile_fit$contract$terms
  expect_true(.msr_needs_fixed_metadata(old_terms))
  expect_equal(.msr_profile_prediction_design(legacy, nd, f$inputs), correct, tolerance = 1e-10)
  expect_identical(legacy$profile_fit$contract$terms, old_terms)
  constant <- nd; constant$z <- .2
  expect_equal(.msr_profile_prediction_design(legacy, constant, f$inputs),
    .msr_profile_prediction_design(object, constant, f$inputs), tolerance = 1e-10)
  thin <- object; thin$opt_mod <- list(family = native$family)
  expect_false(.msr_needs_fixed_metadata(contract$terms))
  expect_equal(.msr_profile_prediction_design(thin, nd, f$inputs), correct, tolerance = 1e-10)
  expect_false(.msr_needs_fixed_metadata(terms(y ~ base::scale(z, center = FALSE, scale = FALSE))))
  expect_false(.msr_needs_fixed_metadata(terms(y ~ stats::poly(z, 2, raw = TRUE))))
})

# Fixed, bounded software conformance. These tests do not estimate ecological
# scales, assess recovery, or qualify intervals beyond native conditional SEs.
native_nested_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- refit_spatial_fixture()
    cache
  }
})

test_that("native frames and unmarked curve bases freeze nested standard learners", {
  train <- data.frame(y = seq_len(24), z = sin(seq_len(24)/3) + seq_len(24)/20)
  next_rows <- train[c(18, 3, 18, 9), , drop = FALSE]
  next_rows$z <- next_rows$z + .15
  expressions <- c(scale = "scale(z)", poly = "poly(z, 2)",
    ns = "splines::ns(z, df = 3)", bs = "splines::bs(z, df = 3)")
  reference <- function(name, data) switch(name,
    scale = (data$z - mean(train$z))/sd(train$z),
    poly = predict(poly(train$z, 2), data$z),
    ns = predict(splines::ns(train$z, df = 3), data$z),
    bs = predict(splines::bs(train$z, df = 3), data$z))
  for (provider in c("native", "unmarked_curve")) for (wrapper in c("bare", "I")) {
    for (name in names(expressions)) {
      ex <- expressions[[name]]
      if (wrapper == "I") ex <- paste0("I(", ex, ")")
      form <- as.formula(paste(if (provider == "native") "y ~" else "~", ex))
      frame <- model.frame(form, train)
      original <- terms(frame)
      tt <- if (provider == "native") .msr_fixed_prediction_terms(original, frame, nested = TRUE) else
        .msr_unmarked_curve_basis(form, train)$terms
      expect_false(.msr_needs_nested_metadata(tt))
      expect_identical(attr(tt, "variables"), attr(original, "variables"))
      expect_identical(attr(tt, "term.labels"), attr(original, "term.labels"))
      for (nd in list(next_rows, next_rows[1, , drop = FALSE], transform(next_rows, z = .2))) {
        mf <- model.frame(delete.response(tt), nd)
        actual <- model.matrix(delete.response(tt), mf)
        expected <- cbind(1, reference(name, nd))
        expect_equal(unname(actual), unname(expected), tolerance = 1e-12, ignore_attr = TRUE)
      }
    }
  }
  form <- y ~ offset(log1p(abs(scale(z))))
  frame <- model.frame(form, train)
  for (tt in list(.msr_fixed_prediction_terms(terms(frame), frame, nested = TRUE),
      .msr_unmarked_curve_basis(~offset(log1p(abs(scale(z)))), train)$terms)) {
    for (nd in list(next_rows, next_rows[1, , drop = FALSE])) {
      actual <- model.offset(model.frame(delete.response(tt), nd))
      expect_equal(as.numeric(actual), log1p(abs(reference("scale", nd))), tolerance = 1e-12)
      expect_identical(dim(actual), c(nrow(nd), 1L))
    }
  }
})

test_that("native nested recovery respects ambiguity, bindings and saved metadata", {
  dat <- data.frame(y = 1:24, z = sin((1:24)/3), w = cos((1:24)/5))
  expect_true(.msr_standard_prediction_function(base::scale, "scale", "base"))
  expect_true(.msr_standard_prediction_function(getExportedValue("terra", "scale"), "scale", "base"))
  expect_false(.msr_standard_prediction_function(function(x) x + 2, "scale", "base"))
  custom_generic <- getExportedValue("terra", "scale")
  method <- methods::slot(custom_generic, "default"); body(method) <- quote(x + 2)
  methods::slot(custom_generic, "default") <- method
  expect_false(.msr_standard_prediction_function(custom_generic, "scale", "base"))
  custom_generic <- getExportedValue("terra", "scale")
  scope <- new.env(parent = environment(custom_generic))
  scope$standardGeneric <- function(...) 0
  environment(custom_generic) <- scope
  expect_true(methods::is(custom_generic, "standardGeneric"))
  expect_false(.msr_standard_prediction_function(custom_generic, "scale", "base"))
  env <- new.env(parent = baseenv()); env$scale <- custom_generic
  custom <- as.formula("y ~ I(scale(z))", env = env)
  expect_error(.msr_cv_prediction_terms(terms(custom), dat, qualify = TRUE), "Unqualified custom learner")
  fixed <- function(form) {
    frame <- model.frame(form, dat)
    .msr_fixed_prediction_terms(terms(frame), frame, nested = TRUE)
  }
  tt <- fixed(y ~ I(base::scale(z, cen = FALSE, sc = TRUE)))
  mf <- model.frame(delete.response(tt), dat[3, , drop = FALSE])
  expected <- dat$z[3]/sqrt(sum(dat$z^2)/(nrow(dat)-1))
  expect_equal(as.numeric(mf[[1]]), expected, tolerance = 1e-12)
  tt <- fixed(y ~ base::scale(z, cen = FALSE, sc = TRUE))
  expect_false(.msr_needs_nested_metadata(tt))
  expect_equal(as.numeric(model.frame(delete.response(tt), dat[3, , drop = FALSE])[[1]]),
    expected, tolerance = 1e-12)
  expect_error(fixed(y ~ I(scale(z) + scale(w))), "unambiguous original training metadata")
  expect_error(fixed(y ~ I(scale(rank(z)))), "supported pointwise")
  frame <- model.frame(y ~ I(scale(z)), dat)
  frame[[2]] <- I(unclass(frame[[2]])); attr(frame[[2]], "scaled:center") <- NULL
  expect_error(.msr_fixed_prediction_terms(terms(frame), frame, nested = TRUE), "Retained attributes")
  original <- terms(model.frame(y ~ I(scale(z)), dat))
  expect_true(.msr_needs_nested_metadata(original))
  expect_identical(.msr_fixed_prediction_terms(original, model.frame(y ~ I(scale(z)), dat)), original)
  expect_false(.msr_needs_nested_metadata(terms(I(scale(y)) ~ z)))
  env <- new.env(parent = baseenv()); env$scale <- function(x) x + 2
  custom <- as.formula("y ~ I(scale(z))", env = env)
  frame <- model.frame(custom, dat)
  expect_identical(.msr_fixed_prediction_terms(terms(frame), frame, nested = TRUE), terms(frame))
  expect_identical(as.numeric(model.frame(delete.response(terms(frame)), dat[3, , drop = FALSE])[[1]]),
    dat$z[3] + 2)
  complete <- fixed(y ~ I(scale(z)))
  ready <- list(terms = complete, model = NULL)
  class(ready) <- "lm"
  expect_identical(.msr_native_fixed_metadata(ready, nested = TRUE), ready)
  missing <- ready; missing$terms <- original
  expect_error(.msr_native_fixed_metadata(missing, nested = TRUE), "retained evaluated training frame")
  complete <- terms(y ~ base::scale(sc = 1.8, x = z, cen = .4))
  attr(complete, "predvars") <- attr(complete, "variables")
  ready <- structure(list(terms = complete, model = NULL), class = "lm")
  expect_identical(.msr_native_fixed_metadata(ready, nested = TRUE), ready)
  expect_error(.msr_native_fixed_metadata(ready), "retained evaluated training frame")
  thin <- list(profile_fit = list(contract = list(terms = complete)))
  expect_identical(.msr_profile_prediction_terms(thin), complete)
  # Raw unmarked covariates can identify several learners independently.
  basis <- .msr_unmarked_curve_basis(~I(scale(z) + scale(w)), dat)
  actual <- model.matrix(basis$terms, model.frame(basis$terms, dat[3, , drop = FALSE]))
  expect_equal(as.numeric(actual), c(1, as.numeric(scale(dat$z)[3] + scale(dat$w)[3])), tolerance = 1e-12)
})

test_that("reordered learned arguments keep their original training constants", {
  train <- data.frame(y = 1:24, z = sin((1:24)/3) + (1:24)/20)
  nd <- train[c(18, 3, 18, 9), , drop = FALSE]; nd$z <- nd$z + .15
  calls <- c("base::scale(sc=TRUE,x=z,cen=FALSE)",
    "base::scale(sc=1.8,x=z,cen=.4)", "base::scale(cen=TRUE,x=z,sc=FALSE)")
  centers <- c(0, .4, mean(train$z))
  spreads <- c(sqrt(sum(train$z^2)/(nrow(train)-1)), 1.8, 1)
  for (i in seq_along(calls)) for (wrapper in c("bare", "I", "offset")) {
    ex <- switch(wrapper, bare = calls[i], I = paste0("I(", calls[i], ")"),
      offset = paste0("offset(log1p(abs(", calls[i], ")))"))
    frame <- model.frame(as.formula(paste("y ~", ex)), train)
    tt <- .msr_fixed_prediction_terms(terms(frame), frame, nested = TRUE)
    expect_false(.msr_needs_nested_metadata(tt))
    for (rows in list(nd, nd[1, , drop = FALSE], transform(nd, z = .2))) {
      expected <- (rows$z - centers[i])/spreads[i]
      if (wrapper == "offset") expected <- log1p(abs(expected))
      actual <- model.frame(delete.response(tt), rows)[[1]]
      expect_equal(as.numeric(actual), expected, tolerance = 1e-12)
    }
  }
})

test_that("native nested predictions and conditional SEs agree across row and map routes", {
  f <- native_nested_fixture(); dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  make <- function(kind, nested) {
    expression <- if (nested) "I(scale(z))" else "scale(z)"
    rhs <- paste("mean + rough * group + fixed +", expression, "+ offset(log(exposure))")
    if (kind == "gam") rhs <- paste("rough * group + fixed +", expression,
      "+ s(mean, k=4) + offset(log(exposure))")
    if (kind %in% c("lmer", "glmmTMB")) rhs <- paste(rhs, "+ (1 | site)")
    form <- as.formula(paste("y ~", rhs),
      env = if (kind == "gam") asNamespace("mgcv") else environment())
    switch(kind, lm = lm(form, data = dat), glm = glm(form, data = dat, family = poisson()),
      gam = mgcv::gam(form, data = dat, method = "ML"),
      lmer = lme4::lmer(form, data = dat, REML = FALSE),
      glmmTMB = glmmTMB::glmmTMB(form, data = dat, family = gaussian()))
  }
  for (kind in c("lm", "glm", "gam", "lmer", "glmmTMB")) {
    pkg <- switch(kind, gam = "mgcv", lmer = "lme4", glmmTMB = "glmmTMB", "stats")
    if (!requireNamespace(pkg, quietly = TRUE)) next
    bare <- make(kind, FALSE); nested <- make(kind, TRUE); original <- nested
    repaired <- .msr_native_fixed_metadata(nested, nested = TRUE)
    expect_identical(nested, original)
    expect_identical(coef(repaired), coef(nested))
    expect_identical(fitted(repaired), fitted(nested))
    expect_identical(logLik(repaired), logLik(nested))
    expect_equal(vcov(repaired), vcov(nested), tolerance = 0)
    expect_identical(.msr_native_fixed_metadata(repaired, nested = TRUE), repaired)
    if (kind == "glmmTMB") expect_identical(repaired$sdr$pdHess, nested$sdr$pdHess)
    object <- f$object; object$opt_mod <- nested
    controls <- if (kind %in% c("lmer", "glmmTMB")) list(re.form = NA) else list()
    nd <- dat[c(5, 23, 17), , drop = FALSE]
    native <- function(data, se = FALSE) do.call(stats::predict,
      c(list(object = bare, newdata = data, type = "response"),
        if (se) list(se.fit = TRUE) else list(), controls))
    live <- function(data, se = FALSE) do.call(stats::predict,
      c(list(object = object, newdata = data, kernel_inputs = f$inputs,
        include_fixed = TRUE, type = "response", se.fit = se, uncertainty = "conditional"), controls))
    for (rows in list(nd, nd[1, , drop = FALSE], transform(nd, z = .2))) {
      expect_equal(as.numeric(live(rows)), as.numeric(native(rows)), tolerance = 1e-8)
      if (kind != "lmer") {
        expect_equal(as.numeric(live(rows, TRUE)$fit), as.numeric(native(rows, TRUE)$fit), tolerance = 1e-8)
        expect_equal(as.numeric(live(rows, TRUE)$se.fit), as.numeric(native(rows, TRUE)$se.fit), tolerance = 1e-8)
      }
    }
    snap <- prepare_prediction(object, nd, f$inputs, include_fixed = TRUE, native_args = controls)
    restored <- unserialize(serialize(snap, NULL))
    expect_null(restored$opt_mod)
    expect_equal(as.numeric(predict(restored)), as.numeric(native(nd)), tolerance = 1e-8)
    if (kind != "lmer") expect_equal(as.numeric(predict(restored, se.fit = TRUE,
      uncertainty = "conditional")$se.fit), as.numeric(native(nd, TRUE)$se.fit), tolerance = 1e-8)
    expect_error(predict(object, nd, f$inputs, se.fit = TRUE), "conditional")
    scenario <- list(group = "b", z = .2, exposure = 1.2)
    if (kind %in% c("lmer", "glmmTMB")) scenario$site <- levels(dat$site)[1]
    grid <- kernel_scale.raster(f$rasters, multiScaleR = object, scale_center = TRUE,
      include_fixed = TRUE, match_buffer = TRUE, site_values = scenario, verbose = FALSE)
    map <- do.call(stats::predict, c(list(object = object, newdata = grid), controls))
    reference <- dat; reference$group <- factor("b", levels = levels(dat$group))
    reference$z <- .2; reference$exposure <- 1.2
    expect_equal(terra::extract(map, f$points, ID = FALSE)[[1]], as.numeric(native(reference)), tolerance = 1e-8)
    expect_identical(object$scale_est, f$object$scale_est)
    if (kind == "glmmTMB") {
      # One declared extra fit. Both aliases lead their calls, and both learned
      # predictors share a provider prefix. The kernel mean is already on its
      # training mean/SD scale, so the bare model is an independent row oracle.
      aliased <- glmmTMB::glmmTMB(y ~ base::scale(sc=TRUE,x=mean,cen=TRUE) + rough * group +
        fixed + base::scale(sc=TRUE,x=z,cen=TRUE) +
        offset(log(exposure)) + (1 | site), data = dat, family = gaussian())
      expect_true(aliased$sdr$pdHess)
      expect_equal(as.numeric(logLik(aliased)), as.numeric(logLik(bare)), tolerance = 1e-8)
      alias_object <- object; alias_object$opt_mod <- aliased
      expect_equal(.msr_selection_parameter_count(alias_object), .msr_selection_parameter_count(object))
      repaired <- .msr_native_fixed_metadata(aliased, nested = TRUE)
      expect_identical(coef(repaired), coef(aliased))
      expect_identical(repaired$sdr$pdHess, aliased$sdr$pdHess)
      expect_equal(vcov(repaired), vcov(aliased), tolerance = 0)
      for (rows in list(nd, nd[1, , drop = FALSE], transform(nd, z = .2))) {
        actual <- predict(alias_object, rows, f$inputs, include_fixed = TRUE,
          se.fit = TRUE, uncertainty = "conditional", re.form = NA)
        expect_equal(as.numeric(actual$fit), as.numeric(native(rows, TRUE)$fit), tolerance = 1e-8)
        expect_equal(as.numeric(actual$se.fit), as.numeric(native(rows, TRUE)$se.fit), tolerance = 1e-8)
      }
      snap <- prepare_prediction(alias_object, nd, f$inputs, include_fixed = TRUE,
        native_args = list(re.form = NA))
      restored <- unserialize(serialize(snap, NULL))
      expect_equal(as.numeric(predict(restored)), as.numeric(native(nd)), tolerance = 1e-8)
      expect_equal(as.numeric(predict(restored, se.fit = TRUE, uncertainty = "conditional")$se.fit),
        as.numeric(native(nd, TRUE)$se.fit), tolerance = 1e-8)
      map <- predict(alias_object, grid, re.form = NA, se.fit = TRUE, uncertainty = "conditional")
      expected <- native(reference, TRUE)
      expect_equal(unname(as.matrix(terra::extract(map, f$points, ID = FALSE))),
        unname(cbind(expected$fit, expected$se.fit)), tolerance = 1e-8)
    }
  }
})

test_that("nested native curves retain the trained basis and native interval target", {
  f <- native_nested_fixture(); dat <- f$data
  bare <- glm(y ~ mean + rough * group + fixed + scale(z) + z:fixed,
    data = dat, family = poisson())
  nested <- glm(y ~ mean + rough * group + fixed + I(scale(z)) + z:fixed,
    data = dat, family = poisson())
  object <- f$object; object$opt_mod <- nested
  row.names(dat) <- row.names(nested$model)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(object, length.out = 7, reference_data = dat)
  vars <- c("mean", "rough", "fixed", "group", "z")
  defaults <- lapply(dat[vars], function(v) if (is.numeric(v)) mean(v) else v[1L])
  nd <- as.data.frame(lapply(defaults, rep, times = 7))
  nd$z <- seq(min(dat$z), max(dat$z), length.out = 7)
  expected <- predict(bare, nd, type = "link", se.fit = TRUE)
  expect_equal(plots$z$data$x, nd$z, tolerance = 1e-10)
  expect_equal(plots$z$data$fit, unname(exp(expected$fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$lwr, unname(exp(expected$fit + qnorm(.025)*expected$se.fit)), tolerance = 1e-8)
  expect_equal(plots$z$data$upr, unname(exp(expected$fit + qnorm(.975)*expected$se.fit)), tolerance = 1e-8)
})

test_that("nested conditional choice predictions retain Cox targets and guards", {
  skip_if_not_installed("survival")
  f <- native_nested_fixture(); dat <- f$data
  dat$step <- factor(rep(seq_len(12), each = 3))
  dat$case <- as.integer(rep(1:3, 12) == rep(rep(1:3, length.out = 12), each = 3))
  coxph <- survival::coxph
  strata <- survival::strata
  bare <- survival::clogit(case ~ mean + rough + scale(exposure) + poly(z, 2) + strata(step),
    dat, method = "exact", model = TRUE, x = TRUE, y = TRUE)
  nested <- survival::clogit(case ~ mean + rough + I(scale(exposure)) + I(poly(z, 2)) + strata(step),
    dat, method = "exact", model = TRUE, x = TRUE, y = TRUE)
  repaired <- .msr_native_fixed_metadata(nested, nested = TRUE)
  for (name in c("coefficients", "var", "linear.predictors", "loglik", "method", "strata", "iter", "y"))
    expect_identical(repaired[[name]], nested[[name]])
  object <- f$object; object$opt_mod <- nested
  nd <- dat[c(5, 23, 17), , drop = FALSE]
  for (rows in list(nd, nd[1, , drop = FALSE], transform(nd, exposure = 1.2))) {
    expected <- predict(bare, rows, type = "lp", se.fit = TRUE, reference = "zero")
    actual <- predict(object, rows, f$inputs, include_fixed = TRUE, type = "lp",
      se.fit = TRUE, uncertainty = "conditional", reference = "zero")
    expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-8)
    expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-8)
  }
  unsupported <- nested; class(unsupported) <- c("coxph.penal", class(unsupported))
  expect_error(.msr_native_fixed_metadata(unsupported, nested = TRUE), "time-fixed, unpenalized")
  expect_error(predict(object, nd, f$inputs, type = "response"), "type")
})

test_that("legacy analytic point and raster terms share the original nested basis", {
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had) seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", seed, .GlobalEnv) else
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  f <- profile_glm_fixture()
  model <- glm(y ~ hab + I(scale(z)) + offset(log(exposure)), data = f$data, family = poisson())
  contract <- .msr_profile_contract(model, f$inputs)
  solved <- contract$objective$evaluate(log(.3))
  native <- .msr_profile_conditional_model(contract, solved, model, log(.3))
  object <- structure(list(backend = "profile_glm", opt_mod = native, kernel_inputs = f$inputs,
    profile_fit = list(contract = contract, par = log(.3))), class = "multiScaleR")
  original <- contract$terms
  tt <- .msr_profile_prediction_terms(object)
  expect_false(.msr_needs_nested_metadata(tt))
  expect_identical(object$profile_fit$contract$terms, original)
  nd <- f$data[c(8, 2, 19), , drop = FALSE]
  for (rows in list(nd, nd[1, , drop = FALSE], transform(nd, z = .2))) {
    design <- .msr_profile_prediction_design(object, rows, f$inputs)
    expect_equal(as.numeric(design$X[, "I(scale(z))"]),
      (rows$z - mean(f$data$z))/sd(f$data$z), tolerance = 1e-12)
    rows$hab <- design$X[, "hab"]
    # Analytic raster chunks evaluate these exact frozen terms on standardized
    # smoothed covariates, without reentering the native training frame.
    chunk <- model.frame(delete.response(tt), rows)
    chunk_X <- model.matrix(delete.response(tt), chunk)[, colnames(design$X), drop = FALSE]
    expect_equal(unname(chunk_X), unname(design$X), tolerance = 1e-12)
    expect_equal(as.numeric(model.offset(chunk)), as.numeric(design$offset), tolerance = 1e-12)
  }
})

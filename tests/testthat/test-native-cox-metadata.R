# Two declared native conditional models at fixed physical ecological scales.
# These qualify software prediction routes, not recovery or interval calibration.
native_cox_metadata_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    coxph <- survival::coxph
    strata <- survival::strata
    f <- refit_spatial_fixture(); dat <- f$data
    dat$step <- factor(rep(seq_len(12), each = 3))
    dat$case <- as.integer(rep(1:3, 12) == rep(rep(1:3, length.out = 12), each = 3))
    forms <- list(bare = case ~ mean + rough + scale(exposure) + poly(z, 2) + strata(step),
      qualified = case ~ mean + rough + base::scale(exposure) + stats::poly(z, 2) + strata(step))
    models <- lapply(forms, function(form) survival::clogit(form, data = dat,
      method = "exact", model = TRUE, x = TRUE, y = TRUE, na.action = stats::na.fail))
    for (model in models) {
      stopifnot(length(coef(model)) == 5L, all(is.finite(coef(model))), is.null(model$fail))
      X <- model$x; score <- numeric(ncol(X)); H <- matrix(0, ncol(X), ncol(X))
      for (rows in split(seq_len(nrow(dat)), dat$step)) {
        x <- X[rows, , drop = FALSE]; eta <- drop(x %*% coef(model))
        p <- exp(eta - max(eta)); p <- p / sum(p)
        score <- score + colSums(x * (dat$case[rows] - p))
        H <- H + crossprod(x, x * p) - tcrossprod(drop(crossprod(x, p)))
      }
      stopifnot(max(abs(score)) < 1e-4,
        min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) > 0,
        min(eigen(model$var, symmetric = TRUE, only.values = TRUE)$values) > 0)
    }
    cache <<- list(f = f, data = dat, models = models)
    cache
  }
})

test_that("ordinary native Cox rows retain trained transformations and targets", {
  skip_if_not_installed("survival")
  f <- native_cox_metadata_fixture(); model <- f$models$qualified; original <- model
  fixed <- .msr_native_fixed_metadata(model)
  expect_identical(model, original)
  for (name in c("coefficients", "var", "linear.predictors", "loglik", "method", "strata", "iter", "y"))
    expect_identical(fixed[[name]], model[[name]])
  expect_identical(.msr_native_fixed_metadata(fixed), fixed)
  object <- f$f$object; object$opt_mod <- model
  expect_equal(.msr_selection_parameter_count(object), as.numeric(attr(logLik(model), "df")) + 2L)
  expect_equal(.msr_selection_nobs(model), 12L)
  for (type in c("lp", "risk", "terms")) {
    full <- predict(model, f$data, type = type, reference = "zero", se.fit = TRUE)
    for (rows in list(seq_len(nrow(f$data)), c(1L, 17L, 24L), 17L)) {
      nd <- f$data[rows, , drop = FALSE]
      slice <- function(x) if (is.matrix(x)) x[rows, , drop = FALSE] else x[rows]
      expected <- lapply(full, slice)
      actual <- predict(object, nd, f$f$inputs, include_fixed = TRUE, type = type,
        reference = "zero", se.fit = TRUE, uncertainty = "conditional")
      for (name in c("fit", "se.fit")) {
        expect_identical(dim(actual[[name]]), dim(expected[[name]]))
        expect_identical(dimnames(actual[[name]]), dimnames(expected[[name]]))
        expect_identical(names(actual[[name]]), names(expected[[name]]))
        expect_equal(as.numeric(actual[[name]]), as.numeric(expected[[name]]), tolerance = 1e-8)
        expect_true(all(is.finite(actual[[name]])))
      }
      snapshot <- prepare_prediction(object, nd, f$f$inputs, include_fixed = TRUE,
        component = type, native_args = list(reference = "zero"))
      snapshot <- unserialize(serialize(snapshot, NULL))
      expect_null(snapshot$opt_mod)
      saved <- predict(snapshot, type = type, se.fit = TRUE, uncertainty = "conditional")
      expect_equal(saved, actual, tolerance = 1e-8)
    }
  }
  one <- f$data[17L, , drop = FALSE]
  expect_error(predict(object, one, f$f$inputs, type = "response"), "no generic response probability")
  expect_error(predict(object, one, f$f$inputs, se.fit = TRUE), "conditional")
  expect_error(predict(object, one, f$f$inputs, se.fit = TRUE, uncertainty = "joint"), "conditional")
})

test_that("native Cox checked maps and raw-reference curves keep native scales", {
  skip_if_not_installed("survival")
  f <- native_cox_metadata_fixture(); dat <- f$data; model <- f$models$qualified
  object <- f$f$object; object$opt_mod <- model
  grid <- kernel_scale.raster(f$f$rasters, multiScaleR = object, include_fixed = TRUE,
    scale_center = TRUE, match_buffer = TRUE,
    site_values = list(step = levels(dat$step)[1L], exposure = 1.25, z = .2), verbose = FALSE)
  nd <- dat; nd$step <- factor(levels(dat$step)[1L], levels = levels(dat$step))
  nd$exposure <- 1.25; nd$z <- .2
  for (type in c("lp", "risk")) {
    # The audit separately establishes that native bare calls keep trained metadata.
    expected <- predict(f$models$bare, nd, type = type, reference = "zero", se.fit = TRUE)
    map <- predict(object, grid, type = type, reference = "zero", se.fit = TRUE, uncertainty = "conditional")
    expect_equal(unname(as.matrix(terra::extract(map, f$f$points, ID = FALSE))),
      unname(cbind(expected$fit, expected$se.fit)), tolerance = 1e-8)
  }
  model$call$data <- quote(expired_native_cox_raw_data)
  attr(model$terms, ".Environment") <- asNamespace("survival")
  attr(attr(model$model, "terms"), ".Environment") <- asNamespace("survival")
  object$opt_mod <- model; original <- model
  reference <- dat; row.names(reference) <- row.names(model$model)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (risk in c(FALSE, TRUE)) {
    curves <- plot_marginal_effects(object, length.out = 7, link = risk,
      reference_data = reference[rev(seq_len(nrow(reference))), , drop = FALSE])
    expect_false("step" %in% names(curves))
    expect_true(all(c("mean", "rough", "exposure", "z") %in% names(curves)))
    nd <- data.frame(mean = mean(dat$mean), rough = mean(dat$rough),
      exposure = mean(dat$exposure), z = seq(min(dat$z), max(dat$z), length.out = 7),
      step = factor(levels(dat$step)[1L], levels = levels(dat$step)))
    expected <- predict(f$models$bare, nd, type = "lp", reference = "zero", se.fit = TRUE)
    inverse <- if (risk) exp else identity
    expect_equal(curves$z$data$fit, unname(inverse(expected$fit)), tolerance = 1e-8)
    expect_equal(curves$z$data$lwr, unname(inverse(expected$fit + qnorm(.025) * expected$se.fit)), tolerance = 1e-8)
    expect_equal(curves$z$data$upr, unname(inverse(expected$fit + qnorm(.975) * expected$se.fit)), tolerance = 1e-8)
    expect_true(attr(curves, "reference_data_contract")$native_design_and_offsets_verified)
  }
  expect_identical(object$opt_mod, original)
  wrong <- reference; wrong$exposure[1L] <- wrong$exposure[1L] + .1
  expect_error(plot_marginal_effects(object, reference_data = wrong), "fitted feature")
})

test_that("Cox transformation recovery retains unsupported-structure boundaries", {
  skip_if_not_installed("survival")
  model <- native_cox_metadata_fixture()$models$qualified
  missing <- model; missing$model <- NULL
  expect_error(.msr_native_fixed_metadata(missing), "retained evaluated training frame")
  ready <- .msr_native_fixed_metadata(model); ready$model <- NULL
  expect_identical(.msr_native_fixed_metadata(ready), ready)
  # Dispatch/terms metadata probes, not fitted penalized, multistate or tt models.
  for (kind in c("coxph.penal", "coxphms")) {
    unsupported <- model; class(unsupported) <- c(kind, class(model))
    expect_error(.msr_native_fixed_metadata(unsupported), "time-fixed, unpenalized")
  }
  unsupported <- model; specials <- attr(unsupported$terms, "specials")
  specials$tt <- 1L; attr(unsupported$terms, "specials") <- specials
  expect_error(.msr_native_fixed_metadata(unsupported), "time-fixed, unpenalized")
})

test_that("Cox replay checks retained rows before restoring a missing training basis", {
  skip_if_not_installed("survival")
  f <- native_cox_metadata_fixture()
  model <- f$models$qualified
  column <- "base::scale(exposure)"
  # Emulate the retained column on providers that omit matrix attributes.
  attr(model$model[[column]], "scaled:center") <- NULL
  attr(model$model[[column]], "scaled:scale") <- NULL
  expect_true(.msr_needs_fixed_metadata(terms(model)))
  attr(model$terms, ".Environment") <- asNamespace("survival")
  attr(attr(model$model, "terms"), ".Environment") <- asNamespace("survival")
  model$call$data <- quote(unavailable_original_cox_training)
  original <- model
  damaged <- f$data
  damaged$exposure[1L] <- damaged$exposure[1L] + .1
  misaligned <- f$data[rev(seq_len(nrow(f$data))), , drop = FALSE]
  for (training in list(NULL, damaged, misaligned)) {
    candidate <- .msr_native_fixed_metadata(model, training = training)
    expect_true(.msr_needs_fixed_metadata(terms(candidate)))
  }
  restored <- .msr_native_fixed_metadata(model, training = f$data)
  expect_false(.msr_needs_fixed_metadata(terms(restored)))
  expect_identical(restored$x, model$x)
  expect_identical(coef(restored), coef(model))
  expect_identical(restored$var, model$var)
  expect_identical(restored$loglik, model$loglik)
  expect_identical(model, original)
  reference <- f$data
  row.names(reference) <- row.names(model$model)
  # These different raw units produce the same scaled exposure values.
  reference$exposure <- 10 + 2 * reference$exposure
  object <- f$f$object
  object$opt_mod <- model
  expect_error(plot_marginal_effects(object, reference_data = reference),
    "original nested training learner")
  nd <- f$data[17L, , drop = FALSE]
  expect_equal(as.numeric(predict(restored, nd, reference = "zero")),
    as.numeric(predict(f$models$bare, nd, reference = "zero")), tolerance = 1e-10)
})

test_that("a complete retained Cox recipe supports nested conditional choice prediction", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture()
  dat <- f$data
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- factor(rep(seq_len(nrow(dat) / 3), each = 3))
  formula <- as.formula("case_ ~ hab + I(base::scale(z)) + strata(choice)",
    env = asNamespace("survival"))
  model <- survival::clogit(formula, dat, model = TRUE, x = TRUE, y = TRUE,
    na.action = stats::na.fail)
  contract <- .msr_clogit_contract(model, f$inputs)
  fit <- .msr_clogit_native(contract, log(.3))
  nd <- dat[c(4, 12), , drop = FALSE]
  nd$z <- c(.1, .2)
  nd$hab <- c(-.3, .5)
  expected <- coef(fit)[1L] * nd$hab + coef(fit)[2L] * (nd$z - mean(dat$z)) / sd(dat$z)
  expect_equal(as.numeric(predict(fit, nd, reference = "zero")), expected, tolerance = 1e-10)
})

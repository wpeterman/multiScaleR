# Three fixed native software specifications, without an outer scale search.
zip_response_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- refit_spatial_fixture(); dat <- f$data
    specifications <- list(
      plain = list(y ~ mean + rough + z, ~1),
      offset = list(y ~ mean + rough + z + offset(log(exposure)),
                    ~1 + offset(log(exposure))),
      factor = list(y ~ mean + rough + group + base::scale(z), ~1))
    models <- lapply(names(specifications), function(kind) {
      model <- mgcv::gam(specifications[[kind]], data = dat,
        family = mgcv::ziplss(), method = "REML")
      model <- .msr_native_fixed_metadata(model)
      X <- stats::predict(model, type = "lpmatrix")
      offsets <- if (kind == "offset") rep(list(log(dat$exposure)), 2) else list(0, 0)
      differential <- model$family$ll(model$y, X, coef(model), rep(1, nrow(dat)),
        model$family, offset = offsets, deriv = 1)
      expect_lt(max(abs(differential$lb)), 1e-4)
      expect_gt(min(eigen(-differential$lbb, symmetric = TRUE, only.values = TRUE)$values), 0)
      expect_gt(min(eigen(model$Vp, symmetric = TRUE, only.values = TRUE)$values), 0)
      if (!is.null(model$converged)) expect_true(model$converged)
      list(model = model, native_df = as.numeric(attr(logLik(model), "df")))
    })
    names(models) <- names(specifications)
    cache <<- list(f = f, data = dat, models = models)
    cache
  }
})

# Native full-context means supply an independent coefficient derivative,
# including when a requested slice contains only one fitted factor level.
zip_native_mean <- function(model, data, training, controls = list()) {
  full <- rbind(training, data)
  value <- do.call(stats::predict, c(list(object = model, newdata = full,
    type = "response"), controls))
  stopifnot(is.numeric(value), length(value) == nrow(full),
    is.null(dim(value)) || length(dim(value)) == 1L ||
      (length(dim(value)) == 2L && ncol(value) == 1L))
  as.numeric(value)[nrow(training) + seq_len(nrow(data))]
}
zip_coefficient_gradient <- function(model, data, training, multiplier = 1, controls = list()) {
  beta <- coef(model); G <- matrix(0, nrow(data), length(beta))
  for (j in seq_along(beta)) {
    h <- multiplier * 1e-5 * max(1, abs(beta[j]))
    plus <- minus <- model
    plus$coefficients[j] <- beta[j] + h
    minus$coefficients[j] <- beta[j] - h
    G[, j] <- (zip_native_mean(plus, data, training, controls) -
      zip_native_mean(minus, data, training, controls)) / (2 * h)
  }
  G
}

test_that("ziplss response uncertainty uses both predictors and signed covariance", {
  skip_if_not_installed("mgcv")
  f <- zip_response_fixture(); dat <- f$data
  native <- utils::getS3method("predict", "gam", envir = asNamespace("mgcv"))
  saved_body <- body(native); saved_environment <- environment(native)
  for (kind in names(f$models)) {
    case <- f$models[[kind]]; model <- case$model; original <- model
    object <- f$f$object; object$opt_mod <- model
    expect_equal(.msr_selection_parameter_count(object), case$native_df + 2L, tolerance = 1e-8)
    covariances <- list(fitted = model$Vp)
    if (kind == "plain") {
      negative <- diag(c(.3, .2, .4, .1, .5))
      negative[1, 5] <- negative[5, 1] <- -.05
      positive <- negative; positive[1, 5] <- positive[5, 1] <- .05
      # Explicit algebraic inputs, not fitted scientific covariance estimates.
      covariances$negative_cross <- negative; covariances$positive_cross <- positive
    }
    for (rows in list(seq_len(nrow(dat)), c(1L, 17L, 24L), 17L)) {
      nd <- dat[rows, , drop = FALSE]
      G <- zip_coefficient_gradient(model, nd, dat)
      half <- zip_coefficient_gradient(model, nd, dat, .5)
      expect_lt(max(abs(G - half)), 1e-6)
      expected_mean <- zip_native_mean(model, nd, dat)
      for (V in covariances) {
        variant <- model; variant$Vp <- V
        object$opt_mod <- variant
        expected_se <- sqrt(rowSums((half %*% V) * half))
        value <- predict(object, nd, f$f$inputs, include_fixed = TRUE,
          type = "response", se.fit = TRUE, uncertainty = "conditional",
          native_args = list(block.size = 1))
        expect_equal(as.numeric(value$fit), expected_mean, tolerance = 1e-8)
        expect_equal(as.numeric(value$se.fit), expected_se, tolerance = 1e-6)
        expect_lt(max(abs(as.numeric(value$fit) - expected_mean)), 1e-8)
        expect_lt(max(abs(as.numeric(value$se.fit) - expected_se)), 1e-6)
        expect_lt(max(abs(as.numeric(value$se.fit)^2 - expected_se^2)), 1e-6)
        expect_true(all(is.finite(value$se.fit)) && all(value$se.fit >= 0))
        contract <- attr(value, "prediction_contract")
        expect_identical(contract$native_components$response_outputs, 1L)
        expect_identical(contract$native_components$response_quantities, "overall_expected_count")
        expect_identical(contract$uncertainty$response_standard_error_method,
          "delta_full_native_coefficient_covariance")
        expect_false(contract$uncertainty$native_provider_response_SE_used)
        expect_false(contract$uncertainty$observation_noise)
        expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
      }
    }
    # Native controls select the covariance before the family hook sees it.
    variant <- model; variant$Vc <- 1.25 * model$Vp
    nd <- dat[c(1L, 17L, 24L), , drop = FALSE]
    G <- zip_coefficient_gradient(model, nd, dat, .5, list(exclude = "mean"))
    value <- .msr_predict_gam(variant, newdata = nd, type = "response", se.fit = TRUE,
      unconditional = TRUE, exclude = "mean", block.size = 1)
    expect_equal(as.numeric(value$fit), zip_native_mean(model, nd, dat,
      list(exclude = "mean")), tolerance = 1e-8)
    expect_equal(as.numeric(value$se.fit), sqrt(rowSums((G %*% variant$Vc) * G)), tolerance = 1e-6)
    expect_identical(model, original)
  }
  expect_identical(body(native), saved_body)
  expect_identical(environment(native), saved_environment)
})

test_that("scalar count means survive snapshots, checked rasters and raw-reference curves", {
  skip_if_not_installed("mgcv")
  f <- zip_response_fixture(); dat <- f$data; model <- f$models$plain$model
  object <- f$f$object; object$opt_mod <- model
  nd <- dat[c(1L, 17L, 24L), , drop = FALSE]
  G <- zip_coefficient_gradient(model, nd, dat, .5)
  expected_se <- sqrt(rowSums((G %*% model$Vp) * G))
  snapshot <- prepare_prediction(object, nd, f$f$inputs, include_fixed = TRUE, component = "response")
  snapshot <- unserialize(serialize(snapshot, NULL))
  expect_null(snapshot$opt_mod)
  value <- predict(snapshot, se.fit = TRUE, uncertainty = "conditional")
  expect_equal(as.numeric(value$fit), zip_native_mean(model, nd, dat), tolerance = 1e-8)
  expect_equal(as.numeric(value$se.fit), expected_se, tolerance = 1e-6)
  grid <- kernel_scale.raster(f$f$rasters, multiScaleR = object, include_fixed = TRUE,
    scale_center = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  map_data <- dat; map_data$z <- .2
  G <- zip_coefficient_gradient(model, map_data, dat, .5)
  map <- predict(object, grid, type = "response", se.fit = TRUE, uncertainty = "conditional")
  expected <- cbind(zip_native_mean(model, map_data, dat), sqrt(rowSums((G %*% model$Vp) * G)))
  expect_equal(unname(as.matrix(terra::extract(map, f$f$points, ID = FALSE))),
    unname(expected), tolerance = 1e-6)
  reference <- dat; row.names(reference) <- row.names(model$model)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot_marginal_effects(object, reference_data = reference, link = TRUE,
    component = 1, length.out = 5)
  curve_data <- dat[rep(1L, 5), , drop = FALSE]
  for (name in names(curve_data)) if (is.numeric(dat[[name]])) curve_data[[name]] <- mean(dat[[name]])
  curve_data$z <- seq(min(dat$z), max(dat$z), length.out = 5)
  G <- zip_coefficient_gradient(model, curve_data, dat, .5)
  fit <- zip_native_mean(model, curve_data, dat)
  se <- sqrt(rowSums((G %*% model$Vp) * G))
  expect_equal(curves$z$data$fit, fit, tolerance = 1e-8)
  expect_equal(curves$z$data$lwr, fit + qnorm(.025) * se, tolerance = 1e-6)
  expect_equal(curves$z$data$upr, fit + qnorm(.975) * se, tolerance = 1e-6)
  expect_match(curves$z$labels$caption, "both predictors")
  expect_identical(attr(curves, "prediction_contract")$native_output,
    list(type = "response", column = 1))
  expect_error(plot_marginal_effects(object, reference_data = reference,
    link = TRUE, component = 2), "one overall expected-count")
  expect_error(plot_marginal_effects(object, reference_data = reference, link = TRUE),
    "explicit native output")
  link <- plot_marginal_effects(object, reference_data = reference,
    link = FALSE, component = 1, length.out = 5)
  expect_true(length(link) > 0)
  native_link <- stats::predict(model, dat, type = "link")
  value <- predict(object, dat, f$f$inputs, include_fixed = TRUE, type = "link")
  expect_equal(as.numeric(value), as.numeric(native_link), tolerance = 1e-8)
  expect_identical(dim(value), c(nrow(dat), 2L))
})

test_that("delta algebra supports shared coefficients and guards native cancellation", {
  skip_if_not_installed("mgcv")
  family <- mgcv::ziplss(); native <- family$predict
  corrected <- .msr_gam_ziplss_predictor(native)
  X <- rbind(c(1, .1, -.2, .4, .2), c(1, -.3, .2, -.1, .4), c(1, .2, .5, .3, -.2))
  attr(X, "lpi") <- list(1:3, c(1L, 4L, 5L))
  beta <- c(.2, .3, -.1, .1, -.2); off <- list(c(-.2, .1, .3), c(.1, .3, -.4))
  G <- matrix(0, 3, 5)
  for (j in 1:5) {
    plus <- minus <- beta; plus[j] <- plus[j] + 1e-5; minus[j] <- minus[j] - 1e-5
    G[, j] <- (native(family, X = X, beta = plus, off = off)[[1]] -
      native(family, X = X, beta = minus, off = off)[[1]]) / 2e-5
  }
  for (cross in c(-.05, .05)) {
    V <- diag(c(.3, .2, .4, .1, .5)); V[1, 5] <- V[5, 1] <- cross
    value <- corrected(family, se = TRUE, X = X, beta = beta, off = off, Vb = V)
    expect_identical(value$fit, native(family, X = X, beta = beta, off = off)[[1]])
    expect_equal(as.numeric(value$se.fit), sqrt(rowSums((G %*% V) * G)), tolerance = 1e-6)
  }
  X <- matrix(1, 1, 2); attr(X, "lpi") <- list(1L, 2L)
  for (offset in list(list(-40, 0), list(40, 0), list(0, Inf), list(0, -Inf))) {
    value <- corrected(family, se = TRUE, X = X, beta = c(0, 0), off = offset, Vb = diag(c(.2, .3)))
    expect_identical(value$fit, native(family, X = X, beta = c(0, 0), off = offset)[[1]])
    expect_true(all(is.finite(value$se.fit)) && all(value$se.fit >= 0))
    if (identical(offset, list(-40, 0)))
      expect_equal(as.numeric(value$se.fit), sqrt(.3) * exp(-1), tolerance = 1e-12)
    if (identical(offset, list(40, 0))) {
      limit <- exp(40) * sqrt(.2 * (1 - exp(-1))^2 + .3 * exp(-2))
      expect_lt(abs(as.numeric(value$se.fit) / limit - 1), 1e-12)
    }
    if (identical(offset, list(0, -Inf))) expect_equal(as.numeric(value$se.fit), 0)
  }
  expect_error(corrected(family, se = TRUE, X = X, beta = c(0, 0),
    off = list(40, -40), Vb = diag(2)), "response mean")
  expect_error(corrected(family, se = TRUE, X = X, beta = c(0, 0),
    off = list(720, 0), Vb = diag(2)), "response mean")
  expect_error(corrected(family, se = TRUE, X = X, beta = c(0, 0),
    Vb = -diag(2)), "nonnegative variance")
  model <- zip_response_fixture()$models$plain$model
  model$family$predict <- function(...) list(0)
  expect_error(.msr_predict_gam(model, newdata = zip_response_fixture()$data,
    type = "response", se.fit = TRUE), "custom ziplss predictor")
})

test_that("response covariance provenance stays within the qualified provider", {
  skip_if_not_installed("mgcv")
  f <- zip_response_fixture(); model <- f$models$plain$model
  expect_true(.msr_gam_ziplss_model(model))
  inherited <- model
  class(inherited) <- c("msr_unregistered_gam_fixture", class(model))
  expect_true(.msr_gam_ziplss_model(inherited))
  # Dispatch-only metadata probe, not a fitted BAM or BAM qualification.
  distinct <- model; class(distinct) <- c("bam", class(model))
  expect_false(.msr_gam_ziplss_model(distinct))
  object <- f$f$object; object$opt_mod <- distinct
  contract <- .msr_native_prediction_contract(object, "gam", "response", list(), se.fit = FALSE)
  expect_null(contract$uncertainty$response_standard_error_method)
  expect_null(contract$uncertainty$native_provider_response_SE_used)
  expect_null(contract$native_components$response_outputs)
})

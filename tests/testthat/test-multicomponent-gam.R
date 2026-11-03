# Native family outputs and basis geometry, using the established spatial data.
# No recovery study or outer ecological-scale optimization is run here.
test_that("multiple GAM predictors keep native fixed bases, labels and uncertainty", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture(); dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  dat$category <- (seq_len(nrow(dat)) %/% 3L) %% 3L
  make <- function(kind, qualified) {
    sc <- if (qualified) "base::scale(z)" else "scale(z)"
    pol <- if (qualified) "stats::poly(exposure,2)" else "poly(exposure,2)"
    first <- as.formula(paste(if (kind == "gaulss") "y ~" else "category ~",
      "rough * group +", sc, "+ s(mean,k=4)"), env = asNamespace("mgcv"))
    second <- as.formula(paste("~ fixed +", sc, "+", pol), env = asNamespace("mgcv"))
    mgcv::gam(list(first, second), data = dat,
      family = if (kind == "gaulss") mgcv::gaulss() else mgcv::multinom(K = 2),
      method = "REML")
  }
  compare <- function(actual, expected) {
    if (is.list(expected)) {
      expect_identical(names(actual), names(expected))
      for (nm in names(expected)) compare(actual[[nm]], expected[[nm]])
    } else {
      expect_identical(dim(actual), dim(expected))
      expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-8)
      expect_true(all(is.finite(actual)))
    }
  }
  for (kind in c("gaulss", "multinom")) {
    bare <- make(kind, FALSE); model <- make(kind, TRUE); original <- model
    fixed <- .msr_native_fixed_metadata(model)
    expect_identical(model, original)
    expect_identical(coef(fixed), coef(model))
    expect_identical(fitted(fixed), fitted(model))
    expect_identical(logLik(fixed), logLik(model))
    expect_identical(fixed$Vp, model$Vp)
    expect_identical(fixed$Vc, model$Vc)
    expect_identical(fixed$nsdf, model$nsdf)
    expect_identical(fixed$outer.info, model$outer.info)
    expect_identical(fixed$outer.info$conv, "full convergence")
    expect_true(all(eigen(fixed$outer.info$hess, symmetric = TRUE)$values > 0))
    expect_identical(attributes(fixed$pterms), attributes(model$pterms))
    expect_identical(.msr_native_fixed_metadata(fixed), fixed)
    expect_setequal(.model_predictors(model), c("rough", "group", "z", "mean", "fixed", "exposure"))
    object <- f$object; object$opt_mod <- model
    expect_equal(.msr_selection_parameter_count(object),
      as.numeric(attr(logLik(model), "df")) + 2L)
    rows <- dat[c(5, 23, 17), , drop = FALSE]
    grids <- list(dat, rows, rows[1, , drop = FALSE], transform(rows, z = .2, exposure = 1.25))
    for (nd in grids) for (type in c("response", "link", "terms", "iterms", "lpmatrix")) {
      # mgcv explicitly warns and uses terms for multi-predictor iterms.
      native <- function(se = FALSE, ...) {
        result <- function() predict(bare, nd, type = type, se.fit = se, ...)
        if (type == "iterms") {
          value <- NULL
          expect_warning(value <- result(), "iterms not available")
          value
        } else result()
      }
      wrapper <- function(se = FALSE, ...) {
        result <- function() predict(object, nd, f$inputs, include_fixed = TRUE,
          type = type, se.fit = se, uncertainty = "conditional", ...)
        if (type == "iterms") {
          value <- NULL
          expect_warning(value <- result(), "iterms not available")
          value
        } else result()
      }
      compare(wrapper(), native())
      if (type != "lpmatrix") compare(wrapper(TRUE), native(TRUE))
      if (type == "link") compare(wrapper(TRUE, unconditional = TRUE),
        native(TRUE, unconditional = TRUE))
      if (type == "lpmatrix") {
        expect_identical(attr(wrapper(), "lpi"), attr(predict(model, dat, type = "lpmatrix"), "lpi"))
        expect_identical(colnames(wrapper()), names(coef(model)))
      }
    }
    expect_error(predict(object, rows, f$inputs, se.fit = TRUE), "conditional")
    for (variable in c("group", "exposure")) {
      absent <- rows; absent[[variable]] <- NULL
      expect_error(predict(object, absent, f$inputs, include_fixed = TRUE), "formula variables")
    }
    selected <- attr(model$pterms, "term.labels")[1L]
    for (control in list(list(terms = selected), list(exclude = selected),
        list(terms = selected, unconditional = TRUE, block.size = 1))) {
      for (type in c("terms", "link")) {
        expected <- do.call(predict, c(list(object = bare, newdata = rows,
          type = type, se.fit = TRUE), control))
        actual <- predict(object, rows, f$inputs, include_fixed = TRUE, type = type,
          se.fit = TRUE, uncertainty = "conditional", native_args = control)
        compare(actual, expected)
      }
    }
    for (type in c("response", "link", "terms", "lpmatrix")) {
      controls <- if (type == "terms") list(terms = selected, unconditional = TRUE) else list()
      snap <- prepare_prediction(object, rows, f$inputs, include_fixed = TRUE,
        component = type, native_args = controls)
      restored <- unserialize(serialize(snap, NULL))
      expect_null(restored$opt_mod)
      compare(predict(restored), do.call(predict, c(list(object = bare,
        newdata = rows, type = type), controls)))
      contract <- restored$prediction_contract[[type]]
      expect_identical(contract$native_components$linear_predictors, 2L)
      expect_identical(contract$native_components$family, model$family$family)
      expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
      if (type == "response") expect_identical(contract$quantity, "native_component")
      if (type != "lpmatrix") compare(predict(restored, se.fit = TRUE, uncertainty = "conditional"),
        do.call(predict, c(list(object = bare, newdata = rows, type = type, se.fit = TRUE), controls)))
      else expect_error(predict(restored, se.fit = TRUE, uncertainty = "conditional"), "se.fit")
    }
    grid <- kernel_scale.raster(f$rasters, multiScaleR = object,
      scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE,
      site_values = list(group = "b", z = .2, exposure = 1.25), verbose = FALSE)
    reference <- dat; reference$group <- factor("b", levels = levels(dat$group))
    reference$z <- .2; reference$exposure <- 1.25
    for (type in c("response", "link")) {
      expected <- predict(bare, reference, type = type, se.fit = TRUE)
      map <- predict(object, grid, type = type, se.fit = TRUE, uncertainty = "conditional")
      expect_equal(unname(as.matrix(terra::extract(map, f$points, ID = FALSE))),
        unname(cbind(expected$fit, expected$se.fit)), tolerance = 1e-8)
      expect_equal(terra::nlyr(map), if (type == "response" && kind == "multinom") 6L else 4L)
    }
    expect_error(plot_marginal_effects(object), "explicit native output selection")
    missing_frame <- model; missing_frame$model <- NULL
    expect_error(.msr_native_fixed_metadata(missing_frame), "retained evaluated training frame")
    bad <- model; attr(bad$pterms, "term.labels") <- NULL
    expect_error(.msr_native_fixed_metadata(bad), "align native GAM component")
  }
})

test_that("native two-predictor GAM refits discover both components without changing the penalty", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture(); dat <- f$data
  dat$y <- dat$y + as.integer(dat$site)^2
  model <- mgcv::gam(list(y ~ mean + rough + group + s(z, k = 4), ~ fixed + exposure),
    data = dat, family = mgcv::gaulss(), method = "REML")
  context <- build_opt_context(model, f$inputs$raw_cov,
    scale_vars = f$inputs$scale_vars, unit_conv = f$inputs$unit_conv,
    resolution = f$inputs$resolution, n_cols = f$inputs$n_cols,
    binned = f$inputs$binned, d_list = f$inputs$d_list)
  expect_setequal(context$covs, c("mean", "rough"))
  expect_setequal(.model_predictors(model), c("mean", "rough", "group", "z", "fixed", "exposure"))
  actual <- kernel_scale_fn(f$object$optim_results$par, f$inputs$d_list,
    f$inputs$raw_cov, f$inputs$kernel, model, opt_context = context, mod_return = TRUE)
  expect_equal(as.numeric(logLik(actual$mod)), as.numeric(logLik(model)), tolerance = 1e-8)
  expect_equal(as.numeric(coef(actual$mod)), as.numeric(coef(model)), tolerance = 1e-8)
  expect_identical(attr(actual$mod$pterms, "term.labels"), attr(model$pterms, "term.labels"))
})

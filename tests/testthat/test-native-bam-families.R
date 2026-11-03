# Existing36-site weighted-SQ fixtures, actual native binomial/Poisson modes.
# Cache failures as well as successes; no retry or changed fitting settings.
native_bam_family_fixture <- local({
  cache <- NULL
  function() {
    if (inherits(cache, "error")) stop(cache)
    if (!is.null(cache)) return(cache)
    cache <<- tryCatch({
      old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
      cases <- list()
      for (kernel in c("gaussian", "exp", "expow", "fixed")) {
        f <- refit_spatial_fixture(kernel = kernel)
        stopifnot(isTRUE(f$object$opt_mod$converged))
        for (family in c("binomial", "poisson")) for (mode in c("ordinary", "discrete")) {
          dat <- f$data
          if (family == "binomial") {
            dat$successes <- round(20 * plogis(.15 + .35 * dat$mean - .2 * dat$rough + .1 * dat$z))
            dat$failures <- 20 - dat$successes
            lhs <- "cbind(successes,failures)"; native_family <- binomial()
          } else {
            dat$y <- round(exp(1 + .35 * dat$mean - .2 * dat$rough + .1 * dat$z))
            lhs <- "y"; native_family <- poisson()
          }
          form <- as.formula(paste(lhs, "~ rough+z+z:group+I(z^2)+s(mean,k=4)"), env = asNamespace("mgcv"))
          model <- mgcv::bam(form, data = dat, family = native_family, method = "fREML",
            discrete = mode == "discrete", nthreads = 1L, cluster = NULL)
          stopifnot(isTRUE(model$converged), all(is.finite(coef(model))), all(is.finite(model$sp)),
            length(model$sp) == 1L, length(model$smooth) == 1L, length(model$smooth[[1L]]$S) == 1L,
            identical(model$family$family, family), isTRUE(all.equal(as.numeric(model$sig2), 1, tolerance = 0)))
          conv <- model$mgcv.conv
          if (is.logical(conv)) stopifnot(isTRUE(conv))
          if (is.list(conv) && !is.null(conv$message)) stopifnot(identical(conv$message, "full convergence"))
          if (is.list(conv)) for (flag in intersect(c("fully.converged", "hess.pos.def"), names(conv))) stopifnot(isTRUE(conv[[flag]]))
          if (!is.null(model$outer.info$conv)) stopifnot(identical(model$outer.info$conv, "full convergence"))
          if (!is.null(model$outer.info$hess)) stopifnot(min(eigen(model$outer.info$hess, symmetric = TRUE, only.values = TRUE)$values) > 0)
          X <- predict(model, dat, type = "lpmatrix", discrete = FALSE, n.threads = 1L)
          S <- matrix(0, ncol(X), ncol(X)); sm <- model$smooth[[1L]]; ix <- sm$first.para:sm$last.para
          S[ix, ix] <- model$sp[[1L]] * sm$S[[1L]]
          mu <- as.numeric(model$fitted.values); eta <- as.numeric(model$linear.predictors)
          expected_weight <- if (family == "binomial") 20 else 1
          stopifnot(length(model$y) == nrow(dat), length(model$prior.weights) == nrow(dat),
            all(is.finite(model$y)), all(model$prior.weights == expected_weight))
          derivative <- model$family$mu.eta(eta); variance <- model$family$variance(mu)
          weights <- model$prior.weights * derivative^2 / variance
          score <- drop(crossprod(X, model$prior.weights * (model$y - mu) * derivative / variance) - S %*% coef(model))
          stopifnot(all(is.finite(weights)), all(weights > 0), all(is.finite(score)), max(abs(score)) < 1e-4,
            min(eigen(crossprod(X, X * weights) + S, symmetric = TRUE, only.values = TRUE)$values) > 0,
            min(eigen(model$Vp, symmetric = TRUE, only.values = TRUE)$values) > 0)
          object <- f$object; object$opt_mod <- model
          cases[[paste(kernel, family, mode, sep = "_")]] <- list(model = model, object = object, data = dat, f = f)
        }
      }
      stopifnot(length(cases) == 16L)
      cases
    }, error = identity)
    if (inherits(cache, "error")) stop(cache)
    cache
  }
})

# Shared prediction helpers are loaded from helper-native-bam.R.

test_that("binomial and Poisson BAM points and snapshots preserve native formulas", {
  skip_if_not_installed("mgcv")
  cases <- native_bam_family_fixture()
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  for (entry in cases) {
    model <- entry$model; original <- model
    for (type in c("link", "response", "terms", "iterms", "lpmatrix")) for (rows in
        list(seq_len(nrow(entry$data)), c(1L, 17L, 24L), 17L)) {
      nd <- entry$data[rows, , drop = FALSE]
      actual <- predict(entry$object, nd, entry$f$inputs, include_fixed = TRUE, type = type,
        se.fit = type != "lpmatrix", uncertainty = "conditional", discrete = FALSE, n.threads = 1L)
      expect_equal(bam_family_clean(actual), bam_family_clean(bam_family_oracle(entry, nd, type)), tolerance = 1e-8)
      expect_identical(attr(actual, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
    }
    nd <- entry$data[c(1L, 17L, 24L), , drop = FALSE]
    for (type in c("response", "iterms", "lpmatrix")) {
      snapshot <- prepare_prediction(entry$object, nd, entry$f$inputs, include_fixed = TRUE,
        component = type, native_args = list(discrete = FALSE, exclude = "s(mean)"))
      actual <- predict(unserialize(serialize(snapshot, NULL)), type = type,
        se.fit = type != "lpmatrix", uncertainty = "conditional")
      expect_equal(bam_family_clean(actual), bam_family_clean(bam_family_oracle(entry, nd, type,
        list(exclude = "s(mean)"))), tolerance = 1e-8)
    }
    expect_error(predict(entry$object, nd, entry$f$inputs, se.fit = TRUE, uncertainty = "joint"), "conditional")
    expect_identical(model, original)
  }
})

test_that("binomial and Poisson BAM maps and curves retain response semantics", {
  skip_if_not_installed("mgcv")
  cases <- native_bam_family_fixture()
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  for (entry in cases) {
    model <- entry$model; original <- model; dat <- entry$data; object <- entry$object
    vars <- intersect(all.vars(delete.response(terms(model))), names(dat)); source <- dat[vars]
    grid <- kernel_scale.raster(entry$f$rasters, multiScaleR = object, include_fixed = TRUE,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(group = "b", z = .2), verbose = FALSE)
    nd <- source; nd$group <- factor("b", levels = levels(dat$group)); nd$z <- .2
    for (type in c("response", "iterms")) {
      target <- bam_family_oracle(entry, nd, type)
      target <- as.matrix(data.frame(fit = target$fit, se.fit = target$se.fit))
      map <- predict(object, grid, type = type, se.fit = TRUE, uncertainty = "conditional",
        discrete = FALSE, raster_args = list(cores = 1L))
      sampled <- as.matrix(terra::extract(map, entry$f$points, ID = FALSE))
      expect_equal(bam_family_clean(sampled), bam_family_clean(target), tolerance = 1e-8)
    }
    reference <- dat; row.names(reference) <- row.names(model$model)
    view <- object; view$opt_mod$call$data <- quote(expired_native_bam_family_reference)
    # Single-predictor BAM inherits glm's automatic response transformation.
    # Its existing link argument differs from the multi-component GAM switch.
    for (link in c(FALSE, TRUE)) {
      plots <- suppressMessages(plot_marginal_effects(view, length.out = 7L,
        reference_data = reference, link = link))
      expect_true(attr(plots, "reference_data_contract")$native_design_and_offsets_verified)
      for (var in c("mean", "rough", "z")) {
        nd <- as.data.frame(lapply(source, function(x) if (is.numeric(x)) rep(mean(x), 7L) else x[rep(1L, 7L)]), check.names = FALSE)
        xs <- seq(min(dat[[var]]), max(dat[[var]]), length.out = 7L); nd[[var]] <- xs
        native <- bam_family_oracle(entry, nd, "link")
        target <- cbind(fit = as.numeric(native$fit), lwr = as.numeric(native$fit + qnorm(.025) * native$se.fit),
          upr = as.numeric(native$fit + qnorm(.975) * native$se.fit))
        target[] <- model$family$linkinv(target)
        axis <- if (var %in% names(object$scl_params$mean)) xs * object$scl_params$sd[[var]] + object$scl_params$mean[[var]] else xs
        expect_equal(bam_family_clean(as.matrix(plots[[var]]$data[c("fit", "lwr", "upr")])),
          bam_family_clean(target), tolerance = 1e-8)
        expect_equal(plots[[var]]$data$x, axis, tolerance = 1e-10, ignore_attr = TRUE)
      }
    }
    expect_identical(model, original)
  }
})

test_that("binomial and Poisson BAM tables preserve family and native penalty counts", {
  skip_if_not_installed("mgcv")
  for (entry in native_bam_family_fixture()) {
    model <- entry$model; models <- list(entry$object, model); labels <- c("wrapped", "native")
    quiet <- capture.output(aic <- aic_tab(models, AICc = FALSE, mod_names = labels, verbose = FALSE))
    quiet2 <- capture.output(bic <- bic_tab(models, mod_names = labels, verbose = FALSE))
    auto <- aic_tab(models, AICc = FALSE, verbose = FALSE)
    extra <- if (entry$f$inputs$kernel == "expow") 4L else 2L
    count <- as.numeric(attr(logLik(model), "df")) + c(extra, 0)
    ll <- as.numeric(logLik(model)); n <- as.numeric(nobs(model))
    ia <- match(labels, aic$Modnames); ib <- match(labels, bic$Modnames)
    expect_equal(aic$K[ia], count); expect_equal(bic$K[ib], count)
    expect_equal(aic$AIC[ia], -2 * ll + 2 * count, tolerance = 1e-8)
    expect_equal(bic$BIC[ib], -2 * ll + log(n) * count, tolerance = 1e-8)
    expect_true(all(grepl(paste0("family=", model$family$family), auto$Modnames, fixed = TRUE)))
    expect_identical(c(quiet, quiet2), character()); expect_equal(n, 36)
  }
})

test_that("non-Gaussian BAM workers retain native parallel shapes and serial values", {
  skip_if_not_installed("mgcv")
  cases <- native_bam_family_fixture()
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  pool <- parallel::makePSOCKcluster(2L); on.exit(parallel::stopCluster(pool), add = TRUE)
  configure <- function(coding) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    options(contrasts = c(coding, "contr.poly")); TRUE
  }
  environment(configure) <- baseenv()
  flatten <- function(x) if (is.list(x)) unlist(lapply(x, flatten), use.names = FALSE) else as.numeric(x)
  for (entry in cases) {
    nd <- entry$data[rep(seq_len(nrow(entry$data)), 6L), , drop = FALSE]
    nd$case_site <- rep(seq_len(nrow(entry$data)), 6L); row.names(nd) <- paste0("bam-family-", seq_len(nrow(nd)))
    keys <- data.frame(case_site = seq_len(nrow(entry$data)))
    vars <- intersect(all.vars(delete.response(terms(entry$model))), names(nd))
    for (type in c("link", "response", "terms", "iterms", "lpmatrix")) {
      serial <- bam_family_oracle(entry, nd, type)
      if (type == "iterms") target <- serial else {
        invisible(parallel::clusterCall(pool, configure, "contr.sum"))
        target <- predict(entry$model, nd[vars], type = type, se.fit = type != "lpmatrix",
          discrete = FALSE, cluster = pool, n.threads = 1L)
      }
      invisible(parallel::clusterCall(pool, configure, "contr.treatment"))
      options(contrasts = c("contr.treatment", "contr.poly"))
      actual <- predict(entry$object, nd, entry$f$inputs, join_by = keys, include_fixed = TRUE,
        type = type, se.fit = type != "lpmatrix", uncertainty = "conditional", cluster = pool, discrete = FALSE)
      options(contrasts = c("contr.sum", "contr.poly"))
      expect_equal(bam_family_clean(actual), bam_family_clean(target), tolerance = 1e-8)
      expect_equal(flatten(bam_family_clean(actual)), flatten(serial), tolerance = 1e-8)
    }
  }
})

test_that("optimized binomial and Poisson BAM predictions retain exact native batches", {
  skip_if_not_installed("mgcv")
  # Reuse the existing failure-retaining fixture, including every fit gate.
  cases <- native_bam_family_fixture()
  cases <- cases[grepl("_discrete$", names(cases))]
  expect_length(cases, 8L)
  old <- options(contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old), add = TRUE)
  coding <- function(expr) {
    before <- options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(before), add = TRUE)
    force(expr)
  }
  capture <- function(expr) {
    conditions <- list()
    value <- withCallingHandlers(force(expr), warning = function(w) {
      conditions[[length(conditions) + 1L]] <<- list(message = conditionMessage(w), classes = class(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = conditions)
  }
  clean <- function(x) {
    attr(x, "prediction_contract") <- NULL
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.list(x)) {
      result <- lapply(x, clean); attributes(result) <- attributes(x)
      return(result)
    }
    if (is.matrix(x)) rownames(x) <- NULL else names(x) <- NULL
    x
  }
  types <- c("link", "response", "terms", "iterms", "lpmatrix")
  native_callback <- .msr_native_raster_eval
  for (entry in cases) {
    model <- entry$model; before <- model; object <- entry$object
    variables <- intersect(all.vars(delete.response(terms(model))), names(entry$data))
    for (type in types) for (rows in list(seq_len(nrow(entry$data)), c(1L, 17L, 24L), 17L)) {
      nd <- entry$data[rows, , drop = FALSE]
      expected <- coding(capture(predict(model, nd[variables], type = type,
        se.fit = type != "lpmatrix", discrete = TRUE, n.threads = 1L)))
      actual <- capture(predict(object, nd, entry$f$inputs, include_fixed = TRUE,
        type = type, se.fit = type != "lpmatrix", uncertainty = "conditional",
        discrete = TRUE, n.threads = 1L))
      expect_equal(clean(actual$value), clean(expected$value), tolerance = 1e-8)
      expect_identical(actual$warnings, expected$warnings)
    }
    nd <- entry$data[c(1L, 17L, 24L), , drop = FALSE]
    for (type in types) {
      snapshot <- prepare_prediction(object, nd, entry$f$inputs, include_fixed = TRUE,
        component = type, native_args = list(discrete = TRUE))
      expected <- coding(capture(predict(model, nd[variables], type = type,
        se.fit = type != "lpmatrix", discrete = TRUE, n.threads = 1L)))
      actual <- capture(predict(unserialize(serialize(snapshot, NULL)), type = type,
        se.fit = type != "lpmatrix", uncertainty = "conditional"))
      expect_equal(clean(actual$value), clean(expected$value), tolerance = 1e-8)
      expect_identical(actual$warnings, expected$warnings)
    }
    grid <- kernel_scale.raster(entry$f$rasters, multiScaleR = object, include_fixed = TRUE,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(group = "b", z = .2), verbose = FALSE)
    batches <- list()
    checked_callback <- function(model, data) {
      seen <- list()
      value <- withCallingHandlers(native_callback(model, data), warning = function(w)
        seen[[length(seen) + 1L]] <<- list(message = conditionMessage(w), classes = class(w)))
      # Discrete predictions can depend on their batch. Compare the exact rows
      # received by the original callback, rather than a substituted grid.
      expected <- coding(capture(do.call(stats::predict, c(list(object = model$object$opt_mod,
        newdata = data, type = model$type, se.fit = model$se.fit), model$native_args))))
      target <- expected$value
      if (is.list(target) && !is.data.frame(target) && !is.null(target$fit))
        target <- data.frame(fit = target$fit, se.fit = target$se.fit)
      if (is.matrix(target) && is.object(target)) target <- unclass(target)
      batches[[length(batches) + 1L]] <<- list(rows = nrow(data),
        equal = isTRUE(all.equal(clean(value), clean(target), tolerance = 1e-8)),
        warnings_equal = identical(seen, expected$warnings))
      value
    }
    for (type in types) {
      batches <- list()
      checked <- with_mocked_bindings(capture(predict(object, grid, type = type,
        se.fit = type != "lpmatrix", uncertainty = "conditional", discrete = TRUE,
        raster_args = list(cores = 1L))), .msr_native_raster_eval = checked_callback)
      actual <- capture(predict(object, grid, type = type, se.fit = type != "lpmatrix",
        uncertainty = "conditional", discrete = TRUE, raster_args = list(cores = 1L)))
      expect_true(length(batches) > 0L)
      expect_true(all(vapply(batches, function(x) x$rows > 0L && x$equal, logical(1))))
      expect_true(all(vapply(batches, `[[`, logical(1), "warnings_equal")))
      expect_equal(terra::values(actual$value), terra::values(checked$value), tolerance = 1e-8)
      expect_true(isTRUE(terra::compareGeom(actual$value, checked$value)))
      expect_identical(actual$warnings, checked$warnings)
    }
    expect_identical(.msr_native_raster_eval, native_callback)
    expect_identical(model, before)
  }
})

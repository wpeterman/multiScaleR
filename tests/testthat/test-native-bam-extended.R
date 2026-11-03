# Fixed36-site software cases. Quantiles are declared once, with no random
# draws, retry, seed search or ecological-scale optimization.
native_bam_extended_fixture_impl <- function() {
  cases <- list(); fit_calls <- preparation_calls <- 0L
  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    f <- refit_spatial_fixture(kernel = kernel)
    preparation_calls <- preparation_calls + 1L
    stopifnot(isTRUE(f$object$opt_mod$converged))
    kinds <- c("Gamma", "nb_estimated", if (kernel == "gaussian") "nb_fixed")
    for (kind in kinds) for (mode in c("ordinary", "discrete")) {
      dat <- f$data; n <- nrow(dat)
      probability <- ((seq_len(n) * 11L) %% 37L) / 37
      mean <- exp(1 + .35 * dat$mean - .2 * dat$rough + .1 * dat$z)
      if (kind == "Gamma") {
        dat$y <- mean * qgamma(probability, shape = 3, rate = 3)
        family <- Gamma(link = "log")
      } else {
        dat$y <- qnbinom(probability, size = 2, mu = mean)
        family <- mgcv::nb(theta = if (kind == "nb_fixed") 2 else NULL)
      }
      form <- as.formula("y ~ rough+z+z:group+I(z^2)+s(mean,k=4)", env = asNamespace("mgcv"))
      key <- paste(kernel, kind, mode, sep = "_")
      fit_calls <- fit_calls + 1L
      cases[[key]] <- tryCatch({
        model <- mgcv::bam(form, data = dat, family = family, method = "fREML",
          discrete = mode == "discrete", nthreads = 1L, cluster = NULL,
          control = if (mode == "ordinary" && kind %in% c("Gamma", "nb_fixed"))
            mgcv::gam.control(epsilon = 1e-10) else mgcv::gam.control())
        stopifnot(isTRUE(model$converged), all(is.finite(coef(model))), all(is.finite(model$sp)),
          length(model$sp) == 1L, length(model$smooth) == 1L, length(model$smooth[[1L]]$S) == 1L)
        conv <- model$mgcv.conv
        if (is.logical(conv)) stopifnot(isTRUE(conv))
        if (is.list(conv) && !is.null(conv$message)) stopifnot(identical(conv$message, "full convergence"))
        if (is.list(conv)) for (flag in intersect(c("fully.converged", "hess.pos.def"), names(conv))) stopifnot(isTRUE(conv[[flag]]))
        if (!is.null(model$outer.info$conv)) stopifnot(identical(model$outer.info$conv, "full convergence"))
        if (!is.null(model$outer.info$hess)) stopifnot(min(eigen(model$outer.info$hess, symmetric = TRUE, only.values = TRUE)$values) > 0)
        phi <- as.numeric(model$sig2)
        stopifnot(length(phi) == 1L, is.finite(phi), phi > 0)
        if (kind == "Gamma") stopifnot(isTRUE(model$scale.estimated)) else
          stopifnot(identical(phi, 1))
        theta <- if (kind == "Gamma") NULL else model$family$getTheta(TRUE)
        if (!is.null(theta)) {
          stopifnot(length(theta) == 1L, is.finite(theta), theta > 0,
            model$family$n.theta == if (kind == "nb_fixed") 0 else 1)
          if (kind == "nb_fixed") stopifnot(identical(as.numeric(theta), 2))
        }
        X <- predict(model, dat, type = "lpmatrix", discrete = FALSE, n.threads = 1L)
        S <- matrix(0, ncol(X), ncol(X)); sm <- model$smooth[[1L]]; ix <- sm$first.para:sm$last.para
        S[ix, ix] <- model$sp[[1L]] * sm$S[[1L]]
        mu <- as.numeric(model$fitted.values); eta <- as.numeric(model$linear.predictors)
        stopifnot(length(model$y) == n, all(is.finite(model$y)), all(model$prior.weights == 1))
        derivative <- model$family$mu.eta(eta); variance <- model$family$variance(mu)
        weights <- model$prior.weights * derivative^2 / variance
        score <- drop(crossprod(X, model$prior.weights * (model$y - mu) * derivative / variance) - S %*% coef(model))
        information <- crossprod(X, X * weights) + S
        minimum_information <- min(eigen(information, symmetric = TRUE, only.values = TRUE)$values)
        minimum_covariance <- min(eigen(model$Vp, symmetric = TRUE, only.values = TRUE)$values)
        stopifnot(all(is.finite(weights)), all(weights > 0), all(is.finite(score)), max(abs(score)) < 1e-4,
          is.finite(minimum_information), minimum_information > 0, is.finite(minimum_covariance), minimum_covariance > 0)
        # Match the captured native logLik.gam count, including corrected edf,
        # estimated scale, the coefficient ceiling and estimated family theta.
        scale_count <- as.numeric(model$scale.estimated)
        edf <- if (is.null(model$edf2)) model$edf else model$edf2
        expected_df <- min(sum(edf) + scale_count, length(coef(model)) + scale_count)
        if (inherits(model$family, "extended.family") && !is.null(model$family$n.theta))
          expected_df <- expected_df + model$family$n.theta
        stopifnot(isTRUE(all.equal(as.numeric(attr(logLik(model), "df")), expected_df, tolerance = 0)))
        object <- f$object; object$opt_mod <- model
        list(model = model, object = object, data = dat, f = f, kind = kind, mode = mode,
          diagnostics = list(converged = TRUE, penalized_coefficient_score = max(abs(score)),
            minimum_information = minimum_information, minimum_covariance = minimum_covariance,
            phi = phi, theta = theta, scale_parameter_count = scale_count,
            theta_parameter_count = if (is.null(theta)) 0L else model$family$n.theta,
            native_parameter_count = expected_df,
            actual_nuisance_parameter_count = if (kind == "Gamma") 1L else model$family$n.theta,
            actual_parameter_count = min(sum(edf), length(coef(model))) +
              if (kind == "Gamma") 1L else model$family$n.theta))
      }, error = function(error) list(error = conditionMessage(error), kind = kind, mode = mode))
    }
  }
  stopifnot(fit_calls == 18L, preparation_calls == 4L, length(cases) == 18L)
  attr(cases, "fit_calls") <- fit_calls
  attr(cases, "preparation_calls") <- preparation_calls
  cases
}

# Cache the declared lifecycle, including failures, across all tests in this file.
native_bam_extended_fixture <- local({
  cache <- NULL
  function() {
    if (inherits(cache, "error")) stop(cache)
    if (!is.null(cache)) return(cache)
    cache <<- tryCatch({
      cases <- native_bam_extended_fixture_impl()
      failed <- vapply(cases, function(x) !is.null(x$error), logical(1))
      if (any(failed)) stop("Unqualified native family fixture: ", paste(names(cases)[failed], collapse = ", "))
      cases
    }, error = identity)
    if (inherits(cache, "error")) stop(cache)
    cache
  }
})

bam_extended_capture <- function(expr) {
  conditions <- list()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    conditions[[length(conditions) + 1L]] <<- list(message = conditionMessage(w), classes = class(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = conditions)
}

bam_extended_native <- function(entry, data, type, optimized, controls = list()) {
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  if (!optimized) return(bam_extended_capture(bam_family_oracle(entry, data, type, controls)))
  variables <- intersect(all.vars(delete.response(terms(entry$model))), names(entry$data))
  bam_extended_capture(do.call(stats::predict, c(list(object = entry$model, newdata = data[variables],
    type = type, se.fit = type != "lpmatrix", discrete = TRUE, n.threads = 1L), controls)))
}

test_that("Gamma and NB BAM points and snapshots preserve native targets and warnings", {
  skip_if_not_installed("mgcv")
  cases <- native_bam_extended_fixture()
  expect_length(cases, 18L)
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  for (entry in cases) {
    original <- entry$model
    expect_lt(entry$diagnostics$penalized_coefficient_score, 1e-4)
    expect_gt(entry$diagnostics$minimum_information, 0)
    expect_gt(entry$diagnostics$minimum_covariance, 0)
    for (optimized in c(FALSE, if (entry$mode == "discrete") TRUE)) {
      for (type in c("link", "response", "terms", "iterms", "lpmatrix")) {
        for (rows in list(seq_len(nrow(entry$data)), c(1L, 17L, 24L), 17L)) {
          nd <- entry$data[rows, , drop = FALSE]
          expected <- bam_extended_native(entry, nd, type, optimized)
          actual <- bam_extended_capture(predict(entry$object, nd, entry$f$inputs,
            include_fixed = TRUE, type = type, se.fit = type != "lpmatrix",
            uncertainty = "conditional", discrete = optimized, n.threads = 1L))
          expect_equal(bam_family_clean(actual$value), bam_family_clean(expected$value), tolerance = 1e-8)
          expect_identical(actual$warnings, expected$warnings)
          expect_identical(attr(actual$value, "prediction_contract")$uncertainty$kernel_parameters,
            "fixed_at_estimate")
        }
        nd <- entry$data[c(1L, 17L, 24L), , drop = FALSE]
        controls <- list(discrete = optimized, exclude = "s(mean)")
        snapshot <- prepare_prediction(entry$object, nd, entry$f$inputs, include_fixed = TRUE,
          component = type, native_args = controls)
        expected <- bam_extended_native(entry, nd, type, optimized, list(exclude = "s(mean)"))
        actual <- bam_extended_capture(predict(unserialize(serialize(snapshot, NULL)), type = type,
          se.fit = type != "lpmatrix", uncertainty = "conditional"))
        expect_equal(bam_family_clean(actual$value), bam_family_clean(expected$value), tolerance = 1e-8)
        expect_identical(actual$warnings, expected$warnings)
      }
    }
    expect_identical(entry$model, original)
  }
})

test_that("Gamma and NB BAM maps use exact native callback batches", {
  skip_if_not_installed("mgcv")
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  callback <- .msr_native_raster_eval
  for (entry in native_bam_extended_fixture()) {
    grid <- kernel_scale.raster(entry$f$rasters, multiScaleR = entry$object, include_fixed = TRUE,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(group = "b", z = .2), verbose = FALSE)
    batches <- list()
    checked <- function(model, data) {
      seen <- list()
      value <- withCallingHandlers(callback(model, data), warning = function(w)
        seen[[length(seen) + 1L]] <<- list(message = conditionMessage(w), classes = class(w)))
      before <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(before), add = TRUE)
      native <- bam_extended_capture(do.call(stats::predict, c(list(object = model$object$opt_mod,
        newdata = data, type = model$type, se.fit = model$se.fit), model$native_args)))
      target <- native$value
      if (is.list(target) && !is.data.frame(target) && !is.null(target$fit))
        target <- data.frame(fit = target$fit, se.fit = target$se.fit)
      if (is.matrix(target) && is.object(target)) target <- unclass(target)
      batches[[length(batches) + 1L]] <<- list(
        equal = isTRUE(all.equal(bam_family_clean(value), bam_family_clean(target), tolerance = 1e-8)),
        warnings_equal = identical(seen, native$warnings))
      value
    }
    for (optimized in c(FALSE, if (entry$mode == "discrete") TRUE)) {
      for (type in c("link", "response", "terms", "iterms", "lpmatrix")) {
        batches <- list()
        observed <- with_mocked_bindings(bam_extended_capture(predict(entry$object, grid,
          type = type, se.fit = type != "lpmatrix", uncertainty = "conditional", discrete = optimized,
          raster_args = list(cores = 1L))), .msr_native_raster_eval = checked)
        plain <- bam_extended_capture(predict(entry$object, grid, type = type,
          se.fit = type != "lpmatrix", uncertainty = "conditional", discrete = optimized,
          raster_args = list(cores = 1L)))
        expect_true(length(batches) > 0L)
        expect_true(all(vapply(batches, `[[`, logical(1), "equal")))
        expect_true(all(vapply(batches, `[[`, logical(1), "warnings_equal")))
        expect_equal(terra::values(plain$value), terra::values(observed$value), tolerance = 1e-8)
        expect_true(isTRUE(terra::compareGeom(plain$value, observed$value)))
        expect_identical(plain$warnings, observed$warnings)
      }
    }
  }
  expect_identical(.msr_native_raster_eval, callback)
})

test_that("Gamma and NB BAM curves keep conditional means and raw reference axes", {
  skip_if_not_installed("mgcv")
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  for (entry in native_bam_extended_fixture()) {
    model <- entry$model; before <- model; object <- entry$object
    reference <- entry$data; row.names(reference) <- row.names(model$model)
    object$opt_mod$call$data <- quote(expired_extended_bam_reference)
    variables <- intersect(all.vars(delete.response(terms(model))), names(entry$data))
    source <- entry$data[variables]
    for (link in c(FALSE, TRUE)) {
      plots <- suppressMessages(plot_marginal_effects(object, length.out = 7L,
        reference_data = reference, link = link))
      expect_true(attr(plots, "reference_data_contract")$native_design_and_offsets_verified)
      # Single-component native curves retain their existing caption contract;
      # structured output-column contracts belong to multi-component curves.
      expect_null(attr(plots, "prediction_contract"))
      for (variable in c("mean", "rough", "z")) {
        nd <- as.data.frame(lapply(source, function(x) if (is.numeric(x)) rep(mean(x), 7L) else x[rep(1L, 7L)]))
        xs <- seq(min(entry$data[[variable]]), max(entry$data[[variable]]), length.out = 7L)
        nd[[variable]] <- xs
        expected <- bam_extended_native(entry, nd, "link", FALSE)$value
        target <- cbind(fit = as.numeric(expected$fit),
          lwr = as.numeric(expected$fit + qnorm(.025) * expected$se.fit),
          upr = as.numeric(expected$fit + qnorm(.975) * expected$se.fit))
        target[] <- model$family$linkinv(target)
        expect_equal(bam_family_clean(as.matrix(plots[[variable]]$data[c("fit", "lwr", "upr")])),
          bam_family_clean(target), tolerance = 1e-8)
        axis <- if (variable %in% names(object$scl_params$mean))
          xs * object$scl_params$sd[[variable]] + object$scl_params$mean[[variable]] else xs
        expect_equal(plots[[variable]]$data$x, axis, tolerance = 1e-10, ignore_attr = TRUE)
        expect_match(plots[[variable]]$labels$caption, "condition on the selected spatial scales")
      }
    }
    expect_identical(entry$model, before)
  }
})

test_that("NB BAM selection counts actual scale and theta parameters exactly once", {
  skip_if_not_installed("mgcv")
  for (entry in native_bam_extended_fixture()) {
    model <- entry$model; before <- model; labels <- c("wrapped", "native")
    quiet <- capture.output(aic <- aic_tab(list(entry$object, model), AICc = FALSE,
      mod_names = labels, verbose = FALSE))
    quiet2 <- capture.output(bic <- bic_tab(list(entry$object, model), mod_names = labels, verbose = FALSE))
    corrected <- entry$kind != "Gamma"
    ll <- if (corrected) sum(model$prior.weights * dnbinom(model$y,
      size = model$family$getTheta(TRUE), mu = model$fitted.values, log = TRUE)) else as.numeric(logLik(model))
    count <- entry$diagnostics$actual_parameter_count
    extra <- if (entry$f$inputs$kernel == "expow") 4L else 2L
    expected_k <- count + c(extra, 0)
    ia <- match(labels, aic$Modnames); ib <- match(labels, bic$Modnames); n <- as.numeric(nobs(model))
    expect_equal(aic$K[ia], expected_k, tolerance = 1e-10)
    expect_equal(bic$K[ib], expected_k, tolerance = 1e-10)
    expect_equal(aic$LL[ia], rep(ll, 2), tolerance = 1e-10)
    expect_equal(bic$LL[ib], rep(ll, 2), tolerance = 1e-10)
    expect_equal(aic$AIC[ia], -2 * ll + 2 * expected_k, tolerance = 1e-8)
    expect_equal(bic$BIC[ib], -2 * ll + log(n) * expected_k, tolerance = 1e-8)
    expect_equal(.msr_parameter_count(model), count, tolerance = 1e-10)
    expect_identical(c(quiet, quiet2), character())
    auto <- aic_tab(list(entry$object, model), AICc = FALSE, verbose = FALSE)
    expect_true(all(grepl(paste0("family=", model$family$family), auto$Modnames, fixed = TRUE)))
    adjustments <- attr(aic, "likelihood_adjustments")
    if (corrected && isTRUE(model$scale.estimated)) {
      expect_length(adjustments, 2L)
      expect_identical(adjustments[[1L]]$input_position, 1L)
      expect_identical(adjustments[[2L]]$input_position, 2L)
      expect_equal(adjustments[[1L]]$estimated_theta_parameters, if (entry$kind == "nb_fixed") 0 else 1)
      expect_equal(adjustments[[1L]]$estimated_scale_parameters, 0)
      expect_equal(adjustments[[1L]]$K + extra, aic$K[ia][1L], tolerance = 1e-10)
      expect_equal(adjustments[[1L]]$native_K - adjustments[[1L]]$K, 1, tolerance = 1e-10)
      expect_equal(adjustments[[1L]]$native_logLik - adjustments[[1L]]$logLik, 1, tolerance = 1e-10)
    } else expect_null(adjustments)
    if (corrected) {
      healthy <- model; healthy$scale.estimated <- FALSE
      check <- .msr_bam_nb_selection(healthy)
      expect_null(check$adjustment)
      expect_equal(check$K, count, tolerance = 1e-10)
      expect_equal(check$loglik, ll, tolerance = 1e-10)
      bad <- model; bad$aic <- bad$aic + .5
      expect_error(.msr_bam_nb_selection(bad), "likelihood metadata do not match")
      bad <- model; bad$sig2 <- 2
      expect_error(.msr_bam_nb_selection(bad), "likelihood metadata do not match")
      bad <- model; bad$family$variance <- function(mu) mu
      expect_error(.msr_bam_nb_selection(bad), "family implementation changed")
    } else expect_null(.msr_bam_nb_selection(model))
    expect_identical(entry$model, before)
  }
})

test_that("Gamma and NB BAM workers retain native matrices and conditional covariance", {
  skip_if_not_installed("mgcv")
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  pool <- parallel::makePSOCKcluster(2L); on.exit(parallel::stopCluster(pool), add = TRUE)
  configure <- function(coding) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    options(contrasts = c(coding, "contr.poly")); TRUE
  }
  environment(configure) <- baseenv()
  for (entry in native_bam_extended_fixture()) {
    before <- entry$model
    nd <- entry$data[rep(seq_len(nrow(entry$data)), 6L), , drop = FALSE]
    nd$case_site <- rep(seq_len(nrow(entry$data)), 6L); row.names(nd) <- paste0("extended-", seq_len(nrow(nd)))
    keys <- data.frame(case_site = seq_len(nrow(entry$data)))
    variables <- intersect(all.vars(delete.response(terms(entry$model))), names(entry$data))
    for (type in c("link", "response", "terms", "iterms", "lpmatrix")) {
      invisible(parallel::clusterCall(pool, configure, "contr.sum"))
      options(contrasts = c("contr.sum", "contr.poly"))
      native_parallel <- bam_extended_capture(predict(entry$model, nd[variables], type = type,
        se.fit = type != "lpmatrix", discrete = FALSE, cluster = pool, n.threads = 1L))
      # Native BAM flattens parallel iterms matrices. Use its serial numerical
      # target, while comparing parent warnings to the same native parallel call.
      # This also preserves pkgload serialization diagnostics in local test runs.
      target <- if (type == "iterms") bam_extended_native(entry, nd, type, FALSE) else native_parallel
      invisible(parallel::clusterCall(pool, configure, "contr.treatment"))
      options(contrasts = c("contr.treatment", "contr.poly"))
      actual <- bam_extended_capture(predict(entry$object, nd, entry$f$inputs, join_by = keys,
        include_fixed = TRUE, type = type, se.fit = type != "lpmatrix", uncertainty = "conditional",
        discrete = FALSE, cluster = pool, n.threads = 1L))
      expect_equal(bam_family_clean(actual$value), bam_family_clean(target$value), tolerance = 1e-8)
      expect_identical(actual$warnings, native_parallel$warnings,
        info = paste(entry$kind, entry$mode, type))
    }
    expect_identical(entry$model, before)
  }
})

# Eight deterministic native software cases, sharing one preparation. This
# checks adapter behavior at selected scales, not ecological-scale recovery.
weighted_bam_fixture <- local({
  cache <- NULL
  function() {
    if (inherits(cache, "error")) stop(cache)
    if (!is.null(cache)) return(cache)
    cache <<- tryCatch({
      old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
      f <- refit_spatial_fixture("gaussian"); dat <- f$data; n <- nrow(dat)
      p <- ((seq_len(n) * 11L) %% 37L) / 37
      mu <- exp(1 + .35 * dat$mean - .2 * dat$rough + .1 * dat$z)
      dat$y <- qnbinom(p, size = 2, mu = mu)
      form <- as.formula("y ~ rough+z+z:group+I(z^2)+s(mean,k=4)", env = asNamespace("mgcv"))
      cases <- list(); calls <- 0L
      for (kind in c("legacy", "extended_fixed", "extended_estimated")) {
        for (weighted in if (kind == "legacy") c(FALSE, TRUE) else TRUE) {
          for (mode in c("ordinary", "discrete")) {
            key <- paste(kind, weighted, mode, sep = "_")
            calls <- calls + 1L
            cases[[key]] <- tryCatch({
              dat$case_weight <- if (weighted) rep(c(0, .5, 1, 1.5, 2, 1), length.out = n) else rep(1, n)
              family <- if (kind == "legacy") mgcv::negbin(2) else
                mgcv::nb(theta = if (kind == "extended_fixed") 2 else NULL)
              model <- mgcv::bam(form, data = dat, weights = case_weight, family = family,
                scale = if (kind == "legacy") 1 else 0, method = "fREML", discrete = mode == "discrete",
                nthreads = 1L, cluster = NULL, control = mgcv::gam.control(epsilon = 1e-10))
              stopifnot(isTRUE(model$converged), all(is.finite(coef(model))),
                all(is.finite(model$sp)), length(model$sp) == 1L, length(model$smooth) == 1L)
              conv <- model$mgcv.conv
              if (is.logical(conv)) stopifnot(isTRUE(conv))
              if (is.list(conv) && !is.null(conv$message)) stopifnot(identical(conv$message, "full convergence"))
              if (is.list(conv)) for (flag in intersect(c("fully.converged", "hess.pos.def"), names(conv))) stopifnot(isTRUE(conv[[flag]]))
              if (!is.null(model$outer.info$conv)) stopifnot(identical(model$outer.info$conv, "full convergence"))
              if (!is.null(model$outer.info$hess)) stopifnot(min(eigen(model$outer.info$hess, symmetric = TRUE, only.values = TRUE)$values) > 0)
              theta <- if (kind == "legacy") model$family$getTheta() else model$family$getTheta(TRUE)
              stopifnot(length(theta) == 1L, is.finite(theta), theta > 0,
                identical(as.numeric(model$sig2), 1), identical(as.numeric(model$prior.weights), dat$case_weight),
                nobs(model) == sum(dat$case_weight > 0))
              if (kind != "extended_estimated") stopifnot(identical(as.numeric(theta), 2))
              X <- predict(model, dat, type = "lpmatrix", discrete = FALSE, n.threads = 1L)
              S <- matrix(0, ncol(X), ncol(X)); sm <- model$smooth[[1L]]; ix <- sm$first.para:sm$last.para
              S[ix, ix] <- model$sp[[1L]] * sm$S[[1L]]
              derivative <- model$family$mu.eta(model$linear.predictors)
              variance <- model$family$variance(model$fitted.values)
              W <- model$prior.weights * derivative^2 / variance
              score <- drop(crossprod(X, model$prior.weights * (model$y - model$fitted.values) * derivative / variance) - S %*% coef(model))
              info <- min(eigen(crossprod(X, X * W) + S, symmetric = TRUE, only.values = TRUE)$values)
              covariance <- min(eigen(model$Vp, symmetric = TRUE, only.values = TRUE)$values)
              stopifnot(all(is.finite(W)), all(W >= 0), any(W > 0), all(is.finite(score)),
                max(abs(score)) < 1e-4, info > 0, covariance > 0)
              edf <- if (is.null(model$edf2)) model$edf else model$edf2
              count <- min(sum(edf), length(coef(model))) + as.integer(kind == "extended_estimated")
              object <- f$object; object$opt_mod <- model
              list(model = model, object = object, data = dat, f = f, kind = kind, mode = mode,
                weighted = weighted, theta = theta, count = count, score = max(abs(score)),
                information = info, covariance = covariance)
            }, error = function(e) list(error = conditionMessage(e)))
          }
        }
      }
      stopifnot(calls == 8L, length(cases) == 8L)
      failed <- vapply(cases, function(x) !is.null(x$error), logical(1))
      if (any(failed)) stop(paste(names(cases)[failed], vapply(cases[failed], `[[`, character(1), "error"), collapse = "; "))
      cases
    }, error = identity)
    if (inherits(cache, "error")) stop(cache)
    cache
  }
})

weighted_bam_capture <- function(expr) {
  seen <- list()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    seen[[length(seen) + 1L]] <<- list(message = conditionMessage(w), classes = class(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = seen)
}

weighted_bam_native <- function(entry, data, type, optimized, exclude = NULL) {
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  if (!optimized) return(weighted_bam_capture(bam_family_oracle(entry, data, type, list(exclude = exclude))))
  vars <- intersect(all.vars(delete.response(terms(entry$model))), names(entry$data))
  weighted_bam_capture(predict(entry$model, data[vars], type = type, se.fit = type != "lpmatrix",
    discrete = TRUE, n.threads = 1L, exclude = exclude))
}

test_that("weighted and legacy NB BAM points and snapshots retain native targets", {
  skip_if_not_installed("mgcv")
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  cases <- weighted_bam_fixture(); expect_length(cases, 8L)
  for (entry in cases) {
    before <- entry$model; transform <- entry$object$scl_params
    expect_lt(entry$score, 1e-4); expect_gt(entry$information, 0); expect_gt(entry$covariance, 0)
    for (optimized in c(FALSE, if (entry$mode == "discrete") TRUE)) {
      for (type in c("link", "response", "terms", "iterms", "lpmatrix")) {
        for (rows in list(seq_len(nrow(entry$data)), c(1L, 17L, 24L), 17L)) {
          nd <- entry$data[rows, , drop = FALSE]
          target <- weighted_bam_native(entry, nd, type, optimized)
          actual <- weighted_bam_capture(predict(entry$object, nd, entry$f$inputs, include_fixed = TRUE,
            type = type, se.fit = type != "lpmatrix", uncertainty = "conditional", discrete = optimized, n.threads = 1L))
          expect_equal(bam_family_clean(actual$value), bam_family_clean(target$value), tolerance = 1e-8)
          expect_identical(actual$warnings, target$warnings)
          expect_identical(attr(actual$value, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
        }
        nd <- entry$data[c(1L, 17L, 24L), , drop = FALSE]
        snapshot <- prepare_prediction(entry$object, nd, entry$f$inputs, include_fixed = TRUE,
          component = type, native_args = list(discrete = optimized, exclude = "s(mean)"))
        target <- weighted_bam_native(entry, nd, type, optimized, "s(mean)")
        actual <- weighted_bam_capture(predict(unserialize(serialize(snapshot, NULL)), type = type,
          se.fit = type != "lpmatrix", uncertainty = "conditional"))
        expect_equal(bam_family_clean(actual$value), bam_family_clean(target$value), tolerance = 1e-8)
        expect_identical(actual$warnings, target$warnings)
      }
    }
    expect_identical(entry$model, before); expect_identical(entry$object$scl_params, transform)
  }
})

test_that("weighted and legacy NB BAM maps use exact native callback batches", {
  skip_if_not_installed("mgcv")
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  callback <- .msr_native_raster_eval
  for (entry in weighted_bam_fixture()) {
    grid <- kernel_scale.raster(entry$f$rasters, multiScaleR = entry$object, include_fixed = TRUE,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(group = "b", z = .2), verbose = FALSE)
    batches <- list()
    checked <- function(model, data) {
      seen <- list()
      value <- withCallingHandlers(callback(model, data), warning = function(w)
        seen[[length(seen) + 1L]] <<- list(message = conditionMessage(w), classes = class(w)))
      before <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(before), add = TRUE)
      native <- weighted_bam_capture(do.call(stats::predict, c(list(object = model$object$opt_mod,
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
        observed <- with_mocked_bindings(weighted_bam_capture(predict(entry$object, grid,
          type = type, se.fit = type != "lpmatrix", uncertainty = "conditional", discrete = optimized,
          raster_args = list(cores = 1L))), .msr_native_raster_eval = checked)
        plain <- weighted_bam_capture(predict(entry$object, grid, type = type,
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

test_that("weighted and legacy NB BAM curves keep conditional means and raw reference axes", {
  skip_if_not_installed("mgcv")
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  for (entry in weighted_bam_fixture()) {
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
        expected <- weighted_bam_native(entry, nd, "link", FALSE)$value
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

test_that("weighted NB BAM selection counts actual nuisance and spatial parameters", {
  skip_if_not_installed("mgcv")
  for (entry in weighted_bam_fixture()) {
    before <- entry$model; model <- entry$model
    direct <- sum(model$prior.weights * dnbinom(model$y, size = entry$theta, mu = model$fitted.values, log = TRUE))
    mods <- list(native = model, spatial = entry$object)
    expect_output(a <- aic_tab(mods, AICc = FALSE, mod_names = names(mods), verbose = FALSE), NA)
    expect_output(b <- bic_tab(mods, mod_names = names(mods), verbose = FALSE), NA)
    ia <- match(c("native", "spatial"), a$Modnames); ib <- match(c("native", "spatial"), b$Modnames)
    expect_equal(a$K[ia], entry$count + c(0, 2), tolerance = 1e-10)
    expect_equal(b$K[ib], entry$count + c(0, 2), tolerance = 1e-10)
    expect_equal(a$LL[ia], rep(direct, 2), tolerance = 1e-10)
    expect_equal(b$LL[ib], rep(direct, 2), tolerance = 1e-10)
    expect_equal(.msr_selection_nobs(model), sum(model$prior.weights > 0))
    if (entry$kind == "legacy") {
      expect_null(.msr_bam_nb_selection(model)); expect_null(attr(a, "likelihood_adjustments"))
      expect_false(model$scale.estimated)
    } else {
      changes <- attr(a, "likelihood_adjustments")
      expect_length(changes, 2L)
      expect_equal(changes[[1L]]$estimated_theta_parameters, as.integer(entry$kind == "extended_estimated"))
      expect_equal(changes[[1L]]$estimated_scale_parameters, 0)
    }
    expect_identical(entry$model, before)
  }
})

test_that("large-weight NB bookkeeping verifies density without cancellation or altered metadata", {
  skip_if_not_installed("mgcv")
  entry <- weighted_bam_fixture()[["extended_fixed_TRUE_ordinary"]]
  original <- entry$model
  for (factor in c(1, 100, 1e4, 1e6, 1e8, 1e12)) {
    # Arithmetic copies only. Coefficients are not refitted and these are not
    # new qualified inference cases or ecological recovery evidence.
    model <- original; model$prior.weights <- original$prior.weights * factor
    model$aic <- model$family$aic(model$y, model$fitted.values, wt = model$prior.weights, dev = 0) + 2 * sum(model$edf)
    direct <- sum(model$prior.weights * dnbinom(model$y, size = entry$theta, mu = model$fitted.values, log = TRUE))
    native <- logLik(model); checked <- .msr_bam_nb_selection(model)
    expect_equal(checked$loglik, direct, tolerance = 0)
    expect_equal(checked$K, entry$count, tolerance = 1e-10)
    expect_equal(checked$adjustment$estimated_theta_parameters, 0)
    if (abs(as.numeric(native) - direct - 1) > 1e-8)
      expect_identical(checked$adjustment$verification, "native_aic_and_rescaled_density")
    bad <- model; bad$aic <- bad$aic + .5
    expect_error(.msr_bam_nb_selection(bad), "likelihood metadata do not match")
    known <- model; known$scale.estimated <- FALSE
    value <- .msr_bam_nb_selection(known)
    expect_null(value$adjustment); expect_equal(value$K, entry$count, tolerance = 1e-10)
    expect_equal(value$loglik, as.numeric(logLik(known)), tolerance = 0)
  }
  expect_identical(entry$model, original)
})

test_that("weighted NB BAM native workers preserve geometry and diagnostics", {
  skip_if_not_installed("mgcv")
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  pool <- parallel::makePSOCKcluster(2L); on.exit(parallel::stopCluster(pool), add = TRUE)
  configure <- function(coding) { options(contrasts = c(coding, "contr.poly")); TRUE }
  environment(configure) <- baseenv()
  for (entry in weighted_bam_fixture()) {
    before <- entry$model
    nd <- entry$data[rep(seq_len(nrow(entry$data)), 6L), , drop = FALSE]
    nd$case_site <- rep(seq_len(nrow(entry$data)), 6L); rownames(nd) <- paste0("weighted-", seq_len(nrow(nd)))
    keys <- data.frame(case_site = seq_len(nrow(entry$data)))
    vars <- intersect(all.vars(delete.response(terms(entry$model))), names(entry$data))
    for (type in c("link", "response", "terms", "iterms", "lpmatrix")) {
      invisible(parallel::clusterCall(pool, configure, "contr.sum"))
      options(contrasts = c("contr.sum", "contr.poly"))
      native <- weighted_bam_capture(predict(entry$model, nd[vars], type = type, se.fit = type != "lpmatrix",
        discrete = FALSE, cluster = pool, n.threads = 1L))
      target <- if (type == "iterms") weighted_bam_native(entry, nd, type, FALSE) else native
      invisible(parallel::clusterCall(pool, configure, "contr.treatment"))
      options(contrasts = c("contr.treatment", "contr.poly"))
      actual <- weighted_bam_capture(predict(entry$object, nd, entry$f$inputs, include_fixed = TRUE,
        join_by = keys, type = type, se.fit = type != "lpmatrix", uncertainty = "conditional",
        discrete = FALSE, cluster = pool, n.threads = 1L))
      expect_equal(bam_family_clean(actual$value), bam_family_clean(target$value), tolerance = 1e-8)
      expect_identical(actual$warnings, native$warnings)
    }
    expect_identical(entry$model, before)
  }
})

test_that("NB roundoff verification rejects inconsistent theta state", {
  skip_if_not_installed("mgcv")
  # Pure extraction arithmetic, not a fitted or qualified inference model.
  family <- mgcv::nb(theta = 2)
  model <- structure(list(family = family, sig2 = 1, y = 0:3, fitted.values = 1:4,
    prior.weights = rep(1e8, 4), coefficients = c(a = 0, b = 0),
    edf = c(a = .75, b = .8), edf2 = c(a = .75, b = .8), scale.estimated = TRUE),
    class = c("bam", "gam", "glm", "lm"))
  model$aic <- family$aic(model$y, model$fitted.values, wt = model$prior.weights, dev = 0) + 2 * sum(model$edf)
  before <- model
  direct <- sum(model$prior.weights * dnbinom(model$y, size = 2, mu = model$fitted.values, log = TRUE))
  expect_equal(.msr_bam_nb_selection(model)$loglik, direct, tolerance = 0)
  expect_equal(.msr_bam_nb_selection(model)$K, 1.55, tolerance = 1e-10)
  same <- model; same$family$getTheta <- mgcv::nb(theta = 2)$getTheta
  expect_equal(.msr_bam_nb_selection(same)$loglik, direct, tolerance = 0)
  bad <- model; bad$family$getTheta <- mgcv::nb(theta = 2 * (1 + 1e-9))$getTheta
  expect_error(.msr_bam_nb_selection(bad), "likelihood metadata do not match")
  bad <- model; bad$family$aic <- mgcv::nb(theta = 2 * (1 + 1e-9))$aic
  bad$aic <- bad$family$aic(bad$y, bad$fitted.values, wt = bad$prior.weights, dev = 0) + 2 * sum(bad$edf)
  expect_error(.msr_bam_nb_selection(bad), "likelihood metadata do not match")
  expect_identical(model, before)
})

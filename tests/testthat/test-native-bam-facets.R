# One failure-retaining software cache per declared kernel. Original fixture
# definitions and all fitting diagnostics are reused without changing settings.
native_bam_facets_fixture <- local({
  cache <- new.env(parent = emptyenv())
  function(kernel) {
    if (exists(kernel, cache, inherits = FALSE)) {
      value <- get(kernel, cache)
      if (inherits(value, "error")) stop(value)
      return(value)
    }
    value <- tryCatch({
      e <- new.env(parent = parent.frame())
      original <- get("refit_spatial_fixture", envir = parent.frame(), inherits = TRUE)
      e$refit_spatial_fixture <- local({
        core <- original; selected <- kernel
        function(...) core(..., kernel = selected)
      })
      eval(parse(testthat::test_path("test-native-bam-coding.R"))[[1L]], e)
      definitions <- parse(testthat::test_path("test-native-bam-terms.R"))
      allowed <- c("native_bam_terms_fixture", "bam_terms_value", "bam_terms_capture",
        "bam_terms_native_capture", "bam_terms_partition_warnings")
      found <- character()
      for (expr in definitions) if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
          is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% allowed) {
        eval(expr, e); found <- c(found, as.character(expr[[2L]]))
      }
      stopifnot(setequal(found, allowed))
      f <- evalq(native_bam_terms_fixture(), e)
      stopifnot(length(f$models) == 6L, identical(f$f$inputs$kernel, kernel))
      f
    }, error = identity)
    assign(kernel, value, cache)
    if (inherits(value, "error")) stop(value)
    value
  }
})

test_that("Gaussian BAM maps and raw-reference curves agree across four kernels", {
  skip_if_not_installed("mgcv")
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    f <- native_bam_facets_fixture(kernel)
    for (kind in names(f$models)) {
      model <- f$models[[kind]]$model; original <- model
      object <- f$f$object; object$opt_mod <- model
      vars <- intersect(all.vars(delete.response(terms(model))), names(f$data))
      scenario <- NULL
      if ("group" %in% vars) scenario$group <- "b"
      if ("z" %in% vars) scenario$z <- .2
      grid <- kernel_scale.raster(f$f$rasters, multiScaleR = object, include_fixed = TRUE,
        scale_center = TRUE, match_buffer = TRUE, site_values = scenario, verbose = FALSE)
      source <- f$data[vars]; map <- source
      if ("group" %in% vars) map$group <- factor("b", levels = levels(f$data$group))
      if ("z" %in% vars) map$z <- .2
      for (type in c("response", "iterms")) {
        expected <- predict(model, rbind(source, map), type = type, se.fit = TRUE,
          discrete = FALSE, n.threads = 1L)
        expected <- f$slice(expected, nrow(source) + seq_len(nrow(map)))
        expected <- as.matrix(data.frame(fit = expected$fit, se.fit = expected$se.fit))
        options(contrasts = c("contr.treatment", "contr.poly"))
        actual <- predict(object, grid, type = type, se.fit = TRUE, uncertainty = "conditional",
          discrete = FALSE, raster_args = list(cores = 1L))
        options(contrasts = c("contr.sum", "contr.poly"))
        sampled <- as.matrix(terra::extract(actual, f$f$points, ID = FALSE))
        expect_equal(unname(sampled), unname(expected), tolerance = 1e-8)
        expect_identical(colnames(sampled), colnames(expected))
      }
      reference <- f$data; row.names(reference) <- row.names(model$model)
      object$opt_mod$call$data <- quote(expired_bam_facet_reference)
      options(contrasts = c("contr.treatment", "contr.poly"))
      plots <- suppressMessages(plot_marginal_effects(object, length.out = 7L, link = TRUE,
        reference_data = reference))
      options(contrasts = c("contr.sum", "contr.poly"))
      expect_true(attr(plots, "reference_data_contract")$native_design_and_offsets_verified)
      numeric_vars <- vars[vapply(f$data[vars], is.numeric, logical(1))]
      expect_true(all(numeric_vars %in% names(plots)))
      for (var in numeric_vars) {
        nd <- as.data.frame(lapply(f$data[vars], function(x)
          if (is.numeric(x)) rep(mean(x), 7L) else x[rep(1L, 7L)]), check.names = FALSE)
        xs <- seq(min(f$data[[var]]), max(f$data[[var]]), length.out = 7L)
        nd[[var]] <- xs
        expected <- predict(model, rbind(f$data[vars], nd), type = "link", se.fit = TRUE, n.threads = 1L)
        expected <- f$slice(expected, nrow(f$data) + seq_len(7L))
        target <- cbind(fit = as.numeric(expected$fit),
          lwr = as.numeric(expected$fit + qnorm(.025) * expected$se.fit),
          upr = as.numeric(expected$fit + qnorm(.975) * expected$se.fit))
        expect_equal(unname(as.matrix(plots[[var]]$data[c("fit", "lwr", "upr")])),
          unname(target), tolerance = 1e-8)
        display <- if (var %in% names(object$scl_params$mean))
          xs * object$scl_params$sd[[var]] + object$scl_params$mean[[var]] else xs
        expect_equal(plots[[var]]$data$x, display, tolerance = 1e-10, ignore_attr = TRUE)
      }
      expect_identical(model, original)
    }
  }
})

test_that("Gaussian BAM tables retain native EDF and selected scale/shape penalties", {
  skip_if_not_installed("mgcv")
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    f <- native_bam_facets_fixture(kernel)
    for (entry in f$models) {
      model <- entry$model; object <- f$f$object; object$opt_mod <- model
      models <- list(object, model); labels <- c("wrapped", "native")
      output <- capture.output(aic <- aic_tab(models, AICc = FALSE, mod_names = labels, verbose = FALSE))
      output2 <- capture.output(bic <- bic_tab(models, mod_names = labels, verbose = FALSE))
      auto <- aic_tab(models, AICc = FALSE, verbose = FALSE)
      extra <- if (kernel == "expow") 4L else 2L
      count <- as.numeric(attr(logLik(model), "df")) + c(extra, 0)
      ll <- as.numeric(logLik(model)); n <- as.numeric(nobs(model))
      ia <- match(labels, aic$Modnames); ib <- match(labels, bic$Modnames)
      expect_equal(aic$K[ia], count)
      expect_equal(bic$K[ib], count)
      expect_equal(aic$AIC[ia], -2 * ll + 2 * count, tolerance = 1e-8)
      expect_equal(bic$BIC[ib], -2 * ll + log(n) * count, tolerance = 1e-8)
      expect_true(all(grepl("family=gaussian", auto$Modnames, fixed = TRUE)))
      expect_identical(c(output, output2), character())
    }
  }
})

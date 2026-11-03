# Two native Gaussian modes at the original fixed physical ecological scales.
# Aggregate native equalities cover values, shapes and labels, not count targets.
native_bam_coding_fixture <- local({
  cache <- NULL
  function() {
    if (inherits(cache, "error")) stop(cache)
    if (!is.null(cache)) return(cache)
    cache <<- tryCatch({
      f <- refit_spatial_fixture(); dat <- f$data
      dat$y <- dat$y + as.integer(dat$site)^2
      old <- options(contrasts = c("contr.sum", "contr.poly"))
      on.exit(options(old), add = TRUE)
      form <- as.formula("y ~ rough + z + z:group + s(mean,k=4)", env = asNamespace("mgcv"))
      types <- c("link", "response", "terms", "lpmatrix")
      controls <- list(default = list(), excluded = list(exclude = "s(mean)"),
        blocks = list(block.size = 1L))
      slice <- function(value, rows) {
        if (is.list(value)) { result <- lapply(value, slice, rows = rows); attributes(result) <- attributes(value); return(result) }
        if (is.matrix(value)) result <- value[rows, , drop = FALSE] else if (length(dim(value)) == 1L) {
          labels <- dimnames(value)
          result <- array(as.numeric(value)[rows], length(rows),
            dimnames = if (is.null(labels)) NULL else list(labels[[1L]][rows]))
        } else result <- value[rows]
        for (name in setdiff(names(attributes(value)), c("names", "dim", "dimnames")))
          attr(result, name) <- attr(value, name)
        result
      }
      models <- list()
      for (kind in c("ordinary", "discrete")) {
        model <- mgcv::bam(form, data = dat, family = gaussian(), method = "fREML",
          discrete = kind == "discrete", nthreads = 1, cluster = NULL)
        stopifnot(isTRUE(model$converged), all(is.finite(coef(model))), length(model$sp) == 1L,
          length(model$smooth) == 1L, length(model$smooth[[1L]]$S) == 1L)
        conv <- model$mgcv.conv
        if (is.logical(conv)) stopifnot(isTRUE(conv))
        if (is.list(conv) && !is.null(conv$message)) stopifnot(identical(conv$message, "full convergence"))
        if (is.list(conv)) for (flag in intersect(c("fully.converged", "hess.pos.def"), names(conv)))
          stopifnot(isTRUE(conv[[flag]]))
        if (!is.null(model$outer.info$conv)) stopifnot(identical(model$outer.info$conv, "full convergence"))
        if (!is.null(model$outer.info$hess))
          stopifnot(min(eigen(model$outer.info$hess, symmetric = TRUE, only.values = TRUE)$values) > 0)
        X <- predict(model, dat, type = "lpmatrix", n.threads = 1)
        S <- matrix(0, ncol(X), ncol(X)); smooth <- model$smooth[[1L]]
        ix <- smooth$first.para:smooth$last.para; S[ix, ix] <- model$sp[[1L]] * smooth$S[[1L]]
        score <- drop(crossprod(X, dat$y - drop(X %*% coef(model))) - S %*% coef(model))
        stopifnot(max(abs(score)) < 1e-4,
          min(eigen(crossprod(X) + S, symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(model$Vp, symmetric = TRUE, only.values = TRUE)$values) > 0)
        source <- dat[c("rough", "z", "group", "mean")]
        outputs <- setNames(lapply(types, function(type) lapply(controls, function(control)
          do.call(stats::predict, c(list(object = model, newdata = source, type = type,
            se.fit = type != "lpmatrix", n.threads = 1), control)))), types)
        map_data <- source; map_data$group <- factor("b", levels = levels(dat$group)); map_data$z <- .2
        map <- predict(model, rbind(source, map_data), type = "response", se.fit = TRUE, n.threads = 1)
        curve_data <- data.frame(rough = mean(dat$rough), z = seq(min(dat$z), max(dat$z), length.out = 7),
          group = factor(levels(dat$group)[1L], levels = levels(dat$group)), mean = mean(dat$mean))
        curve <- predict(model, rbind(source, curve_data), type = "link", se.fit = TRUE, n.threads = 1)
        models[[kind]] <- list(model = model, outputs = outputs,
          map = slice(map, nrow(source) + seq_len(nrow(map_data))),
          curve = slice(curve, nrow(source) + seq_len(nrow(curve_data))))
      }
      list(f = f, data = dat, models = models, controls = controls, types = types, slice = slice)
    }, error = identity)
    # A failed fixture is retained in this test lifecycle, never refit for a pass.
    if (inherits(cache, "error")) stop(cache)
    cache
  }
})

bam_native_value <- function(x) { attr(x, "prediction_contract") <- NULL; x }

test_that("serial native BAM routes retain fitted coding and numerical outputs", {
  skip_if_not_installed("mgcv")
  f <- native_bam_coding_fixture()
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  method <- utils::getS3method("predict", "bam", envir = asNamespace("mgcv"))
  for (kind in names(f$models)) {
    entry <- f$models[[kind]]; model <- entry$model; original <- model
    object <- f$f$object; object$opt_mod <- model
    expect_equal(.msr_selection_parameter_count(object), as.numeric(attr(logLik(model), "df")) + 2L)
    for (type in f$types) for (control in names(f$controls)) for (rows in
        list(seq_len(nrow(f$data)), c(1L, 17L, 24L), 17L)) {
      nd <- f$data[rows, , drop = FALSE]
      expected <- f$slice(entry$outputs[[type]][[control]], rows)
      actual <- do.call(stats::predict, c(list(object = object, newdata = nd,
        kernel_inputs = f$f$inputs, include_fixed = TRUE, type = type,
        se.fit = type != "lpmatrix", uncertainty = "conditional", n.threads = 1), f$controls[[control]]))
      expect_equal(bam_native_value(actual), expected, tolerance = 1e-8)
    }
    nd <- f$data[c(1L, 17L, 24L), , drop = FALSE]
    for (type in f$types) for (control in c("excluded", "blocks")) {
      snapshot <- prepare_prediction(object, nd, f$f$inputs, include_fixed = TRUE,
        component = type, native_args = f$controls[[control]])
      restored <- unserialize(serialize(snapshot, NULL))
      actual <- predict(restored, type = type, se.fit = type != "lpmatrix", uncertainty = "conditional")
      expect_equal(bam_native_value(actual), f$slice(entry$outputs[[type]][[control]], c(1L, 17L, 24L)), tolerance = 1e-8)
    }
    expect_identical(model, original)
  }
  expect_identical(utils::getS3method("predict", "bam", envir = asNamespace("mgcv")), method)
})

test_that("native BAM checked maps and raw curves retain conditional means", {
  skip_if_not_installed("mgcv")
  f <- native_bam_coding_fixture()
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (kind in names(f$models)) {
    entry <- f$models[[kind]]; object <- f$f$object; object$opt_mod <- entry$model
    grid <- kernel_scale.raster(f$f$rasters, multiScaleR = object, include_fixed = TRUE,
      scale_center = TRUE, match_buffer = TRUE, site_values = list(group = "b", z = .2), verbose = FALSE)
    map <- predict(object, grid, type = "response", se.fit = TRUE, uncertainty = "conditional", n.threads = 1)
    expect_equal(unname(as.matrix(terra::extract(map, f$f$points, ID = FALSE))),
      unname(cbind(entry$map$fit, entry$map$se.fit)), tolerance = 1e-8)
    reference <- f$data; row.names(reference) <- row.names(entry$model$model)
    object$opt_mod$call$data <- quote(expired_native_bam_raw_data)
    curves <- plot_marginal_effects(object, length.out = 7, reference_data = reference)
    expected <- entry$curve
    expect_equal(as.matrix(curves$z$data[c("fit", "lwr", "upr")]),
      unname(cbind(as.numeric(expected$fit), as.numeric(expected$fit + qnorm(.025) * expected$se.fit),
        as.numeric(expected$fit + qnorm(.975) * expected$se.fit))), tolerance = 1e-8, ignore_attr = TRUE)
    expect_true(attr(curves, "reference_data_contract")$native_design_and_offsets_verified)
    point <- predict(object, f$data[17L, , drop = FALSE], f$f$inputs,
      se.fit = TRUE, uncertainty = "conditional")
    expect_identical(attr(point, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
  }
})

test_that("BAM snapshots keep external pools and joint uncertainty guarded", {
  skip_if_not_installed("mgcv")
  f <- native_bam_coding_fixture(); object <- f$f$object; object$opt_mod <- f$models$ordinary$model
  nd <- f$data[c(1L, 17L, 24L), , drop = FALSE]
  expect_error(predict(object, nd, f$f$inputs, se.fit = TRUE, uncertainty = "joint"), "conditional")
  expect_error(prepare_prediction(object, nd, f$f$inputs,
    native_args = list(cluster = structure(list(), class = "cluster"))), "supported native prediction controls")
})

test_that("native BAM prediction workers retain coding and explicit exclusions", {
  skip_if_not_installed("mgcv")
  f <- native_bam_coding_fixture()
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  # Repeated prediction rows, not new observations or an additional fit. 216
  # exceeds the native two-worker prediction threshold of 100 rows per worker.
  long <- f$data[rep(seq_len(nrow(f$data)), 6L), , drop = FALSE]
  long$case_site <- rep(seq_len(nrow(f$data)), 6L)
  row.names(long) <- paste0("bam-prediction-", seq_len(nrow(long)))
  keys <- data.frame(case_site = seq_len(nrow(f$data)))
  pool <- parallel::makePSOCKcluster(2L)
  on.exit(parallel::stopCluster(pool), add = TRUE)
  configure <- function(coding) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    options(contrasts = c(coding, "contr.poly")); TRUE
  }
  environment(configure) <- baseenv()
  invisible(parallel::clusterCall(pool, configure, "contr.sum"))
  oracles <- excluded <- list()
  for (kind in names(f$models)) {
    model <- f$models[[kind]]$model
    oracles[[kind]] <- setNames(lapply(f$types, function(type)
      predict(model, long, type = type, se.fit = type != "lpmatrix", cluster = pool,
        discrete = FALSE, block.size = 1L)), f$types)
    excluded[[kind]] <- setNames(lapply(c("link", "response"), function(type)
      predict(model, long, type = type, se.fit = TRUE, discrete = FALSE, exclude = "s(mean)")),
      c("link", "response"))
  }
  options(contrasts = c("contr.treatment", "contr.poly"))
  invisible(parallel::clusterCall(pool, configure, "contr.treatment"))
  for (kind in names(f$models)) {
    object <- f$f$object; object$opt_mod <- f$models[[kind]]$model
    for (type in f$types) {
      actual <- predict(object, long, f$f$inputs, join_by = keys, include_fixed = TRUE,
        type = type, se.fit = type != "lpmatrix", uncertainty = "conditional", cluster = pool,
        discrete = FALSE, block.size = 1L)
      expect_equal(bam_native_value(actual), oracles[[kind]][[type]], tolerance = 1e-8)
    }
    for (type in c("link", "response")) {
      # Preserve native parallel geometry, using the independent serial native
      # target for exclusions. Some provider workers omit that explicit control.
      expected <- oracles[[kind]][[type]]
      expected$fit[] <- as.numeric(excluded[[kind]][[type]]$fit)
      expected$se.fit[] <- as.numeric(excluded[[kind]][[type]]$se.fit)
      actual <- predict(object, long, f$f$inputs, join_by = keys, include_fixed = TRUE,
        type = type, se.fit = TRUE, uncertainty = "conditional", cluster = pool,
        discrete = FALSE, block.size = 1L, exclude = "s(mean)")
      expect_equal(bam_native_value(actual), expected, tolerance = 1e-8)
    }
  }
})

test_that("BAM worker provider changes fail before coefficient prediction", {
  skip_if_not_installed("mgcv")
  native <- get("predict.gam", envir = asNamespace("mgcv"))
  worker <- .msr_bam_saved_worker(native)
  environment(worker)$.expected_version <- "unqualified-version-fixture"
  expect_error(worker(list()), "verify the native GAM provider")
})

# Actual native Gaussian fixtures from audit169. Cache failures as well as
# successes; changing a seed, formula or fitting setting is not a repair.
native_bam_terms_fixture <- local({
  cache <- NULL
  function() {
    if (inherits(cache, "error")) stop(cache)
    if (!is.null(cache)) return(cache)
    cache <<- tryCatch({
      original <- get0("native_bam_coding_fixture", envir = parent.frame(),
        mode = "function", inherits = TRUE)
      if (is.null(original)) {
        e <- new.env(parent = parent.frame())
        expr <- parse(testthat::test_path("test-native-bam-coding.R"))[[1L]]
        stopifnot(identical(expr[[2L]], as.name("native_bam_coding_fixture")))
        eval(expr, e); original <- e$native_bam_coding_fixture
      }
      f <- original()
      old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
      formulas <- c(numeric = "y ~ rough + z + s(mean,k=4)",
        smooth_only = "y ~ 0 + s(mean,k=4) + s(rough,k=4)")
      for (name in names(formulas)) for (kind in c("ordinary", "discrete")) {
        model <- mgcv::bam(as.formula(formulas[[name]], env = asNamespace("mgcv")),
          data = f$data, family = gaussian(), method = "fREML",
          discrete = kind == "discrete", nthreads = 1L, cluster = NULL)
        stopifnot(isTRUE(model$converged), all(is.finite(coef(model))), all(is.finite(model$sp)),
          length(model$sp) == length(model$smooth), length(model$contrasts) == 0L)
        conv <- model$mgcv.conv
        if (is.logical(conv)) stopifnot(isTRUE(conv))
        if (is.list(conv) && !is.null(conv$message)) stopifnot(identical(conv$message, "full convergence"))
        if (is.list(conv)) for (flag in intersect(c("fully.converged", "hess.pos.def"), names(conv))) stopifnot(isTRUE(conv[[flag]]))
        if (!is.null(model$outer.info$conv)) stopifnot(identical(model$outer.info$conv, "full convergence"))
        if (!is.null(model$outer.info$hess)) stopifnot(min(eigen(model$outer.info$hess,
          symmetric = TRUE, only.values = TRUE)$values) > 0)
        X <- predict(model, f$data, type = "lpmatrix", n.threads = 1L)
        S <- matrix(0, ncol(X), ncol(X))
        for (j in seq_along(model$smooth)) {
          sm <- model$smooth[[j]]; stopifnot(length(sm$S) == 1L)
          ix <- sm$first.para:sm$last.para
          S[ix, ix] <- S[ix, ix] + model$sp[[j]] * sm$S[[1L]]
        }
        score <- drop(crossprod(X, f$data$y - drop(X %*% coef(model))) - S %*% coef(model))
        stopifnot(all(is.finite(score)), max(abs(score)) < 1e-4,
          min(eigen(crossprod(X) + S, symmetric = TRUE, only.values = TRUE)$values) > 0,
          min(eigen(model$Vp, symmetric = TRUE, only.values = TRUE)$values) > 0,
          name != "smooth_only" || identical(as.integer(model$nsdf), 0L))
        f$models[[paste(name, kind, sep = "_")]] <- list(model = model)
      }
      for (kind in names(f$models)) {
        model <- f$models[[kind]]$model
        labels <- colnames(predict(model, f$data, type = "iterms", discrete = FALSE))
        controls <- list(default = list(), type2 = list(iterms.type = 2L),
          selected = list(terms = "s(mean)"), excluded = list(exclude = "s(mean)"),
          empty = list(exclude = labels), blocks = list(block.size = 1L),
          type2_selected = list(iterms.type = 2L, terms = "s(mean)", block.size = 1L))
        f$models[[kind]]$term_controls <- controls
        f$models[[kind]]$terms <- lapply(controls, function(control) setNames(lapply(c(FALSE, TRUE), function(se)
          bam_terms_native_capture(do.call(stats::predict, c(list(object = model, newdata = f$data, type = "iterms",
            se.fit = se, discrete = FALSE, n.threads = 1L), control)), model)$value), c("FALSE", "TRUE")))
        f$models[[kind]]$row_warnings <- setNames(lapply(list(c(1L, 17L, 24L), 17L), function(rows)
          lapply(controls[c("type2_selected", "excluded")], function(control) setNames(lapply(c(FALSE, TRUE), function(se)
            bam_terms_native_capture(do.call(stats::predict, c(list(object = model,
              newdata = f$data[rows, , drop = FALSE], type = "iterms", se.fit = se,
              discrete = FALSE, n.threads = 1L), control)), model)), c("FALSE", "TRUE")))), c("1,17,24", "17"))
      }
      f
    }, error = identity)
    if (inherits(cache, "error")) stop(cache)
    cache
  }
})

bam_terms_value <- function(x) { attr(x, "prediction_contract") <- NULL; x }
bam_terms_capture <- function(expr) {
  warnings <- list()
  value <- withCallingHandlers(force(expr), warning = function(w) {
    warnings[[length(warnings) + 1L]] <<- w; invokeRestart("muffleWarning")
  })
  list(value = value, messages = vapply(warnings, conditionMessage, ""),
    classes = lapply(warnings, class))
}

# Native mgcv may inspect min(assign) on an empty parametric frame when
# selecting terms. Preserve that exact provider warning, without treating it
# as a failed fit or accepting any unrecognized diagnostic.
bam_terms_native_capture <- function(expr, model) {
  value <- bam_terms_capture(expr)
  allowed <- "no non-missing arguments to min; returning Inf"
  if (length(value$messages) && (any(value$messages != allowed) ||
      !identical(as.integer(model$nsdf), 0L)))
    stop("Unexpected native prediction warning: ", paste(value$messages, collapse = "; "))
  value
}
bam_terms_partition_warnings <- function(model, newdata, type, se, control) {
  parts <- list(seq_len(nrow(newdata) / 2L), seq.int(nrow(newdata) / 2L + 1L, nrow(newdata)))
  values <- lapply(parts, function(rows) bam_terms_native_capture(do.call(stats::predict,
    c(list(object = model, newdata = newdata[rows, , drop = FALSE], type = type,
      se.fit = se, discrete = FALSE, n.threads = 1L), control)), model))
  list(messages = as.character(unlist(lapply(values, `[[`, "messages"), use.names = FALSE)),
    classes = Reduce(c, lapply(values, `[[`, "classes"), init = list()))
}

test_that("native BAM ordinary term rows and numerical snapshots retain their targets", {
  skip_if_not_installed("mgcv")
  f <- native_bam_terms_fixture()
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old), add = TRUE)
  for (kind in names(f$models)) {
    entry <- f$models[[kind]]; model <- entry$model; original <- model
    object <- f$f$object; object$opt_mod <- model
    expect_equal(.msr_selection_parameter_count(object), as.numeric(attr(logLik(model), "df")) + 2L)
    for (control in c("type2_selected", "excluded")) for (se in c(FALSE, TRUE)) for (rows in list(c(1L, 17L, 24L), 17L)) {
      nd <- f$data[rows, , drop = FALSE]
      actual <- bam_terms_capture(do.call(stats::predict, c(list(object = object, newdata = nd, kernel_inputs = f$f$inputs,
        include_fixed = TRUE, type = "iterms", se.fit = se, uncertainty = "conditional",
        discrete = FALSE), entry$term_controls[[control]])))
      expect_equal(bam_terms_value(actual$value), f$slice(entry$terms[[control]][[as.character(se)]], rows), tolerance = 1e-8)
      warning_target <- entry$row_warnings[[paste(rows, collapse = ",")]][[control]][[as.character(se)]]
      expect_identical(actual$messages, warning_target$messages)
      expect_identical(actual$classes, warning_target$classes)
    }
    for (rows in list(c(1L, 17L, 24L), 17L)) {
      nd <- f$data[rows, , drop = FALSE]
      snapshot <- prepare_prediction(object, nd, f$f$inputs, include_fixed = TRUE,
        component = "iterms", native_args = list(disc = FALSE, iterms.t = 2L,
          terms = "s(mean)", block.size = 1L))
      restored <- unserialize(serialize(snapshot, NULL))
      actual <- bam_terms_capture(predict(restored, type = "iterms", se.fit = TRUE, uncertainty = "conditional"))
      expect_equal(bam_terms_value(actual$value), f$slice(entry$terms$type2_selected[["TRUE"]], rows), tolerance = 1e-8)
      warning_target <- entry$row_warnings[[paste(rows, collapse = ",")]]$type2_selected[["TRUE"]]
      expect_identical(actual$messages, warning_target$messages)
      expect_identical(actual$classes, warning_target$classes)
      expect_identical(restored$provenance$native_args$discrete, FALSE)
      expect_identical(attr(actual$value, "prediction_contract")$term_selection$discrete_requested, FALSE)
    }
    expect_identical(model, original)
  }
  expect_error(prepare_prediction(f$f$object, f$data[17L, , drop = FALSE], f$f$inputs,
    component = "response", native_args = list(discrete = FALSE)), "supported native prediction controls")
})

test_that("native BAM workers preserve term geometry, controls and covariance warnings", {
  skip_if_not_installed("mgcv")
  f <- native_bam_terms_fixture()
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  long <- f$data[rep(seq_len(nrow(f$data)), 6L), , drop = FALSE]
  long$case_site <- rep(seq_len(nrow(f$data)), 6L)
  row.names(long) <- paste0("bam-term-prediction-", seq_len(nrow(long)))
  keys <- data.frame(case_site = seq_len(nrow(f$data)))
  pool <- parallel::makePSOCKcluster(2L); on.exit(parallel::stopCluster(pool), add = TRUE)
  configure <- function(coding) {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    options(contrasts = c(coding, "contr.poly")); TRUE
  }
  environment(configure) <- baseenv(); invisible(parallel::clusterCall(pool, configure, "contr.sum"))
  targets <- list()
  for (kind in names(f$models)) {
    entry <- f$models[[kind]]; model <- entry$model
    term_targets <- lapply(entry$term_controls, function(control) setNames(lapply(c(FALSE, TRUE), function(se)
      bam_terms_native_capture(do.call(stats::predict, c(list(object = model, newdata = long, type = "iterms", se.fit = se,
        discrete = FALSE, n.threads = 1L), control)), model)$value), c("FALSE", "TRUE")))
    term_warnings <- lapply(entry$term_controls, function(control) setNames(lapply(c(FALSE, TRUE), function(se)
      bam_terms_partition_warnings(model, long, "iterms", se, control)), c("FALSE", "TRUE")))
    scalar_targets <- setNames(lapply(c("link", "response"), function(type) {
      # Native scalar parallel geometry is preserved. Exclusion values come
      # from the independent serial provider because native workers omit it.
      template <- predict(model, long, type = type, se.fit = TRUE, discrete = FALSE, cluster = pool)
      serial <- bam_terms_native_capture(predict(model, long, type = type, se.fit = TRUE,
        discrete = FALSE, exclude = "s(mean)"), model)$value
      excluded <- template; excluded$fit[] <- as.numeric(serial$fit); excluded$se.fit[] <- as.numeric(serial$se.fit)
      list(default = template, excluded = excluded, warnings = bam_terms_partition_warnings(
        model, long, type, TRUE, list(exclude = "s(mean)")))
    }), c("link", "response"))
    targets[[kind]] <- list(terms = term_targets, term_warnings = term_warnings, scalar = scalar_targets,
      optimized = if (!is.null(model$dinfo)) predict(model, long, type = "terms", se.fit = TRUE) else NULL)
  }
  options(contrasts = c("contr.treatment", "contr.poly"))
  invisible(parallel::clusterCall(pool, configure, "contr.treatment"))
  missing_covariance <- "Smoothness uncertainty corrected covariance not available"
  for (kind in names(f$models)) {
    entry <- f$models[[kind]]; model <- entry$model; original <- model
    object <- f$f$object; object$opt_mod <- model
    for (control in names(entry$term_controls)) for (se in c(FALSE, TRUE)) {
      actual <- bam_terms_capture(do.call(stats::predict, c(list(object = object, newdata = long, kernel_inputs = f$f$inputs,
        join_by = keys, include_fixed = TRUE, type = "iterms", se.fit = se,
        uncertainty = "conditional", cluster = pool, discrete = FALSE), entry$term_controls[[control]])))
      expect_equal(bam_terms_value(actual$value), targets[[kind]]$terms[[control]][[as.character(se)]], tolerance = 1e-8)
      warning_target <- targets[[kind]]$term_warnings[[control]][[as.character(se)]]
      expect_identical(actual$messages, warning_target$messages)
      expect_identical(actual$classes, warning_target$classes)
    }
    if (!length(model$contrasts)) for (type in c("link", "response")) {
      actual <- bam_terms_capture(predict(object, long, f$f$inputs, join_by = keys, include_fixed = TRUE,
        type = type, se.fit = TRUE, uncertainty = "conditional", cluster = pool,
        discrete = FALSE, exclude = "s(mean)"))
      expect_equal(bam_terms_value(actual$value), targets[[kind]]$scalar[[type]]$excluded, tolerance = 1e-8)
      expect_identical(actual$messages, targets[[kind]]$scalar[[type]]$warnings$messages)
      expect_identical(actual$classes, targets[[kind]]$scalar[[type]]$warnings$classes)
    }
    for (copy in c("actual", "software_only_Vc_1.25Vp")) for (type in c("link", "iterms")) {
      probe <- model
      if (copy != "actual") probe$Vc <- 1.25 * probe$Vp
      object$opt_mod <- probe
      actual <- bam_terms_capture(predict(object, long, f$f$inputs, join_by = keys,
        include_fixed = TRUE, type = type, se.fit = TRUE, uncertainty = "conditional",
        cluster = pool, discrete = FALSE, uncond = TRUE))
      expected <- if (type == "iterms") targets[[kind]]$terms$default[["TRUE"]] else targets[[kind]]$scalar$link$default
      expect_equal(bam_terms_value(actual$value), expected, tolerance = 1e-8)
      expect_identical(unname(actual$messages), rep(missing_covariance, 2L))
      expect_true(all(vapply(actual$classes, function(x) "warning" %in% x, logical(1))))
      contract <- attr(actual$value, "prediction_contract")
      expect_identical(contract$uncertainty$smoothing_parameters, "native_unconditional_requested")
      expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
    }
    object$opt_mod <- model
    if (!is.null(model$dinfo)) {
      actual <- bam_terms_capture(predict(object, long, f$f$inputs, join_by = keys,
        include_fixed = TRUE, type = "iterms", se.fit = TRUE, uncertainty = "conditional"))
      expect_equal(bam_terms_value(actual$value), targets[[kind]]$optimized, tolerance = 1e-8)
      expect_identical(unname(actual$messages), "iterms reset to terms")
    }
    expect_identical(model, original)
  }
})

test_that("native BAM worker controls are forwarded without changing provider admission", {
  skip_if_not_installed("mgcv")
  native <- get("predict.gam", envir = asNamespace("mgcv"))
  worker <- .msr_bam_saved_worker(native, controls = list(unconditional = TRUE, iterms.type = 2L))
  environment(worker)$.msr_predict_gam_saved <- function(object, args, native) args
  value <- worker(list(object = list(), newdata = data.frame(x = 1), type = "link", se.fit = FALSE))
  expect_identical(value$unconditional, TRUE)
  expect_identical(value$iterms.type, 2L)
  expect_error(worker(list(object = list(), iterms.type = 1L)), "align the native BAM worker controls")
  environment(worker)$.msr_predict_gam_saved <- function(object, args, native)
    stop("retained native worker failure", call. = FALSE)
  expect_error(worker(list(object = list())), "retained native worker failure")
  nocoding <- .msr_bam_saved_worker(native, coding = FALSE)
  environment(nocoding)$.expected_version <- "unqualified-version-fixture"
  expect_error(nocoding(list()), "verify the native GAM provider")
})

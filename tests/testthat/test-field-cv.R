test_that("native Gaussian field CV estimates ranges and retains exact joint scores", {
  skip_if_not_installed("glmmTMB")
  for (covariance in c("exponential", "gaussian")) {
    f <- field_cv_fixture(covariance); x <- f$result
    expect_s3_class(x, "multiScaleR_field_cv")
    expect_true(x$complete)
    expect_identical(x$settings$workers, 1L)
    expect_identical(x$settings$parameter_uncertainty, "not_integrated")
    expect_true(x$settings$observation_noise_included)
    for (name in names(x$folds)) {
      reference <- field_cv_dense_reference(f, name)
      p <- x$predictions[x$predictions$fold == name, ]
      expect_identical(p$row_id, x$folds[[name]]$assessment)
      expect_equal(p$prediction, unname(reference$prediction), tolerance = 1e-8)
      expect_equal(p$observation_predictive_SE^2, unname(diag(reference$covariance)), tolerance = 1e-8)
      expect_equal(x$joint_scores[[name]]$neg_log_score, reference$score, tolerance = 1e-8)
      expect_equal(p$fitted_training_field_variance + p$new_location_field_variance +
        p$observation_noise_variance, unname(diag(reference$covariance)), tolerance = 1e-8)
      marginal_score <- -sum(dnorm(p$observed, p$prediction, p$observation_predictive_SE, log = TRUE))
      expect_gt(abs(marginal_score - reference$score), 1e-6)
      model <- x$fits[[name]]
      nd <- f$data[reference$rows, , drop = FALSE]
      nd$.msr_field_position <- glmmTMB::numFactor(unname(f$xy[reference$rows, ]))
      nd$.msr_field_group <- factor("shared")
      native <- predict(model, nd, type = "response", allow.new.levels = TRUE)
      expect_equal(p$prediction, as.numeric(native), tolerance = 1e-7)
      expect_equal(x$selection_records[[name]]$native_K, 7L)
      expect_equal(x$training[[name]]$parameters$estimated_range_parameters, 1L)
      expect_equal(x$training[[name]]$parameters$integrated_state_penalty, 0L)
      expect_true(all(vapply(x$diagnostics[[name]]$attempts, function(a) isTRUE(a$eligible), logical(1))))
      expect_identical(x$prediction_contracts[[name]]$schema_version, 1L)
      expect_true(x$prediction_contracts[[name]]$uncertainty$observation_noise)
      expect_false(x$prediction_contracts[[name]]$metadata_establishes_eligibility)
      train <- f$data[reference$training, , drop = FALSE]
      permutation <- rev(seq_len(nrow(nd)))
      reordered <- .msr_field_cv_predict(model, y ~ mean + rough + z + offset(log(exposure)),
        train, f$data[reference$rows[permutation], , drop = FALSE],
        unname(f$xy[reference$training, ]), unname(f$xy[reference$rows[permutation], ]), covariance, 512)
      expect_equal(unname(reordered$prediction), rev(p$prediction), tolerance = 1e-8)
      expect_equal(reordered$joint$neg_log_score, reference$score, tolerance = 1e-8)
    }
    expect_output(print(x), "2 of 2 folds scored")
    expect_equal(cv_tab(list(field = x), "neg_log_score")$Score,
      sum(vapply(x$joint_scores, `[[`, 0, "neg_log_score"))/nrow(x$predictions), tolerance = 0)
  }
  cat("GAUSSIAN_FIELD_CV families=1 covariances=2 estimated_range=1 joint_score=1 native_counts=1 parameter_uncertainty=0\n")
})

test_that("field points survive covariance work without losing final native SEs or calls", {
  f <- field_cv_fixture(); model <- f$result$fits[[1L]]
  points <- list.files(f$audit, pattern = "-point[12][.]rds$", recursive = TRUE, full.names = TRUE)
  expect_length(points, 4L)
  packets <- lapply(points, readRDS)
  expect_true(all(vapply(packets, inherits, logical(1), "msr_field_cv_point")))
  expect_false(any(vapply(packets, `[[`, logical(1), "inference_eligible")))
  expect_true(all(file.info(points)$size < 20000))
  expect_equal(unname(sqrt(diag(vcov(model)$cond))),
    unname(summary(model)$coefficients$cond[, "Std. Error"]), tolerance = 1e-12)
  replay <- update(model)
  expect_equal(as.numeric(logLik(replay)), as.numeric(logLik(model)), tolerance = 1e-7)
  expect_equal(vcov(replay)$cond, vcov(model)$cond, tolerance = 1e-6)
})

test_that("held-out and buffer responses cannot train field parameters or predictions", {
  f <- field_cv_fixture(); name <- names(f$folds)[1L]; split <- f$folds[[name]]
  changed <- f$data
  outside <- setdiff(rownames(changed), split$analysis)
  changed$y[match(outside, rownames(changed))] <- changed$y[match(outside, rownames(changed))] + 100
  altered <- cross_validate_field(y ~ mean + rough + z + offset(log(exposure)), changed,
    f$xy, f$folds[name], crs = f$crs, keep_fits = TRUE,
    audit_dir = gdr_test_audit_path("field-cv-leakage-"))
  expect_true(altered$complete)
  expect_equal(altered$training[[name]]$parameters, f$result$training[[name]]$parameters, tolerance = 1e-10)
  expect_equal(altered$training[[name]]$coefficients, f$result$training[[name]]$coefficients, tolerance = 1e-10)
  expect_equal(altered$predictions$prediction,
    f$result$predictions$prediction[f$result$predictions$fold == name], tolerance = 1e-10)
  expect_gt(altered$joint_scores[[name]]$neg_log_score, f$result$joint_scores[[name]]$neg_log_score)
  expect_true(length(f$excluded) > 0)
})

test_that("field CV compares joint targets and selects only from training criteria", {
  exp <- field_cv_fixture()$result; gau <- field_cv_fixture("gaussian")$result
  candidates <- list(exp = exp, gau = gau)
  expect_equal(nrow(cv_tab(candidates, "neg_log_score", by_fold = TRUE)), 4L)
  chosen <- select_cv(candidates, criterion = "AIC", keep_fits = TRUE)
  expect_true(chosen$complete)
  for (name in names(exp$folds)) {
    native <- aic_tab(list(exp = exp$fits[[name]], gau = gau$fits[[name]]), AICc = FALSE, verbose = FALSE, mod_names = names(candidates))
    winner <- native$Modnames[which.min(native$AIC)]
    expect_identical(chosen$selection$chosen[[name]], winner)
    expect_equal(chosen$joint_scores[[name]], candidates[[winner]]$joint_scores[[name]], tolerance = 0)
  }
  tied <- select_cv(list(first = exp, second = exp), criterion = "AIC")
  expect_false(tied$complete)
  expect_true(all(vapply(tied$errors, function(e) grepl("tie", e$message), logical(1))))
  resolved <- select_cv(list(first = exp, second = exp), criterion = "BIC", ties = "first")
  expect_true(resolved$complete)
  expect_true(all(resolved$selection$chosen == "first"))
  changed <- exp; changed$joint_scores[[1L]]$quadratic <- changed$joint_scores[[1L]]$quadratic + .1
  expect_error(cv_tab(list(a = exp, b = changed), "neg_log_score"), "metadata changed")
  expect_error(cv_tab(list(a = exp, b = list()), "rmse"), "matching spatial-field CV")
})

test_that("field CV preserves learned products, factors, offsets and native columns", {
  f <- field_cv_fixture(); d <- f$data
  d$group <- factor(rep(c("a", "b"), length.out = nrow(d)))
  form <- y ~ I(scale(mean, , TRUE)) * group + I(rough^2) + splines::ns(z, df = 3) + offset(log(exposure))
  name <- names(f$folds)[1L]
  x <- cross_validate_field(form, d, f$xy, f$folds[name], crs = f$crs,
    keep_fits = TRUE, audit_dir = gdr_test_audit_path("field-cv-formula-"))
  expect_true(x$complete)
  s <- f$folds[[name]]; tr <- d[match(s$analysis, rownames(d)), ]; te <- d[match(s$assessment, rownames(d)), ]
  basis <- .msr_field_cv_basis(form, tr, x$fits[[name]])
  actual <- .msr_unmarked_fixed_design(basis, te)$X
  center <- mean(tr$mean); spread <- sd(tr$mean)
  ns <- splines::ns(tr$z, df = 3)
  explicit <- ~ I((mean - center)/spread) * group + I(rough^2) +
    splines::ns(z, knots = attr(ns, "knots"), Boundary.knots = attr(ns, "Boundary.knots")) + offset(log(exposure))
  expected <- model.matrix(explicit, te)
  expect_equal(unname(actual), unname(expected), tolerance = 1e-12)
  expect_equal(x$selection_records[[name]]$native_K, ncol(actual) + 3L)
})

test_that("field folds refit actual mean and weighted-SQ scales using training only", {
  f <- field_cv_fixture(kernel = TRUE); x <- f$result
  expect_true(x$complete)
  for (name in names(x$folds)) {
    fit <- x$fits[[name]]
    expect_s3_class(fit, "multiScaleR")
    expect_identical(rownames(fit$kernel_inputs$kernel_dat), x$folds[[name]]$analysis)
    expect_equal(x$selection_records[[name]]$ecological_K, 2L)
    expect_equal(x$selection_records[[name]]$native_K, 7L)
    expect_true(all(is.finite(x$training[[name]]$scale_est$Mean)))
    expect_lte(x$training[[name]]$refits, 200L)
    expect_false(identical(fit$scl_params, f$original$object$scl_params))
  }
  cat("FIELD_CV_KERNELS mean=1 weighted_SQ=1 training_refit=1 heldout_normalization=frozen range_refit=1\n")
})

test_that("invalid geometry, budgets and failed field folds remain explicit", {
  f <- field_cv_fixture()
  run <- function(...) cross_validate_field(y ~ mean + rough, f$data, f$xy,
    f$folds, audit_dir = gdr_test_audit_path("field-cv-guard-"), crs = f$crs, ...)
  expect_error(run(max_memory_mb = .001), "max_memory_mb")
  expect_error(run(max_refits = 1.5), "max_refits")
  expect_error(run(control = list(iter.max = 1.5)), "integers")
  expect_error(cross_validate_field(y ~ mean, f$data, f$xy, f$folds,
    crs = 4326, audit_dir = gdr_test_audit_path("field-cv-crs-")), "projected")
  bad <- list(one = list(analysis = rownames(f$data)[-1L], assessment = rownames(f$data)[1L]))
  expect_error(cross_validate_field(y ~ mean, f$data, f$xy, bad,
    crs = f$crs, audit_dir = gdr_test_audit_path("field-cv-location-")), "complete locations")
  failed <- with_mocked_bindings(run(), .msr_field_cv_fit = function(...) stop("declared fit failure"))
  expect_false(failed$complete)
  expect_true(all(failed$fold_results$status == "failed"))
  expect_true(all(vapply(failed$errors, function(e) grepl("declared fit failure", e$message), logical(1))))
  expect_error(cv_tab(list(failed = failed), "neg_log_score"), "failed folds")
})

test_that("named raw rows override numeric-label guesses without changing fallback", {
  dat <- data.frame(x = c(3, 1, 4, 2), y = c(1, 0, 2, 1), row.names = c("3", "1", "4", "2"))
  ids <- c("2", "4", "1", "3")
  expect_identical(.model_complete_indices(dat, 4L, ids)$complete_idx, c(4L, 3L, 2L, 1L))
  expect_identical(.model_complete_indices(dat, 4L)$complete_idx, c(3L, 1L, 4L, 2L))
  missing <- dat; missing$x[2L] <- NA_real_
  expect_identical(.model_complete_indices(missing, 4L, ids)$complete_idx, c(4L, 2L, 1L))
  raw <- stats::setNames(lapply(1:4, function(i) matrix(c(i, i + 1), ncol = 1L,
    dimnames = list(NULL, "x"))), ids)
  fit <- glm(y ~ x, dat, family = poisson())
  context <- build_opt_context(fit, raw)
  expect_identical(context$complete_idx, c(4L, 3L, 2L, 1L))
})

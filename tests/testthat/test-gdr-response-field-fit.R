test_that("canonical three-bin GDR targets retain native joint derivatives", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fit_fixture()
  for (K in c(400L, max(rowSums(f$data$yd)))) {
    f$data$K <- K; f$native$K <- as.integer(K)
    context <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, 50000)
    program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
    native_data <- f$native
    L <- t(chol(exp(-as.matrix(dist(f$locations)) / 50000)))
    native_data$Z_lambda <- Matrix::Matrix(native_data$Z_lambda %*% L, sparse = TRUE)
    actual <- RTMB::MakeADFun(program$evaluate, f$pars, random = NULL, silent = TRUE)
    native <- TMB::MakeADFun(native_data, f$pars, random = NULL,
      DLL = "unmarked_TMBExports", silent = TRUE)
    key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
    index <- match(key(names(actual$par)), key(names(native$par)))
    expect_false(anyNA(index))
    expect_lt(abs(actual$fn(actual$par) - native$fn(native$par)), 1e-8)
    expect_lt(max(abs(actual$gr(actual$par) - native$gr(native$par)[index])), 1e-7)
    expect_lt(max(abs(actual$he(actual$par) - native$he(native$par)[index, index])), 1e-7)
    expect_identical(program$fixed_coordinate_count, 6L)
    expect_identical(program$latent_coordinate_count, 12L)
    expect_false(program$fit_admission_established)
    rm(actual, native); gc(verbose = FALSE)
  }
  bad <- f$data; bad$breaks <- c(0, 20, 50, 75)
  expect_error(.msr_gdr_response_field_context(bad, f$pars, f$coords, f$crs, 50000),
    "canonical")
})

test_that("optimized spatial GDR reference preserves starts and full covariance", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fit_fixture()
  context <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, 50000)
  key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
  results <- list()
  for (cold in c(FALSE, TRUE)) for (engine in c("TMB", "RTMB")) {
    fit <- gdr_response_field_fit_attempt(f, context, engine, cold)
    expect_identical(as.integer(fit$record$optimizer$convergence), 0L)
    expect_length(fit$record$warnings, 0L)
    expect_true(fit$record$pdHess)
    expect_true(fit$record$fixed_covariance_positive)
    expect_true(fit$record$joint_precision_positive)
    expect_true(fit$record$current_latent_precision_positive)
    expect_lt(fit$record$maximum_score, 1e-3)
    expect_identical(fit$record$coordinate_count, 6L)
    expect_identical(fit$record$latent_count, 12L)
    expect_false(fit$record$public_field_admission)
    results[[paste(engine, cold)]] <- fit
  }
  reference <- results[[1L]]
  for (fit in results[-1L]) {
    index <- match(key(names(reference$fixed)), key(names(fit$fixed)))
    joint_index <- match(key(names(reference$point)), key(names(fit$point)))
    latent_index <- match(key(names(reference$point)[reference$random]),
      key(names(fit$point)[fit$random]))
    expect_false(anyNA(c(index, joint_index, latent_index)))
    expect_lt(abs(reference$record$value - fit$record$value), 1e-6)
    expect_lt(max(abs(reference$fixed - fit$fixed[index])), 1e-3)
    expect_lt(max(abs(reference$point - fit$point[joint_index])), 1e-3)
    expect_equal(reference$covariance, fit$covariance[index, index], tolerance = 1e-3,
      ignore_attr = TRUE)
    expect_equal(reference$precision, fit$precision[latent_index, latent_index],
      tolerance = 1e-3, ignore_attr = TRUE)
    expect_equal(reference$joint_covariance, fit$joint_covariance[joint_index, joint_index],
      tolerance = 1e-3, ignore_attr = TRUE)
    # Two retained location contrasts use complete fixed/latent cross-terms.
    projection <- context$project(f$locations[c(2, 7), ], c("left", "right"), f$crs)
    gradient <- matrix(0, 2L, length(reference$point))
    gradient[, names(reference$point) == "beta_lambda"] <- rbind(c(1, -.5), c(1, .5))
    gradient[, names(reference$point) == "b_lambda"] <- projection$loading
    expected <- gradient %*% reference$point
    actual <- gradient %*% fit$point[joint_index]
    expect_equal(expected, actual, tolerance = 1e-3)
    expected_variance <- gradient %*% reference$joint_covariance %*% t(gradient)
    actual_variance <- gradient %*% fit$joint_covariance[joint_index, joint_index] %*% t(gradient)
    expect_equal(expected_variance, actual_variance, tolerance = 1e-3)
    expect_true(all(diag(expected_variance) > 0))
    expect_identical(projection$conditional_unit_variance, c(0, 0))
  }
  expect_false(context$fit_admission_established)
  expect_false(context$prediction_recipe_qualified)
  cat("GDR_RESPONSE_FIELD_OPTIMIZED starts=2 native_reference=1 fixed=6 latent=12 full_joint=1 public_admission=0\n")
})

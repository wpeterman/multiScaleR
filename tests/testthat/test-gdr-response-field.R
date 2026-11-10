test_that("fixed-range GDR spatial fields retain native finite-K likelihood derivatives", {
  gdr_gaussian_require_RTMB()
  for (K in c(2L, 40L)) for (range in c(10, 50, 150)) {
    f <- gdr_response_field_fixture(K)
    context <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, range)
    program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
    R <- exp(-as.matrix(dist(f$locations)) / range); L <- t(chol(R))
    native_data <- f$captured$data
    native_data$Z_lambda <- Matrix::Matrix(native_data$Z_lambda %*% L, sparse = TRUE)
    for (cold in c(FALSE, TRUE)) {
      pars <- f$pars
      if (cold) { pars$b_lambda[] <- 0; pars$b_dist[] <- 0 }
      actual <- RTMB::MakeADFun(program$evaluate, pars, random = NULL, silent = TRUE)
      native <- TMB::MakeADFun(data = c(list(model = f$captured$model), native_data),
        parameters = pars, random = NULL, DLL = "unmarked_TMBExports", silent = TRUE)
      key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
      index <- match(key(names(actual$par)), key(names(native$par)))
      expect_false(anyNA(index))
      expect_equal(actual$fn(actual$par), native$fn(native$par), tolerance = 1e-9)
      expect_lt(max(abs(actual$gr(actual$par) - native$gr(native$par)[index])), 1e-8)
      expect_lt(max(abs(actual$he(actual$par) - native$he(native$par)[index, index])), 1e-8)
      # Independent original-margin score with explicit Gaussian normalization.
      physical <- pars; physical$b_lambda <- drop(L %*% pars$b_lambda)
      expected <- gdr_gaussian_numeric(f$data, physical) +
        sum(dnorm(physical$b_lambda, 0, exp(pars$lsigma_lambda), log = TRUE)) -
        sum(dnorm(pars$b_lambda, 0, exp(pars$lsigma_lambda), log = TRUE))
      expect_equal(program$evaluate(pars), expected, tolerance = 1e-9)
      sd <- exp(pars$lsigma_lambda)
      physical_prior <- 3 * log(sd) + sum(log(diag(L))) + 3 * log(2 * pi) / 2 +
        sum(solve(L, physical$b_lambda)^2) / (2 * sd^2)
      expect_equal(physical_prior - sum(log(diag(L))),
        -sum(dnorm(pars$b_lambda, 0, sd, log = TRUE)), tolerance = 1e-12)
      rm(actual, native); gc(verbose = FALSE)
    }
    expect_identical(program$fixed_coordinate_count, 15L)
    expect_identical(program$latent_coordinate_count, 5L)
    expect_identical(context$fixed_coordinate_count_change, 0L)
    expect_identical(program$response_fingerprint, context$response_fingerprint)
    expect_identical(context$annual_field_coordinates, 0L)
    expect_true(context$requires_serial_R_callbacks)
    expect_false(context$portable_objective)
    expect_false(context$field$range_estimated)
    expect_false(program$fit_admission_established)
    expect_false(program$uncertainty_evaluated)
    expect_false(program$prediction_recipe_qualified)
    expect_identical(program$field$state_coordinates, "correlation_whitened_abundance")
  }
  cat("GDR_RESPONSE_FIELD fixed_points=6 native_joint=1 cold=1 year=0 public_admission=0\n")
})

test_that("GDR field independence and zero-SD limits retain normalized targets", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fixture()
  base <- .msr_gdr_gaussian_program(f$data, f$pars)
  independent <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, 1e-4)
  program <- independent$prepare(f$data$Xlambda, f$data$row_ids$survey)
  expect_equal(program$evaluate(f$pars), base$evaluate(f$pars), tolerance = 1e-12)
  context <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, 50)
  program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
  difference <- numeric(2)
  for (i in 1:2) {
    sd <- c(1e-4, 1e-7)[i]; p <- f$pars
    p$lsigma_lambda <- log(sd); p$b_lambda <- sd * c(.1, -.2, .3)
    difference[i] <- abs(program$evaluate(p) - base$evaluate(p))
    expect_true(is.finite(program$evaluate(p)))
  }
  expect_lt(difference[2], difference[1])
  expect_lt(difference[2], 1e-5)
  original <- program$evaluate(f$pars)
  X <- f$data$Xlambda; X[, 2L] <- X[, 2L] + seq_len(8) / 100
  changed <- context$prepare(X, f$data$row_ids$survey)
  expect_false(identical(changed$fingerprint, program$fingerprint))
  expect_identical(changed$field_id, program$field_id)
  expect_identical(changed$response_fingerprint, program$response_fingerprint)
  f$data$yd[] <- 0; f$coords[] <- 0
  expect_identical(program$evaluate(f$pars), original)
  expect_error(context$prepare(X[8:1, ], rev(f$data$row_ids$survey)), "exact survey")
})

test_that("GDR spatial Laplace targets retain cold modes and current latent curvature", {
  gdr_gaussian_require_RTMB()
  key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
  for (K in c(2L, 40L)) for (range in c(10, 50, 150)) {
    f <- gdr_response_field_fixture(K)
    context <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, range)
    program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
    L <- t(chol(exp(-as.matrix(dist(f$locations)) / range)))
    native_data <- f$captured$data
    native_data$Z_lambda <- Matrix::Matrix(native_data$Z_lambda %*% L, sparse = TRUE)
    attempts <- list()
    for (cold in c(FALSE, TRUE)) {
      pars <- f$pars
      if (cold) { pars$b_lambda[] <- 0; pars$b_dist[] <- 0 }
      actual <- RTMB::MakeADFun(program$evaluate, pars, random = f$captured$random, silent = TRUE)
      native <- TMB::MakeADFun(data = c(list(model = f$captured$model), native_data),
        parameters = pars, random = f$captured$random, DLL = "unmarked_TMBExports", silent = TRUE)
      value <- actual$fn(actual$par); reference <- native$fn(native$par)
      score <- as.vector(actual$gr(actual$par)); reference_score <- as.vector(native$gr(native$par))
      index <- match(key(names(actual$par)), key(names(native$par)))
      expect_false(anyNA(index))
      expect_lt(abs(value - reference), 1e-8)
      expect_lt(max(abs(score - reference_score[index])), 1e-7)
      # fn/gr above explicitly restore the current conditional modes. Never
      # evaluate latent curvature at the tape's original state vector.
      mode <- actual$env$last.par[actual$env$random]
      reference_mode <- native$env$last.par[native$env$random]
      latent_index <- match(key(names(mode)), key(names(reference_mode)))
      expect_false(anyNA(latent_index))
      expect_lt(max(abs(as.numeric(mode) - as.numeric(reference_mode)[latent_index])), 1e-6)
      precision <- as.matrix(actual$env$spHess(actual$env$last.par, random = TRUE))
      reference_precision <- as.matrix(native$env$spHess(native$env$last.par, random = TRUE))
      expect_silent(chol(precision))
      expect_silent(chol(reference_precision))
      expect_lt(max(abs(precision - reference_precision[latent_index, latent_index])), 1e-6)
      attempts[[length(attempts) + 1L]] <- list(value = value, mode = as.numeric(mode))
      expect_false(program$fit_admission_established)
      expect_false(program$uncertainty_evaluated)
      rm(actual, native); gc(verbose = FALSE)
    }
    expect_lt(abs(attempts[[1]]$value - attempts[[2]]$value), 1e-8)
    expect_lt(max(abs(attempts[[1]]$mode - attempts[[2]]$mode)), 1e-6)
  }
  cat("GDR_RESPONSE_FIELD_LAPLACE points=6 starts=2 current_latent_precision=1 outer_fit_admission=0\n")
})

test_that("new-coordinate GDR field variance stays separate and row aligned", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fixture(); range <- 50
  context <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, range)
  xy <- rbind(f$locations[2L, ], c(430015, 4400020), f$locations[2L, ], c(430030, 4400040))
  ids <- c("known", "new-a", "known", "new-b")
  projection <- context$project(xy, ids, f$crs)
  R <- exp(-as.matrix(dist(f$locations)) / range); L <- t(chol(R))
  distance <- sqrt(outer(xy[, 1L], f$locations[, 1L], "-")^2 +
    outer(xy[, 2L], f$locations[, 2L], "-")^2)
  K <- exp(-distance / range)
  expected <- K %*% solve(R, L)
  variance <- 1 - rowSums(K * t(solve(R, t(K))))
  expect_equal(projection$loading, expected, tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(projection$conditional_unit_variance, pmax(0, variance), tolerance = 1e-12)
  expect_identical(projection$row_ids, ids)
  expect_equal(projection$conditional_unit_variance[c(1, 3)], c(0, 0), tolerance = 0)
  expect_true(all(projection$conditional_unit_variance[c(2, 4)] > 0))
  expect_false(projection$full_prediction_covariance_retained)
  expect_false(projection$fit_admission_established)
  expect_false(projection$public_prediction_admission)
  expect_identical(projection$variance_scale, "multiply_by_exp(2*lsigma_lambda)")
  expect_true(all(c("parameter_uncertainty", "training_field_uncertainty",
    "ecological_scale_uncertainty", "observation_noise", "cross_row_field_covariance") %in%
    projection$omitted))
  permuted <- context$project(xy[4:1, ], rev(ids), f$crs)
  expect_equal(permuted$loading, projection$loading[4:1, ], tolerance = 0)
  expect_equal(permuted$conditional_unit_variance, rev(projection$conditional_unit_variance), tolerance = 0)
  empty <- context$project(matrix(numeric(), 0, 2), character(), f$crs)
  expect_identical(dim(empty$loading), c(0L, 3L))
  expect_identical(empty$omitted, projection$omitted)
  expect_true(context$numeric_projection_portable)
  expect_false(exists("base", environment(context$project), inherits = FALSE))
  expect_false(exists("data", environment(context$project), inherits = FALSE))
  expect_lt(length(serialize(context$project, NULL)), 128 * 1024)
  cluster <- parallel::makePSOCKcluster(1L)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  installed <- find.package("multiScaleR")
  remote <- parallel::clusterCall(cluster, function(project, xy, ids, crs, installed) {
    stopifnot(identical(normalizePath(find.package("multiScaleR"), winslash = "/"),
      normalizePath(installed, winslash = "/")))
    list(projection = project(xy, ids, crs),
      package_version = as.character(utils::packageVersion("multiScaleR")))
  }, context$project, xy, ids, f$crs, installed)[[1L]]
  expect_identical(remote$package_version, as.character(utils::packageVersion("multiScaleR")))
  expect_equal(remote$projection$loading, projection$loading, tolerance = 0)
  expect_equal(remote$projection$conditional_unit_variance, projection$conditional_unit_variance, tolerance = 0)
  expect_error(context$project(xy, ids, 4326), "known projected CRS")
  expect_error(context$project(xy, ids, 32618), "retain the declared CRS")
  expect_error(context$project(xy, ids[-1L], f$crs), "aligned row identities")
  small <- .msr_gdr_response_field_context(f$data, f$pars, f$coords, f$crs, range, max_memory_mb = .01)
  expect_error(small$project(matrix(rep(xy[1, ], each = 100), 100, 2), as.character(1:100), f$crs),
    "exceeds max_memory_mb")
})

test_that("private GDR fields reject ambiguous geometry without changing guards", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fixture()
  prepare <- function(coords = f$coords, data = f$data, range = 50, budget = 512)
    .msr_gdr_response_field_context(data, f$pars, coords, f$crs, range, budget)
  expect_error(prepare(f$coords[8:1, ]), "exact named GDR survey order")
  bad <- f$coords; bad[1L, 1L] <- bad[1L, 1L] + 1
  expect_error(prepare(bad), "one distinct consistent")
  bad <- f$coords; rows <- which(f$data$block == 3L)
  bad[rows, ] <- f$locations[rep(2L, length(rows)), , drop = FALSE]
  expect_error(prepare(bad), "one distinct consistent")
  unused <- f$data; unused$block[unused$block == 3L] <- 2L
  expect_error(prepare(data = unused), "unused states")
  expect_error(prepare(range = 0), "range")
  expect_error(prepare(range = Inf), "finite")
  expect_error(prepare(budget = .001), "exceeds max_memory_mb")
  expect_error(multiScale_optim(NULL, backend = "joint_mixed"))
})

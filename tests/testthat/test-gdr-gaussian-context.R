test_that("one observation context refreshes independent abundance designs", {
  gdr_gaussian_require_RTMB()
  f <- gdr_gaussian_fixture(TRUE, TRUE, 2L)
  counts <- c(cdf = 0L, annulus = 0L)
  original_cdf <- .msr_gdr_poisson_logcdf_program
  original_annulus <- .msr_gdr_annulus_program
  local_mocked_bindings(
    .msr_gdr_poisson_logcdf_program = function(upper) {
      counts['cdf'] <<- counts['cdf'] + 1L
      original_cdf(upper)
    },
    .msr_gdr_annulus_program = function() {
      counts['annulus'] <<- counts['annulus'] + 1L
      original_annulus()
    }, .package = "multiScaleR")
  context <- .msr_gdr_gaussian_context(f$data, f$pars)
  first <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
  X <- f$data$Xlambda
  X[, 2] <- X[, 2] + seq_len(nrow(X)) * .02
  second <- context$prepare(X, f$data$row_ids$survey)
  third <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
  expect_identical(counts, c(cdf = 1L, annulus = 1L))
  expect_identical(first$response_fingerprint, second$response_fingerprint)
  expect_false(identical(first$fingerprint, second$fingerprint))
  expect_identical(first$fingerprint, third$fingerprint)
  expect_identical(first$fixed_coordinate_count, 15L)
  expect_identical(first$latent_coordinate_count, 5L)
  expect_true(first$shared_observation_preparation)
  expect_identical(first$format_version, 2L)
  expect_false(first$fit_admission_established)
  expect_false(first$uncertainty_evaluated)
  expect_false(first$prediction_recipe_qualified)
  original_value <- first$evaluate(f$pars)
  for (program in list(first, second)) {
    data <- f$data
    if (identical(program$fingerprint, second$fingerprint)) data$Xlambda <- X
    expect_lt(abs(program$evaluate(f$pars) - gdr_gaussian_numeric(data, f$pars)), 1e-9)
    fresh <- .msr_gdr_gaussian_program(data, f$pars)
    expect_identical(program$fingerprint, fresh$fingerprint)
    actual <- RTMB::MakeADFun(program$evaluate, f$pars, random = NULL, silent = TRUE)
    native_data <- f$captured$data
    native_data$X_lambda <- data$Xlambda
    independent <- TMB::MakeADFun(data = c(list(model = f$captured$model), native_data),
      parameters = f$pars, random = NULL, DLL = "unmarked_TMBExports", silent = TRUE)
    key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")
    index <- match(key(names(actual$par)), key(names(independent$par)))
    expect_false(anyNA(index))
    expect_lt(abs(actual$fn(actual$par) - independent$fn(independent$par)), 1e-9)
    expect_lt(max(abs(as.vector(actual$gr(actual$par)) -
      as.vector(independent$gr(independent$par))[index])), 1e-8)
    expect_lt(max(abs(actual$he(actual$par) - independent$he(independent$par)[index, index])), 1e-8)
    rm(actual, independent); gc(verbose = FALSE)
  }
  # Subsequent caller edits and preparations cannot alter a previous program.
  X[] <- 0
  f$data$yd[] <- 0
  expect_identical(first$evaluate(f$pars), original_value)
  expect_true(first$fingerprint_is_not_full_provider_provenance)
})

test_that("refreshed numeric contexts require exact original survey and design layouts", {
  gdr_gaussian_require_RTMB()
  f <- gdr_gaussian_fixture(TRUE, TRUE, 2L)
  context <- .msr_gdr_gaussian_context(f$data, f$pars)
  X <- f$data$Xlambda
  ids <- f$data$row_ids$survey
  expect_error(context$prepare(X), "survey identities")
  expect_error(context$prepare(X, rev(ids)), "survey identities")
  expect_error(context$prepare(X, rep(ids[1L], length(ids))), "survey identities")
  for (changed in list(X[-1L, , drop = FALSE], X[, -1L, drop = FALSE],
      X[, rev(seq_len(ncol(X))), drop = FALSE], as.data.frame(X)))
    expect_error(context$prepare(changed, ids), "original abundance columns")
  changed <- X; changed[1L, 1L] <- NA_real_
  expect_error(context$prepare(changed, ids), "original abundance columns")
  changed <- X; changed[1L, 1L] <- Inf
  expect_error(context$prepare(changed, ids), "original abundance columns")
  changed <- X; rownames(changed) <- paste0("changed-row", seq_len(nrow(X)))
  expect_error(context$prepare(changed, ids), "row labels")
  changed <- X; colnames(changed)[1L] <- "changed"
  expect_error(context$prepare(changed, ids), "original abundance columns")
  changed_data <- f$data
  changed_data$yr[2L, ] <- rev(changed_data$yr[2L, ])
  new_context <- .msr_gdr_gaussian_context(changed_data, f$pars)
  expect_false(identical(context$response_fingerprint, new_context$response_fingerprint))
  expect_false(context$fit_admission_established)
  expect_false(context$prediction_recipe_qualified)
  cat("\nGDR_REUSABLE_CONTEXT observations=frozen abundance=refreshed fitting_admission=0\n")
})

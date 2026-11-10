# Existing declared outcomes only. Nonzero physical field coordinates exercise
# range/SD/state cross-derivatives rather than a zero-state special case.
gdr_range_point <- function(f, cold = FALSE) {
  p <- f$pars
  if (cold) for (i in seq_len(nrow(f$constructor_map)))
    p[[f$constructor_map$block[i]]][f$constructor_map$within_block[i]] <- unname(f$cold[i])
  p$b_lambda <- seq(-.3, .3, length.out = length(p$b_lambda))
  p$b_dist <- seq(-.15, .15, length.out = length(p$b_dist))
  p
}

gdr_range_key <- function(x) paste(x, ave(seq_along(x), x, FUN = seq_along), sep = "#")

test_that("physical range joint density and derivatives retain independent references", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fit_fixture()
  distance <- as.matrix(dist(f$locations)); step <- distance[1L, 2L]
  for (cold in c(FALSE, TRUE)) for (range in c(25000, 50000, 100000)) {
    p <- gdr_range_point(f, cold)
    context <- .msr_gdr_range_context(f$data, p, f$coords, f$crs, range)
    program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
    parameters <- context$parameters
    base <- .msr_gdr_gaussian_program(f$data, p)$evaluate
    # Equally spaced collinear nodes have an independent stationary AR1
    # transition factorization. This includes every state/SD/range normalizer.
    reference <- function(q) {
      sd <- exp(q$lsigma_lambda); rho <- exp(-step / exp(q$lrange_lambda))
      u <- q$b_lambda
      ar <- RTMB::dnorm(u[1L], 0, sd, log = TRUE) +
        sum(RTMB::dnorm(u[-1L], rho * u[-length(u)], sd * sqrt(1 - rho^2), log = TRUE))
      base(q[names(p)]) + sum(RTMB::dnorm(u, 0, sd, log = TRUE)) - ar
    }
    actual <- RTMB::MakeADFun(program$evaluate, parameters, random = NULL, silent = TRUE)
    independent <- RTMB::MakeADFun(reference, parameters, random = NULL, silent = TRUE)
    expect_equal(actual$fn(actual$par), independent$fn(independent$par), tolerance = 1e-9)
    expect_equal(actual$gr(actual$par), independent$gr(independent$par), tolerance = 1e-8)
    expect_equal(actual$he(actual$par), independent$he(independent$par), tolerance = 1e-8)
    # An independent compiled native joint section uses whitened coordinates.
    # Apply the exact state Jacobian once to compare physical densities.
    lower <- t(chol(exp(-distance/range)))
    white <- p; white$b_lambda <- drop(forwardsolve(lower, p$b_lambda))
    native_data <- f$native
    native_data$Z_lambda <- Matrix::Matrix(native_data$Z_lambda %*% lower, sparse = TRUE)
    native <- TMB::MakeADFun(native_data, white, random = NULL,
      DLL = "unmarked_TMBExports", silent = TRUE)
    retained <- names(actual$par) != "lrange_lambda"
    index <- match(gdr_range_key(names(actual$par)[retained]), gdr_range_key(names(native$par)))
    expect_false(anyNA(index))
    expect_lt(abs(actual$fn(actual$par) - native$fn(native$par) - sum(log(diag(lower)))), 1e-8)
    J <- diag(sum(retained)); states <- which(names(actual$par)[retained] == "b_lambda")
    J[states, states] <- solve(lower)
    expect_equal(drop(actual$gr(actual$par)[retained]),
      drop(t(J) %*% native$gr(native$par)[index]), tolerance = 1e-7)
    expect_equal(actual$he(actual$par)[retained, retained],
      t(J) %*% native$he(native$par)[index, index] %*% J, tolerance = 1e-7,
      ignore_attr = TRUE)
    expect_identical(program$fixed_coordinate_count, 7L)
    expect_identical(program$latent_coordinate_count, 12L)
    expect_identical(program$range_coordinate_count, 1L)
    expect_false(program$added_range_prior)
    expect_false(program$hyperparameter_jacobian_included)
    expect_false(program$fit_admission_established)
    expect_false(program$CV_admission_established)
    rm(actual, independent, native); gc(verbose = FALSE)
  }
  cat("GDR_RANGE_JOINT ranges=3 starts=2 physical_Jacobian=1 AR1_reference=1 native_reference=1 public_admission=0\n")
})

test_that("physical range Laplace sections agree with the native whitened measure", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fit_fixture()
  for (range in c(25000, 50000, 100000)) {
    p <- gdr_range_point(f)
    context <- .msr_gdr_range_context(f$data, p, f$coords, f$crs, range)
    program <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
    actual <- RTMB::MakeADFun(program$evaluate, context$parameters,
      random = f$random, silent = TRUE)
    lower <- t(chol(exp(-as.matrix(dist(f$locations))/range)))
    white <- p; white$b_lambda[] <- 0; white$b_dist[] <- 0
    native_data <- f$native
    native_data$Z_lambda <- Matrix::Matrix(native_data$Z_lambda %*% lower, sparse = TRUE)
    native <- TMB::MakeADFun(native_data, white, random = f$random,
      DLL = "unmarked_TMBExports", silent = TRUE)
    retained <- names(actual$par) != "lrange_lambda"
    index <- match(gdr_range_key(names(actual$par)[retained]), gdr_range_key(names(native$par)))
    expect_false(anyNA(index))
    # Integrating physical rather than whitened states cancels the Jacobian.
    expect_lt(abs(actual$fn(actual$par) - native$fn(native$par)), 1e-6)
    expect_equal(drop(actual$gr(actual$par)[retained]),
      drop(native$gr(native$par)[index]), tolerance = 1e-6)
    expect_true(all(is.finite(actual$gr(actual$par))))
    actual$gr(actual$par)
    precision <- as.matrix(actual$env$spHess(actual$env$last.par, random = TRUE))
    expect_true(isTRUE(tryCatch({chol(precision); TRUE}, error = function(e) FALSE)))
    expect_false(context$fit_admission_established)
    rm(actual, native); gc(verbose = FALSE)
  }
  cat("GDR_RANGE_LAPLACE ranges=3 native_reference=1 current_latent_precision=1 public_admission=0\n")
})

test_that("range contexts retain geometry, counts and exact observation guards", {
  gdr_gaussian_require_RTMB()
  f <- gdr_response_field_fit_fixture()
  context <- .msr_gdr_range_context(f$data, f$pars, f$coords, f$crs, 50000)
  expect_identical(context$fixed_coordinate_count_change, 1L)
  expect_identical(context$annual_field_coordinates, 0L)
  expect_identical(context$field$state_coordinates, "physical_abundance_field")
  expect_equal(context$parameters$lrange_lambda, log(50000), tolerance = 0)
  expect_false(context$fit_admission_established)
  expect_false(context$prediction_recipe_qualified)
  expect_false(context$CV_admission_established)
  expect_error(.msr_gdr_range_context(f$data, f$pars, f$coords, f$crs, 0), "range")
  expect_error(.msr_gdr_range_context(f$data, f$pars, f$coords, f$crs, Inf), "finite")
  expect_error(.msr_gdr_range_context(f$data, f$pars, f$coords, f$crs, 50000,
    max_memory_mb = .001), "exceeds")
  expect_error(.msr_gdr_range_context(f$data, f$pars, f$coords, 4326, 50000), "projected")
  expect_error(.msr_gdr_range_context(f$data, f$pars, f$coords[nrow(f$coords):1, ], f$crs, 50000), "survey order")
  expect_error(context$prepare(f$data$Xlambda, rev(f$data$row_ids$survey)), "survey identities")
  changed <- f$data$Xlambda; changed[, 2L] <- .8 * changed[, 2L] + .1
  before <- context$prepare(f$data$Xlambda, f$data$row_ids$survey)
  after <- context$prepare(changed, f$data$row_ids$survey)
  expect_false(identical(before$fingerprint, after$fingerprint))
  expect_identical(before$response_fingerprint, after$response_fingerprint)
  expect_false(before$prediction_recipe_qualified)
  expect_true(before$coordinate_count_is_not_AIC_admission)
})

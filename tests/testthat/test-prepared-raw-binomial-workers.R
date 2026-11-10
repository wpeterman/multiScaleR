test_that("a fresh installed worker evaluates frozen raw contexts without native models", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  f <- raw_binomial_fixture()
  expected <- .msr_evaluate_raw_binomial(f$context, f$reference,
    c(sigma_forest = .2, sigma_developed = .4))
  pool <- parallel::makePSOCKcluster(1L); on.exit(parallel::stopCluster(pool))
  configure <- function(paths) { .libPaths(paths); normalizePath(find.package("multiScaleR")) }
  environment(configure) <- baseenv()
  installed <- parallel::clusterCall(pool, configure, .libPaths())[[1L]]
  expect_identical(installed, normalizePath(find.package("multiScaleR")))
  evaluate <- function(context, features) {
    fun <- get(".msr_evaluate_raw_binomial", asNamespace("multiScaleR"))
    point <- fun(context, features, c(sigma_forest = .2, sigma_developed = .4))
    list(value = point$value, coefficients = point$coefficients, target_id = point$target_id,
      accepted = point$coefficient_fit_accepted, finite = point$finite_mle$status,
      selection = point$selection_eligible, uncertainty = point$uncertainty_available)
  }
  environment(evaluate) <- baseenv()
  observed <- parallel::clusterCall(pool, evaluate, f$context, f$reference)[[1L]]
  expect_equal(observed$value, expected$value, tolerance = 1e-12)
  expect_equal(observed$coefficients, expected$coefficients, tolerance = 1e-10)
  expect_identical(observed$target_id, expected$target_id)
  expect_true(observed$accepted)
  expect_identical(observed$finite, "finite")
  expect_false(observed$selection)
  expect_false(observed$uncertainty)
})

test_that("native finite-K loss is separated from the same-bin infinite sum", {
  for (mixture in c("P", "NB")) {
    model <- gdr_guard_fixture()[[mixture]]
    original <- gdr_guard_state(model)
    # Native GDR stores row-major response vectors, not the original matrices.
    expect_identical(environment(model@nllFun)$gd$yDist, as.vector(t(model@data@yDistance)))
    expect_identical(environment(model@nllFun)$gd$yRem, as.vector(t(model@data@yRemoval)))
    actual <- audit_gdr_integration(model, c(5, 400))
    expect_true(actual$stored_point$agrees)
    expect_equal(actual$cutoffs$cutoff, c(40, 5, 400))
    expect_equal(actual$cutoffs$original_cutoff, c(TRUE, FALSE, FALSE))
    expect_gt(actual$cutoffs$truncation_loss[2], 1e-5)
    expect_false(actual$cutoffs$truncation_qualified[2])
    expect_true(actual$cutoffs$truncation_qualified[3])
    expect_true(actual$distance$analytic_equivalence_qualified)
    expect_identical(actual$parameter_count, .msr_parameter_count(model))
    expect_identical(actual$coefficients, unmarked::coef(model))
    expect_false(actual$fitting_admission_established)
    expect_false(actual$density_identification_established)
    expect_identical(actual$uncertainty, "none_read_only_numerical_audit")
    expect_identical(gdr_guard_state(model), original)
    expect_match(actual$provider$native_DLL_sha256, "^[a-f0-9]{64}$")
    expect_match(actual$provider$audit_adapter_sha256, "^[a-f0-9]{64}$")
  }
})

test_that("native units, inactive distance and interval predictors keep their target", {
  base <- gdr_guard_fixture()$P@data
  configurations <- list(
    list(key = "exp", mixture = "P", unit = "m", output = "ha"),
    list(key = "exp", mixture = "NB", unit = "km", output = "kmsq"),
    list(key = "uniform", mixture = "P", unit = "m", output = "ha"),
    list(key = "uniform", mixture = "NB", unit = "km", output = "kmsq"),
    list(key = "halfnorm", mixture = "P", unit = "km", output = "ha"))
  for (configuration in configurations) {
    unit <- configuration$unit
    conversion <- if (unit == "m") 1 else .001
    area <- pi * .075^2 * if (configuration$output == "ha") 100 else 1
    frame <- unmarked::unmarkedFrameGDR(base@yDistance, base@yRemoval,
      siteCovs = unmarked::siteCovs(base), unitsIn = unit,
      obsCovs = data.frame(interval = rep(1:3, 5)),
      dist.breaks = base@dist.breaks * conversion, period.lengths = c(1, 2, 1))
    starts <- c(log(8 / area), .2, if (configuration$mixture == "NB") log(2),
      if (configuration$key != "uniform") log(40 * conversion),
      if (configuration$key == "halfnorm") .1, stats::qlogis(.3), .4)
    model <- unmarked::gdistremoval(lambdaformula = ~I(hab^2), removalformula = ~I(interval == 1),
      distanceformula = if (configuration$key == "halfnorm") ~hab else ~1,
      data = frame, mixture = configuration$mixture, keyfun = configuration$key,
      output = "density", unitsOut = configuration$output, engine = "C", threads = 1,
      K = 400, starts = starts, se = FALSE, control = list(maxit = 0))
    result <- audit_gdr_integration(model)
    expect_true(result$stored_point$agrees)
    expect_true(result$cutoffs$truncation_qualified)
    expect_equal(result$references$native_bins$nll, -as.numeric(unmarked::logLik(model)), tolerance = 1e-9)
    expect_identical(result$units$input, unit)
    expect_identical(result$units$density, configuration$output)
    expect_identical(result$parameter_count, .msr_parameter_count(model))
    if (configuration$key == "exp") {
      expect_match(result$references$native_bins$method, "100_trapezoids", fixed = TRUE)
      expect_gt(abs(result$distance$native_minus_analytic_nll), 1e-9)
      expect_false(result$distance$analytic_equivalence_qualified)
    } else {
      expect_true(result$distance$analytic_equivalence_qualified)
    }
  }
})

test_that("integration audit rejects unsupported layouts and invalid work requests", {
  model <- gdr_guard_fixture()$P
  expect_error(audit_gdr_integration(gdr_guard_fixture()$ZIP), "single-primary")
  for (bad in list(0, -1, 1.5, Inf, NA_real_, c(5, 5), 100001,
      factor("5"), TRUE, matrix(5), "5", 4:40))
    expect_error(audit_gdr_integration(model, bad), "cutoff|work bound")
  original <- gdr_guard_state(model)
  stale <- model
  stale@data@siteCovs$hab <- 2 * stale@data@siteCovs$hab
  expect_error(audit_gdr_integration(stale), "do not match")
  expect_identical(gdr_guard_state(model), original)
})

test_that("a failed native-bin reference retains completed cutoffs and stable reference", {
  local_mocked_bindings(.msr_gdr_audit_native_bins = function(...) stop("native-bin test failure"))
  result <- audit_gdr_integration(gdr_guard_fixture()$P, 400)
  expect_identical(result$references$native_bins$error, "native-bin test failure")
  expect_null(result$references$stable_analytic$error)
  expect_true(is.finite(result$references$stable_analytic$nll))
  expect_true(all(is.finite(result$cutoffs$nll)))
  expect_true(all(is.na(result$cutoffs$truncation_loss)))
  expect_false(any(result$cutoffs$truncation_qualified))
  expect_false(result$distance$analytic_equivalence_qualified)
  expect_true(result$stored_point$agrees)
  expect_false(result$fitting_admission_established)
})

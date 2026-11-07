test_that("GDR selection checks both paired histories and full native penalties", {
  models <- gdr_guard_fixture()
  signatures <- lapply(models, .msr_selection_signature)
  expect_null(.msr_assert_selection_signatures(signatures))
  expect_equal(dim(signatures$P$response), c(5L, 5L))
  expect_equal(signatures$P$gdr_protocol$numPrimary, 1)
  for (model in models)
    expect_equal(.msr_parameter_count(model), length(unmarked::coef(model)))
  model <- models$P
  model@data@yDistance[3, ] <- c(1, 1) # Same total and identical removal history.
  expect_equal(rowSums(model@data@yDistance), rowSums(models$P@data@yDistance))
  expect_error(aic_tab(list(models$P, model), AICc = FALSE), "different retained response")
  expect_error(bic_tab(list(models$P, model)), "different retained response")
  model <- models$P
  model@data@yRemoval[3, ] <- c(0, 1, 1)
  model@data@y <- model@data@yRemoval
  expect_error(aic_tab(list(models$P, model), AICc = FALSE), "different retained response")
  model@data@y[1, 1] <- 2
  expect_error(.msr_selection_signature(model), "inherited removal mirror")
})

test_that("GDR comparison protocols retain bins, units, durations and layouts", {
  model <- gdr_guard_fixture()$P
  cases <- list(breaks = function(x) { x@data@dist.breaks[2] <- 30; x },
    units = function(x) { x@data@unitsIn <- "km"; x },
    durations = function(x) { x@data@period.lengths <- c(1, 1, 2); x },
    survey = function(x) { x@data@survey <- "line"; x })
  for (change in cases)
    expect_error(aic_tab(list(model, change(model)), AICc = FALSE), "different GDR observation protocols")
  signature <- .msr_selection_signature(model)
  other <- signature
  other$gdr_protocol$numPrimary <- 2
  expect_error(.msr_assert_selection_signatures(list(signature, other)), "different GDR observation protocols")
  other <- signature
  other$gdr_protocol$distance_columns <- 1L
  expect_error(.msr_assert_selection_signatures(list(signature, other)), "different GDR observation protocols")
  other$gdr_protocol <- NULL
  expect_error(.msr_assert_selection_signatures(list(signature, other)), "different GDR observation protocols")
  # Key-normalized row ordering is preserved for both response histories.
  reordered <- model
  rows <- c(5, 3, 1, 4, 2)
  reordered@data@yDistance <- model@data@yDistance[rows, , drop = FALSE]
  reordered@data@yRemoval <- model@data@yRemoval[rows, , drop = FALSE]
  reordered@data@y <- reordered@data@yRemoval
  reordered@data@siteCovs <- model@data@siteCovs[rows, , drop = FALSE]
  expect_identical(.msr_selection_signature(reordered), signature)
})

test_that("inactive GDR formulas cannot introduce optimized spatial parameters", {
  model <- gdr_guard_fixture()$P
  model@formlist$phi <- ~annual
  expect_identical(.unmarked_model_predictors(model), "hab")
  expect_error(build_opt_context(model, list(data.frame(annual = 0))), "No raster covariates")
  model@data@numPrimary <- 2
  expect_setequal(.unmarked_model_predictors(model), c("hab", "annual"))
  model <- gdr_guard_fixture()$P
  model@keyfun <- "uniform"
  model@formlist$lambda <- ~1
  model@formlist$dist <- ~hab
  expect_length(.unmarked_model_predictors(model), 0L)
  expect_error(build_opt_context(model, list(data.frame(hab = 0))), "No raster covariates")
  # A native model with no covariate tables still has inspectable ~1 formulas.
  model@data@siteCovs <- NULL
  expect_null(.msr_gdr_native_guard(model))
})

test_that("GDR dot and equation inputs use the active component's training names", {
  model <- gdr_guard_fixture()$P
  model@formlist$lambda <- ~.
  expect_setequal(.unmarked_model_predictors(model), c("hab", "exposure"))
  model@data@yearlySiteCovs <- data.frame(annual = seq_len(5))
  model@formlist$phi <- ~.
  expect_setequal(.unmarked_model_predictors(model), c("hab", "exposure"))
  model@data@numPrimary <- 2
  expect_setequal(.unmarked_model_predictors(model), c("hab", "exposure", "annual"))
  model <- gdr_guard_fixture()$P
  model@formlist$lambda <- ~hab * exposure + I(hab^2) + poly(exposure, 2)
  expect_setequal(.unmarked_model_predictors(model), c("hab", "exposure"))
})

test_that("optimized GDR site columns cannot be shadowed by other native slots", {
  model <- gdr_guard_fixture()$P
  expect_null(.msr_gdr_spatial_guard(model, "hab"))
  for (slot in c("yearlySiteCovs", "obsCovs")) {
    altered <- model
    n <- if (slot == "yearlySiteCovs") 5L else 15L
    methods::slot(altered@data, slot) <- data.frame(hab = seq_len(n))
    expect_error(.msr_gdr_spatial_guard(altered, "hab"), "Optimized GDR site covariates are shadowed")
    expect_null(.msr_gdr_spatial_guard(altered, "other"))
    expect_error(build_opt_context(altered, rep(list(data.frame(hab = 0)), 5)), "Optimized GDR site covariates are shadowed")
    fit <- structure(list(opt_mod = altered, scale_est = data.frame(Mean = 30,
      row.names = "hab")), class = "multiScaleR")
    expect_error(prepare_unmarked_data(fit, data.frame(hab = 1), NULL,
      join_by = data.frame(id = 1)), "Optimized GDR site covariates are shadowed")
  }
})

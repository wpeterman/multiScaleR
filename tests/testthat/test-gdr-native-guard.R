test_that("qualified GDR paths preserve native likelihood and full penalties", {
  fixtures <- gdr_guard_fixture()
  for (mixture in names(fixtures)) {
    model <- fixtures[[mixture]]
    before <- gdr_guard_state(model)
    expect_null(.msr_gdr_native_guard(model))
    expect_identical(.msr_unmarked_prediction_function(model), unmarked::predict)
    expect_equal(.msr_model_loglik(model), as.numeric(unmarked::logLik(model)))
    expect_equal(.msr_parameter_count(model), 4L + as.integer(mixture != "P"))
    expect_identical(gdr_guard_state(model), before)
    # A nonempty native TMB slot identifies its engine. This metadata fixture
    # never represents a TMB fit and is used only to exercise early admission.
    if (mixture != "ZIP") {
      model@TMB <- list(guard_metadata_fixture = TRUE)
      expect_null(.msr_gdr_native_guard(model))
    }
  }
  tab <- aic_tab(fixtures, AICc = FALSE, mod_names = names(fixtures))
  expect_equal(stats::setNames(tab$K, tab$Modnames)[names(fixtures)],
    c(P = 4, NB = 5, ZIP = 5))
})

test_that("GDR TMB ZIP fails before selection, optimization and prediction", {
  model <- gdr_guard_fixture()$ZIP
  model@TMB <- list(guard_metadata_fixture = TRUE)
  before <- gdr_guard_state(model)
  message <- "GDR TMB ZIP models are not qualified"
  expect_error(.msr_gdr_native_guard(model), message)
  expect_error(.msr_model_loglik(model), message)
  expect_error(.msr_parameter_count(model), message)
  expect_error(aic_tab(list(model), AICc = FALSE), message)
  expect_error(bic_tab(list(model)), message)
  expect_error(build_opt_context(model, NULL), message)
  expect_error(.msr_unmarked_prediction_function(model), message)
  fit <- structure(list(opt_mod = model, scale_est = data.frame(Mean = 30,
    row.names = "hab"), scl_params = list(mean = c(hab = 0), sd = c(hab = 1))),
    class = "multiScaleR")
  expect_error(prepare_unmarked_data(fit, data.frame(hab = 1), NULL, NULL), message)
  expect_error(prepare_prediction(fit, component = "lambda"), message)
  expect_error(plot_marginal_effects(fit), message)
  expect_identical(gdr_guard_state(model), before)
})

test_that("GDR active offsets are checked without evaluating formula inputs", {
  base <- gdr_guard_fixture()$P
  for (component in c("lambda", "dist", "rem")) {
    model <- base
    model@formlist[[component]] <- ~hab + offset(exposure)
    message <- paste("Active GDR", component, "offsets are not qualified")
    expect_error(.msr_gdr_native_guard(model), message)
    expect_error(.msr_model_loglik(model), message)
    expect_error(build_opt_context(model, NULL), message)
    expect_error(.msr_unmarked_prediction_function(model), message)
    model@TMB <- list(guard_metadata_fixture = TRUE)
    expect_error(.msr_gdr_native_guard(model), message)
  }
  model <- base
  model@formlist$lambda <- ~. + offset(exposure)
  expect_error(.msr_gdr_native_guard(model), "Active GDR lambda offsets")
  model@formlist$lambda <- ~stop("must not evaluate") + offset(exposure)
  expect_error(.msr_gdr_native_guard(model), "Active GDR lambda offsets")
  model@formlist$lambda <- ~hab * exposure + I(offset(exposure)) + `offset(hab)`
  expect_null(.msr_gdr_native_guard(model))
  model@formlist$lambda <- ~hab + (1 | group) + offset(exposure)
  expect_error(.msr_gdr_native_guard(model), "Active GDR lambda offsets")
  # Inactive native components do not acquire coefficients or restrictions.
  model <- base
  model@formlist$phi <- ~offset(exposure)
  expect_null(.msr_gdr_native_guard(model))
  model@data@numPrimary <- 2
  expect_error(.msr_gdr_native_guard(model), "Active GDR phi offsets")
  model <- base
  model@keyfun <- "uniform"
  model@formlist$dist <- ~offset(exposure)
  expect_null(.msr_gdr_native_guard(model))
  model@formlist$rem <- NULL
  expect_error(.msr_gdr_native_guard(model), "Could not inspect the active GDR rem formula")
})

test_that("GDR guards leave other model families' legitimate offsets alone", {
  model <- stats::glm(y ~ x + offset(z), family = poisson(),
    data = data.frame(y = c(0, 2, 1, 4, 2, 5), x = 1:6, z = rep(0.4, 6)))
  # Formula environments can include this test's changing local variables.
  # Compare the actual fitted data, coefficients, offsets and uncertainty.
  fields <- c("coefficients", "fitted.values", "deviance", "aic", "offset",
    "model", "weights", "prior.weights")
  before <- model[fields]
  covariance <- vcov(model)
  expect_null(.msr_gdr_native_guard(model))
  expect_equal(.msr_model_loglik(model), as.numeric(logLik(model)))
  expect_identical(.msr_parameter_count(model), 2L)
  expect_identical(model[fields], before)
  expect_identical(vcov(model), covariance)
})

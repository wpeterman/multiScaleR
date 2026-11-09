# Two fixed native likelihood constructors, no optimization or recovery panel.
gdr_design_audit_fixture <- local({
  cache <- list()
  function(adjusted = FALSE) {
    key <- as.character(adjusted)
    if (!is.null(cache[[key]])) return(cache[[key]])
    reference <- gdr_guard_fixture()$P@data
    rows <- rep(seq_len(nrow(reference@yDistance)), length.out = 30)
    dat <- data.frame(hab = c(rep(0, 12), seq(.02, .98, length.out = 18)),
      group = factor(rep(c("a", "b", "c"), 10)))
    frame <- unmarked::unmarkedFrameGDR(reference@yDistance[rows, ], reference@yRemoval[rows, ],
      dist.breaks = reference@dist.breaks, unitsIn = "m", period.lengths = reference@period.lengths,
      siteCovs = dat)
    formula <- stats::as.formula(if (adjusted)
      "~splines::ns(hab, df=4, Boundary.knots=c(0,1))" else
      "~splines::ns(hab, knots=c(.25,.6), Boundary.knots=c(0,1))*group", env = baseenv())
    warnings <- character()
    capture <- function(expr) withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    columns <- ncol(capture(stats::model.matrix(formula, dat)))
    model <- capture(unmarked::gdistremoval(lambdaformula = formula, data = frame, output = "abund",
      mixture = "P", K = 40, engine = "C", threads = 1, se = FALSE,
      starts = c(log(8), rep(0, columns - 1L), log(40), qlogis(.3)), control = list(maxit = 0)))
    cache[[key]] <<- list(model = model, data = dat, formula = formula, warnings = warnings)
    cache[[key]]
  }
})

test_that("opt-in native GDR design audits preserve default numerics and all fitted slots", {
  skip_if_not_installed("unmarked")
  skip_if(as.character(packageVersion("unmarked")) != "1.5.2", "Native GDR provider needs qualification")
  f <- gdr_design_audit_fixture(); model <- f$model
  before <- gdr_guard_state(model)
  original <- audit_gdr_integration(model)
  scenarios <- data.frame(hab = c(-.2, 0, .2, 1.2), group = c("c", "a", "b", "a"))
  audited <- audit_gdr_integration(model, design = TRUE, design_data = scenarios)
  expect_identical(audited[setdiff(names(audited), "design")], original)
  expect_identical(gdr_guard_state(model), before)
  expect_false("design" %in% names(original))
  expect_true(audited$design$state$native_matrix_verified)
  expect_identical(audited$design$state$status, "exact_native_matrix_verified")
  expect_false(audited$design$fitting_admission_established)
  expect_false(audited$design$prediction_context_established)
  expect_false(audited$design$numerical_gates_changed)
  gd <- get("gd", environment(model@nllFun), inherits = FALSE)
  for (component in c("lambda", "dist", "rem")) {
    actual <- audited$design$components[[component]]
    native <- gd[[paste0("X_", component)]]
    expect_equal(actual$rows, nrow(native))
    expect_equal(actual$columns, ncol(native))
    expect_equal(actual$rank, qr(native)$rank)
    expect_equal(actual$condition_number, kappa(native, exact = TRUE))
  }
  recipe <- audited$design$state$splines[[1L]]
  expect_identical(recipe$basis, "ns")
  expect_equal(unname(recipe$knots), c(.25, .6))
  expect_equal(recipe$Boundary.knots, c(0, 1))
  expect_equal(recipe$training$zeros, 12)
  expect_identical(audited$design$projection$status, "retained_state_basis_evaluated")
  support <- audited$design$projection$summary$splines[[1L]]
  expect_equal(support$below_boundary, 1)
  expect_equal(support$above_boundary, 1)
  expect_equal(support$outside_training_range, 2)
  # Literal fixed knots provide an independent native scenario matrix oracle.
  native_scenarios <- scenarios
  native_scenarios$group <- factor(native_scenarios$group, levels = levels(f$data$group))
  expected <- stats::model.matrix(f$formula, native_scenarios)
  expected_hash <- .msr_gdr_hash(list(dimensions = dim(expected), columns = colnames(expected), values = as.numeric(expected)))
  expect_identical(audited$design$projection$summary$matrix$matrix_sha256, expected_hash)
  expect_equal(audited$parameter_count, length(unmarked::coef(model)))
  expect_identical(audited$actual_cutoff, original$actual_cutoff)
  invalid_level <- scenarios; invalid_level$group[1L] <- "unseen"
  invalid <- audit_gdr_integration(model, design = TRUE, design_data = invalid_level)
  expect_identical(invalid$design$projection$status, "failed")
  expect_match(invalid$design$projection$error, "new level")
  expect_true(invalid$design$state$native_matrix_verified)
})

test_that("native adjusted knots, failed reconstructions and scenario errors remain explicit", {
  skip_if_not_installed("unmarked")
  skip_if(as.character(packageVersion("unmarked")) != "1.5.2", "Native GDR provider needs qualification")
  f <- gdr_design_audit_fixture(TRUE)
  expect_true(length(f$warnings) > 0)
  expect_true(all(grepl("shoving", f$warnings)))
  warnings <- character()
  audited <- withCallingHandlers(audit_gdr_integration(f$model, design = TRUE,
    design_data = data.frame(hab = c(.1, .9))), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_identical(audited$design$state$warnings, warnings)
  expect_true(length(warnings) > 0)
  expect_true(audited$design$state$native_matrix_verified)
  recipe <- audited$design$state$splines[[1L]]
  expect_match(recipe$requested_call, "df = 4")
  expect_true(all(recipe$knots > 0 & recipe$knots < 1))
  expect_match(recipe$retained_prediction_call, "knots")
  # Separate batches retain the same adjusted training recipe.
  gd <- get("gd", environment(f$model@nllFun), inherits = FALSE)
  one <- suppressWarnings(.msr_gdr_design_audit(gd, c("lambda", "dist", "rem"), f$formula,
    f$data, data.frame(hab = .1)))
  expect_identical(one$state$splines, audited$design$state$splines)
  expect_equal(one$projection$summary$splines[[1L]]$outside_training_range, 0)
  changed <- .msr_gdr_design_audit(gd, c("lambda", "dist", "rem"), ~hab, f$data)
  expect_false(changed$state$native_matrix_verified)
  expect_match(changed$state$error, "exactly reproduce")
  expect_identical(changed$components$lambda, audited$design$components$lambda)
  bad <- suppressWarnings(.msr_gdr_design_audit(gd, c("lambda", "dist", "rem"), f$formula,
    f$data, data.frame(hab = NA_real_)))
  expect_identical(bad$projection$status, "failed")
  expect_type(bad$projection$error, "character")
  expect_true(nzchar(bad$projection$error))
  expect_error(audit_gdr_integration(NULL, design = NA), "design")
  expect_error(audit_gdr_integration(NULL, design_data = f$data), "requires design")
})

test_that("basis diagnostics retain rank failures and native bs extrapolation warnings", {
  data <- data.frame(x = seq(0, 1, length.out = 12))
  formula <- ~splines::bs(x, knots = .4, Boundary.knots = c(0, 1)) + I(x*0)
  native <- stats::model.matrix(formula, data)
  warnings <- character()
  audit <- withCallingHandlers(.msr_gdr_design_audit(list(X_lambda = native), "lambda", formula, data,
    data.frame(x = c(-.2, 1.2))), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  expect_true(audit$components$lambda$column_rank_deficient)
  expect_identical(audit$components$lambda$condition_number, Inf)
  expect_identical(audit$projection$warnings, warnings)
  expect_true(any(grepl("boundary knots", warnings)))
  expect_equal(audit$projection$summary$splines[[1L]]$below_boundary, 1)
  expect_equal(audit$projection$summary$splines[[1L]]$above_boundary, 1)
  transformed <- ~splines::ns(I(x^2), df = 3)
  audit <- .msr_gdr_design_audit(list(X_lambda = stats::model.matrix(transformed, data)), "lambda", transformed, data)
  expect_true(audit$state$native_matrix_verified)
  expect_identical(audit$state$splines[[1L]]$training$status, "nonliteral_input_not_recomputed")
  intercept <- .msr_gdr_design_audit(list(X_lambda = matrix(1, 5, 1,
    dimnames = list(NULL, "(Intercept)"))), "lambda", ~1, NULL)
  expect_true(intercept$state$native_matrix_verified)
  nested <- ~splines::ns(base::scale(x), df = 3)
  native_nested <- stats::model.matrix(nested, data)
  scenarios <- data.frame(x = c(-1, .25, 2))
  nested_audit <- .msr_gdr_design_audit(list(X_lambda = native_nested), "lambda", nested, data, scenarios)
  trained <- splines::ns(as.numeric(base::scale(data$x)), df = 3)
  oracle <- cbind(1, stats::predict(trained, (scenarios$x - mean(data$x))/stats::sd(data$x)))
  colnames(oracle) <- colnames(native_nested)
  expect_true(nested_audit$state$native_matrix_verified)
  expect_identical(nested_audit$projection$status, "retained_state_basis_evaluated")
  expect_identical(nested_audit$projection$summary$matrix$matrix_sha256,
    .msr_gdr_hash(list(dimensions = dim(oracle), columns = colnames(oracle), values = as.numeric(oracle))))
  data$group <- rep(c("a", "b", "c"), 4)
  factored <- ~splines::ns(x, df = 3)*factor(group)
  native_factor <- stats::model.matrix(factored, data)
  factor_audit <- .msr_gdr_design_audit(list(X_lambda = native_factor), "lambda", factored, data,
    data.frame(x = c(.2, .8), group = c("c", "a")))
  expect_identical(factor_audit$projection$status, "retained_state_basis_evaluated")
  expect_identical(factor_audit$projection$summary$matrix$column_names, colnames(native_factor))
  contrast_data <- data; contrast_data$group <- factor(contrast_data$group)
  contrast_formula <- ~splines::ns(x, df = 3)*group
  contrast_native <- stats::model.matrix(contrast_formula, contrast_data,
    contrasts.arg = list(group = stats::contr.sum(3)))
  contrast_audit <- .msr_gdr_design_audit(list(X_lambda = contrast_native), "lambda",
    contrast_formula, contrast_data, data.frame(x = c(.2, .8), group = c("c", "a")))
  expect_true(contrast_audit$state$native_matrix_verified)
  expect_identical(contrast_audit$projection$summary$matrix$column_names, colnames(contrast_native))
  expect_false(any(vapply(audit$components, is.matrix, logical(1))))
  limited <- .msr_gdr_design_summary(matrix(1, 2, 257))
  expect_identical(limited$status, "work_limit_exceeded")
  expect_false(any(c("rank", "matrix_sha256") %in% names(limited)))
  cohort <- function(x) x - mean(x)
  custom <- ~splines::ns(x, df = 3) + cohort(x)
  guarded <- .msr_gdr_design_audit(list(X_lambda = stats::model.matrix(custom, data)),
    "lambda", custom, data, data.frame(x = c(.1, .9)))
  expect_identical(guarded$projection$status, "failed")
  expect_match(guarded$projection$error, "qualified fixed training")
  cat("\nGDR_DESIGN_AUDIT native_matrix=verified knots=retained admission=0\n")
})

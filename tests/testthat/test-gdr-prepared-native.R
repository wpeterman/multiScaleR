test_that("prepared GDR binds trained designs and rebuilds the same native objective", {
  f <- gdr_trained_fixture()
  settings <- list(method = "BFGS", se = TRUE, threads = 1L,
    control = list(maxit = 2500L, reltol = 1e-10))
  for (kind in c("TMB", "mixed")) {
    model <- f$models[[kind]]
    raw <- unmarked:::get_covariates(model, "dist")
    keys <- data.frame(survey = raw$survey)
    prepared <- .msr_gdr_prepare_native(model, keys, "survey", "minutes", settings)
    input <- prepared$payload
    expect_identical(input$event_keys, keys)
    expect_identical(input$native, .msr_gdr_TMB_preparation(model))
    expect_identical(input$native$K, model@K)
    expect_identical(input$random_index, model@TMB$env$random)
    expect_identical(input$rows$rem$event, rep(seq_len(nrow(keys)), each = 3L))
    expect_identical(input$rows$rem$interval, rep(1:3, nrow(keys)))
    expect_identical(names(input$recipes), names(.msr_gdr_active_formulas(model)))
    for (type in names(input$recipes)) {
      raw <- unmarked:::get_covariates(model, type)
      expect_identical(environment(input$recipes[[type]]$terms), baseenv())
      X <- .msr_unmarked_fixed_design(input$recipes[[type]], raw)$X
      expect_equal(X, input$native[[paste0("X_", type)]], tolerance = 1e-12, ignore_attr = TRUE)
      group <- input$groups[[type]]
      if (!is.null(group))
        expect_identical(group$ids[group$row_coordinate], as.character(raw[[group$variable]]))
    }
    if (kind == "mixed") expect_identical(input$random, c("b_lambda", "b_dist"))
    restored <- unserialize(serialize(prepared, NULL, version = 2L))
    rebuilt <- .msr_gdr_rebuild_native(restored)
    expect_identical(rebuilt$input_fingerprint, prepared$fingerprint)
    expect_identical(rebuilt$parameter_count, .msr_parameter_count(model))
    expect_false(rebuilt$inference_eligible)
    expect_identical(rebuilt$uncertainty, "none_objective_reconstruction_only")
    actual <- rebuilt$objective
    expect_false(identical(actual$env, model@TMB$env))
    expect_identical(actual$env$data, model@TMB$env$data)
    native <- TMB::MakeADFun(model@TMB$env$data, model@TMB$env$parameters,
      random = unique(names(model@TMB$env$last.par)[model@TMB$env$random]),
      DLL = "unmarked_TMBExports", silent = TRUE)
    at <- model@opt$par
    expect_lt(abs(actual$fn(at) - model@opt$value), 1e-7)
    expect_lt(max(abs(actual$gr(at) - model@TMB$gr(at))), 1e-7)
    nearby <- at
    nearby[1L] <- nearby[1L] + .01
    for (point in list(at, nearby)) {
      expect_equal(actual$fn(point), native$fn(point), tolerance = 1e-12)
      expect_equal(actual$gr(point), native$gr(point), tolerance = 1e-12)
      if (length(input$random)) {
        expect_equal(actual$env$spHess(actual$env$last.par, random = TRUE),
          native$env$spHess(native$env$last.par, random = TRUE), tolerance = 1e-12)
      }
    }
    expect_identical(model@data, f$snapshots[[kind]]$data)
    expect_identical(model@estimates, f$snapshots[[kind]]$estimates)
    expect_identical(model@AIC, f$snapshots[[kind]]$AIC)

    expect_error(.msr_gdr_prepare_native(model, keys[rev(seq_len(nrow(keys))), , drop = FALSE],
      "survey", "minutes", settings), "verified survey order")
    fake <- data.frame(guessed = seq_len(nrow(keys)))
    expect_error(.msr_gdr_prepare_native(model, fake, "guessed", "minutes", settings), "retained in the original")
    expect_error(.msr_gdr_prepare_native(model, keys, "survey", "minutes",
      c(settings, list(random.start = "zero"))), "original inner-mode policy")
    changed <- restored
    changed$payload$native$K <- if (is.integer(input$native$K)) as.double(input$native$K) else as.integer(input$native$K)
    expect_error(.msr_gdr_rebuild_native(changed), "input identity changed")
    changed <- restored
    changed$payload$recipes$lambda$columns <- rev(changed$payload$recipes$lambda$columns)
    expect_error(.msr_gdr_rebuild_native(changed), "input identity changed")
    if (kind == "mixed") {
      changed <- restored
      changed$payload$groups$lambda$ids <- rev(changed$payload$groups$lambda$ids)
      expect_error(.msr_gdr_rebuild_native(changed), "input identity changed")
    }
  }
})

test_that("group identities follow verified native columns and membership", {
  raw <- data.frame(block = factor(c("B", "A", "B"), levels = c("B", "unused", "A")))
  Z <- Matrix::sparseMatrix(i = 1:3, j = c(2, 1, 2), x = 1, dims = c(3, 2),
    dimnames = list(NULL, c("A", "B")))
  group <- .msr_gdr_prepared_groups(~1 + (1 | block), raw, Z)
  expect_identical(group$ids, c("A", "B"))
  expect_identical(group$row_coordinate, c(2L, 1L, 2L))
  colnames(Z) <- rev(colnames(Z))
  expect_error(.msr_gdr_prepared_groups(~1 + (1 | block), raw, Z), "do not match")
  expect_error(.msr_gdr_prepared_groups(~1 + (block | block), raw, Z), "bare random intercept")
})

test_that("transport compacts only recognized native GAM warning call data", {
  skip_if_not_installed("mgcv")
  native <- get("predict.gam", asNamespace("mgcv"))
  model <- structure(list(payload = seq_len(1000L)), class = "gam")
  rows <- data.frame(z = seq_len(20L))
  call <- as.call(c(list(native), list(object = model, newdata = rows,
    type = "link", se.fit = TRUE, unconditional = TRUE, iterms.type = 2L)))
  original <- simpleWarning("Native warning context fixture", call = call)
  original$code <- 17L; class(original) <- c("contextFixtureWarning", class(original))
  result <- .msr_bam_warning_context(original, native, as.character(packageVersion("mgcv")))
  expect_identical(conditionMessage(result), conditionMessage(original))
  expect_identical(class(result), class(original)); expect_identical(result$code, original$code)
  expect_identical(conditionCall(original), call)
  compact <- conditionCall(result)
  expect_identical(compact[[1L]], quote(mgcv::predict.gam))
  expect_identical(compact[["object"]], quote(.native_BAM_model))
  expect_identical(compact[["newdata"]], quote(.native_prediction_rows))
  for (name in c("type", "se.fit", "unconditional", "iterms.type"))
    expect_identical(compact[[name]], call[[name]])
  provenance <- attr(result, "multiScaleR_native_warning")
  expect_identical(provenance$provider, "mgcv")
  expect_identical(provenance$method, "predict.gam")
  expect_identical(provenance$version, as.character(packageVersion("mgcv")))
  expect_identical(provenance$symbolic_data_arguments, c("object", "newdata"))
  expect_identical(.msr_bam_warning_context(result, native), result)
  for (other in list(NULL, quote(min(assign)), as.call(c(list(function(...) NULL),
      list(object = model, newdata = rows))))) {
    condition <- simpleWarning("Unchanged diagnostic", call = other)
    expect_identical(.msr_bam_warning_context(condition, native), condition)
  }
})

test_that("actual BAM worker warnings retain diagnostics without transporting fitted data", {
  skip_if_not_installed("mgcv")
  # This cache is normally shared with the facet test. Isolated-file execution
  # loads its definition only, not the tests, with the same bounded fixtures.
  fixture <- get0("native_bam_facets_fixture", envir = parent.frame(), mode = "function", inherits = TRUE)
  if (is.null(fixture)) {
    e <- new.env(parent = parent.frame())
    eval(parse(testthat::test_path("test-native-bam-facets.R"))[[1L]], e)
    fixture <- e$native_bam_facets_fixture
  }
  f <- fixture("gaussian"); native <- get("predict.gam", asNamespace("mgcv"))
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old), add = TRUE)
  for (kind in c("ordinary", "numeric_ordinary")) {
    model <- f$models[[kind]]$model
    for (name in c("Sl", "qrx", "R", "F", "Ve", "Vc", "G", "residuals", "fitted.values", "linear.predictors"))
      model[[name]] <- NULL
    model$model <- model$model[1:2, , drop = FALSE]
    worker <- .msr_bam_saved_worker(native, list(unconditional = TRUE),
      coding = length(model$contrasts) > 0L)
    sizes <- integer(); markers <- list()
    for (copies in c(1L, 6L)) {
      rows <- f$data[rep(seq_len(nrow(f$data)), copies), , drop = FALSE]
      options(contrasts = c("contr.sum", "contr.poly"))
      expected <- native(model, rows, type = "link", se.fit = TRUE, unconditional = FALSE)
      options(contrasts = c("contr.treatment", "contr.poly"))
      observed <- list()
      actual <- withCallingHandlers(worker(list(object = model, newdata = rows, type = "link", se.fit = TRUE)),
        warning = function(w) { observed[[length(observed) + 1L]] <<- w; invokeRestart("muffleWarning") })
      marker <- attr(actual$fit, ".msr_bam_delegate")
      expect_length(observed, 1L); expect_length(marker$warnings, 1L)
      expect_match(conditionMessage(observed[[1L]]), "covariance not available", fixed = TRUE)
      expect_identical(conditionMessage(marker$warnings[[1L]]), conditionMessage(observed[[1L]]))
      expect_identical(class(marker$warnings[[1L]]), class(observed[[1L]]))
      original_call <- conditionCall(observed[[1L]])
      expect_true(inherits(original_call[["object"]], "gam"))
      expect_identical(original_call[["newdata"]], rows)
      compact <- conditionCall(marker$warnings[[1L]])
      expect_identical(compact[[1L]], quote(mgcv::predict.gam))
      expect_identical(compact[["object"]], quote(.native_BAM_model))
      expect_identical(compact[["newdata"]], quote(.native_prediction_rows))
      expect_true(compact[["unconditional"]])
      sizes <- c(sizes, length(serialize(marker, NULL)))
      markers[[length(markers) + 1L]] <- marker
      expect_lt(tail(sizes, 1L), 16384L)
      for (name in c("fit", "se.fit")) attr(actual[[name]], ".msr_bam_delegate") <- NULL
      expect_equal(actual, expected, tolerance = 1e-8)
    }
    # Equal row IDs can use different R serialization representations. Verify
    # semantic context identity, allowing only the declared origin row count.
    # Each actual payload still has the independent16KiB bound above.
    for (i in seq_along(markers)) markers[[i]]$origin$rows <- 0L
    expect_identical(markers[[1L]], markers[[2L]])
  }
})

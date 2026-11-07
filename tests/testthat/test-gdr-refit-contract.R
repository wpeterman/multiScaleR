test_that("C settings do not re-evaluate mutable call variables", {
  model <- gdr_guard_fixture()$P
  settings <- .msr_gdr_C_settings(model)
  expect_identical(settings$engine, "C")
  expect_identical(settings$method, "BFGS")
  expect_false(settings$se)
  expect_equal(settings$K, model@K)
  expect_identical(settings$coefficient_names, names(unmarked::coef(model)))
  changed_call <- model
  changed_call@call$method <- as.name("changed_caller_method")
  changed_call@call$K <- 999
  changed_call@call$se <- TRUE
  unchanged <- .msr_gdr_C_settings(changed_call)
  expect_identical(unchanged$method, "BFGS")
  expect_equal(unchanged$K, model@K)
  expect_false(unchanged$se)
  altered <- model
  altered@K <- model@K + 1
  expect_error(.msr_gdr_C_settings(altered), "differs from its likelihood")
  altered <- model
  altered@TMB <- list(guard_metadata_fixture = TRUE)
  expect_error(.msr_gdr_C_settings(altered), "native C engine")
})

test_that("invalid starts and existing audit evidence fail without new fits", {
  model <- gdr_contract_model()
  cold <- unmarked::coef(model)
  cold[] <- 0
  expect_error(make_gdr_refitter(gdr_guard_fixture()$P, list(cold = cold), tempfile()), "positive integer maxit")
  expect_error(make_gdr_refitter(model, list(cold = unmarked::coef(model)), tempfile()), "independent")
  bad <- cold
  names(bad)[1] <- "wrong"
  expect_error(make_gdr_refitter(model, list(cold = bad), tempfile()), "exactly the finite named")
  expect_error(make_gdr_refitter(model, list(a = cold, b = cold), tempfile()), "distinct from one another")
  path <- tempfile("gdr-existing-")
  dir.create(path)
  saveRDS(list(retain = TRUE), file.path(path, "previous.rds"))
  expect_error(make_gdr_refitter(model, list(cold = cold), path), "new or empty")
  expect_true(readRDS(file.path(path, "previous.rds"))$retain)
  original <- tempfile("gdr-exact-record-", fileext = ".rds")
  saveRDS(list(retain = TRUE), original)
  expect_error(.msr_gdr_atomic_audit(list(retain = FALSE), original), class = "msr_gdr_audit_error")
  expect_true(readRDS(original)$retain)
})

test_that("nonconverged highest likelihood and invalid data retain complete audits", {
  model <- gdr_contract_model() # One iteration deliberately cannot fit this case.
  cold <- unmarked::coef(model)
  cold[] <- 0
  refit <- make_gdr_refitter(model, list(cold = cold), gdr_test_audit_path("gdr-failure-"))
  error <- tryCatch(refit(model, model@data), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_match(conditionMessage(error), "highest declared native likelihood")
  expect_true(file.exists(error$audit_file))
  expect_length(error$audit$attempts, 2L)
  expect_false(any(vapply(error$audit$attempts, `[[`, logical(1), "converged")))
  expect_true(all(file.exists(vapply(error$audit$attempts, `[[`, character(1), "audit_file"))))
  invalid <- model@data
  invalid@yDistance[1, 1] <- 1
  error <- tryCatch(refit(model, invalid), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_match(conditionMessage(error), "changed histories")
  expect_length(error$audit$attempts, 0L)
  expect_true(file.exists(error$audit_file))
  error <- tryCatch(refit(model, model@data, context = FALSE), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_match(conditionMessage(error), "context must contain")
})

test_that("bad requested covariance cannot be labeled available", {
  model <- gdr_guard_fixture()$P
  expect_identical(.msr_gdr_variance_status(model, FALSE)$status, "not_requested")
  model@opt$hessian <- matrix(NA_real_, 4, 4)
  expect_identical(.msr_gdr_variance_status(model, TRUE)$status, "unusable")
})

test_that("restored provider mismatches fail with an audit before fitting", {
  model <- gdr_contract_model()
  cold <- unmarked::coef(model)
  cold[] <- 0
  refit <- make_gdr_refitter(model, list(cold = cold), gdr_test_audit_path("gdr-provider-code-"))
  environment(refit)$source$native_function_sha256 <- "changed-contract-fixture"
  error <- tryCatch(refit(model, model@data), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_match(conditionMessage(error), "qualified provider or adapter")
  expect_length(error$audit$attempts, 0L)
  expect_match(error$audit$provider_check$native_function_sha256, "^[[:xdigit:]]{64}$")
  expect_true(file.exists(error$audit_file))
  refit <- make_gdr_refitter(model, list(cold = cold), gdr_test_audit_path("gdr-provider-files-"))
  # Exercise the restored-process branch without altering any provider binding.
  environment(refit)$counter$provider_pid <- 0L
  environment(refit)$source$native_DLL_sha256 <- "changed-contract-fixture"
  error <- tryCatch(refit(model, model@data), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_match(conditionMessage(error), "native files differ")
  expect_length(error$audit$attempts, 0L)
  expect_true(error$audit$provider_check$files_rechecked)
  expect_true(file.exists(error$audit_file))
})

test_that("native mixture labels are distinct without changing likelihood penalties", {
  models <- gdr_guard_fixture()
  labels <- vapply(models, .msr_selection_family_label, character(1))
  expect_length(unique(labels), 3L)
  expect_match(labels[["P"]], "mixture=P", fixed = TRUE)
  expect_match(labels[["NB"]], "mixture=NB", fixed = TRUE)
  expect_match(labels[["ZIP"]], "mixture=ZIP", fixed = TRUE)
  table <- aic_tab(models, AICc = FALSE)
  expect_length(unique(table$Modnames), 3L)
  expect_equal(sort(table$K), c(4, 5, 5))
  model <- models$P
  model@formlist$lambda <- hab ~ exposure
  expect_identical(.unmarked_model_predictors(model), "exposure")
})

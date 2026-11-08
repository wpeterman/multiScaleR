test_that("optimizer checkpoints retain native diagnostics without executable objects", {
  # Pure metadata fixture, not an optimizer result or a qualified fit.
  native <- list(par = c(a = Inf, b = NaN), value = Inf,
    counts = c("function" = 12L, gradient = 4L), convergence = 1L,
    message = "retained failed optimizer", hessian = function() stop("do not call"),
    large_unused = new.env(parent = emptyenv()))
  record <- .msr_gdr_optimizer_checkpoint(native)
  expect_identical(record, native[c("par", "value", "counts", "convergence", "message")])
  expect_identical(unserialize(serialize(record, NULL, version = 2)), record)
  expect_false(any(vapply(record, function(x) is.function(x) || is.environment(x), logical(1))))
  expect_error(.msr_gdr_optimizer_checkpoint(list(par = 1)), "incomplete")
  for (field in c("par", "value", "counts", "convergence", "message")) {
    invalid <- native
    invalid[[field]] <- function() stop("do not call")
    expect_error(.msr_gdr_optimizer_checkpoint(invalid), "plain numeric")
  }
  invalid <- native
  invalid$par <- matrix(1:4, 2)
  expect_error(.msr_gdr_optimizer_checkpoint(invalid), "plain numeric")
  invalid <- native
  attr(invalid$par, "payload") <- new.env(parent = emptyenv())
  expect_error(.msr_gdr_optimizer_checkpoint(invalid), "plain numeric")
  invalid <- native
  names(invalid)[2] <- "par"
  expect_error(.msr_gdr_optimizer_checkpoint(invalid), "incomplete")
})

test_that("native C refit checkpoints bind each attempt to its data and settings", {
  for (mixture in c("P", "NB", "ZIP")) {
    fixture <- gdr_refit_fixture(mixture)
    before <- gdr_guard_state(fixture$model)
    path <- gdr_test_audit_path(paste0("gdr-checkpoint-", mixture, "-"))
    callback <- make_gdr_refitter(fixture$model, fixture$starts, path,
      provenance = list(scope = "original native C regression fixtures"))
    contract <- readRDS(file.path(path, "contract.rds"))
    checksum <- contract$contract_sha256
    contract$contract_sha256 <- NULL
    expect_identical(.msr_gdr_hash(contract), checksum)
    expect_identical(contract$original_data_sha256, .msr_gdr_hash(fixture$model@data))
    expect_identical(contract$original_design_sha256,
      .msr_gdr_design_hash(fixture$model@data, fixture$model@formlist))
    expect_identical(contract$original_coordinates, gdr_start_values(fixture$model))
    expect_identical(contract$original_optimizer, .msr_gdr_optimizer_checkpoint(fixture$model@opt))
    expect_identical(contract$controls$dots$control, list(maxit = 2500, reltol = 1e-10))
    expect_true(contract$controls$se)
    changed <- contract
    changed$controls$se <- FALSE
    expect_false(identical(.msr_gdr_hash(changed), checksum))
    changed <- contract
    changed$original_design_sha256 <- "different design"
    expect_false(identical(.msr_gdr_hash(changed), checksum))
    final <- callback(fixture$model, fixture$data, list(return_model = TRUE, covs = "hab"))
    audit <- attr(final, "msr_gdr_refit_audit")
    expect_identical(audit$contract_sha256, checksum)
    expect_identical(audit$input_sha256, .msr_gdr_hash(fixture$data))
    expect_identical(audit$input_design_sha256,
      .msr_gdr_design_hash(fixture$data, fixture$model@formlist))
    expect_false(identical(audit$input_design_sha256, contract$original_design_sha256))
    expect_identical(audit$input_context_sha256, .msr_gdr_hash(list(
      contract_sha256 = checksum, input_sha256 = audit$input_sha256,
      input_design_sha256 = audit$input_design_sha256, spatial_columns = "hab")))
    input_path <- sub("[.]rds$", "-input.rds", audit$path)
    input <- readRDS(input_path)
    expect_identical(input$input_context_sha256, audit$input_context_sha256)
    expect_length(input$attempts, 0L)
    expect_identical(audit$contract_file_sha256,
      digest::digest(file = audit$contract_file, algo = "sha256"))
    expect_length(audit$attempts, 4L)
    for (attempt in audit$attempts) {
      request <- readRDS(attempt$request_file)
      completion <- readRDS(attempt$audit_file)
      expect_identical(request$status, "started")
      expect_null(request$optimizer)
      expect_false(request$eligible)
      expect_identical(request$start, attempt$start)
      expect_identical(request$input_context_sha256, audit$input_context_sha256)
      expect_identical(request$contract_sha256, checksum)
      expect_identical(completion, attempt)
      expect_identical(attempt$status, "returned")
      expect_identical(as.numeric(attempt$optimizer$par), as.numeric(attempt$coefficients))
      expect_identical(attempt$optimizer$convergence, attempt$convergence)
      expect_equal(-attempt$optimizer$value, attempt$logLik, tolerance = 0)
    }
    last <- audit$attempts[[4L]]
    expect_identical(last$optimizer, .msr_gdr_optimizer_checkpoint(final@opt))
    expect_true(last$SE_requested)
    expect_identical(last$covariance$status, "available")
    expect_true(is.matrix(final@opt$hessian))
    expect_identical(.msr_parameter_count(final), .msr_parameter_count(fixture$model))
    expect_identical(gdr_guard_state(fixture$model), before)
    expect_identical(unserialize(serialize(audit, NULL, version = 2)), audit)
    expect_lt(length(serialize(audit, NULL, version = 2)), 65536L)
  }
  cat("GDR_REFIT_CHECKPOINT mixtures=3 preentry_requests=12 final_SE=original\n")
})

test_that("interrupted requests survive before native entry", {
  fixture <- gdr_refit_fixture("P")
  path <- gdr_test_audit_path("gdr-checkpoint-interrupted-")
  writer <- .msr_gdr_atomic_audit
  # Audit-writer stub only. Stop after persisting the request, before any fit.
  interrupt <- function(value, path) {
    result <- writer(value, path)
    if (grepl("-attempt-1-request[.]rds$", path)) stop("declared pre-entry interruption")
    result
  }
  testthat::with_mocked_bindings({
    callback <- make_gdr_refitter(fixture$model, fixture$starts, path)
    expect_error(callback(fixture$model, fixture$data, list(covs = "hab")),
      "declared pre-entry interruption")
  }, .msr_gdr_atomic_audit = interrupt, .package = "multiScaleR")
  expect_identical(.msr_gdr_atomic_audit, writer)
  files <- list.files(path, full.names = TRUE)
  request_path <- files[grepl("-attempt-1-request[.]rds$", files)]
  expect_length(request_path, 1L)
  request <- readRDS(request_path)
  expect_identical(request$status, "started")
  expect_null(request$optimizer)
  expect_false(file.exists(request$audit_file))
  expect_match(request$input_context_sha256, "^[[:xdigit:]]{64}$")
  expect_match(request$contract_sha256, "^[[:xdigit:]]{64}$")
})

test_that("changed retained settings fail with an audit before native entry", {
  fixture <- gdr_refit_fixture("P")
  path <- gdr_test_audit_path("gdr-checkpoint-settings-")
  callback <- make_gdr_refitter(fixture$model, fixture$starts, path)
  environment(callback)$settings$se <- FALSE
  error <- tryCatch(callback(fixture$model, fixture$data, list(covs = "hab")), error = identity)
  expect_s3_class(error, "msr_gdr_refit_error")
  expect_match(conditionMessage(error), "retained GDR contract or evaluated settings changed")
  expect_length(error$audit$attempts, 0L)
  expect_true(file.exists(error$audit_file))
  expect_match(error$audit$contract_sha256, "^[[:xdigit:]]{64}$")
  expect_length(list.files(path, pattern = "-attempt-"), 0L)
})

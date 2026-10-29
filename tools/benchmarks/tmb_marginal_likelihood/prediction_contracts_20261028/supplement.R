library(testthat)
library(multiScaleR)
task_study <- Sys.getenv("MSR_CONTRACT_STUDY")
task_expected_package <- normalizePath(Sys.getenv("MSR_CONTRACT_INSTALL"), winslash = "/")
stopifnot(identical(normalizePath(find.package("multiScaleR"), winslash = "/"), task_expected_package))
task_output <- file.path(task_study, "supplement05.json")
stopifnot(!file.exists(task_output), identical(Sys.getenv("NOT_CRAN"), "true"))
task_archive_manifest <- jsonlite::read_json("archive-test-manifest.json", simplifyVector = TRUE)
task_before <- tools::md5sum(task_archive_manifest$file)
task_result <- testthat::test_check("multiScaleR", reporter = "summary",
  filter = "^(coverage-expansion|internals-and-helpers|joint-nb2-pool|optimization-workflow)$",
  stop_on_failure = FALSE)
task_table <- as.data.frame(task_result)
task_table$file <- basename(task_table$file)
task_pool <- "reusable workers preserve all starts, predictions and uncertainty"
task_parallel <- c("PSOCK workers load the same multiScaleR namespace as the main session",
  "screened starts run on real PSOCK workers",
  "PSOCK optimization works for unqualified MASS model calls",
  "PSOCK optimization works for nested clogit wrappers")
task_selected <- task_table$file %in% c("test-coverage-expansion.R", "test-internals-and-helpers.R") |
  (task_table$file == "test-joint-nb2-pool.R" & task_table$test == task_pool) |
  (task_table$file == "test-optimization-workflow.R" & task_table$test %in% task_parallel)
task_missing <- task_table[task_selected, , drop = FALSE]
task_stable <- identical(task_before, tools::md5sum(task_archive_manifest$file))
task_passed <- !any(task_table$failed > 0 | task_table$error) &&
  !any(task_missing$skipped > 0) && all(c(task_pool, task_parallel) %in% task_missing$test) &&
  all(c("test-coverage-expansion.R", "test-internals-and-helpers.R") %in% task_missing$file) && task_stable
task_record <- list(passed = task_passed, installed_package = task_expected_package,
  archive_tests_unchanged = task_stable,
  results = task_table[, c("file", "test", "passed", "failed", "error", "warning", "skipped")],
  previously_skipped = task_missing[, c("file", "test", "passed", "failed", "error", "warning", "skipped")],
  additional_unique_assertions = sum(task_missing$passed),
  supplemental_assertions = sum(task_table$passed),
  failures = sum(task_table$failed), errors = sum(task_table$error),
  warnings = sum(task_table$warning), skips = sum(task_table$skipped))
writeLines(jsonlite::toJSON(task_record, auto_unbox = TRUE, null = "null"), paste0(task_output, ".tmp"))
stopifnot(file.rename(paste0(task_output, ".tmp"), task_output), task_passed)

# Compact test metadata only. Keep the existing count ledger and assertions
# unchanged, and never serialize test results that might contain fitted models.
msr_test_timings <- function(d, elapsed, setup = NA_real_) {
  required <- c("file", "test", "passed", "failed", "error", "warning", "skipped")
  if (!is.data.frame(d) || !all(required %in% names(d)))
    stop("Unsupported testthat result columns; retain the count ledger.")
  if (length(elapsed) != 1L || !is.finite(elapsed) || elapsed < 0)
    stop("Test-run elapsed time must be one nonnegative finite value.")
  times <- if (is.numeric(d$real)) as.numeric(d$real) else rep(NA_real_, nrow(d))
  times[!is.finite(times) | times < 0] <- NA_real_
  blocks <- lapply(seq_len(nrow(d)), function(i) list(
    index = i, file = basename(d$file[i]), test = as.character(d$test[i]),
    elapsed_seconds = times[i], passed = d$passed[i], failed = d$failed[i],
    errors = d$error[i], warnings = d$warning[i], skips = d$skipped[i]))
  files <- lapply(split(seq_len(nrow(d)), basename(d$file)), function(i) list(
    blocks = length(i), timing_available = all(!is.na(times[i])),
    timed_block_seconds = if (all(!is.na(times[i]))) sum(times[i]) else NA_real_,
    passed = sum(d$passed[i]), failed = sum(d$failed[i]), errors = sum(d$error[i]),
    warnings = sum(d$warning[i]), skips = sum(d$skipped[i])))
  complete <- all(!is.na(times))
  total <- if (complete) sum(times) else NA_real_
  list(schema_version = 1L, testthat_version = as.character(utils::packageVersion("testthat")),
    R_version = R.version.string, timing_available = complete,
    test_run_elapsed_seconds = elapsed, setup_elapsed_seconds = setup,
    timed_block_seconds = total,
    unassigned_elapsed_seconds = if (complete && total <= elapsed) elapsed - total else NA_real_,
    overlapping_or_incomplete_block_times = !complete || total > elapsed,
    scope = "Descriptive block elapsed times; file/helper setup and child CPU are not separately measured.",
    controlled_speed_comparison_eligible = FALSE, by_file = files, blocks = blocks)
}

msr_write_test_timings <- function(d, result_path, elapsed, setup = NA_real_) {
  path <- paste0(result_path, ".timings.json")
  temporary <- paste0(path, ".tmp")
  if (file.exists(path) || file.exists(temporary)) stop("Retain the existing timing record.")
  value <- msr_test_timings(d, elapsed, setup)
  jsonlite::write_json(value, temporary, auto_unbox = TRUE, pretty = TRUE,
    digits = 16, na = "null", null = "null")
  if (!file.rename(temporary, path)) stop("Could not publish the atomic timing record.")
  invisible(value)
}

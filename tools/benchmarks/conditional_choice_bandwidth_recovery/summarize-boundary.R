args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else "tools/benchmarks/conditional_choice_bandwidth_recovery/results"
primary <- read.csv(file.path(root, "flagged-combined.csv"), stringsAsFactors = FALSE)
if (!nrow(primary)) quit(status = 0L)
refit <- read.csv(file.path(root, "boundary-starts", "scenarios.csv"), stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(refit$id), setequal(refit$id, primary$id))
refit <- refit[match(primary$id, refit$id), ]
stopifnot(identical(refit$id, primary$id))
out <- cbind(primary[c("id", "kernel", "n_analysis", "signal", "replicate")],
  primary_status = primary$fitted_status, primary_bandwidth = primary$fitted_bandwidth,
  primary_objective_gap = primary$objective_gap,
  diagnostic_status = refit$fitted_status, diagnostic_bandwidth = refit$fitted_bandwidth,
  diagnostic_objective_gap = refit$fitted_nll - refit$independent_nll,
  diagnostic_eligible_starts = refit$eligible_starts,
  primary_test_nll_per_set = primary$fitted_assessment_nll/160,
  diagnostic_test_nll_per_set = refit$fitted_assessment_nll/160,
  reference_bandwidth = refit$independent_bandwidth)
stopifnot(isTRUE(all.equal(refit$independent_nll, primary$independent_nll, tolerance = 1e-12)))
write.csv(out, file.path(root, "boundary-start-comparison.csv"), row.names = FALSE)
print(out, row.names = FALSE)

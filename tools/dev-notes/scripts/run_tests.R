# Run selected test files against the worktree source and print one summary row per file.
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(devtools))
load_all(".", quiet = TRUE)
res <- lapply(args, function(f) {
  r <- as.data.frame(testthat::test_file(file.path("tests/testthat", f), reporter = "silent"))
  bad <- r[r$failed > 0 | r$error, c("test", "failed", "error")]
  if (nrow(bad)) { cat("PROBLEMS in", f, "\n"); print(bad) }
  data.frame(file = f, tests = nrow(r), expectations = sum(r$nb), failed = sum(r$failed),
             errors = sum(r$error), skipped = sum(r$skipped), warnings = sum(r$warning))
})
cat("SUMMARY\n"); print(do.call(rbind, res), row.names = FALSE)

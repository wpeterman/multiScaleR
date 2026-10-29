study <- "tools/benchmarks/tmb_marginal_likelihood/spatial_reference_20261028"
preceding <- "tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028"
previous_study <- "tools/benchmarks/tmb_marginal_likelihood/latent_gaussian_20261028"
previous_study <- normalizePath(previous_study, mustWork = TRUE)
mode <- commandArgs(trailingOnly = TRUE)[1L]
stopifnot(grepl("^unit[0-9]+$", mode), !file.exists(file.path(study, paste0(mode, ".json"))))
qualified <- jsonlite::read_json(file.path(preceding, "check01.json"), simplifyVector = TRUE)
stopifnot(isTRUE(qualified$passed))
snapshot <- jsonlite::read_json(file.path(preceding, "snapshot01-input.json"), simplifyVector = TRUE)
libraries <- jsonlite::read_json(file.path(preceding, "compatibility01.json"), simplifyVector = TRUE)$libraries
.libPaths(c(file.path(qualified$workspace, "multiScaleR.Rcheck"), libraries, .libPaths()))
library(multiScaleR)
source(file.path(snapshot$package, "tests/testthat/helper-bayesian-scales.R"))
source(file.path(study, "adapter.R"))
critical <- file.path(study, c("adapter.R", "test-adapter.R", "PROTOCOL.md", "run-fit.R"))
hash_before <- tools::md5sum(critical)
stopifnot(length(parse(file.path(study, "run-fit.R"))) > 0L)
result <- testthat::test_file(file.path(study, "test-adapter.R"), reporter = "summary",
  env = new.env(parent = globalenv()))
tab <- as.data.frame(result)
stable <- identical(hash_before, tools::md5sum(critical))
record <- list(passed = all(tab$failed == 0 & !tab$error) && stable,
  assertions = sum(tab$passed), failures = sum(tab$failed), errors = sum(tab$error),
  warnings = sum(tab$warning), skips = sum(tab$skipped), source_stable = stable,
  source = as.list(hash_before), package_path = find.package("multiScaleR"),
  version = as.character(utils::packageVersion("multiScaleR")), R = R.version.string,
  qualification = "Deterministic adapter checks only; no sampler fitted.")
sr_atomic(record, file.path(study, paste0(mode, ".json")))
stopifnot(record$passed)

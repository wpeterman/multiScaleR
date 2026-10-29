.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/spatial_source_contract_20261029"
mode <- commandArgs(trailingOnly = TRUE)[1L]
path <- file.path(study, paste0(mode, ".json"))
stopifnot(!file.exists(path))
atomic <- function(value, path) {
  temporary <- paste0(path, ".tmp")
  writeLines(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null"), temporary)
  stopifnot(file.rename(temporary, path))
}
critical <- c(list.files("R", full.names = TRUE), list.files("man", full.names = TRUE),
  list.files("tests/testthat", "[.]R$", full.names = TRUE),
  "NAMESPACE", "DESCRIPTION", "vignettes/bayesian-computation.Rmd")
source_hash_before <- tools::md5sum(critical)
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
for (helper in list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE)) source(helper)
files <- if (grepl("^unit", mode)) "test-spatial-source-report.R" else
  if (grepl("^compatibility", mode)) c("test-spatial-scale-prototype.R",
    "test-spatial-scale-integration.R", "test-scale-model-preparation.R",
    "test-kernel-prep.R", "test-year-specific-example.R", "test-prediction-contract.R") else stop("Unknown mode")
results <- lapply(files, function(file) {
  result <- testthat::test_file(file.path("tests/testthat", file), reporter = "summary",
    env = new.env(parent = globalenv()))
  tab <- as.data.frame(result)
  record <- list(file = file, assertions = sum(tab$passed), failures = sum(tab$failed),
    errors = sum(tab$error), warnings = sum(tab$warning), skips = sum(tab$skipped))
  atomic(record, file.path(study, paste0(mode, "-", sub("[.]R$", ".json", file))))
  record
})
passed <- all(vapply(results, function(x) x$failures == 0 && x$errors == 0, TRUE))
stable <- identical(source_hash_before, tools::md5sum(critical))
atomic(list(passed = passed && stable, results = results, R = R.version.string,
  source_stable = stable, source = as.list(source_hash_before), libraries = .libPaths()), path)
stopifnot(passed, stable)



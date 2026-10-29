.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/annual_sd_profile_20261029"
stopifnot(!file.exists(file.path(study, "documentation01.json")))
before <- tools::md5sum("NAMESPACE")
roxygen2::roxygenise(".", load_code = function(path) {
  pkgload::load_all(path, helpers = FALSE, compile = FALSE)
  asNamespace("multiScaleR")
})
stopifnot(identical(before, tools::md5sum("NAMESPACE")))
help <- "man/dot-msr_factorized_scale_sd_profile.Rd"
stopifnot(!length(capture.output(tools::checkRd(help))))
record <- list(passed = TRUE, namespace_unchanged = TRUE,
  help = as.list(tools::md5sum(help)),
  vignette = "Unchanged; full candidate check will rebuild all vignettes and examples.")
temporary <- file.path(study, "documentation01.json.tmp")
writeLines(jsonlite::toJSON(record, auto_unbox = TRUE), temporary)
stopifnot(file.rename(temporary, file.path(study, "documentation01.json")))

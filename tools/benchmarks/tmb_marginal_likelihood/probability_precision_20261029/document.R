.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/probability_precision_20261029"
stopifnot(!file.exists(file.path(study, "documentation01.json")))
before <- tools::md5sum("NAMESPACE")
roxygen2::roxygenise(".", load_code = function(path) {
  pkgload::load_all(path, helpers = FALSE, compile = FALSE)
  asNamespace("multiScaleR")
})
stopifnot(identical(before, tools::md5sum("NAMESPACE")))
help <- "man/diagnose_bayes_probabilities.Rd"
stopifnot(!length(capture.output(tools::checkRd(help))))
output <- file.path(study, "rendered")
dir.create(output, showWarnings = FALSE)
rendered <- rmarkdown::render("vignettes/bayesian-computation.Rmd", output_dir = normalizePath(output),
  intermediates_dir = tempfile("msr-probability-availability-"), quiet = TRUE)
record <- list(passed = TRUE, namespace_unchanged = TRUE, help = as.list(tools::md5sum(help)),
  vignette = rendered, rendered_bytes = file.info(rendered)$size,
  visual_preview = "Not claimed: existing local-file browser policy rejection is preserved; no alternate bypass attempted.")
temporary <- file.path(study, "documentation01.json.tmp")
writeLines(jsonlite::toJSON(record, auto_unbox = TRUE), temporary)
stopifnot(file.rename(temporary, file.path(study, "documentation01.json")))

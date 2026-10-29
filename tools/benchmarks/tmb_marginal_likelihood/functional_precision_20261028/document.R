.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028"
stopifnot(!file.exists(file.path(study, "documentation01.json")))
namespace_before <- readLines("NAMESPACE")
roxygen2::roxygenise(".", load_code = function(path) {
  pkgload::load_all(path, helpers = FALSE, compile = FALSE)
  asNamespace("multiScaleR")
})
namespace_after <- readLines("NAMESPACE")
stopifnot(identical(setdiff(namespace_after, namespace_before), "export(diagnose_bayes_probabilities)"),
  !length(setdiff(namespace_before, namespace_after)))
help <- "man/diagnose_bayes_probabilities.Rd"
stopifnot(file.exists(help), !length(capture.output(tools::checkRd(help))))
output <- file.path(study, "rendered")
dir.create(output, showWarnings = FALSE)
rendered <- rmarkdown::render("vignettes/bayesian-computation.Rmd", output_dir = normalizePath(output),
  intermediates_dir = tempfile("msr-probability-vignette-"), quiet = TRUE)
stopifnot(file.exists(rendered))
result <- list(passed = TRUE, namespace_only_added = "diagnose_bayes_probabilities",
  help = as.list(tools::md5sum(help)), vignette = rendered, rendered_bytes = file.info(rendered)$size,
  visual_preview = "Pending: local-file browser preview was policy-blocked in preceding batch; no bypass attempted.")
writeLines(jsonlite::toJSON(result, auto_unbox = TRUE), file.path(study, "documentation01.json.tmp"))
stopifnot(file.rename(file.path(study, "documentation01.json.tmp"), file.path(study, "documentation01.json")))

.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/prediction_contracts_20261028"
stopifnot(!file.exists(file.path(study, "documentation01.json")))
before <- tools::md5sum("NAMESPACE")
roxygen2::roxygenise(".", load_code = function(path) {
  pkgload::load_all(path, helpers = FALSE, compile = FALSE)
  asNamespace("multiScaleR")
})
stopifnot(identical(before, tools::md5sum("NAMESPACE")))
help <- file.path("man", c("prediction_contract.Rd", "prepare_prediction.Rd",
  "predict.multiScaleR_bayes.Rd", "predict_scales.Rd", "simulate.multiScaleR_bayes.Rd"))
stopifnot(all(file.exists(help)))
for (path in help) {
  problem <- capture.output(tools::checkRd(path))
  stopifnot(!length(problem))
}
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
output <- file.path(study, "rendered")
dir.create(output, showWarnings = FALSE)
rendered <- rmarkdown::render("vignettes/bayesian-prediction.Rmd", output_dir = normalizePath(output),
  intermediates_dir = tempfile("msr-contract-vignette-"), quiet = TRUE)
stopifnot(file.exists(rendered))
result <- list(passed = TRUE, namespace_unchanged = TRUE, help = as.list(tools::md5sum(help)),
               vignette = rendered, rendered_bytes = file.info(rendered)$size)
writeLines(jsonlite::toJSON(result, auto_unbox = TRUE), file.path(study, "documentation01.json.tmp"))
stopifnot(file.rename(file.path(study, "documentation01.json.tmp"), file.path(study, "documentation01.json")))

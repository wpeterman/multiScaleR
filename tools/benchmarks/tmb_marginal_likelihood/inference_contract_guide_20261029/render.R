.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
study <- "tools/benchmarks/tmb_marginal_likelihood/inference_contract_guide_20261029"
mode <- commandArgs(trailingOnly = TRUE)[1L]
stopifnot(grepl("^render[0-9]+$", mode), !file.exists(file.path(study, paste0(mode, ".json"))))
critical <- c(list.files("R", full.names = TRUE), list.files("src", full.names = TRUE),
  list.files("man", full.names = TRUE), list.files("tests/testthat", recursive = TRUE, full.names = TRUE),
  list.files("inst/extdata", full.names = TRUE), list.files("vignettes", full.names = TRUE),
  "NAMESPACE", "DESCRIPTION", "NEWS")
critical <- critical[file.info(critical)$isdir %in% FALSE]
source_before <- tools::md5sum(critical)
pkgload::load_all(".", quiet = TRUE, compile = FALSE)
help_topics <- c("multiScale_bayes", "multiScale_optim", "prediction_contract", "prepare_prediction")
help_valid <- vapply(help_topics, function(x) {
  path <- file.path("man", paste0(x, ".Rd"))
  file.exists(path) && any(grepl(paste0("\\alias{", x, "}"), readLines(path), fixed = TRUE)) &&
    !length(capture.output(tools::checkRd(path)))
}, TRUE)
vignette_links <- c("quickstart", "bayesian-scales", "bayesian-computation", "bayesian-prediction", "year-specific-rasters")
stopifnot(all(help_valid), all(file.exists(file.path("vignettes", paste0(vignette_links, ".Rmd")))))
output <- file.path(study, paste0(mode, "-rendered"))
dir.create(output, showWarnings = FALSE)
environment <- new.env(parent = globalenv())
rendered <- rmarkdown::render("vignettes/inference-contracts.Rmd", output_dir = normalizePath(output),
  intermediates_dir = tempfile("msr-inference-guide-render-"), envir = environment, quiet = TRUE)
stopifnot(file.info(rendered)$size > 10000,
  identical(environment$info$target, "posterior"),
  identical(environment$info$population_priors, "proper"),
  identical(environment$info$uncertainty$population_parameters, "not_computed"),
  identical(environment$contract$quantity, "expected_response"),
  identical(environment$contract$uncertainty$observation_noise, FALSE),
  identical(environment$contract$metadata_establishes_eligibility, FALSE),
  nrow(environment$expected) == 3L)
stable <- identical(source_before, tools::md5sum(critical))
stopifnot(stable)
record <- list(passed = TRUE, source_stable = stable, help_topics = as.list(help_valid),
  vignette_sources_exist = vignette_links, libraries = .libPaths(),
  rendered = rendered, rendered_bytes = file.info(rendered)$size,
  prepared_target = environment$info[c("target", "method", "population_priors", "uncertainty")],
  prediction_contract = environment$contract, prediction_rows = nrow(environment$expected),
  no_new_fits_chains_data_or_replicates = TRUE,
  visual_preview = "Not claimed: prior browser policy rejection remains; no alternate bypass.")
path <- file.path(study, paste0(mode, ".json"))
writeLines(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null"), paste0(path, ".tmp"))
stopifnot(file.rename(paste0(path, ".tmp"), path))

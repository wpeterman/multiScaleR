# Supplemental exact-target sampler comparison, with no new NUTS fits.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[1] %in% c("fit", "assess"))
  stop("Usage: Rscript run_slice.R fit|assess")
root <- "tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112"
directory <- file.path(root, "slice_agreement")
source(file.path(root, "simulate.R"))
cases <- data.frame(case = c("gaussian-annual-01", "gaussian-annual-02",
  "nbinom2-repeated-01", "nbinom2-repeated-02"),
  family = c("gaussian", "gaussian", "nbinom2", "nbinom2"),
  design = c("annual", "annual", "repeated", "repeated"),
  seed = c(910020L, 910040L, 910180L, 910200L))
files <- c(file.path(root, c("simulate.R", "run_slice.R")), "R/bayes_scale_slice.R",
  "R/bayes_scale_likelihood.R", "R/bayes_scale_diagnostics.R", "R/scale_response_likelihood.R",
  "R/scale_model_prepare.R")
fingerprint <- vapply(files, function(x) digest::digest(file = x, algo = "sha256"), character(1))
fit_one <- function(case) {
  input <- simulate_case(case)
  stage <- "fit"; warnings <- character()
  start <- proc.time()[["elapsed"]]
  ans <- tryCatch(withCallingHandlers({
    fit <- multiScale_bayes(input$formula, input$inputs, input$data[input$analysis, ],
      scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
      method = "slice", chains = 4L, iter = 3000L, warmup = 1500L, seed = case$seed + 9L)
    stage <- "assessment"
    derived <- diagnose_bayes_scales(fit)
    eligible <- isTRUE(fit$diagnostics$passed) && isTRUE(derived$passed)
    probs <- diagnose_bayes_probabilities(fit, list(habitat = c(150, 350)))$probabilities
    list(status = if (eligible) "eligible" else "rejected",
      raw = fit$diagnostics, derived = derived$parameters, derived_passed = derived$passed,
      probability = probs[, c("row", "threshold", "estimate", "mcse")])
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(status = "failed", stage = stage, error = conditionMessage(e)))
  c(as.list(case), list(hash = fingerprint, elapsed = proc.time()[["elapsed"]] - start,
    warnings = unique(warnings), session = sessionInfo()), ans)
}
read_case <- function(case) {
  file <- file.path(directory, paste0(case$case, ".slice.rds"))
  if (!file.exists(file)) return(c(as.list(case), list(status = "unrun")))
  r <- readRDS(file)
  stopifnot(identical(r$hash, fingerprint))
  r
}
if (args[1] == "fit") {
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    if (read_case(case)$status != "unrun") next
    message(case$case, " slice")
    result <- fit_one(case)
    tmp <- tempfile(tmpdir = directory)
    saveRDS(result, tmp)
    if (!file.rename(tmp, file.path(directory, paste0(case$case, ".slice.rds"))))
      stop("Cannot publish slice fit record.")
    message("  ", result$status, if (!is.null(result$error)) paste0(" (", result$error, ")"))
  }
} else {
  result <- lapply(seq_len(nrow(cases)), function(i) read_case(cases[i, ]))
  fits <- do.call(rbind, lapply(result, function(r) data.frame(case = r$case,
    family = r$family, design = r$design, status = r$status,
    elapsed = if (is.null(r$elapsed)) NA_real_ else r$elapsed,
    error = if (is.null(r$error)) NA_character_ else r$error)))
  write.csv(fits, file.path(directory, "fits.csv"), row.names = FALSE)
  pairs <- do.call(rbind, lapply(result, function(r) {
    ref <- readRDS(file.path(root, "results", paste0(r$case, ".nuts.rds")))
    if (is.null(r$probability) || is.null(ref$probability)) return(NULL)
    p <- merge(r$probability, ref$probability, by = c("row", "threshold"), suffixes = c("_slice", "_nuts"))
    p$case <- r$case
    p$status_slice <- r$status
    p$status_nuts <- ref$status
    p$absolute_difference <- abs(p$estimate_slice - p$estimate_nuts)
    p$combined_mcse <- sqrt(p$mcse_slice^2 + p$mcse_nuts^2)
    p$material <- p$absolute_difference > pmax(.05, 5 * p$combined_mcse)
    p
  }))
  if (!is.null(pairs)) write.csv(pairs, file.path(directory, "paired_probabilities.csv"), row.names = FALSE)
  print(fits)
}

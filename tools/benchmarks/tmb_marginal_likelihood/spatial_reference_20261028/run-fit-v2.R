study <- normalizePath("tools/benchmarks/tmb_marginal_likelihood/spatial_reference_20261028")
preceding <- "tools/benchmarks/tmb_marginal_likelihood/functional_precision_20261028"
run_id <- commandArgs(trailingOnly = TRUE)[1L]
attempt_id <- paste0(run_id, "-v2")
plan <- jsonlite::read_json(file.path(study, "RUNLIST.json"), simplifyVector = TRUE)
run <- plan$runs[plan$runs$id == run_id, , drop = FALSE]
stopifnot(nrow(run) == 1L, !file.exists(file.path(study, paste0(attempt_id, ".json"))),
  !file.exists(file.path(study, paste0(attempt_id, "-started.json"))))
qualified <- jsonlite::read_json(file.path(preceding, "check01.json"), simplifyVector = TRUE)
snapshot <- jsonlite::read_json(file.path(preceding, "snapshot01-input.json"), simplifyVector = TRUE)
libraries <- jsonlite::read_json(file.path(preceding, "compatibility01.json"), simplifyVector = TRUE)$libraries
.libPaths(c(file.path(qualified$workspace, "multiScaleR.Rcheck"), libraries, .libPaths()))
library(multiScaleR)
source(file.path(snapshot$package, "tests/testthat/helper-bayesian-scales.R"))
source(file.path(study, "adapter.R"))
stopifnot(isTRUE(qualified$passed), isTRUE(jsonlite::read_json(file.path(study, "unit04.json"))$passed))
stopifnot(isTRUE(jsonlite::read_json(file.path(study, "repair-unit01.json"))$passed))
producer_files <- file.path(study, c("adapter.R", "run-fit-v2.R", "PROTOCOL.md", "RUNLIST.json"))
producer_before <- vapply(producer_files, function(path) digest::digest(file = path, algo = "sha256"), "")
stopifnot(identical(unname(producer_before[3L]), tolower(plan$protocol_sha256)))
source_files <- names(snapshot$source_sha256)
source_before <- vapply(source_files, function(path) digest::digest(file = path, algo = "sha256"), "")
stopifnot(identical(unname(source_before), unname(unlist(snapshot$source_sha256))))
record <- list(run = run_id, attempt = attempt_id, driver_revision = 2L,
  case = run[["case"]], sampler_seed = run$seed,
  settings = plan$settings, protocol_sha256 = plan$protocol_sha256,
  source_package = find.package("multiScaleR"), archive_SHA256 = qualified$archive_SHA256,
  producer_sha256 = as.list(producer_before), source_matches_311_file_gate = TRUE,
  started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), completed = FALSE,
  reference_qualified = FALSE, scope = "Selected-case computation; no calibration or public spatial inference.")
warnings <- character(); warning_count <- 0L
started <- proc.time()[[3L]]
result <- tryCatch(withCallingHandlers({
  case <- sr_case(run[["case"]])
  stopifnot(identical(case$input_hash, run$input_hash))
  record$input_hash <- case$input_hash
  sr_atomic(record, file.path(study, paste0(attempt_id, "-started.json")))
  objective <- multiScaleR:::.msr_spatial_scale_adfun(case$object)
  # RTMB updates last.par.best during fn/gr evaluation. Keep its block names so
  # tmbstan can map the unchanged numeric initialization into its Stan vector.
  named_initial <- lapply(seq_len(4L), function(i) objective$par + .05 * sin(seq_along(objective$par) + i))
  initial <- lapply(named_initial, unname)
  stopifnot(all(vapply(named_initial, function(p) is.finite(objective$fn(p)) && all(is.finite(objective$gr(p))), TRUE)),
    identical(names(objective$env$last.par.best), names(objective$par)))
  clock <- proc.time()[[3L]]
  sampler <- tmbstan::tmbstan(obj = objective, chains = 4L, iter = 16000L, warmup = 4000L,
    control = list(adapt_delta = .999, max_treedepth = 14L, metric = "diag_e"),
    cores = 1L, refresh = 0L, init = initial, seed = as.integer(run$seed), laplace = FALSE, silent = TRUE)
  record$sampling_seconds <- proc.time()[[3L]] - clock
  record$sampling_completed <- TRUE
  record$raw_diagnostics <- multiScaleR:::.msr_bayes_diagnostics(sampler, list(max_treedepth = 14L), 4L)
  sr_atomic(record, file.path(study, paste0(attempt_id, "-raw.json")))
  raw <- rstan::extract(sampler, permuted = FALSE, inc_warmup = FALSE)
  stopifnot(identical(dim(raw)[1:2], c(12000L, 4L)))
  derived <- sr_derived_array(raw, case)
  table <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(derived),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail))
  derived_pass <- all(is.finite(table$rhat) & table$rhat <= 1.01 &
    is.finite(table$ess_bulk) & table$ess_bulk >= 400 & is.finite(table$ess_tail) & table$ess_tail >= 400)
  old <- jsonlite::read_json(file.path("tools/benchmarks/tmb_marginal_likelihood/latent_gaussian_20261028",
    paste0("pilot-", run[["case"]], ".json")), simplifyVector = TRUE)
  fixed <- sr_functional(derived, old$functional[c("variable", "quantile", "threshold")])
  quantiles <- sr_quantiles(derived)
  own <- sr_functional(derived, quantiles)
  fixed$historical_laplace_cdf <- old$functional$laplace_cdf
  fixed$historical_laplace_mcse <- old$functional$laplace_mcse
  fixed$combined_mcse <- sqrt(fixed$mcse^2 + fixed$historical_laplace_mcse^2)
  fixed$historical_difference_flag <- ifelse(is.finite(fixed$combined_mcse) & fixed$combined_mcse > 0,
    abs(fixed$estimate - fixed$historical_laplace_cdf) > 3 * fixed$combined_mcse, NA)
  summaries <- old$comparison
  summaries$nuts_lower <- quantiles$threshold[quantiles$quantile == .025]
  summaries$nuts_median <- quantiles$threshold[quantiles$quantile == .5]
  summaries$nuts_upper <- quantiles$threshold[quantiles$quantile == .975]
  sampler_parameters <- rstan::get_sampler_params(sampler, inc_warmup = FALSE)
  compact <- function(x) if (length(x) && all(is.finite(x)))
    as.list(setNames(unname(quantile(x, c(0, .025, .5, .975, 1))), c("min", "lower", "median", "upper", "max"))) else
      list(unavailable = TRUE, count = length(x))
  unpacked <- sr_unpack(raw, case$object$parameters)
  latent_norm <- matrix(sqrt(rowSums(unpacked$scale_z^2)), 12000L, 4L)
  geometry <- lapply(seq_len(4L), function(i) {
    sp <- sampler_parameters[[i]]
    stopifnot(nrow(sp) == 12000L)
    divergent <- sp[, "divergent__"] > 0
    distributions <- lapply(list(all = rep(TRUE, 12000L), divergent = divergent, nondivergent = !divergent),
      function(index) list(count = sum(index), year_sd = compact(derived[index, i, "year_sd"]),
        field_sd = compact(derived[index, i, "field_sd"]), scale_latent_norm = compact(latent_norm[index, i])))
    list(chain = i, divergences = sum(divergent), depth_hits = sum(sp[, "treedepth__"] >= 14),
      E_BFMI = mean(diff(sp[, "energy__"])^2) / var(sp[, "energy__"]),
      energy = compact(sp[, "energy__"]), accept_stat = compact(sp[, "accept_stat__"]),
      stepsize = compact(sp[, "stepsize__"]), treedepth = compact(sp[, "treedepth__"]),
      leapfrog = compact(sp[, "n_leapfrog__"]), distributions = distributions)
  })
  record$derived_diagnostics <- table; record$derived_pass <- derived_pass
  record$fixed_functional <- fixed; record$fixed_functional_pass <- all(fixed$precision_pass)
  record$quantiles <- quantiles; record$own_functional <- own
  record$own_functional_pass <- all(own$precision_pass)
  record$comparison <- summaries; record$geometry <- geometry
  record$reference_qualified <- isTRUE(record$raw_diagnostics$passed) && derived_pass &&
    record$fixed_functional_pass && record$own_functional_pass
  record$historical_comparison_note <- paste("Laplace CDFs are retained historical outputs; old thresholds were serialized to 12 digits.",
    "They were not recomputed at the parsed rounded values. New independent runs use identical stored thresholds.")
  record$completed <- TRUE
  record
}, warning = function(w) {
  warning_count <<- warning_count + 1L
  warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
}), error = function(e) { record$error <- conditionMessage(e); record$completed <- FALSE; record })
result$warnings <- unique(warnings); result$warning_count <- warning_count
result$elapsed_seconds <- proc.time()[[3L]] - started
result$source_unchanged <- identical(source_before,
  vapply(source_files, function(path) digest::digest(file = path, algo = "sha256"), ""))
result$producer_unchanged <- identical(producer_before,
  vapply(producer_files, function(path) digest::digest(file = path, algo = "sha256"), ""))
result$reference_qualified <- isTRUE(result$reference_qualified) && result$source_unchanged && result$producer_unchanged
result$completed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
sr_atomic(result, file.path(study, paste0(attempt_id, ".json")))
cat(run_id, "completed:", result$completed, "reference qualified:", result$reference_qualified, "\n")
if (!isTRUE(result$completed)) quit(status = 1L)


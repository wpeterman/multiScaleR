# Coverage check of joint_mixed profile-likelihood intervals for kernel sigma.
#
# Run from the source root:
#   Rscript tools/benchmarks/joint_mixed_wald/profile_run.R <output_dir> [workers]
# This reuses the Wald study's design (design.R): the same 900 fixed-seed
# datasets, the same default fit, and the same search range. For each dataset
# it refits the model, calls summary(fit, profile = TRUE), and records the
# profile limits and their labels. Each dataset is saved as
# <output_dir>/profile_records/<id>.rds, so an interrupted run resumes where
# it stopped. Workers load this source tree with pkgload::load_all(). The
# driver writes its own PID and the worker PIDs to <output_dir>/profile_pids.txt.
# profile_summarize.R joins these records with results/metrics.csv from the
# Wald study, which holds the likelihood-ratio decision at the true sigma.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply an output directory.")
output <- normalizePath(args[1L], mustWork = FALSE)
workers <- if (length(args) > 1L) as.integer(args[2L]) else 8L
study <- "tools/benchmarks/joint_mixed_wald"
source(file.path(study, "design.R"))
dir.create(file.path(output, "profile_records"), recursive = TRUE, showWarnings = FALSE)
plan <- jm_plan(150L)
todo <- plan[!file.exists(file.path(output, "profile_records", paste0(plan$id, ".rds"))), ]
cat("Driver PID", Sys.getpid(), "-", nrow(todo), "of", nrow(plan), "datasets to run\n")

# One dataset: the study's default fit, then profile limits from summary().
jm_profile_case <- function(case) {
  sample <- jm_case(case)
  started <- proc.time()[["elapsed"]]
  start <- jm_capture(jm_native(sample$data, case$family))
  fit <- if (inherits(start$value, "error")) start else jm_capture(multiScale_optim(
    start$value, sample$inputs, backend = "joint_mixed", n_cores = 1, verbose = FALSE))
  fit_seconds <- proc.time()[["elapsed"]] - started
  ok <- !inherits(fit$value, "error")
  out <- data.frame(case[c("id", "family", "sigma", "replicate", "seed")],
    true_sigma = case$sigma * 1000, fit_ok = ok, fit_seconds = fit_seconds,
    stringsAsFactors = FALSE)
  if (!ok) return(out)
  m <- fit$value
  started <- proc.time()[["elapsed"]]
  prof <- jm_capture(summary(m, profile = TRUE)$opt_scale)
  profile_seconds <- proc.time()[["elapsed"]] - started
  profile_ok <- !inherits(prof$value, "error")
  tab <- if (profile_ok) prof$value else NULL
  bounded <- if (profile_ok) attr(tab, "bound_limited") else NULL
  lower <- if (profile_ok) tab[1L, "2.5%"] else NA_real_
  upper <- if (profile_ok) tab[1L, "97.5%"] else NA_real_
  truth <- case$sigma * 1000
  cbind(out, data.frame(estimate = m$scale_est$Mean, nll_hat = m$optim_results$value,
    profile_ok = profile_ok,
    profile_error = if (profile_ok) "" else conditionMessage(prof$value),
    profile_warnings = paste(prof$warnings, collapse = " | "),
    method = if (profile_ok) unname(attr(tab, "interval_method")) else NA_character_,
    bound_lower = if (profile_ok) bounded[1L, 1L] else NA,
    bound_upper = if (profile_ok) bounded[1L, 2L] else NA,
    profile_lower = lower, profile_upper = upper,
    profile_cover95 = if (profile_ok) truth >= lower & truth <= upper else NA,
    profile_above = if (profile_ok) truth > upper else NA,
    profile_below = if (profile_ok) truth < lower else NA,
    profile_seconds = profile_seconds, stringsAsFactors = FALSE))
}

root <- normalizePath(".")
cl <- parallel::makePSOCKcluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
pids <- unlist(parallel::clusterCall(cl, Sys.getpid))
writeLines(c(paste("driver", Sys.getpid()), paste("worker", pids)),
           file.path(output, "profile_pids.txt"))
parallel::clusterCall(cl, function(root, study) {
  Sys.setenv(OMP_NUM_THREADS = "1")
  setwd(root)
  suppressMessages(pkgload::load_all(root, compile = FALSE, quiet = TRUE))
  source(file.path(study, "design.R"))
  invisible(NULL)
}, root, study)
parallel::clusterExport(cl, "jm_profile_case")
started <- Sys.time()
done <- parallel::parLapplyLB(cl, split(todo, seq_len(nrow(todo))), function(case, output) {
  destination <- file.path(output, "profile_records", paste0(case$id, ".rds"))
  record <- jm_profile_case(case)
  tmp <- tempfile(tmpdir = dirname(destination))
  saveRDS(record, tmp)
  file.rename(tmp, destination)
  case$id
}, output)
cat("Finished", length(done), "datasets in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")

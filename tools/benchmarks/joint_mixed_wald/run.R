# Coverage of joint_mixed log-scale Wald intervals for kernel sigma.
#
# Run from the source root:
#   Rscript tools/benchmarks/joint_mixed_wald/run.R <output_dir> [workers]
# Each dataset is saved as <output_dir>/records/<id>.rds, so an interrupted run
# resumes where it stopped. Every worker loads this source tree with
# pkgload::load_all(), never the installed package. The driver writes its own
# PID and the worker PIDs to <output_dir>/pids.txt.
#
# The design is fixed before any results are seen:
# * Families: Poisson and NB2 (size 2), log link, one random intercept.
# * n = 200 points in 40 sites (5 points per site), random-intercept SD 0.5.
# * One Gaussian kernel mean, landscape coefficient 0.8, a fixed covariate z
#   with coefficient 0.3, intercept 0.4, and a known offset log(effort).
# * True sigma 120, 240, or 450 map units; search range 50 to 800.
# * 150 fixed-seed replicates per family and sigma (900 datasets).
# * Default joint_mixed starts (six), as a user would run it.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply an output directory.")
output <- normalizePath(args[1L], mustWork = FALSE)
workers <- if (length(args) > 1L) as.integer(args[2L]) else 8L
study <- "tools/benchmarks/joint_mixed_wald"
source(file.path(study, "design.R"))
dir.create(file.path(output, "records"), recursive = TRUE, showWarnings = FALSE)
plan <- jm_plan(150L)
utils::write.csv(plan, file.path(output, "plan.csv"), row.names = FALSE)
todo <- plan[!file.exists(file.path(output, "records", paste0(plan$id, ".rds"))), ]
cat("Driver PID", Sys.getpid(), "-", nrow(todo), "of", nrow(plan), "datasets to run\n")

root <- normalizePath(".")
cl <- parallel::makePSOCKcluster(workers)
on.exit(parallel::stopCluster(cl), add = TRUE)
pids <- unlist(parallel::clusterCall(cl, Sys.getpid))
writeLines(c(paste("driver", Sys.getpid()), paste("worker", pids)),
           file.path(output, "pids.txt"))
parallel::clusterCall(cl, function(root, study) {
  Sys.setenv(OMP_NUM_THREADS = "1")
  setwd(root)
  suppressMessages(pkgload::load_all(root, compile = FALSE, quiet = TRUE))
  source(file.path(study, "design.R"))
  invisible(NULL)
}, root, study)
started <- Sys.time()
done <- parallel::parLapplyLB(cl, split(todo, seq_len(nrow(todo))), function(case, output) {
  destination <- file.path(output, "records", paste0(case$id, ".rds"))
  record <- jm_run_case(case)
  tmp <- tempfile(tmpdir = dirname(destination))
  saveRDS(record, tmp)
  file.rename(tmp, destination)
  case$id
}, output)
cat("Finished", length(done), "datasets in",
    format(round(difftime(Sys.time(), started, units = "mins"), 1)), "\n")

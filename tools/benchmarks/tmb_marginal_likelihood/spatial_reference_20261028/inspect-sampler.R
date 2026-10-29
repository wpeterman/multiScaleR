.libPaths(c("tools/benchmarks/joint_scale_assessment/phase31/library",
  "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library", .libPaths()))
cat("tmbstan", as.character(packageVersion("tmbstan")), "rstan", as.character(packageVersion("rstan")), "\n")
cat(paste(deparse(tmbstan::tmbstan), collapse = "\n"), "\n")

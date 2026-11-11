args <- commandArgs(trailingOnly = TRUE)
path <- args[1L]
setup <- readRDS(file.path(path, "setup.rds"))
.libPaths(setup$paths)
if (setup$source) {
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Install pkgload or install multiScaleR before running parallel fits from source.")
  pkgload::load_all(setup$path, compile = FALSE, attach = FALSE, helpers = FALSE, quiet = TRUE)
} else {
  loadNamespace("multiScaleR", lib.loc = dirname(setup$path))
}
ns <- asNamespace("multiScaleR")
loaded <- normalizePath(getNamespaceInfo(ns, "path"), winslash = "/", mustWork = TRUE)
if (!identical(loaded, setup$path)) stop("Worker loaded a different multiScaleR installation.")
get(".msr_bayes_worker_main", ns, inherits = FALSE)(path, file.path(path, args[2L]))

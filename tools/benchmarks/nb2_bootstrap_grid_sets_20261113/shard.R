# Continue one disjoint subset of the frozen grid plan.
# Rscript shard.R ORIGINAL_GRID SHARD_OUTPUT IDS_TXT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L)
original <- normalizePath(args[1L], mustWork = TRUE)
out <- normalizePath(args[2L], mustWork = FALSE)
wanted <- scan(args[3L], what = "", quiet = TRUE)
stopifnot(length(wanted) > 0L, !anyDuplicated(wanted))
config <- readRDS(file.path(original, "configuration.rds"))
stopifnot(identical(config$panel, "grid"), all(wanted %in% config$plan$id))
runner <- "tools/benchmarks/nb2_bootstrap_grid_sets_20261113/run.R"
stopifnot(identical(config$script, digest::digest(runner, file = TRUE, algo = "sha256")))
code <- parse(runner)
# The frozen runner takes command-line arguments in its first expression and
# invokes main() last. Load its exact definitions with the grid configuration.
stopifnot(identical(code[[1L]][[2L]], quote(args)),
  identical(code[[length(code)]], quote(main())))
env <- new.env(parent = globalenv())
env$args <- c("grid", out)
env$commandArgs <- function(trailingOnly = FALSE)
  if (trailingOnly) env$args else paste0("--file=", runner)
for (i in seq.int(2L, length(code) - 1L)) eval(code[[i]], env)
stopifnot(identical(env$plan, config$plan),
  identical(env$fixed_grid, config$grid),
  identical(env$coarse_nodes, config$coarse),
  identical(env$B, config$B), identical(env$level, config$level),
  identical(env$hash_files(), config$source))

main <- function() {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out, "cases")
  dir.create(path, showWarnings = FALSE)
  lock <- file.path(out, "run.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Shard output already in use.")
  stopifnot(identical(dirname(normalizePath(lock, winslash = "/")),
    normalizePath(out, winslash = "/")))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  pkgload::load_all(path = env$source_root, compile = FALSE, quiet = TRUE)
  stopifnot(identical(env$hash_files(), config$source))
  manifest <- list(wanted = wanted, source = config$source, script = config$script)
  meta <- file.path(out, "shard.rds")
  if (file.exists(meta)) stopifnot(identical(readRDS(meta), manifest)) else {
    env$atomic(manifest, meta)
    writeLines(wanted, file.path(out, "ids.txt"))
  }
  for (id in wanted) {
    dest <- file.path(path, paste0(id, ".rds"))
    if (file.exists(dest)) next
    if (file.exists(file.path(out, "STOP"))) break
    if (file.exists(file.path(original, "cases", basename(dest))))
      stop("Shard overlaps an original checkpoint at ", id)
    cat(format(Sys.time()), "grid shard starting", id, "\n")
    flush.console()
    z <- env$fit_case(env$plan[env$plan$id == id, , drop = FALSE])
    env$atomic(z, dest)
    cat(format(Sys.time()), "grid shard finished", id, z$status, "seconds", z$seconds,
      "grid", z$sets$fine99$available, "covered", z$sets$fine99$covers,
      "point", z$point$accept[1L], "\n")
    flush.console()
  }
  stopifnot(identical(env$hash_files(), config$source))
}
main()

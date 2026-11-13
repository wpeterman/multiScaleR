# Verify the archived protocol identities and all production keys independently.
study <- "tools/benchmarks/bayesian_scales/cv_refits_20261113"
source(file.path(study, "design.R"))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% c("frozen", "records"))
dir <- file.path(study, "records")
config <- lapply(c("gaussian", "nbinom2"), function(f)
  readRDS(file.path(dir, paste0("config-", f, ".rds"))))
stopifnot(identical(config[[1L]], config[[2L]]),
  identical(config[[1L]]$plan, plan))
if (args == "frozen") {
  expected <- config[[1L]]$code
  actual <- vapply(names(expected), function(path)
    digest::digest(readLines(path, warn = FALSE), algo = "sha256"), character(1L))
  stopifnot(identical(expected, actual))
}
all_paths <- unlist(lapply(plan$case, function(x)
  file.path(dir, paste(x, config[[1L]]$model, "rds", sep = "."))))
stopifnot(all(file.exists(all_paths)),
  setequal(list.files(dir, pattern = "[.]rds$", full.names = TRUE),
    c(all_paths, file.path(dir, c("config-gaussian.rds", "config-nbinom2.rds")))))
attempts <- read.csv(file.path(study, "attempts.csv"))
refs <- read.csv(file.path(study, "references.csv"))
old_refs <- read.csv(file.path(study, "references-frozen-oracle.csv"))
aug <- read.csv(file.path(study, "augmentation.csv"))
pits <- read.csv(file.path(study, "pits.csv"))
old_pits <- read.csv(file.path(study, "pits-frozen-oracle.csv"))
repair_config <- readRDS(file.path(study, "oracle-repair-configuration.rds"))
stopifnot(identical(repair_config$code, vapply(names(repair_config$code), function(path)
  digest::digest(readLines(path, warn = FALSE), algo = "sha256"), character(1L))))
n_fold <- length(config[[1L]]$folds)
stopifnot(nrow(attempts) == n_fold * length(all_paths),
  !any(attempts$status == "pending"),
  !anyDuplicated(attempts[c("case", "model", "fold")]),
  identical(refs[setdiff(names(refs), c("oracle_nll", "oracle_refinement"))],
    old_refs[setdiff(names(old_refs), c("oracle_nll", "oracle_refinement"))]),
  identical(pits[setdiff(names(pits), c("oracle_pit", "oracle_refinement"))],
    old_pits[setdiff(names(old_pits), c("oracle_pit", "oracle_refinement"))]),
  all(is.finite(refs$oracle_nll)), all(is.finite(pits$oracle_pit)))
close <- function(a, b) isTRUE(all.equal(a, b, tolerance = 1e-12, check.attributes = FALSE))
for (path in all_paths) {
  a <- readRDS(path)
  key <- which(attempts$case == a$case$case & attempts$model == a$model)
  q <- attempts[key, ]
  stopifnot(nrow(q) == n_fold, identical(unname(a$status[q$fold]), q$status),
    !anyNA(a$status))
  for (j in seq_len(n_fold)) {
    r <- a$result[[q$fold[j]]]
    keep <- function(d) d[d$case == a$case$case & d$model == a$model & d$fold == q$fold[j], ]
    x <- keep(refs)
    x_old <- keep(old_refs)
    y <- keep(aug)
    z <- keep(pits)
    z_old <- keep(old_pits)
    if (is.null(r$reference)) {
      stopifnot(!nrow(x), !nrow(y), !nrow(z))
      next
    }
    stopifnot(nrow(x) == 2L, nrow(y) == 2L * nrow(r$streams),
      nrow(z) == nrow(r$outcomes), identical(r$outcomes$visit, z$visit),
      close(r$outcomes$observed, z$observed), close(r$reference$pit, z$fit_pit),
      close(r$oracle$score$pit, z_old$oracle_pit), close(r$score, q$package_nll[j]),
      close(r$package_replay_diff, q$replay_error[j]))
    for (event in c("observed", "stress")) {
      k <- x[x$event == event, ]
      k_old <- x_old[x_old$event == event, ]
      s <- y[y$event == event, ]
      stopifnot(close(k$reference_nll, r$reference[[event]]),
        close(k_old$oracle_nll, r$oracle$score[[event]]),
        close(k$package_seed_nll, r$original[[event]]),
        close(k$refinement, r$reference_refine[[event]]),
        close(k_old$oracle_refinement, r$oracle$refine[[event]]),
        identical(s$score_seed, r$streams$seed), identical(s$budget, r$streams$budget),
        close(s$nll, r$streams[[event]]),
        close(s$delta, r$streams[[event]] - r$reference[[event]]),
        max(abs(s$delta - (s$nll - k$reference_nll))) < 1e-10)
    }
  }
}
old_joint <- is.finite(old_refs$oracle_nll)
old_pit <- is.finite(old_pits$oracle_pit)
stopifnot(max(abs(refs$oracle_nll[old_joint] - old_refs$oracle_nll[old_joint])) < 1e-8,
  max(abs(pits$oracle_pit[old_pit] - old_pits$oracle_pit[old_pit])) < 1e-8)
cat(length(all_paths), "fit-call records and", nrow(attempts), "fold attempts verified\n")
cat(nrow(refs), "paired joint references,", nrow(aug), "augmentation values and",
  nrow(pits), "marginal probabilities verified against records\n")

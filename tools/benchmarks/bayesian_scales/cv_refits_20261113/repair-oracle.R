# Reconstruct only the oracle with stable distance exponents. No refits.
study <- "tools/benchmarks/bayesian_scales/cv_refits_20261113"
source(file.path(study, "design.R"))
source(file.path(study, "reference.R"))
pkgload::load_all(compile = FALSE, quiet = TRUE)
dir <- file.path(study, "records")
csv <- function(nm) file.path(study, paste0(nm, ".csv"))
refs <- read.csv(csv("references"))
pits <- read.csv(csv("pits"))
# The raw fitted records and fitted-reference outputs preserve the frozen run.
# Capture every original oracle value before writing the corrected flat exports.
write.csv(refs, csv("references-frozen-oracle"), row.names = FALSE)
write.csv(pits, csv("pits-frozen-oracle"), row.names = FALSE)

# The oracle kernel originally exponentiated the unshifted squared distances,
# underflowing at some extremely small bandwidths. Subtracting the squared
# nearest distance leaves ratios invariant and guarantees one unit weight.
original_oracle <- oracle_score
stable_oracle <- original_oracle
code <- deparse(body(original_oracle), width.cutoff = 500L)
old <- "input$distances^2))"
new <- "(input$distances^2 - min(input$distances)^2)))"
stopifnot(sum(grepl(old, code, fixed = TRUE)) == 1L)
body(stable_oracle) <- parse(text = sub(old, new, code, fixed = TRUE))[[1L]]
original_joint <- is.finite(refs$oracle_nll)
original_pit <- is.finite(pits$oracle_pit)
audit <- list()
for (case in unique(refs$case)) {
  input <- simulate_cv_case(plan[plan$case == case, , drop = FALSE])
  for (fold in unique(refs$fold[refs$case == case])) {
    rows <- which(pits$case == case & pits$fold == fold)
    single <- rows[!duplicated(pits$visit[rows])]
    match_row <- single[match(pits$visit[rows], pits$visit[single])]
    stopifnot(identical(pits$randomized_u[rows], pits$randomized_u[match_row]),
      identical(pits$stress[rows], pits$stress[match_row]),
      identical(pits$observed[rows], pits$observed[match_row]))
    out <- stable_oracle(input, fold, pits$stress[single], pits$randomized_u[single])
    stopifnot(all(is.finite(unlist(out))))
    matching <- which(refs$case == case & refs$fold == fold)
    values <- unlist(out$score[refs$event[matching]], use.names = FALSE)
    old_values <- refs$oracle_nll[matching]
    ok <- is.finite(old_values)
    stopifnot(!any(ok) || max(abs(old_values[ok] - values[ok])) < 1e-8)
    joint_delta <- if (any(ok)) max(abs(old_values[ok] - values[ok])) else NA_real_
    old_pit <- pits$oracle_pit[rows]
    new_pit <- out$score$pit[match(pits$visit[rows], pits$visit[single])]
    ok <- is.finite(old_pit)
    stopifnot(!any(ok) || max(abs(old_pit[ok] - new_pit[ok])) < 1e-8)
    pit_delta <- if (any(ok)) max(abs(old_pit[ok] - new_pit[ok])) else NA_real_
    audit[[length(audit) + 1L]] <- data.frame(case = case, fold = fold,
      finite_joint_before = sum(original_joint[matching]), finite_joint_after = length(matching),
      finite_pit_before = sum(original_pit[rows]), finite_pit_after = length(rows),
      unchanged_joint_max = joint_delta, unchanged_pit_max = pit_delta)
    refs$oracle_nll[matching] <- values
    refs$oracle_refinement[matching] <- out$refine[refs$event[matching]]
    pits$oracle_pit[rows] <- new_pit
    pits$oracle_refinement[rows] <- out$refine[["pit"]]
  }
}
stopifnot(all(is.finite(refs$oracle_nll)), all(is.finite(pits$oracle_pit)))
write.csv(refs, csv("references"), row.names = FALSE)
write.csv(pits, csv("pits"), row.names = FALSE)
write.csv(do.call(rbind, audit), csv("oracle-repair-audit"), row.names = FALSE)
paths <- file.path(study, c("design.R", "reference.R", "repair-oracle.R"))
saveRDS(list(code = setNames(vapply(paths, function(path) digest::digest(
  readLines(path, warn = FALSE), algo = "sha256"), character(1L)), paths),
  R = R.version.string, statmod = as.character(packageVersion("statmod"))),
  file.path(study, "oracle-repair-configuration.rds"))

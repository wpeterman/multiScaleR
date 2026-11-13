# Read complete exported records and summarize paired contrasts within design cells.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(!length(args))
study <- "tools/benchmarks/bayesian_scales/cv_refits_20261113"
source(file.path(study, "design.R"))
read <- function(name) read.csv(file.path(study, paste0(name, ".csv")))
a <- read("attempts")
r <- read("references")
u <- read("augmentation")
p <- read("pits")
# The frozen fit runner appends both the data and score seeds. Give the two
# columns distinct names in the final CSV without altering any fit record.
if (all(c("seed", "seed.1") %in% names(u))) {
  names(u)[match(c("seed", "seed.1"), names(u))] <- c("data_seed", "score_seed")
  write.csv(u, file.path(study, "augmentation.csv"), row.names = FALSE)
}
stopifnot(nrow(a) == 4L * nrow(plan),
  !anyDuplicated(a[c("case", "model", "fold")]),
  all(a$status %in% c("scored", "diagnostic_failed", "failed", "pending")),
  !anyDuplicated(r[c("case", "model", "fold", "event")]),
  !anyDuplicated(u[c("case", "model", "fold", "event", "budget", "score_seed")]),
  !anyDuplicated(p[c("case", "model", "fold", "visit")]),
  all(is.na(a$replay_passed) | a$replay_passed),
  all(a$status != "scored" | is.finite(a$package_nll)),
  all(a$status != "scored" | a$n_analysis > 0),
  all(is.finite(r$reference_nll)), all(is.finite(r$oracle_nll)),
  all(is.finite(u$delta)),
  all(p$fit_pit >= 0 & p$fit_pit <= 1), all(p$oracle_pit >= 0 & p$oracle_pit <= 1))
key <- c("case", "family", "generator", "replicate", "model", "fold")
a <- merge(a, unique(r[c(key, "uses_latent")]), by = key, all.x = TRUE)
make_splits <- function(x, names) {
  labels <- interaction(x[names], drop = TRUE, lex.order = TRUE)
  split(x, labels)
}
aggregate_rows <- function(x, key, summarize) do.call(rbind,
  lapply(make_splits(x, key), function(cell)
    cbind(cell[1L, key, drop = FALSE], summarize(cell), row.names = NULL)))
q <- function(x, p) unname(quantile(x, p, names = FALSE))
replicate_id <- c("case", "family", "generator", "replicate", "fold")

availability <- aggregate_rows(a, c("family", "generator", "model", "fold"), function(x)
  data.frame(planned = nrow(x), finite = sum(is.finite(x$package_nll)),
    scored = sum(x$status == "scored"), rejected = sum(x$status == "diagnostic_failed"),
    failed = sum(x$status == "failed"), pending = sum(x$status == "pending"),
    replayed = sum(x$replay_passed %in% TRUE), with_reference = sum(!is.na(x$uses_latent))))

main <- r[r$event == "observed", ]
main$gap_vs_oracle <- (main$reference_nll - main$oracle_nll) / main$n
main$package_gap <- main$package_seed_nll - main$reference_nll
eligible <- main[main$status == "scored", ]
original_precision <- aggregate_rows(eligible[eligible$uses_latent, ], c("family", "fold"), function(x)
  data.frame(datasets = length(unique(x$case)), blocks = nrow(x),
    mean_error = mean(x$package_gap), sd_error = sd(x$package_gap),
    q95_abs = q(abs(x$package_gap), .95), max_abs = max(abs(x$package_gap))))
oracles <- aggregate_rows(eligible, c("family", "generator", "model", "fold"), function(x)
  data.frame(datasets = nrow(x), mean_gap = mean(x$gap_vs_oracle),
    sd_gap = sd(x$gap_vs_oracle), minimum_gap = min(x$gap_vs_oracle), maximum_gap = max(x$gap_vs_oracle)))

shared <- eligible[eligible$model == "pooled", c(replicate_id, "reference_nll", "package_seed_nll", "oracle_nll", "n")]
annual <- eligible[eligible$model == "annual", c(replicate_id, "reference_nll", "package_seed_nll", "oracle_nll", "n")]
pairs <- merge(shared, annual, by = replicate_id, suffixes = c("_pooled", "_annual"))
stopifnot(all(pairs$oracle_nll_pooled == pairs$oracle_nll_annual),
  all(pairs$n_pooled == pairs$n_annual))
pairs$annual_minus_pooled <- (pairs$reference_nll_annual - pairs$reference_nll_pooled) / pairs$n_annual
pairs$package_annual_minus_pooled <- (pairs$package_seed_nll_annual - pairs$package_seed_nll_pooled) / pairs$n_annual
pairs$ranking_changed <- sign(pairs$annual_minus_pooled) != sign(pairs$package_annual_minus_pooled)
pairs_by_cell <- aggregate_rows(pairs, c("family", "generator", "fold"), function(x) {
  n <- nrow(x)
  average <- mean(x$annual_minus_pooled)
  spread <- sd(x$annual_minus_pooled)
  margin <- if (n < 2L) NA_real_ else qt(.975, n - 1L) * spread / sqrt(n)
  data.frame(datasets = n, mean_delta = average, sd_delta = spread,
    min_delta = min(x$annual_minus_pooled), max_delta = max(x$annual_minus_pooled),
    se_delta = spread / sqrt(n), descriptive_t_lower = average - margin,
    descriptive_t_upper = average + margin, annual_better = sum(x$annual_minus_pooled < 0),
    ranking_changed = sum(x$ranking_changed))
})

ee <- merge(u[u$uses_latent, ], r[c(key, "event", "status", "refinement", "density_draws", "n")],
  by = c(key, "event"), all.x = TRUE)
stopifnot(!anyNA(ee$status))
ee <- ee[ee$status == "scored", ]
# Streams and cases overlap across candidates and designs. These are descriptive
# summaries of numerical perturbations, with an explicit number of independent data sets.
prec <- aggregate_rows(ee, c("family", "event", "budget"), function(x)
  data.frame(datasets = length(unique(x$case)), blocks = nrow(unique(x[c("case", "model", "fold")])),
    streams = nrow(x), signed_mean = mean(x$delta), mae = mean(abs(x$delta)),
    rmse = sqrt(mean(x$delta^2)), q95_abs = q(abs(x$delta), .95),
    max_abs = max(abs(x$delta)), beyond_point_one = sum(abs(x$delta) > .1),
    beyond_point_one_per_visit = sum(abs(x$delta) / x$n > .1)))
worst <- ee[order(-abs(ee$delta)), c(key, "event", "budget", "score_seed", "nll", "delta", "density_draws", "refinement")]
budget_pairs <- merge(ee[ee$event == "observed" & ee$model == "pooled", c(replicate_id, "budget", "score_seed", "nll")],
  ee[ee$event == "observed" & ee$model == "annual", c(replicate_id, "budget", "score_seed", "nll")],
  by = c(replicate_id, "budget", "score_seed"), suffixes = c("_pooled", "_annual"))
# In a whole-year fold the pooled candidate has no new effect and therefore
# exactly the same score at every augmentation budget and score seed.
whole <- ee[ee$fold == "whole_year" & ee$event == "observed" & ee$model == "annual", ]
whole <- merge(whole[c(replicate_id, "budget", "score_seed", "nll")],
  shared[c(replicate_id, "reference_nll")], by = replicate_id)
names(whole)[match(c("nll", "reference_nll"), names(whole))] <- c("nll_annual", "nll_pooled")
budget_pairs <- rbind(budget_pairs, whole[names(budget_pairs)])
budget_pairs <- merge(budget_pairs, pairs[c(replicate_id, "annual_minus_pooled", "n_annual")],
  by = replicate_id)
budget_pairs$delta <- (budget_pairs$nll_annual - budget_pairs$nll_pooled) / budget_pairs$n_annual
budget_pairs$ranking_changed <- sign(budget_pairs$delta) != sign(budget_pairs$annual_minus_pooled)
budget_ranking <- aggregate_rows(budget_pairs, c("family", "generator", "fold", "budget"), function(x)
  data.frame(datasets = length(unique(x$case)), streams = nrow(x),
    ranking_changes = sum(x$ranking_changed),
    cases_with_change = length(unique(x$case[x$ranking_changed])),
    max_abs_delta_error = max(abs(x$delta - x$annual_minus_pooled))))
refinement <- aggregate_rows(r, c("family", "event"), function(x)
  data.frame(refined_blocks = nrow(x), max_fitted = max(x$refinement),
    max_oracle = max(x$oracle_refinement), min_density_draws = min(x$density_draws),
    below_100_draws = sum(x$density_draws < 100)))
p <- p[p$status == "scored", ]
p$fit_central90 <- p$fit_pit >= .05 & p$fit_pit <= .95
p$oracle_central90 <- p$oracle_pit >= .05 & p$oracle_pit <= .95
pit_cases <- aggregate_rows(p, c("family", "generator", "case", "model", "fold"), function(x)
  data.frame(n = nrow(x), central90 = mean(x$fit_central90),
    oracle_central90 = mean(x$oracle_central90), mean_pit = mean(x$fit_pit),
    mean_oracle_pit = mean(x$oracle_pit)))
pit_summary <- aggregate_rows(pit_cases, c("family", "generator", "model", "fold"), function(x)
  data.frame(datasets = nrow(x), visits = sum(x$n), mean_central90 = mean(x$central90),
    sd_central90 = sd(x$central90), mean_oracle_central90 = mean(x$oracle_central90),
    mean_pit = mean(x$mean_pit), mean_oracle_pit = mean(x$mean_oracle_pit)))

results <- list(availability = availability, oracle_gaps = oracles, paired = pairs,
  paired_summary = pairs_by_cell, original_precision = original_precision,
  precision = prec, worst = worst, budget_pairs = budget_pairs,
  budget_ranking = budget_ranking, refinement = refinement,
  pit_cases = pit_cases, pit_summary = pit_summary)
for (name in names(results)) {
  x <- results[[name]]
  write.csv(x, file.path(study, paste0("summary-", name, ".csv")), row.names = FALSE)
  if (name %in% c("availability", "paired_summary", "precision", "refinement")) print(x, row.names = FALSE)
}
if (!any(a$status == "pending")) {
  capture.output(sessionInfo(), file = file.path(study, "summary-session.txt"))
}

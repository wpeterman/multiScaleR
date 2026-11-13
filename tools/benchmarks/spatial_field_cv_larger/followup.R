source("tools/benchmarks/spatial_field_cv_larger/run.R")
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "tools/benchmarks/spatial_field_cv_larger/results"
pkgload::load_all(compile = FALSE, quiet = TRUE)
cases <- read.csv(file.path(path, "scenarios.csv"), stringsAsFactors = FALSE)
cases <- subset(cases, cohort == "larger")
unresolved <- !cases$joint_resolved | cases$means_resolved < cases$n_assessment |
  cases$variances_resolved < cases$n_assessment
nb <- subset(cases, family == "nbinom2" & joint_resolved)
tail_case <- nb$id[which.max(abs(nb$fitted_loss - nb$reference_loss)/nb$n_assessment)]
cases <- subset(cases, unresolved | id == tail_case)
score_rows <- moment_rows <- diagnostics <- list()
for (i in seq_len(nrow(cases))) {
  case <- cases[i, , drop = FALSE]
  d <- simulate_larger(case)
  score <- !case$joint_resolved || case$id == tail_case
  temperatures <- if (score) seq(0, 1, .125) else 0
  draws <- if (score) 3200L else 6400L
  result <- capture(function() {
    states <- sample_temperatures(d, case$seed + 873203L, temperatures, 700L, draws, 3L)
    subset <- monte_carlo_replicates(draws, 3L)
    diag <- do.call(rbind, lapply(states, attr, "diagnostics"))
    list(diag = diag,
      score = if (score) reference_score_mc(states, temperatures, subset) else NULL,
      moments = reference_moments_mc(d, states, subset))
  })
  cat(case$id, if (is.null(result$value)) result$message else "sampled", "\n")
  run <- result$value
  score_rows[[i]] <- data.frame(id = case$id, reason = if (case$id == tail_case)
      "resolved_tail_replication" else if (score) "unresolved_score" else "unresolved_moments",
    draws_per_chain = draws, temperatures = length(temperatures),
    primary_loss = case$reference_loss, primary_mcse = case$reference_loss_mcse,
    primary_resolved = case$joint_resolved, fitted_loss = case$fitted_loss,
    n_assessment = case$n_assessment, loss = if (is.null(run$score)) NA_real_ else run$score$loss,
    loss_mcse = if (is.null(run$score)) NA_real_ else run$score$mcse,
    minimum_overlap = if (is.null(run$score)) NA_real_ else min(run$score$overlap),
    max_rhat = if (is.null(run)) NA_real_ else max(run$diag$max_rhat),
    message = result$message, warnings = result$warnings)
  if (!is.null(run)) {
    diagnostics[[i]] <- cbind(id = case$id, run$diag)
    moment_rows[[i]] <- data.frame(id = case$id, assessment_row = seq_len(case$n_assessment),
      mean = run$moments$point[, "mean"], mean_mcse = run$moments$mcse[, "mean"],
      variance = run$moments$point[, "variance"], variance_mcse = run$moments$mcse[, "variance"],
      analysis_rhat = run$diag$max_rhat[1L])
  }
}
score <- do.call(rbind, score_rows)
score$resolved <- with(score, is.finite(loss_mcse) & loss_mcse/n_assessment <= .015 &
  minimum_overlap >= .02 & max_rhat <= 1.05)
write.csv(score, file.path(path, "followup-scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, diagnostics), file.path(path, "followup-temperatures.csv"), row.names = FALSE)
write.csv(do.call(rbind, moment_rows), file.path(path, "followup-moments.csv"), row.names = FALSE)
print(score, row.names = FALSE)

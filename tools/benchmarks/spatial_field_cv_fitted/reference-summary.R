args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "tools/benchmarks/spatial_field_cv_fitted/results"
dirs <- file.path(path, c("reference1", "reference2"))
read <- function(name) do.call(rbind, lapply(file.path(dirs, paste0("reference-", name, ".csv")),
  read.csv, stringsAsFactors = FALSE))
scores <- read("scenarios")
moments <- read("predictions")
temperatures <- read("temperatures")
stopifnot(nrow(scores) == 24L, !anyDuplicated(scores[c("id", "method")]),
  length(unique(scores$id)) == 12L)
scores$score_error_per_row <- (scores$public_loss - scores$reference_loss)/scores$n_assessment
scores$score_mcse_per_row <- scores$reference_loss_mcse/scores$n_assessment
moments$abs_relative_mean_error <- abs(moments$public_mean/moments$reference_mean - 1)
moments$abs_relative_variance_error <- abs(moments$public_variance/moments$reference_variance - 1)
moments$inside90 <- moments$reference_pit >= .05 & moments$reference_pit <= .95
score_table <- do.call(rbind, lapply(split(scores, list(scores$family, scores$method), drop = TRUE),
  function(x) {
    resolved <- x$reference_status == "resolved" & is.finite(x$public_loss)
    data.frame(family = x$family[1L], method = x$method[1L], attempted = nrow(x),
      no_fit = sum(x$reference_status == "no_fit"), failed = sum(x$reference_status == "failed"),
      unresolved = sum(x$reference_status == "unresolved"), compared = sum(resolved),
      mean_abs_score_error = mean(abs(x$score_error_per_row[resolved])),
      max_abs_score_error = if (any(resolved)) max(abs(x$score_error_per_row[resolved])) else NA_real_,
      max_score_mcse = if (any(resolved)) max(x$score_mcse_per_row[resolved]) else NA_real_)
  }))
moment_table <- do.call(rbind, lapply(split(moments, list(moments$family, moments$method), drop = TRUE),
  function(x) {
    err <- function(value, flag) value[flag & is.finite(value)]
    m <- err(x$abs_relative_mean_error, x$mean_resolved)
    v <- err(x$abs_relative_variance_error, x$variance_resolved)
    pit <- subset(x, pit_resolved)
    per_case <- split(pit, pit$id)
    coverage <- vapply(per_case, function(z) mean(z$inside90), 0)
    data.frame(family = x$family[1L], method = x$method[1L], observed_rows = nrow(x),
      means_compared = length(m), max_mean_error = if (length(m)) max(m) else NA_real_,
      variance_compared = length(v), max_variance_error = if (length(v)) max(v) else NA_real_,
      pit_resolved_rows = nrow(pit), pit_fields = length(coverage),
      pit_inside90 = mean(coverage), pit_max_refinement = max(x$pit_refinement_gap))
  }))
scores <- scores[order(scores$id, scores$method), ]
moments <- moments[order(moments$id, moments$method, moments$row), ]
write.csv(scores, file.path(path, "reference-scores-combined.csv"), row.names = FALSE)
write.csv(moments, file.path(path, "reference-moments-combined.csv"), row.names = FALSE)
write.csv(score_table, file.path(path, "reference-score-summary.csv"), row.names = FALSE)
write.csv(moment_table, file.path(path, "reference-moment-summary.csv"), row.names = FALSE)
print(score_table, row.names = FALSE)
print(moment_table, row.names = FALSE)
cat("Largest temperature R-hat:", max(temperatures$max_rhat), "\n")

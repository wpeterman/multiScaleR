args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else "tools/benchmarks/spatial_field_cv_larger/results"
design <- read.csv(file.path(path, "design.csv"), stringsAsFactors = FALSE)
cases <- read.csv(file.path(path, "scenarios.csv"), stringsAsFactors = FALSE)
predictions <- read.csv(file.path(path, "predictions.csv"), stringsAsFactors = FALSE)
temperatures <- read.csv(file.path(path, "temperatures.csv"), stringsAsFactors = FALSE)
stopifnot(identical(cases$id, design$id), !anyDuplicated(cases$id),
  identical(as.integer(table(factor(predictions$id, levels = cases$id))), cases$n_assessment))
predictions <- cbind(predictions, cases[match(predictions$id, cases$id),
  c("family", "cohort", "field_sd", "setting", "visits")])
cases$joint_abs_error <- abs(cases$fitted_loss - cases$reference_loss) / cases$n_assessment
cases$joint_mcse <- cases$reference_loss_mcse / cases$n_assessment
predictions$mean_abs_error <- abs(predictions$fitted_mean - predictions$reference_mean)
predictions$mean_scaled_error <- predictions$mean_abs_error
log_family <- predictions$family != "binomial"
predictions$mean_scaled_error[log_family] <-
  predictions$mean_abs_error[log_family] / abs(predictions$reference_mean[log_family])
predictions$variance_scaled_error <- abs(predictions$fitted_variance / predictions$reference_variance - 1)

summary <- function(error, eligible, mcse) {
  ok <- eligible & is.finite(error) & is.finite(mcse)
  n <- sum(ok)
  q <- if (n) unname(quantile(error[ok], c(.5, .95, 1))) else rep(NA_real_, 3L)
  data.frame(attempted = length(error), comparable = n, p50 = q[1L], p95 = q[2L], max = q[3L],
    median_reference_mcse = if (n) median(mcse[ok]) else NA_real_,
    discrepancies_over_2_mcse = sum(error[ok] > 2 * mcse[ok]))
}

large <- subset(cases, cohort == "larger")
large_pred <- subset(predictions, cohort == "larger")
score_groups <- c(list(all = large), split(large, large$family),
  split(large, paste(large$family, large$field_sd)), split(large, large$setting))
score_summary <- do.call(rbind, Map(function(name, x) {
  cbind(group = name, metric = "joint_nats_per_observation",
    summary(x$joint_abs_error, x$joint_resolved & x$evaluation_status == "returned", x$joint_mcse))
}, names(score_groups), score_groups))
mean_groups <- c(split(large_pred, large_pred$family),
  split(large_pred, paste(large_pred$family, large_pred$field_sd)),
  split(large_pred, paste(large_pred$family, large_pred$setting)))
mean_summary <- do.call(rbind, Map(function(name, x) {
  binary <- x$family[1L] == "binomial"
  mean_mcse <- if (binary) x$reference_mean_mcse else x$reference_mean_mcse/abs(x$reference_mean)
  rbind(cbind(group = name, metric = if (binary) "mean_absolute_probability" else "mean_relative",
    summary(x$mean_scaled_error, x$mean_resolved, mean_mcse)),
    cbind(group = name, metric = "variance_relative",
      summary(x$variance_scaled_error, x$variance_resolved, x$reference_variance_mcse/x$reference_variance)))
}, names(mean_groups), mean_groups))
summaries <- rbind(score_summary, mean_summary)
rownames(summaries) <- NULL
write.csv(summaries, file.path(path, "error-summary.csv"), row.names = FALSE)
control <- subset(cases, cohort == "quadrature_control")
control$joint_reference_error <- abs(control$quadrature_loss - control$reference_loss)
cp <- subset(predictions, cohort == "quadrature_control")
cp$mean_reference_error <- abs(cp$quadrature_mean - cp$reference_mean)
cp$variance_reference_error <- abs(cp$quadrature_variance - cp$reference_variance)
write.csv(control, file.path(path, "score-controls.csv"), row.names = FALSE)
write.csv(cp, file.path(path, "moment-controls.csv"), row.names = FALSE)
cat("Larger-field score resolutions\n")
print(with(large, table(family, joint_resolved, evaluation_status)))
cat("\nLarger-field response moment resolutions\n")
print(with(large_pred, table(family, mean_resolved, variance_resolved)))
cat("\nJoint loss error at resolved cases\n")
print(score_summary, row.names = FALSE)
cat("\nResponse moments at resolved cases\n")
print(mean_summary[mean_summary$group %in% unique(large$family), ], row.names = FALSE)
cat("\nPositive-control bridge estimates\n")
print(control[c("small_id", "family", "joint_resolved", "reference_loss", "quadrature_loss",
  "joint_reference_error", "reference_loss_mcse")], row.names = FALSE)
cat("\nReference diagnostics\n")
print(summary(temperatures$max_rhat, rep(TRUE, nrow(temperatures)), numeric(nrow(temperatures))))

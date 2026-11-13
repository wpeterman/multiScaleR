# Reproduce the grouped error table after run.R has completed.
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1L] else
  "tools/benchmarks/spatial_field_cv_quadrature/results"
planned <- read.csv(file.path(path, "design.csv"), stringsAsFactors = FALSE)
cases <- read.csv(file.path(path, "scenarios.csv"), stringsAsFactors = FALSE)
means <- read.csv(file.path(path, "predictions.csv"), stringsAsFactors = FALSE)
stopifnot(identical(planned$id, cases$id), !anyDuplicated(cases$id),
  identical(as.integer(table(factor(means$id, levels = cases$id))), cases$n_assessment))
means <- cbind(means, cases[match(means$id, cases$id),
  c("cohort", "family", "field_sd", "visits", "rho", "layout", "covariance")])

quantities <- function(error) {
  x <- error[is.finite(error)]
  percentile <- if (length(x)) unname(quantile(x, c(.5, .9, .95, .99, 1))) else rep(NA_real_, 5L)
  exceeds <- if (length(x)) vapply(c(.01, .05, .1), function(z) sum(x > z), 0L) else rep(NA_integer_, 3L)
  as.list(setNames(c(length(x), percentile, exceeds),
    c("n", "median", "p90", "p95", "p99", "max", "over_0.01", "over_0.05", "over_0.10")))
}
group_results <- function(variable, rows) {
  p <- means[means$id %in% rows$id, , drop = FALSE]
  comparison <- data.frame(group = variable, scenarios = nrow(rows),
    reference_failed = sum(rows$reference_status == "failed"),
    reference_unresolved = sum(rows$reference_status == "unresolved"),
    evaluation_failed = sum(rows$evaluation_status == "failed"),
    comparable = sum(is.finite(rows$loss_error)))
  score <- list(
    abs_joint_nats = abs(rows$loss_error),
    abs_nats_per_observation = abs(rows$loss_error)/rows$n_assessment)
  # Relative count/continuous mean errors and absolute Bernoulli probabilities
  # have different interpretations. Do not pool them across families.
  if (length(unique(rows$family)) == 1L) {
    unit <- if (rows$family[1L] == "binomial") "abs_probability" else
      if (rows$family[1L] == "gaussian") "abs_response" else "relative_mean"
    score[[paste0("row_", unit)]] <- p$scaled_abs_mean_error
    score[[paste0("case_max_", unit)]] <- rows$max_scaled_mean_error
  }
  do.call(rbind, lapply(names(score), function(metric)
    cbind(comparison, data.frame(metric = metric), as.data.frame(quantities(score[[metric]])))))
}

groups <- list(all_Laplace = subset(cases, cohort == "Laplace"),
  zero_field = subset(cases, cohort == "zero_field"),
  exact_gaussian = subset(cases, cohort == "exact_gaussian"))
laplace <- groups$all_Laplace
for (variables in list("family", "field_sd", "visits", "rho", "layout", "covariance",
    c("family", "field_sd"), c("family", "field_sd", "visits"))) {
  names_for_row <- apply(laplace[variables], 1L, function(v)
    paste(paste(variables, v, sep = "="), collapse = " "))
  groups <- c(groups, split(laplace, names_for_row))
}
summary <- do.call(rbind, Map(group_results, names(groups), groups))
rownames(summary) <- NULL
write.csv(summary, file.path(path, "error-summary.csv"), row.names = FALSE)
print(summary[summary$group %in% c("all_Laplace", "zero_field", "exact_gaussian",
  paste0("family=", unique(laplace$family))), ], row.names = FALSE)
cat("\nReference/evaluation status\n")
print(with(cases, table(cohort, reference_status, evaluation_status)))
cat("\nLargest absolute loss errors per assessment observation\n")
top <- laplace[order(-abs(laplace$loss_error)/laplace$n_assessment, na.last = NA), ]
print(head(top[c("id", "family", "field_sd", "rho", "visits", "layout", "covariance",
  "n_assessment", "reference_loss", "evaluated_loss", "loss_error")], 12L), row.names = FALSE)

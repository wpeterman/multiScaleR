# Run after run.R from the package source root.
directory <- "tools/benchmarks/mixed_cv_replicated"
read <- function(name) utils::read.csv(file.path(directory, paste0(name, ".csv")))
write <- function(x, name) utils::write.csv(x, file.path(directory, paste0(name, ".csv")), row.names = FALSE)
plan <- read("plan")
fold <- read("fold-results")
groups <- read("group-results")
rows <- read("row-results")
stopifnot(identical(fold$id, plan$id), !anyDuplicated(fold$id),
  setequal(unique(groups$id), plan$id),
  !anyDuplicated(groups[c("id", "group")]), !anyDuplicated(rows[c("id", "row_id")]),
  all(groups$reference_status %in% c("available", "failed", "not_run_fit_failed")),
  all(fold$fit_status %in% c("available", "failed")),
  identical(unname(as.integer(tapply(groups$n_test, groups$id, sum))), fold$n_assessment))
keys <- c("model", "scenario", "prediction")
groups <- merge(groups, fold[c("id", keys)], by = "id", sort = FALSE)
rows <- merge(rows, fold[c("id", keys)], by = "id", sort = FALSE)
rows <- merge(rows, groups[c("id", "group", "reference_status")], by = c("id", "group"), sort = FALSE)
row_counts <- table(rows$id)
fit_available <- fold$fit_status == "available"
stopifnot(all(fold$n_assessment[fit_available] == unname(row_counts[fold$id[fit_available]])),
  all(groups$reference_status[groups$id %in% fold$id[!fit_available]] == "not_run_fit_failed"),
  all(is.finite(rows$estimate)),
  all(is.finite(rows$reference[rows$reference_status == "available"])),
  all(is.na(rows$reference[rows$reference_status == "failed"])))

counts <- function(f, g) {
  status <- g$reference_status
  data.frame(folds = nrow(f), fitted_and_scored = sum(f$fit_status == "available"),
    failed_folds = sum(f$fit_status == "failed"),
    folds_with_warnings = sum(!is.na(f$warnings) & nzchar(f$warnings)),
    groups = nrow(g), reference_available = sum(status == "available"),
    reference_failed = sum(status == "failed"),
    reference_not_run = sum(status == "not_run_fit_failed"))
}
labels <- unique(fold[keys])
availability <- do.call(rbind, lapply(seq_len(nrow(labels)), function(i) {
  label <- labels[i, , drop = FALSE]
  same <- function(d) Reduce(function(x, y) x & y, lapply(keys, function(k) d[[k]] == label[[k]]))
  cbind(label, counts(fold[same(fold), , drop = FALSE], groups[same(groups), , drop = FALSE]))
}))
write(availability, "availability")
tail_quantiles <- function(values) if (length(values)) as.numeric(stats::quantile(values,
  probs = c(.5, .95, .99, 1), names = FALSE)) else rep(NA_real_, 4L)
factors <- c("all", unique(fold$model))
tail_errors <- do.call(rbind, lapply(factors, function(model) {
  a <- groups[groups$reference_status == "available" & (model == "all" | groups$model == model), ]
  b <- rows[rows$reference_status == "available" & (model == "all" | rows$model == model), ]
  data.frame(model = model, metric = rep(c("joint_log_score_absolute",
    "response_mean_absolute", "response_mean_scaled"), each = 4L),
    percentile = rep(c(50, 95, 99, 100), times = 3L),
    n = rep(c(nrow(a), nrow(b), nrow(b)), each = 4L),
    error = c(tail_quantiles(a$absolute_score_error), tail_quantiles(b$absolute_error),
      tail_quantiles(b$scaled_error)))
}))
write(tail_errors, "error-summary")
print(counts(fold, groups))
cat("Availability by model\n")
print(do.call(rbind, lapply(unique(fold$model), function(model)
  cbind(model = model, counts(fold[fold$model == model, ], groups[groups$model == model, ])))),
  row.names = FALSE)
cat("Availability by scenario and target\n")
sets <- unique(fold[c("scenario", "prediction")])
print(do.call(rbind, lapply(seq_len(nrow(sets)), function(i) {
  same <- function(d) d$scenario == sets$scenario[i] & d$prediction == sets$prediction[i]
  cbind(sets[i, , drop = FALSE], counts(fold[same(fold), ], groups[same(groups), ]))
})), row.names = FALSE)
cat("Tail errors\n")
print(tail_errors, row.names = FALSE)
cat("Failed folds\n")
print(fold[!fit_available, c("id", "model", "scenario", "seed", "prediction", "failure_stage",
  "failure_message", "warnings")], row.names = FALSE)
cat("Failed references\n")
print(groups[groups$reference_status == "failed", c("id", "group", keys, "reference_message")],
  row.names = FALSE)
cat("Worst finite comparisons\n")
print(groups[which.max(groups$absolute_score_error), ], row.names = FALSE)
print(rows[which.max(rows$absolute_error), ], row.names = FALSE)
print(rows[which.max(rows$scaled_error), ], row.names = FALSE)
cat("Largest scaled quadrature refinement discrepancy\n")
print(max(groups$refinement, na.rm = TRUE))

# Summarize the joint_mixed Wald coverage records written by run.R.
#
# Run from the source root:
#   Rscript tools/benchmarks/joint_mixed_wald/summarize.R <output_dir>
# Reads <output_dir>/records/*.rds and writes results/metrics.csv (one row per
# dataset) and results/scenario-summary.csv (one row per family and sigma, plus
# pooled rows) next to this script.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply the run.R output directory.")
study <- "tools/benchmarks/joint_mixed_wald"
files <- list.files(file.path(args[1L], "records"), pattern = "[.]rds$", full.names = TRUE)
records <- lapply(files, readRDS)
columns <- unique(unlist(lapply(records, names)))
metrics <- do.call(rbind, lapply(records, function(r) {
  r[setdiff(columns, names(r))] <- NA
  r[columns]
}))
metrics <- metrics[order(metrics$family, metrics$sigma, metrics$replicate), ]
dir.create(file.path(study, "results"), showWarnings = FALSE)
utils::write.csv(metrics, file.path(study, "results", "metrics.csv"), row.names = FALSE)

wilson <- function(x, n, z = stats::qnorm(.975)) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

summarize_cell <- function(d, label_family, label_sigma) {
  ok <- d[d$fit_ok, , drop = FALSE]
  interior <- ok[!ok$at_lower & !ok$at_upper, , drop = FALSE]
  raw <- interior[is.finite(interior$raw_log_se), , drop = FALSE]
  reported <- ok[is.finite(ok$reported_log_se), , drop = FALSE]
  lr <- ok[is.finite(ok$lr_truth), , drop = FALSE]
  raw_ci <- wilson(sum(raw$raw_cover95), nrow(raw))
  rep_ci <- wilson(sum(reported$reported_cover95), nrow(reported))
  lr_ci <- wilson(sum(lr$lr_cover95), nrow(lr))
  data.frame(family = label_family, true_sigma = label_sigma,
    attempts = nrow(d), fits = nrow(ok), at_lower = sum(ok$at_lower),
    at_upper = sum(ok$at_upper), interior = nrow(interior),
    # Bound-only rule: every interior scale with a finite sdreport SE.
    raw_available = nrow(raw), raw_cover95 = sum(raw$raw_cover95),
    raw_rate95 = mean(raw$raw_cover95), raw_lo95 = raw_ci[1L], raw_hi95 = raw_ci[2L],
    raw_cover90 = sum(raw$raw_cover90), raw_miss_low = sum(raw$raw_above),
    raw_miss_high = sum(raw$raw_below),
    # Shipped rule: SEs that pass every per-scale check.
    reported = nrow(reported), reported_cover95 = sum(reported$reported_cover95),
    reported_rate95 = mean(reported$reported_cover95),
    reported_lo95 = rep_ci[1L], reported_hi95 = rep_ci[2L],
    withheld_se_gt_one = sum(ok$reason == "log-kernel standard error exceeds one"),
    withheld_competing = sum(ok$reason == "competing scale modes within the 95% likelihood threshold"),
    withheld_bound = sum(ok$reason == "parameter at an optimization bound"),
    covered_over_attempts = sum(reported$reported_cover95),
    # Likelihood-ratio reference at the generating scale.
    lr_available = nrow(lr), lr_cover95 = sum(lr$lr_cover95),
    lr_rate95 = mean(lr$lr_cover95), lr_lo95 = lr_ci[1L], lr_hi95 = lr_ci[2L],
    lr_cover90 = sum(lr$lr_cover90), lr_negative = sum(lr$lr_truth < -1e-4),
    # Calibration: mean SE against the spread of log estimates (interior).
    median_log_error = stats::median(raw$log_error),
    sd_log_error = stats::sd(raw$log_error),
    rms_raw_log_se = sqrt(mean(raw$raw_log_se^2)),
    public_matches = all(ok$public_matches),
    median_seconds = stats::median(d$seconds),
    stringsAsFactors = FALSE)
}

cells <- split(metrics, list(metrics$family, metrics$true_sigma), drop = TRUE)
summary_table <- do.call(rbind, lapply(cells, function(d)
  summarize_cell(d, d$family[1L], as.character(d$true_sigma[1L]))))
pooled <- do.call(rbind, c(lapply(split(metrics, metrics$family), function(d)
  summarize_cell(d, d$family[1L], "all")), list(summarize_cell(metrics, "all", "all"))))
summary_table <- rbind(summary_table[order(summary_table$family,
  as.numeric(summary_table$true_sigma)), ], pooled)
rownames(summary_table) <- NULL
utils::write.csv(summary_table, file.path(study, "results", "scenario-summary.csv"), row.names = FALSE)
print(summary_table, digits = 3)
cat("\nFit errors:\n"); print(table(metrics$fit_error[!metrics$fit_ok]))
cat("\nWithholding reasons:\n"); print(table(metrics$family, metrics$reason, useNA = "ifany"))
cat("\nWarnings (nonempty):", sum(nzchar(metrics$warnings)), "\n")
print(utils::head(table(metrics$warnings[nzchar(metrics$warnings)]), 10))
cat("\nTotal fit minutes:", round(sum(metrics$seconds) / 60, 1), "\n")

# Summarize the joint_mixed profile-interval records written by profile_run.R.
#
# Run from the source root after summarize.R has written results/metrics.csv:
#   Rscript tools/benchmarks/joint_mixed_wald/profile_summarize.R <output_dir>
# Reads <output_dir>/profile_records/*.rds, joins them by dataset with the Wald
# study's metrics (likelihood-ratio decision at the true sigma and the Wald
# limits), and writes results/profile-metrics.csv (one row per dataset) and
# results/profile-summary.csv (one row per family and sigma, plus pooled rows).
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply the profile_run.R output directory.")
study <- "tools/benchmarks/joint_mixed_wald"
files <- list.files(file.path(args[1L], "profile_records"), pattern = "[.]rds$",
                    full.names = TRUE)
records <- lapply(files, readRDS)
columns <- unique(unlist(lapply(records, names)))
profile <- do.call(rbind, lapply(records, function(r) {
  r[setdiff(columns, names(r))] <- NA
  r[columns]
}))
wald <- utils::read.csv(file.path(study, "results", "metrics.csv"), stringsAsFactors = FALSE)
keep <- c("id", "estimate", "nll_hat", "raw_lower95", "raw_upper95", "raw_cover95",
          "lr_truth", "lr_cover95")
wald <- wald[keep]
names(wald)[-1L] <- paste0("wald_", names(wald)[-1L])
metrics <- merge(profile, wald, by = "id", all.x = TRUE, sort = FALSE)
metrics <- metrics[order(metrics$family, metrics$sigma, metrics$replicate), ]
# The profile run refits each dataset; its estimate should match the Wald study.
metrics$log_estimate_change <- log(metrics$estimate / metrics$wald_estimate)
metrics$agree_lr <- metrics$profile_cover95 == metrics$wald_lr_cover95
# Interval shape on the log scale: width relative to the Wald interval, and
# upper over lower half-width (1 for a symmetric interval).
metrics$width_ratio <- log(metrics$profile_upper / metrics$profile_lower) /
  log(metrics$wald_raw_upper95 / metrics$wald_raw_lower95)
metrics$asymmetry <- log(metrics$profile_upper / metrics$estimate) /
  log(metrics$estimate / metrics$profile_lower)
utils::write.csv(metrics, file.path(study, "results", "profile-metrics.csv"), row.names = FALSE)

wilson <- function(x, n, z = stats::qnorm(.975)) {
  if (!n) return(c(NA_real_, NA_real_))
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(centre - half, centre + half)
}

summarize_cell <- function(d, label_family, label_sigma) {
  ok <- d[d$fit_ok & d$profile_ok, , drop = FALSE]
  ci <- wilson(sum(ok$profile_cover95), nrow(ok))
  wald_ok <- ok[!is.na(ok$wald_raw_cover95), , drop = FALSE]
  data.frame(family = label_family, true_sigma = label_sigma,
    attempts = nrow(d), fits = sum(d$fit_ok), profiles = nrow(ok),
    profile_errors = sum(d$fit_ok & !d$profile_ok),
    bound_limited = sum(ok$method == "profile-bound-limited"),
    profile_cover95 = sum(ok$profile_cover95), profile_rate95 = mean(ok$profile_cover95),
    profile_lo95 = ci[1L], profile_hi95 = ci[2L],
    miss_truth_above = sum(ok$profile_above), miss_truth_below = sum(ok$profile_below),
    lr_cover95 = sum(ok$wald_lr_cover95), agree_lr = sum(ok$agree_lr),
    wald_cover95 = sum(wald_ok$wald_raw_cover95), wald_available = nrow(wald_ok),
    median_width_ratio = stats::median(ok$width_ratio, na.rm = TRUE),
    median_asymmetry = stats::median(ok$asymmetry, na.rm = TRUE),
    max_abs_log_estimate_change = max(abs(ok$log_estimate_change)),
    median_profile_seconds = stats::median(ok$profile_seconds),
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
utils::write.csv(summary_table, file.path(study, "results", "profile-summary.csv"),
                 row.names = FALSE)
print(summary_table, digits = 3)
cat("\nProfile errors:\n"); print(table(metrics$profile_error[metrics$fit_ok & !metrics$profile_ok]))
cat("\nDisagreements with the likelihood-ratio decision:\n")
print(metrics[metrics$fit_ok & metrics$profile_ok & !metrics$agree_lr,
  c("id", "estimate", "profile_lower", "profile_upper", "true_sigma", "wald_lr_truth", "method")])
cat("\nWarnings (nonempty):", sum(nzchar(metrics$profile_warnings), na.rm = TRUE), "\n")
cat("\nTotal profile minutes:", round(sum(metrics$profile_seconds, na.rm = TRUE) / 60, 1), "\n")

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1L] else "tools/benchmarks/spatial_field_cv_fitted/results"
dirs <- file.path(root, c("part1", "part2"))
read <- function(name) {
  files <- file.path(dirs, paste0(name, ".csv"))
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}
scenarios <- read("scenarios"); predictions <- read("predictions"); starts <- read("starts")
planned <- read("design")
stopifnot(!anyDuplicated(scenarios$id), !anyDuplicated(planned$id), nrow(planned) == 360L,
  setequal(scenarios$id, planned$id))
scenarios <- scenarios[order(scenarios$id), ]
stopifnot(identical(as.integer(table(factor(predictions$id, levels = scenarios$id))), scenarios$n_assessment))
stopifnot(!anyDuplicated(starts[c("id", "attempt")]))
eligible_starts <- rowSums(table(factor(starts$id, levels = scenarios$id), starts$eligible)[, "TRUE", drop = FALSE])
scenarios$eligible_starts <- eligible_starts
scenarios$paired <- with(scenarios, fit_status == "scored" & generating_status == "scored")
scenarios$joint_difference <- with(scenarios, (fitted_loss - generating_loss)/n_assessment)
scenarios$fit_disagreement <- with(scenarios, eligible_starts == 2L &
  (start_value_gap > 1e-5 | start_parameter_gap > 1e-3))
stopifnot(!anyNA(scenarios$fit_disagreement[scenarios$fit_status == "scored" & eligible_starts == 2L]))
groups <- split(scenarios, list(scenarios$family, scenarios$design), drop = TRUE)
keys <- names(groups)
set.seed(903012)
boot_mean <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= 1L) return(c(estimate = if (n) x else NA_real_, lower = NA_real_, upper = NA_real_))
  values <- vapply(seq_len(4000L), function(i) mean(x[sample.int(n, replace = TRUE)]), 0)
  setNames(c(mean(x), quantile(values, c(.025, .975), names = FALSE)), c("estimate", "lower", "upper"))
}
availability <- do.call(rbind, lapply(groups, function(x) {
  n <- nrow(x); scored <- sum(x$fit_status == "scored")
  interval <- binom.test(scored, n)$conf.int
  data.frame(family = x$family[1L], design = x$design[1L], attempted = n,
    fitted_scored = scored, generating_scored = sum(x$generating_status == "scored"),
    paired_scored = sum(x$paired), availability = scored/n, lower = interval[1L], upper = interval[2L],
    failed_initial_fit = sum(x$fit_status == "failed" & x$fit_stage == "initial_field"),
    failed_prediction = sum(x$fit_status == "failed" & x$fit_stage == "joint_prediction"),
    failed_other = sum(x$fit_status == "failed" & !x$fit_stage %in% c("initial_field", "joint_prediction")),
    one_eligible_start = sum(x$fit_status == "scored" & x$eligible_starts == 1L),
    two_eligible_starts = sum(x$fit_status == "scored" & x$eligible_starts == 2L),
    returned_start_disagreement = sum(x$fit_status == "scored" & x$fit_disagreement, na.rm = TRUE))
}))
scores <- do.call(rbind, lapply(groups, function(x) {
  p <- x[x$paired, ]
  cbind(data.frame(family = x$family[1L], design = x$design[1L], attempted = nrow(x),
    paired = nrow(p), fitted_loss_per_row = mean(p$fitted_loss/p$n_assessment),
    generating_loss_per_row_paired = mean(p$generating_loss/p$n_assessment),
    generating_loss_per_row_all = mean(x$generating_loss/x$n_assessment, na.rm = TRUE)),
    as.data.frame(as.list(boot_mean(p$joint_difference))))
}))
recovery <- do.call(rbind, lapply(groups, function(x) {
  p <- x[x$fit_status == "scored", ]
  quantity <- c("intercept", "slope", "field_sd", "range", if (x$family[1L] != "poisson") "dispersion")
  do.call(rbind, lapply(quantity, function(nm) {
    value <- p[[paste0("fitted_", nm)]]
    target <- x[[paste0("true_", nm)]][1L]
    ratio <- if (nm %in% c("intercept", "slope")) value - target else value/target
    q <- if (length(ratio)) unname(quantile(ratio, c(.1, .5, .9))) else rep(NA_real_, 3L)
    data.frame(family = x$family[1L], design = x$design[1L], parameter = nm, true_value = target,
      attempted = nrow(x), scored = nrow(p), scale = if (nm %in% c("intercept", "slope")) "difference" else "ratio",
      p10 = q[1L], median = q[2L], p90 = q[3L],
      outside_third_to_three = if (!nm %in% c("intercept", "slope") && length(ratio))
        mean(ratio < 1/3 | ratio > 3) else NA_real_)
  }))
}))
field_metrics <- function(y, mean, variance, pit) {
  error <- y - mean
  squared <- error^2
  c(bias = mean(error), mse = mean(squared), average_variance = mean(variance),
    standardized_mse = mean(squared/variance), rmse_to_pred_sd = sqrt(sum(squared)/sum(variance)),
    coverage90 = if (all(is.na(pit))) NA_real_ else mean(pit >= .05 & pit <= .95),
    pit_mean = if (all(is.na(pit))) NA_real_ else mean(pit))
}
rows_by_case <- split(predictions, predictions$id)
per_field <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  x <- scenarios[i, ]
  p <- rows_by_case[[x$id]]
  do.call(rbind, lapply(c("fitted", "generating"), function(nm) {
    metric <- field_metrics(p$observed, p[[paste0(nm, "_mean")]], p[[paste0(nm, "_variance")]],
      p[[paste0(nm, "_pit")]])
    cbind(x[c("id", "family", "design", "paired", "fit_status", "generating_status")],
      method = nm, as.data.frame(as.list(metric)))
  }))
}))
metric_names <- c("bias", "mse", "average_variance", "standardized_mse", "rmse_to_pred_sd", "coverage90", "pit_mean")
calibration <- do.call(rbind, lapply(groups, function(x) {
  these <- per_field[per_field$id %in% x$id, ]
  fits <- these[these$method == "fitted" & these$paired, ]
  gen <- these[these$method == "generating" & these$paired, ]
  gen <- gen[match(fits$id, gen$id), ]
  all_gen <- these[these$method == "generating" & these$generating_status == "scored", ]
  cohorts <- list(fitted_paired = fits, generating_paired = gen,
    generating_all = all_gen, fitted_minus_generating_paired = fits)
  do.call(rbind, lapply(names(cohorts), function(label) {
    values <- cohorts[[label]]
    do.call(rbind, lapply(metric_names, function(nm) {
      metric <- values[[nm]]
      if (label == "fitted_minus_generating_paired") metric <- metric - gen[[nm]]
      cbind(data.frame(family = x$family[1L], design = x$design[1L], cohort = label,
        metric = nm, attempted = nrow(x), scored = sum(is.finite(metric))),
        as.data.frame(as.list(boot_mean(metric))))
    }))
  }))
}))
for (name in c("scenarios", "availability", "scores", "recovery", "per_field", "calibration"))
  write.csv(get(name), file.path(root, paste0(name, "-combined.csv")), row.names = FALSE)
write.csv(starts, file.path(root, "starts-combined.csv"), row.names = FALSE)
write.csv(predictions, file.path(root, "predictions-combined.csv"), row.names = FALSE)
print(availability, row.names = FALSE)
print(scores, row.names = FALSE)
print(subset(calibration, metric %in% c("standardized_mse", "coverage90") &
  cohort %in% c("fitted_paired", "generating_paired")), row.names = FALSE)
cat("Native start warning/error frequency\n")
print(sort(table(starts$warnings[nzchar(starts$warnings)]), decreasing = TRUE))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the consolidated metrics CSV and one output directory.")
study <- "tools/benchmarks/analytic_wald_intervals"
source(file.path(study, "reference.R"))
metrics <- utils::read.csv(args[1L], stringsAsFactors = FALSE)
# Sample before looking at agreement. Assess curvature only where the default
# fit returned a Wald interval and the independent full-range set returned.
take <- with(metrics, level == .95 & replicate <= 4L & wald_available & profile_available)
planned <- metrics[metrics$level == .95 & metrics$replicate <= 4L, ]
planned$curvature_status <- ifelse(!planned$wald_available, "wald_unavailable",
  ifelse(!planned$profile_available, "profile_unavailable", "pending"))
rows <- lapply(which(take), function(i) {
  case <- metrics[i, , drop = FALSE]
  design <- wald_case(case)
  oracle <- wald_reference(design, case$family, case$backend)
  fn <- oracle$fn
  estimate <- log(case$estimated_scale / 1000)
  gap <- min(estimate - log(.05), log(.8) - estimate)
  h <- min(0.002, gap / 5)
  result <- wald_capture({
    if (!is.finite(h) || h <= 1e-5) stop("Insufficient distance from the bound for stable central differences.")
    curvature <- function(step) {
      values <- vapply(estimate + step * (-2:2), fn, 0)
      sum(values * c(-1, 16, -30, 16, -1)) / (12 * step^2)
    }
    values <- c(curvature(h), curvature(h / 2))
    if (any(!is.finite(values) | values <= 0)) stop("Independent profile curvature was nonpositive or nonfinite.")
    independent <- sqrt(1 / values[2L])
    discrepancy <- abs(values[1L] - values[2L]) / values[2L]
    if (discrepancy > 1e-3) stop("Independent curvature changed by over 0.1% after halving the step.")
    c(independent = independent, refinement = discrepancy)
  })
  good <- !inherits(result$value, "error")
  data.frame(case[c("id", "family", "scenario", "n", "replicate", "backend")],
    log_scale_step = h, package_log_se = case$wald_log_se, status = if (good) "available" else "failed",
    independent_log_se = if (good) result$value[["independent"]] else NA_real_,
    relative_se_error = if (good) abs(case$wald_log_se / result$value[["independent"]] - 1) else NA_real_,
    relative_refinement_error = if (good) result$value[["refinement"]] else NA_real_,
    message = if (good) "" else conditionMessage(result$value),
    warnings = paste(result$warnings, collapse = " | "))
})
observed <- do.call(rbind, rows)
plan_key <- paste(planned$id, planned$backend)
match_row <- match(paste(observed$id, observed$backend), plan_key)
stopifnot(!anyNA(match_row), !anyDuplicated(match_row))
planned$curvature_status[match_row] <- observed$status
stopifnot(!any(planned$curvature_status == "pending"))
out <- merge(planned[c("id", "family", "scenario", "n", "replicate", "backend", "curvature_status")],
  observed[setdiff(names(observed), c("family", "scenario", "n", "replicate", "status"))],
  by = c("id", "backend"), sort = FALSE, all.x = TRUE)
utils::write.csv(out, file.path(args[2L], "local-curvature.csv"), row.names = FALSE)
print(table(out$curvature_status))
print(summary(out$relative_se_error))

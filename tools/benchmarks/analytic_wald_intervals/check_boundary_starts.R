args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply the consolidated metrics CSV and one output directory.")
study <- "tools/benchmarks/analytic_wald_intervals"
source(file.path(study, "reference.R"))
pkgload::load_all(compile = FALSE, quiet = TRUE)
metrics <- utils::read.csv(args[1L], stringsAsFactors = FALSE)
metrics <- metrics[metrics$level == .95 & is.finite(metrics$global_nll_gap) & metrics$global_nll_gap > 1e-5, ]
metrics <- metrics[order(metrics$id, metrics$backend), , drop = FALSE]
dir.create(args[2L], showWarnings = FALSE, recursive = TRUE)
# This follow-up deliberately targets observed discrepancies. It does not
# replace any primary fit, uncertainty result, or coverage denominator.
fits <- lapply(seq_len(nrow(metrics)), function(i) {
  case <- metrics[i, , drop = FALSE]
  sample <- wald_case(case)
  family <- switch(case$family, gaussian = gaussian(), poisson = poisson(),
    binomial = binomial(), nbinom2 = MASS::negative.binomial(2))
  initial <- wald_capture(glm(sample$formula, data = sample$data, family = family))
  if (inherits(initial$value, "error")) stop(conditionMessage(initial$value))
  # Keep this starting design independent of the oracle winner. Both search
  # endpoints and a fixed set of additional interior points are included.
  starts <- matrix(c(.05, .075, .14, .25, .43, .66, .76, .8), ncol = 1L)
  refit <- wald_capture(multiScale_optim(initial$value, sample$inputs,
    backend = case$backend, n_cores = 1, verbose = FALSE,
    profile_control = list(starts = starts)))
  available <- !inherits(refit$value, "error")
  fit <- if (available) refit$value else NULL
  after <- if (available) fit$optim_results$value else NA_real_
  diagnostics <- data.frame(case[c("id", "family", "scenario", "n", "seed", "backend")],
    default_estimate = case$estimated_scale, default_nll = case$nll,
    oracle_estimate = case$profile_estimate, oracle_nll = case$nll - case$global_nll_gap,
    default_wald = case$wald_available,
    follow_up_returned = available, follow_up_estimate = if (available) fit$scale_est$Mean[1L] else NA_real_,
    follow_up_nll = after, follow_up_oracle_gap = after - (case$nll - case$global_nll_gap),
    follow_up_wald = if (available) isTRUE(fit$uncertainty$available) else NA,
    reasons = if (available) paste(fit$uncertainty$reasons, collapse = " | ") else conditionMessage(refit$value),
    warnings = paste(refit$warnings, collapse = " | "))
  result <- if (available) fit$profile_fit$start_results else refit$value$start_results
  if (!is.null(result)) result <- cbind(id = case$id, backend = case$backend, result)
  list(diagnostics = diagnostics, starts = result)
})
for (item in c("diagnostics", "starts")) {
  rows <- lapply(fits, function(x) x[[item]])
  fields <- unique(unlist(lapply(rows, names), use.names = FALSE))
  out <- do.call(rbind, lapply(rows, function(x) {
    if (is.null(x)) return(NULL)
    for (name in setdiff(fields, names(x))) x[[name]] <- NA
    x[fields]
  }))
  if (is.null(out)) out <- data.frame(id = character(), backend = character())
  utils::write.csv(out, file.path(args[2L], paste0("boundary-start-", item, ".csv")), row.names = FALSE)
}

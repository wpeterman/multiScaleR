#' Score Named Simulation Targets on Their Response Scale
#'
#' Summarize recovery of named parameters or expected quantities while keeping
#' failed and ineligible attempts in the denominator. Estimates and interval
#' endpoints are transformed before bias, variation, or coverage is calculated.
#'
#' @param estimates Data frame with columns `replicate`, `target`, `units`,
#'   `truth`, `estimate`, and logical `eligible`. Each replicate/target pair must
#'   be unique. `truth` is always on the response scale and must be constant
#'   within a target. Use explicit target names, such as `"density_intercept"`,
#'   rather than selecting coefficients by their positions. Optional columns
#'   `lower` and `upper` supply interval endpoints on `input_scale`. Optional
#'   character `reason` records failed-fit diagnostics. Logical `identified`
#'   may mark a target as unidentified; it must be constant within each target.
#'   All supplied columns are retained in the returned attempt table.
#' @param input_scale One of `"response"`, `"log"`, or `"logit"`, describing
#'   estimates and both interval endpoints. This argument is required. Call
#'   separately for targets estimated on different scales.
#' @param level Nominal interval coverage, between zero and one. This labels
#'   the supplied intervals; it does not construct or recalibrate them.
#' @return A list with class `msr_simulation_score` containing `summary`, one
#'   row per named target in first-appearance order, and `attempts`,
#'   every supplied attempt with response-scale estimates and endpoints,
#'   effective eligibility, failure reason, and interval status. `summary`
#'   contains target, units, truth, attempted/eligible/failed counts, eligible
#'   interval count, response-scale bias and RMSE, relative bias and coefficient
#'   of variation in percent, interval coverage, and reasons for unavailable
#'   metrics. Metrics condition on eligible attempts; failure counts describe
#'   the full supplied set. `input_scale` and `level` retain the scoring contract.
#'   Printing displays at most ten summary rows and returns the complete object
#'   invisibly. All target summaries and attempts remain available in the list.
#' @details
#' Relative bias is `100 * (mean(estimate) - truth) / truth`. The coefficient of
#' variation (CV) is `100 * sd(estimate) / mean(estimate)` using the sample
#' standard deviation. CV describes between-replicate variation, not a fitted
#' standard error. Bias and RMSE retain the target's response units. Coverage is
#' the proportion of eligible finite intervals containing the truth.
#'
#' Relative bias is undefined at zero truth; CV is undefined at zero mean or
#' with fewer than two eligible attempts. Invalid intervals do not silently
#' remove an otherwise eligible point estimate. Their status and denominator
#' are reported separately. A nonfinite back-transformed estimate is ineligible.
#'
#' Set `identified = FALSE` for ecological scales simulated under a zero
#' habitat effect. Such a simulation does not identify a scale even if its
#' optimizer returns one. All recovery metrics for that target are then missing.
#' Eligibility must also reflect the fit's original convergence, curvature,
#' boundary, and other diagnostic gates; this helper does not relax them.
#' @examples
#' attempts <- data.frame(replicate = 1:3, target = "density", units = "birds/km2",
#'   truth = 0.4, estimate = log(c(0.2, 0.4, 0.6)), eligible = TRUE)
#' score_simulation(attempts, input_scale = "log")$summary
#' @export
score_simulation <- function(estimates, input_scale, level = 0.95) {
  required <- c("replicate", "target", "units", "truth", "estimate", "eligible")
  if (!is.data.frame(estimates) || !nrow(estimates) ||
      any(!required %in% names(estimates)) || anyDuplicated(names(estimates)))
    stop("estimates needs uniquely named replicate, target, units, truth, estimate, and eligible columns.", call. = FALSE)
  if (missing(input_scale)) stop("Supply input_scale explicitly.", call. = FALSE)
  input_scale <- match.arg(input_scale, c("response", "log", "logit"))
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be a probability strictly between zero and one.", call. = FALSE)
  for (nm in c("replicate", "target", "units")) {
    x <- estimates[[nm]]
    if (!is.atomic(x) || !is.null(dim(x)) || anyNA(x) || any(!nzchar(trimws(as.character(x)))))
      stop(nm, " must contain complete, nonempty identifiers.", call. = FALSE)
  }
  if (anyDuplicated(estimates[c("replicate", "target")]))
    stop("Each replicate/target pair must occur only once.", call. = FALSE)
  for (nm in intersect(c("truth", "estimate", "eligible", "identified", "reason", "lower", "upper"),
                       names(estimates))) {
    if (!is.atomic(estimates[[nm]]) || !is.null(dim(estimates[[nm]])))
      stop(nm, " must be an atomic vector, not a matrix or list column.", call. = FALSE)
  }
  if (!is.numeric(estimates$truth) || any(!is.finite(estimates$truth)) ||
      !is.numeric(estimates$estimate) || !is.logical(estimates$eligible) || anyNA(estimates$eligible))
    stop("truth must be finite numeric; estimate numeric; eligible nonmissing logical.", call. = FALSE)
  identified <- if ("identified" %in% names(estimates)) estimates$identified else rep(TRUE, nrow(estimates))
  if (!is.logical(identified) || anyNA(identified))
    stop("identified must be nonmissing logical.", call. = FALSE)
  original <- estimates
  estimates$target <- as.character(estimates$target)
  estimates$units <- as.character(estimates$units)
  if ("reason" %in% names(estimates) && !is.character(estimates$reason))
    stop("reason must be a character column.", call. = FALSE)
  reason <- if ("reason" %in% names(estimates)) estimates$reason else rep("", nrow(estimates))
  reason[is.na(reason)] <- ""
  reason[!estimates$eligible & !nzchar(trimws(reason))] <- "excluded_by_caller"
  transform <- switch(input_scale, response = identity, log = exp, logit = stats::plogis)
  response <- transform(estimates$estimate)
  reason[estimates$eligible & !is.finite(response)] <- "nonfinite_response_estimate"
  eligible <- estimates$eligible & is.finite(response)
  has_interval <- c("lower", "upper") %in% names(estimates)
  if (any(has_interval) && !all(has_interval))
    stop("Supply both lower and upper, or neither.", call. = FALSE)
  lower <- upper <- rep(NA_real_, nrow(estimates))
  interval_reason <- rep("not_supplied", nrow(estimates))
  if (all(has_interval)) {
    if (!is.numeric(estimates$lower) || !is.numeric(estimates$upper))
      stop("lower and upper must be numeric.", call. = FALSE)
    lower <- transform(estimates$lower); upper <- transform(estimates$upper)
    interval_reason[] <- "nonfinite_endpoints"
    finite <- is.finite(lower) & is.finite(upper)
    interval_reason[finite] <- "reversed_endpoints"
    interval_reason[finite & lower <= upper] <- "eligible"
  }
  interval_reason[!eligible] <- "ineligible_estimate"
  # Group once rather than rescanning every attempt for each target.
  groups <- unname(split(seq_len(nrow(estimates)),
    factor(estimates$target, levels = unique(estimates$target))))
  rows <- lapply(groups, function(i) {
    target <- estimates$target[i[1L]]
    if (length(unique(estimates$truth[i])) != 1L || length(unique(estimates$units[i])) != 1L ||
        length(unique(identified[i])) != 1L)
      stop("Each target needs one response-scale truth, unit, and identification status.", call. = FALSE)
    truth <- estimates$truth[i[1L]]
    values <- response[i[eligible[i]]]
    n <- length(values)
    use_interval <- i[eligible[i] & interval_reason[i] == "eligible"]
    metric_reason <- if (!identified[i[1L]]) "unidentified_target" else if (!n) "no_eligible_attempts" else ""
    bias <- rb <- rmse <- cv <- coverage <- NA_real_
    rb_reason <- cv_reason <- coverage_reason <- metric_reason
    if (!nzchar(metric_reason)) {
      bias <- mean(values) - truth
      rmse <- sqrt(mean((values - truth)^2))
      if (truth == 0) rb_reason <- "zero_truth" else rb <- 100 * bias / truth
      if (n < 2L) cv_reason <- "fewer_than_two_eligible_attempts" else if (mean(values) == 0)
        cv_reason <- "zero_mean" else cv <- 100 * stats::sd(values) / mean(values)
      if (length(use_interval)) coverage <- mean(lower[use_interval] <= truth & upper[use_interval] >= truth)
      else coverage_reason <- "no_eligible_intervals"
    }
    data.frame(target = target, units = estimates$units[i[1L]], truth = truth,
      n_attempted = length(i), n_eligible = n, n_failed = length(i) - n,
      n_interval = length(use_interval), bias = bias, relative_bias_pct = rb,
      rmse = rmse, cv_pct = cv, coverage = coverage, metric_reason = metric_reason,
      relative_bias_reason = rb_reason, cv_reason = cv_reason,
      coverage_reason = coverage_reason, stringsAsFactors = FALSE)
  })
  # Keep supplied estimates/diagnostics intact; append the scoring decisions.
  added <- c("response_estimate", "response_lower", "response_upper", "score_eligible",
             "score_reason", "interval_status")
  if (any(added %in% names(original))) stop("estimates uses reserved scoring column names.", call. = FALSE)
  original$response_estimate <- response; original$response_lower <- lower; original$response_upper <- upper
  original$score_eligible <- eligible; original$score_reason <- reason
  original$interval_status <- interval_reason
  structure(list(summary = do.call(rbind, rows), attempts = original,
    input_scale = input_scale, level = level), class = "msr_simulation_score")
}

#' @rdname score_simulation
#' @param x An object returned by `score_simulation()`.
#' @param ... Further arguments passed to the data-frame print method.
#' @export
print.msr_simulation_score <- function(x, ...) {
  cat("Response-scale simulation scores (", x$input_scale, " input)\n", sep = "")
  print(utils::head(x$summary, 10L), ...)
  if (nrow(x$summary) > 10L)
    cat("Showing the first 10 of", nrow(x$summary), "targets; see $summary for all targets.\n")
  cat("Full attempt records and failure reasons: $attempts\n")
  invisible(x)
}

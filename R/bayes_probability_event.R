# Shared event diagnostics preserve the original scalar probability rules.
# Flattening uses iteration first, then chain, exactly as the retained sampler.
.msr_bayes_probability_event <- function(log_scale, threshold, iterations,
    chains, mcse_target = NULL) {
  metrics <- c(mcse = NA_real_, ess_mean = NA_real_, rhat = NA_real_)
  finite <- all(is.finite(log_scale)); reason <- ""
  indicator <- matrix(NA_real_, iterations, chains)
  if (finite) indicator[] <- as.numeric(log_scale <= threshold)
  constant <- if (finite) length(unique(as.vector(indicator))) == 1L else NA
  if (!finite) reason <- "Nonfinite derived log scales; no iterations were dropped." else
    if (constant) reason <- "Constant indicator; probability precision is nonestimable." else
      if (chains < 2L || iterations < 4L)
        reason <- "At least two chains with four iterations each are required." else {
          metrics <- tryCatch(c(mcse = posterior::mcse_mean(indicator),
            ess_mean = posterior::ess_mean(indicator), rhat = posterior::rhat(indicator)),
            error = function(e) { reason <<- conditionMessage(e); metrics })
          if (any(!is.finite(metrics[c("mcse", "ess_mean")])) ||
              metrics[["mcse"]] <= 0 || metrics[["ess_mean"]] <= 0)
            if (!nzchar(reason)) reason <- "Indicator precision diagnostics are nonestimable; inspect chain mixing."
        }
  usable <- !nzchar(reason)
  list(probabilities = data.frame(
    estimate = if (finite) mean(indicator) else NA_real_,
    mcse = unname(metrics[["mcse"]]), ess_mean = unname(metrics[["ess_mean"]]),
    rhat = unname(metrics[["rhat"]]), constant = constant,
    precision_available = usable, precision_target_met = if (usable && !is.null(mcse_target))
      metrics[["mcse"]] <= mcse_target else NA, reason = reason,
    rhat_available = is.finite(metrics[["rhat"]]),
    rhat_reason = if (is.finite(metrics[["rhat"]])) "" else if (nzchar(reason)) reason else
      if (is.na(metrics[["rhat"]])) "Indicator Rhat is undefined; mean precision is reported separately." else
        "Indicator Rhat is nonfinite; mean precision is reported separately."),
    chain_probabilities = data.frame(chain = seq_len(chains),
      estimate = if (finite) colMeans(indicator) else rep(NA_real_, chains),
      constant = if (finite) apply(indicator, 2L, function(x) length(unique(x)) == 1L) else
        rep(NA, chains)))
}

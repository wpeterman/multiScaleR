#' Check Sampling Precision of Fitted Bayesian Scales
#'
#' Check chain mixing and simulation precision for the quantities you report
#' from a sampled fit, the log scale of each landscape predictor and the
#' between-group SDs. These derived quantities can mix differently from the
#' individual parameters checked by the sampler's own diagnostics.
#'
#' @param object A fitted [multiScale_bayes()] object. New diagnostics need the
#'   original unpermuted post-warmup chains from NUTS (`sampler`) or slice
#'   sampling (`chains`). Prepared, approximate, and compact fits return an
#'   unavailable report.
#' @param probs Numeric vector of distinct probabilities strictly between 0
#'   and 1. The default `c(0.025, 0.5, 0.975)` gives the median and the limits
#'   of a central 95% interval. Changing it changes which quantiles have their
#'   simulation error reported, not any diagnostic threshold.
#'
#' @return A list with these elements.
#'   * `available`, `reason`: whether full chain diagnostics could be computed,
#'     and an explanation when they could not.
#'   * `raw_diagnostics`: the sampler or approximation diagnostics stored in
#'     the fit, unchanged.
#'   * `parameters`: a data frame with one row per derived quantity.
#'     `variable` is an internal label. `quantity` is `"log_distance"`,
#'     `"scale_group_sd"`, or `"response_group_sd"`. `covariate`, `group`, and
#'     a representative observation `row` identify it, and `units` states its
#'     units. `rhat`, `ess_bulk`, and `ess_tail` give the rank-normalized
#'     Rhat and the bulk and tail effective sample sizes.
#'   * `quantiles`: a data frame with the same identifying columns plus
#'     `probability`, the quantile `estimate`, and its Monte Carlo standard
#'     error `mcse`, both in the stated units. Undefined values are `NA`.
#'   * `chains`, `iterations`: the number of retained chains and of
#'     post-warmup iterations per chain, or `NA` when unavailable.
#'   * `derived_passed`: `TRUE` when there are at least two chains and every
#'     derived quantity has Rhat at most 1.01 and bulk and tail ESS of at
#'     least 100 per chain. `NA` when unavailable.
#'   * `passed`: `TRUE` when both the original sampler checks and the derived
#'     checks pass. `NA` when unavailable. For a slice fit sampled without
#'     the posterior package, the original parameter checks are computed here
#'     from the stored chains, while the fit's stored diagnostics stay as
#'     recorded.
#'   * `criteria`, `notes`: the fixed derived checks and their limits.
#'
#' @details
#' # Read the report
#'
#' | Column | Acceptable | If it fails |
#' |:-------|:-----------|:------------|
#' | `rhat` | 1.01 or less | Refit with another `seed` and compare. A weak landscape effect leaves its scale poorly identified, and longer chains do not cure that |
#' | `ess_bulk`, `ess_tail` | At least 100 per chain (400 for four chains) | Increase `iter`, or switch from slice sampling to NUTS |
#' | `mcse` | No fixed cutoff. Small next to the precision your question needs | Increase `iter`, and compare independent runs when an interval limit matters |
#'
#' Landscape scales are diagnosed on the natural log axis in map units, which
#' avoids overflow and the infinite arithmetic distance moments that the
#' hierarchical priors allow (see [posterior_scales()]). Exponentiate a
#' log-distance quantile to express it in map units. Its MCSE is on the log
#' axis, so an MCSE of 0.01 is about 1% of the distance. The MCSE is
#' computational uncertainty, not a posterior interval or a biological SD.
#' Scale-group SDs measure population variation in log distance.
#' Response-group SDs use the response linear-predictor axis.
#'
#' The posterior package computes rank-normalized Rhat, effective sample
#' sizes, and quantile MCSE. A passed screen guarantees neither a precise tail
#' endpoint nor discovery of every mode. Compare the MCSE with the precision
#' your scientific decision needs, inspect independent runs when tails
#' matter, and assess the likelihood and priors separately.
#'
#' # Why compact fits cannot be checked
#'
#' Permuted or thinned stored draws cannot rebuild the original chain order or
#' full-chain precision. Compact examples therefore keep their original
#' diagnostic evidence but cannot provide new chain diagnostics. Laplace
#' draws and weighted importance candidates are not MCMC chains, so use the
#' stored numerical checks for those fits.
#'
#' # Draw chain plots
#'
#' For chain plots in bayesplot, pass `fit$chains` for a slice fit, or
#' `rstan::extract(fit$sampler, permuted = FALSE, inc_warmup = FALSE)` for a
#' NUTS fit. `vignette("v13-bayesian-scales", package = "multiScaleR")` shows
#' trace plots alongside these diagnostics, and
#' `vignette("v15-bayesian-computation", package = "multiScaleR")` works
#' through a fit that fails them.
#'
#' @seealso [multiScale_bayes()], [posterior_scales()],
#'   [plot.multiScaleR_bayes()], [diagnose_bayes_probabilities()],
#'   [multiScaleR_bayes-class]
#' @examples
#' # A stored slice fit that keeps its four original chains
#' shared <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
#'   package = "multiScaleR"))$fit
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   report <- diagnose_bayes_scales(shared)
#'   print(report$passed)   # sampler and derived checks together
#'   # Agreement among chains and effective sample size for each log sigma
#'   print(report$parameters[c("covariate", "rhat", "ess_bulk", "ess_tail")])
#'   # Interval limits of log sigma, with their Monte Carlo standard errors
#'   q <- report$quantiles
#'   print(data.frame(q[c("covariate", "probability")],
#'                    sigma_m = exp(q$estimate), mcse_log = q$mcse))
#' }
#'
#' # A compact fit keeps thinned, permuted draws and cannot be checked
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'                               package = "multiScaleR"))
#' compact <- diagnose_bayes_scales(example$fits$pooled)
#' compact$available  # FALSE: the compact example contains no full sampler
#' compact$reason
#' @export
diagnose_bayes_scales <- function(object, probs = c(0.025, 0.5, 0.975)) {
  if (!inherits(object, "multiScaleR_bayes"))
    stop("object must be a multiScaleR_bayes fit.", call. = FALSE)
  if (!is.numeric(probs) || !is.null(dim(probs)) || !length(probs) ||
      any(!is.finite(probs) | probs <= 0 | probs >= 1) || anyDuplicated(probs))
    stop("probs must contain distinct finite probabilities strictly between 0 and 1.", call. = FALSE)
  metadata <- data.frame(variable = character(), quantity = character(),
    covariate = character(), group = character(), row = integer(), units = character())
  out <- list(available = FALSE, reason = "", raw_diagnostics = object$diagnostics,
    parameters = cbind(metadata, rhat = numeric(), ess_bulk = numeric(), ess_tail = numeric()),
    quantiles = cbind(metadata, probability = numeric(), estimate = numeric(), mcse = numeric()),
    chains = NA_integer_, iterations = NA_integer_, derived_passed = NA, passed = NA,
    criteria = "At least 2 chains; Rhat <= 1.01; bulk/tail ESS >= 100 per chain for every derived quantity.",
    notes = c("Quantile MCSE measures simulation error, not posterior uncertainty.",
      "Passed checks do not guarantee endpoint precision, complete mode exploration, or model adequacy."))
  if (!is.null(object$laplace) || !is.null(object$marginal) || !is.null(object$importance)) {
    out$reason <- "MCMC chain diagnostics are unavailable for approximate draws or importance candidates."
    return(out)
  }
  if (is.null(object$sampler) && is.null(object$chains)) {
    out$reason <- paste("Full unpermuted MCMC chains are not retained.",
      "Stored permuted or thinned draws cannot recover original chain diagnostics.")
    return(out)
  }
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Install posterior to compute derived sampling diagnostics.", call. = FALSE)
  a <- .msr_bayes_fit_chains(object)
  if (!is.numeric(a) || length(dim(a)) != 3L || any(dim(a) < 1L) ||
      is.null(dimnames(a)[[3L]]) || anyDuplicated(dimnames(a)[[3L]]))
    stop("Sampler output must have named iterations-by-chains-by-parameters dimensions.", call. = FALSE)
  out$iterations <- dim(a)[1L]; out$chains <- dim(a)[2L]
  # Flatten in column-major order, retaining iteration order within each chain.
  # Never use the permuted object$draws blocks for chain diagnostics.
  p <- lapply(names(object$parameters), function(nm) {
    n <- length(object$parameters[[nm]])
    columns <- if (n == 1L && nm %in% dimnames(a)[[3L]]) nm else
      paste0(nm, "[", seq_len(n), "]")
    if (n && !all(columns %in% dimnames(a)[[3L]]))
      stop("Sampler output did not preserve the parameter layout.", call. = FALSE)
    if (!n) return(matrix(numeric(), prod(dim(a)[1:2]), 0L))
    matrix(a[, , columns, drop = FALSE], nrow = prod(dim(a)[1:2]))
  })
  names(p) <- names(object$parameters)
  values <- list()
  add <- function(x, quantity, covariate, group, row, units) {
    x <- as.matrix(x)
    ids <- paste0("derived[", nrow(metadata) + seq_len(ncol(x)), "]")
    metadata <<- rbind(metadata, data.frame(variable = ids, quantity = quantity,
      covariate = covariate, group = group, row = row, units = units))
    values[[length(values) + 1L]] <<- x
  }
  for (nm in object$specification$covariates) {
    sm <- object$specification$scales[[nm]]; rows <- sm$representatives
    x <- p$scale_beta[, sm$beta_index, drop = FALSE] %*% t(sm$X[rows, , drop = FALSE])
    for (term in sm$random)
      x <- x + p$scale_z[, term$z_index[term$index[rows]], drop = FALSE] *
        exp(p$log_scale_sd[, term$sd_index])
    add(x + log(object$specification$unit_conv), "log_distance", nm,
      sm$labels, rows, "log(map units)")
    for (group in names(sm$random))
      add(exp(p$log_scale_sd[, sm$random[[group]]$sd_index, drop = FALSE]),
        "scale_group_sd", nm, group, NA_integer_, "log distance")
  }
  for (group in names(object$specification$response$random)) {
    term <- object$specification$response$random[[group]]
    add(exp(p$log_response_sd[, term$sd_index, drop = FALSE]), "response_group_sd",
      NA_character_, group, NA_integer_, "response linear predictor")
  }
  x <- do.call(cbind, values)
  arr <- array(x, c(out$iterations, out$chains, ncol(x)),
    dimnames = list(NULL, NULL, metadata$variable))
  tab <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(arr),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail))
  out$parameters <- cbind(metadata, tab[match(metadata$variable, tab$variable),
    c("rhat", "ess_bulk", "ess_tail"), drop = FALSE])
  out$quantiles <- do.call(rbind, lapply(seq_len(ncol(x)), function(j) {
    chain <- matrix(x[, j], out$iterations, out$chains)
    estimate <- if (all(is.finite(chain))) unname(stats::quantile(chain, probs)) else
      rep(NA_real_, length(probs))
    cbind(metadata[rep(j, length(probs)), , drop = FALSE], probability = unname(probs),
      estimate = estimate, mcse = unname(posterior::mcse_quantile(chain, probs = probs)))
  }))
  rownames(out$parameters) <- rownames(out$quantiles) <- NULL
  d <- out$parameters
  out$derived_passed <- out$chains >= 2L &&
    all(is.finite(d$rhat) & d$rhat <= 1.01) &&
    all(is.finite(d$ess_bulk) & d$ess_bulk >= 100 * out$chains) &&
    all(is.finite(d$ess_tail) & d$ess_tail >= 100 * out$chains)
  raw_passed <- object$diagnostics$passed
  if (!is.null(object$chains) && !isTRUE(raw_passed) && !isFALSE(raw_passed))
    raw_passed <- .msr_bayes_chain_checks(.msr_bayes_chain_summary(a), out$chains)
  out$passed <- isTRUE(raw_passed) && out$derived_passed
  out$available <- TRUE
  out
}

.msr_bayes_fit_chains <- function(object) {
  if (!is.null(object$chains)) return(object$chains)
  .msr_bayes_sampler_array(object$sampler)
}

.msr_bayes_chain_summary <- function(a) {
  as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a)))
}

.msr_bayes_chain_checks <- function(tab, chains) {
  chains >= 2L && all(is.finite(tab$rhat) & tab$rhat <= 1.01) &&
    all(is.finite(tab$ess_bulk) & tab$ess_bulk >= 100 * chains) &&
    all(is.finite(tab$ess_tail) & tab$ess_tail >= 100 * chains)
}

.msr_bayes_chain_criteria <- function() {
  "At least 2 chains; Rhat <= 1.01; bulk/tail ESS >= 100 per chain."
}

.msr_bayes_sampler_array <- function(sampler) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Install rstan to read the retained NUTS sampler.", call. = FALSE)
  rstan::extract(sampler, permuted = FALSE, inc_warmup = FALSE)
}

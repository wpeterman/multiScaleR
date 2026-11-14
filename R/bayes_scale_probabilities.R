#' Check Monte Carlo Precision of Scale Probabilities
#'
#' Estimate the posterior probability that each fitted kernel scale lies at or
#' below thresholds you name, and report how precisely the original sampling
#' chains pin down each probability.
#'
#' @param object A fitted [multiScale_bayes()] object that keeps its original,
#'   unpermuted post-warmup chains from NUTS (`sampler`) or slice sampling
#'   (`chains`). Compact, thinned, and approximate fits return an unavailable
#'   report, because permuted draws cannot recover chain order.
#' @param thresholds A nonempty named list of numeric vectors. Each name
#'   selects a fitted landscape predictor, and each vector holds distinct
#'   finite scale thresholds. Every threshold is applied to every distinct
#'   fitted combination of scale covariates and groups for that predictor. For
#'   example, `list(habitat = c(100, 250))` checks two scales in projected map
#'   units.
#' @param scale `"distance"` (default) means `thresholds` are positive kernel
#'   scales in map units. `"log"` means they are natural logs of those
#'   scales. Thresholds refer to the scale parameters summarized by
#'   [posterior_scales()], not to transformed kernel-weight radii such as
#'   those from [kernel_dist()].
#' @param mcse_target Optional target for the Monte Carlo standard error
#'   (MCSE) on the probability axis, strictly between 0 and 1. For example,
#'   `0.005` asks for an MCSE no larger than half a percentage point. The
#'   default `NULL` reports precision without judging it against a target.
#'
#' @return A list with these elements.
#'   * `available`, `reason`: whether the original chains could be read, and
#'     an explanation when they could not. A single event can lack precision
#'     estimates even when the report is available.
#'   * `raw_diagnostics`: the fit's original diagnostics, unchanged.
#'   * `probabilities`: a data frame with one row per fitted scale combination
#'     and threshold. `covariate`, `group`, `row`, `threshold`, and
#'     `threshold_scale` identify the event. `estimate` is the posterior
#'     probability that the scale is at or below the threshold. `mcse`,
#'     `ess_mean`, and `rhat` are the Monte Carlo standard error, mean
#'     effective sample size, and rank-normalized Rhat of the event indicator.
#'     `constant` flags an indicator that never changed. `precision_available`
#'     is `TRUE` when the MCSE and mean ESS are finite and positive.
#'     `precision_target_met` compares `mcse` with `mcse_target`, and is `NA`
#'     when no target was given. `reason` explains missing precision.
#'     `rhat_available` and `rhat_reason` report separately whether Rhat could
#'     be computed. Undefined values are `NA`.
#'   * `chain_probabilities`: the same event identifiers with `chain`, that
#'     chain's `estimate`, and whether its indicator is `constant`. Large
#'     differences among chains signal poor mixing.
#'   * `chains`, `iterations`: the number of chains and of post-warmup
#'     iterations per chain, or `NA` when unavailable.
#'   * `mcse_target`, `event`, `notes`: the requested target, the event
#'     definition (`"kernel_scale <= threshold"`), and the limits of
#'     interpretation.
#'   * `metadata_establishes_eligibility`: always `FALSE`. The report does not
#'     change the fit's sampling or derived-scale diagnostic gates.
#'
#' @details
#' # Read the precision
#'
#' The MCSE is computational uncertainty in the probability estimate, not
#' posterior uncertainty about the scale and not observation variation.
#' `precision_target_met` only compares that MCSE with your target. It
#' does not override Rhat, ESS, the raw sampler checks, or model checks, and a
#' small MCSE cannot rule out an unvisited tail or mode. Read
#' [diagnose_bayes_scales()] alongside this report, and compare independent
#' runs when a consequential probability remains uncertain. Longer chains
#' reduce the MCSE.
#'
#' # How the events are computed
#'
#' Each event is evaluated on fitted log scales in map units, which avoids
#' overflow from exponentiating distances. The posterior package computes the
#' mean ESS, MCSE, and rank-normalized Rhat of each binary indicator. At least
#' two chains with four retained iterations each are needed for precision
#' diagnostics. A constant indicator has no estimable precision, even when its
#' probability is exactly zero or one. Indicator Rhat can be undefined for a
#' balanced binary event even when the MCSE is finite. It is then reported as
#' `NA` with its own flag and reason, which alone does not make the precision
#' unavailable. A nonfinite derived log scale is an error, and iterations are
#' never dropped.
#'
#' The function does not fit a model, change priors, simulate new groups, or
#' treat importance candidates as MCMC chains. Reading a full NUTS sampler can
#' take substantial memory, but the report keeps only summaries.
#' `vignette("v15-bayesian-computation", package = "multiScaleR")` shows how to
#' report a probability with its MCSE.
#'
#' @seealso [diagnose_bayes_scales()], [posterior_scales()],
#'   [predict_scales()], [diagnose_bayes_laplace()],
#'   [multiScaleR_bayes-class]
#' @examples
#' # A stored slice fit that keeps its four original chains
#' shared <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
#'   package = "multiScaleR"))$fit
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   report <- diagnose_bayes_probabilities(shared,
#'     thresholds = list(land1 = c(200, 300)), mcse_target = 0.01)
#'   # Probability that sigma for land1 is at most 200 m and at most 300 m
#'   print(report$probabilities[c("threshold", "estimate", "mcse",
#'                                "precision_target_met")])
#'   print(report$chain_probabilities[c("threshold", "chain", "estimate")])
#' }
#'
#' # A compact fit keeps thinned, permuted draws and cannot be checked
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' compact <- diagnose_bayes_probabilities(example$fits$pooled,
#'   thresholds = list(habitat = c(100, 250)))
#' compact$available   # FALSE
#' compact$reason
#' @export
diagnose_bayes_probabilities <- function(object, thresholds,
    scale = c("distance", "log"), mcse_target = NULL) {
  if (!inherits(object, "multiScaleR_bayes"))
    stop("object must be a multiScaleR_bayes fit.", call. = FALSE)
  scale <- match.arg(scale)
  .msr_bayes_validate_probability_thresholds(thresholds, object$specification, scale)
  if (!is.null(mcse_target)) {
    if (!is.null(dim(mcse_target)) || is.complex(mcse_target))
      stop("mcse_target must be a numeric scalar strictly between zero and one.", call. = FALSE)
    validate_scalar_numeric(mcse_target, "mcse_target", lower = 0,
      upper = 1, inclusive_lower = FALSE, inclusive_upper = FALSE)
  }
  identifiers <- data.frame(covariate = character(), group = character(), row = integer(),
    threshold = numeric(), threshold_scale = character())
  out <- list(available = FALSE, reason = "", raw_diagnostics = object$diagnostics,
    probabilities = cbind(identifiers, estimate = numeric(), mcse = numeric(), ess_mean = numeric(),
      rhat = numeric(), constant = logical(), precision_available = logical(),
      precision_target_met = logical(), reason = character(),
      rhat_available = logical(), rhat_reason = character()),
    chain_probabilities = cbind(identifiers, chain = integer(), estimate = numeric(), constant = logical()),
    chains = NA_integer_, iterations = NA_integer_, mcse_target = mcse_target,
    event = "kernel_scale <= threshold", metadata_establishes_eligibility = FALSE,
    notes = c("MCSE measures computational precision, not posterior or observation uncertainty.",
      "A target decision does not replace original Rhat, ESS or sampler gates.",
      "Constant indicators cannot certify an absent tail or an unvisited mode."))
  if (!is.null(object$laplace) || !is.null(object$marginal) || !is.null(object$importance)) {
    out$reason <- "MCMC probability diagnostics are unavailable for approximate draws or importance candidates."
    return(out)
  }
  if (is.null(object$sampler) && is.null(object$chains)) {
    out$reason <- "Full unpermuted MCMC chains are not retained; stored permuted or thinned draws cannot recover them."
    return(out)
  }
  if (.msr_has_gaussian_blocks(object$specification)) {
    out$reason <- "Private Gaussian projection blocks require a dedicated diagnostic adapter."
    return(out)
  }
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Install posterior to compute probability diagnostics.", call. = FALSE)
  a <- .msr_bayes_fit_chains(object)
  if (!is.numeric(a) || length(dim(a)) != 3L || any(dim(a) < 1L) ||
      is.null(dimnames(a)[[3L]]) || anyNA(dimnames(a)[[3L]]) ||
      any(!nzchar(dimnames(a)[[3L]])) || anyDuplicated(dimnames(a)[[3L]]))
    stop("Sampler output must have named iterations-by-chains-by-parameters dimensions.", call. = FALSE)
  out$iterations <- dim(a)[1L]; out$chains <- dim(a)[2L]
  p <- list()
  # Validate every original block, but copy only the blocks needed for scales.
  for (nm in names(object$parameters)) {
    n <- length(object$parameters[[nm]])
    bracketed <- paste0(nm, "[", seq_len(n), "]")
    if (n == 1L && nm %in% dimnames(a)[[3L]] && bracketed %in% dimnames(a)[[3L]])
      stop("Sampler output has an ambiguous original parameter layout.", call. = FALSE)
    columns <- if (n == 1L && nm %in% dimnames(a)[[3L]]) nm else bracketed
    if (n && !all(columns %in% dimnames(a)[[3L]]))
      stop("Sampler output did not preserve the original parameter layout.", call. = FALSE)
    if (nm %in% c("scale_beta", "scale_z", "log_scale_sd"))
      p[[nm]] <- if (!n) matrix(numeric(), prod(dim(a)[1:2]), 0L) else
        matrix(a[, , columns, drop = FALSE], nrow = prod(dim(a)[1:2]))
  }
  validate_scalar_numeric(object$specification$unit_conv, "fitted unit_conv", positive = TRUE)
  summaries <- list(); chains <- list()
  for (nm in names(thresholds)) {
    sm <- object$specification$scales[[nm]]
    for (j in seq_along(sm$representatives)) {
      row <- sm$representatives[j]
      log_scale <- drop(p$scale_beta[, sm$beta_index, drop = FALSE] %*% sm$X[row, ])
      for (term in sm$random)
        log_scale <- log_scale + exp(p$log_scale_sd[, term$sd_index]) *
          p$scale_z[, term$z_index[term$index[row]]]
      log_scale <- log_scale + log(object$specification$unit_conv)
      for (threshold in thresholds[[nm]]) {
        id <- data.frame(covariate = nm, group = sm$labels[j], row = as.integer(row),
          threshold = unname(threshold), threshold_scale = scale)
        event <- .msr_bayes_probability_event(log_scale,
          if (scale == "log") threshold else log(threshold), out$iterations, out$chains, mcse_target)
        summaries[[length(summaries) + 1L]] <- cbind(id, event$probabilities)
        chains[[length(chains) + 1L]] <- cbind(id[rep(1L, out$chains), , drop = FALSE],
          event$chain_probabilities)
      }
    }
  }
  out$probabilities <- do.call(rbind, summaries)
  out$chain_probabilities <- do.call(rbind, chains)
  rownames(out$probabilities) <- rownames(out$chain_probabilities) <- NULL
  out$available <- TRUE
  out
}

.msr_bayes_validate_probability_thresholds <- function(thresholds, specification, scale) {
  if (!is.list(thresholds) || !length(thresholds) || is.null(names(thresholds)) ||
      anyNA(names(thresholds)) || any(!nzchar(names(thresholds))) || anyDuplicated(names(thresholds)) ||
      !all(names(thresholds) %in% specification$covariates) ||
      any(vapply(thresholds, function(x) !is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x) ||
        any(!is.finite(x)) || anyDuplicated(x) > 0L || (scale == "distance" && any(x <= 0)), TRUE)))
    stop("thresholds must uniquely name fitted kernel predictors with distinct finite numeric thresholds; distance thresholds must be positive.", call. = FALSE)
}

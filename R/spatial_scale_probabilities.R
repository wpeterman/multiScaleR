#' Diagnose Private Spatial-Only Scale Probabilities
#'
#' Evaluate location-specific scale events on original retained NUTS chains
#' using the frozen spatial basis and the existing probability-precision rules.
#'
#' @param object A private \code{msr_spatial_scale_prototype} with a posterior
#'   target. Original unpermuted post-warmup NUTS chains must be retained in
#'   \code{sampler}. Stored permuted vectors cannot recover their chain order.
#' @param thresholds A nonempty numeric vector of distinct finite thresholds
#'   for the selected spatial kernel predictor. Distance thresholds are positive
#'   kernel scales in projected map units, not kernel-weight radii.
#' @param newdata Optional nonempty data frame of projected coordinates and
#'   required fixed scale covariates. NULL uses every original training row.
#' @param scale \code{"distance"} (default) supplies physical thresholds;
#'   \code{"log"} supplies their natural logarithms.
#' @param mcse_target Optional absolute Monte Carlo standard error target for
#'   the probability estimate, strictly between zero and one. NULL reports
#'   precision without a target decision.
#' @return A list with \code{available} and an unavailability \code{reason};
#'   unchanged \code{raw_diagnostics}; \code{chains} and \code{iterations};
#'   \code{mcse_target}, \code{event}, \code{notes}, original computation
#'   \code{method}, \code{coordinate_columns} and \code{chain_source}.
#'   The \code{spatial} element records the selected predictor, frozen basis,
#'   fixed range and absence of new-location innovation.
#'   \code{metadata_establishes_eligibility} is always FALSE.
#'   \code{probabilities} retains every input location and threshold, identified
#'   by \code{covariate}, input \code{row}, projected \code{x}/\code{y},
#'   \code{threshold} and \code{threshold_scale}. Its remaining fields are
#'   \code{estimate}, indicator \code{mcse}, \code{ess_mean}, \code{rhat},
#'   \code{constant}, \code{precision_available}, \code{precision_target_met},
#'   \code{reason}, \code{rhat_available} and \code{rhat_reason}, with the same
#'   meanings as \code{\link{diagnose_bayes_probabilities}}.
#'   \code{chain_probabilities} retains the same event identifiers plus
#'   \code{chain}, its empirical \code{estimate} and \code{constant} flag.
#'   Unavailable reports have empty event tables. The report retains no sampler
#'   or chain arrays and does not alter the object or random-number state.
#' @details Only the selected frozen spatial field and fixed scale designs
#'   are admitted; ordinary annual or other scale-group effects remain guarded.
#'   Every original parameter block is validated, including response parameters
#'   that do not enter ecological scales. Invalid or nonfinite original blocks
#'   cause an error. Nonfinite derived log scales retain explicit unavailable
#'   event precision, without dropping iterations or locations.
#'
#'   Comparisons use log scale to avoid overflow from exponentiating distances.
#'   Constant events retain descriptive probabilities zero or one but have
#'   unavailable precision, never zero MCSE. Mean precision requires at least
#'   two chains with four iterations, and positive finite mean MCSE and ESS.
#'   Indicator Rhat availability is reported separately. Approximation,
#'   importance and unretained-chain objects return unavailable reports.
#'
#'   MCSE measures computational precision, not ecological-scale posterior
#'   uncertainty or observation noise. A target decision does not override fit
#'   diagnostics or establish interval calibration, an absent tail or an
#'   unvisited mode. Range, prior normalization and training projection remain
#'   fixed; no annual effect or independent new-location innovation is added.
#' @seealso \code{\link{diagnose_bayes_probabilities}},
#'   \code{\link{.msr_spatial_scale_summary}}
#' @keywords internal
.msr_spatial_scale_probabilities <- function(object, thresholds, newdata = NULL,
    scale = c("distance", "log"), mcse_target = NULL) {
  if (!inherits(object, "msr_spatial_scale_prototype"))
    stop("Use a prepared private spatial scale prototype.", call. = FALSE)
  .msr_spatial_scale_objective(object)
  if (!identical(object$target, "posterior"))
    stop("Spatial MCMC probability diagnostics require a posterior target.", call. = FALSE)
  scale <- match.arg(scale)
  if (!is.numeric(thresholds) || is.complex(thresholds) || !is.null(dim(thresholds)) ||
      !length(thresholds) || any(!is.finite(thresholds)) || anyDuplicated(thresholds) ||
      (scale == "distance" && any(thresholds <= 0)))
    stop("thresholds must be distinct finite numeric values; distance thresholds must be positive.", call. = FALSE)
  if (!is.null(mcse_target)) {
    if (!is.null(dim(mcse_target)) || is.complex(mcse_target))
      stop("mcse_target must be a numeric scalar strictly between zero and one.", call. = FALSE)
    validate_scalar_numeric(mcse_target, "mcse_target", lower = 0, upper = 1,
      inclusive_lower = FALSE, inclusive_upper = FALSE)
  }
  nm <- object$prototype$covariate; spec <- object$specification
  term <- spec$scales[[nm]]$random[[".spatial_field"]]
  field <- term$gaussian_block
  allowed <- vapply(names(spec$scales), function(covariate) {
    effects <- spec$scales[[covariate]]$random
    if (identical(covariate, nm))
      identical(names(effects), ".spatial_field") && length(effects) == 1L else !length(effects)
  }, TRUE)
  if (is.null(field$basis) || !all(allowed))
    stop("Spatial-only probabilities require the frozen field and no ordinary scale-group effects.", call. = FALSE)
  if (is.null(newdata)) newdata <- spec$data
  if (!is.data.frame(newdata) || nrow(newdata) < 1L || anyDuplicated(names(newdata)))
    stop("newdata must be a nonempty data frame with distinct column names.", call. = FALSE)
  xy <- .msr_spatial_coordinates(newdata, field$coordinates)
  ids <- data.frame(covariate = character(), row = integer(), x = numeric(), y = numeric(),
    threshold = numeric(), threshold_scale = character())
  out <- list(available = FALSE, reason = "", raw_diagnostics = object$diagnostics,
    probabilities = cbind(ids, estimate = numeric(), mcse = numeric(), ess_mean = numeric(),
      rhat = numeric(), constant = logical(), precision_available = logical(),
      precision_target_met = logical(), reason = character(),
      rhat_available = logical(), rhat_reason = character()),
    chain_probabilities = cbind(ids, chain = integer(), estimate = numeric(), constant = logical()),
    chains = NA_integer_, iterations = NA_integer_, mcse_target = mcse_target,
    event = "kernel_scale <= threshold", method = object$inference,
    coordinate_columns = field$coordinates, chain_source = "original_unpermuted_post_warmup_NUTS",
    spatial = list(covariate = nm, basis = "frozen_training_projection_and_normalization",
      range = "fixed", new_location_innovation = FALSE),
    metadata_establishes_eligibility = FALSE,
    notes = c("MCSE measures computational precision, not posterior or observation uncertainty.",
      "A target decision does not replace original Rhat, ESS or sampler gates.",
      "Constant indicators cannot certify an absent tail or an unvisited mode."))
  if (!is.null(object$laplace) || !is.null(object$marginal) || !is.null(object$importance)) {
    out$reason <- "MCMC probability diagnostics are unavailable for approximate draws or importance candidates."
    return(out)
  }
  if (is.null(object$sampler)) {
    out$reason <- "Full unpermuted NUTS chains are not retained; stored permuted or thinned draws cannot recover them."
    return(out)
  }
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Install posterior to compute probability diagnostics.", call. = FALSE)
  a <- .msr_bayes_sampler_array(object$sampler)
  if (!is.numeric(a) || is.complex(a) || length(dim(a)) != 3L || any(dim(a) < 1L) ||
      is.null(dimnames(a)[[3L]]) || anyNA(dimnames(a)[[3L]]) ||
      any(!nzchar(dimnames(a)[[3L]])) || anyDuplicated(dimnames(a)[[3L]]))
    stop("Sampler output must have named iterations-by-chains-by-parameters dimensions.", call. = FALSE)
  out$iterations <- dim(a)[1L]; out$chains <- dim(a)[2L]
  p <- list()
  # Validate all original blocks, but retain only scale blocks for evaluation.
  for (block in names(object$parameters)) {
    n <- length(object$parameters[[block]])
    bracketed <- paste0(block, "[", seq_len(n), "]")
    if (n == 1L && block %in% dimnames(a)[[3L]] && bracketed %in% dimnames(a)[[3L]])
      stop("Sampler output has an ambiguous original parameter layout.", call. = FALSE)
    columns <- if (n == 1L && block %in% dimnames(a)[[3L]]) block else bracketed
    if (n && !all(columns %in% dimnames(a)[[3L]]))
      stop("Sampler output did not preserve the original parameter layout.", call. = FALSE)
    if (n && any(!is.finite(a[, , columns, drop = FALSE])))
      stop("Sampler output must preserve every complete finite original parameter block.", call. = FALSE)
    if (block %in% c("scale_beta", "scale_z", "log_scale_sd"))
      p[[block]] <- if (!n) matrix(numeric(), prod(dim(a)[1:2]), 0L) else
        matrix(a[, , columns, drop = FALSE], nrow = prod(dim(a)[1:2]))
  }
  # Reuse fixed-design validation and the original training projection, once.
  spec$scales[[nm]]$random[[".spatial_field"]] <- NULL
  context <- .msr_bayes_prediction_scales(spec, newdata, "error", TRUE)
  sm <- context$spec$scales[[nm]]
  design <- .msr_spatial_basis_predict(field$basis, xy, sm$X) %*% field$chol
  validate_scalar_numeric(spec$unit_conv, "fitted unit_conv", positive = TRUE)
  summaries <- chains <- list()
  # Keep only one location's derived iterations, rather than a draws-by-map array.
  for (row in seq_len(nrow(newdata))) {
    log_scale <- drop(p$scale_beta[, sm$beta_index, drop = FALSE] %*% sm$X[row, ]) +
      drop(p$scale_z[, term$z_index, drop = FALSE] %*% design[row, ]) *
        as.numeric(exp(p$log_scale_sd[, term$sd_index])) + log(spec$unit_conv)
    for (threshold in thresholds) {
      id <- data.frame(covariate = nm, row = as.integer(row), x = xy[row, 1L], y = xy[row, 2L],
        threshold = unname(threshold), threshold_scale = scale)
      event <- .msr_bayes_probability_event(log_scale,
        if (scale == "log") threshold else log(threshold), out$iterations, out$chains, mcse_target)
      summaries[[length(summaries) + 1L]] <- cbind(id, event$probabilities)
      chains[[length(chains) + 1L]] <- cbind(id[rep(1L, out$chains), , drop = FALSE],
        event$chain_probabilities)
    }
  }
  out$probabilities <- do.call(rbind, summaries)
  out$chain_probabilities <- do.call(rbind, chains)
  rownames(out$probabilities) <- rownames(out$chain_probabilities) <- NULL
  out$available <- TRUE
  out
}

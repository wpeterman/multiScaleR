#' Check Monte Carlo Precision of Scale Probabilities
#'
#' Estimate probabilities that fitted kernel scales fall below specified
#' thresholds, and diagnose their simulation precision using the original chains.
#'
#' @param object A \code{\link{multiScale_bayes}} fit retaining its original
#'   unpermuted, post-warmup NUTS sampler. Compact and approximate fits return
#'   an unavailable report; permuted draws cannot recover chain order.
#' @param thresholds A nonempty named list of numeric vectors. Each name selects
#'   a fitted kernel predictor, and each vector contains distinct finite scale
#'   thresholds. The thresholds apply to every distinct fitted scale-formula
#'   combination for that predictor. For example,
#'   \code{list(habitat = c(100, 250))} checks two scales in projected map units.
#' @param scale \code{"distance"} (default) supplies positive kernel-scale
#'   parameters in map units. \code{"log"} supplies their natural logarithms.
#'   These are the parameters summarized by \code{\link{posterior_scales}},
#'   not automatically transformed kernel-weight radii.
#' @param mcse_target Optional absolute Monte Carlo standard error target on the
#'   probability axis, strictly between zero and one. For example, \code{0.005}
#'   requests a probability MCSE no greater than half a percentage point.
#'   The default \code{NULL} reports precision without judging a target.
#' @return A list containing:
#'   \describe{
#'     \item{available, reason}{Whether original chains could be read, and an
#'       explanation when unavailable. Individual events can have unavailable
#'       precision even when the chain report is available.}
#'     \item{raw_diagnostics}{The fit's original diagnostics, unchanged.}
#'     \item{probabilities}{A data frame with \code{covariate}, \code{group}, and
#'       representative fitted \code{row}; the supplied \code{threshold} and
#'       \code{threshold_scale}; estimated event probability \code{estimate};
#'       indicator \code{mcse}, \code{ess_mean}, and \code{rhat}; \code{constant},
#'       \code{precision_available}, \code{precision_target_met}, and \code{reason}.
#'       The event includes equality: kernel scale is less than or equal to the
#'       threshold. Probability estimates and MCSE use the probability axis.
#'       Undefined diagnostics and unrequested target decisions remain \code{NA}.}
#'     \item{chain_probabilities}{The same event identifiers plus \code{chain},
#'       its empirical probability \code{estimate}, and whether its indicator
#'       is \code{constant}. This retains chain differences without storing draws.}
#'     \item{chains, iterations}{Numbers of original chains and post-warmup
#'       iterations per chain, or \code{NA} when unavailable.}
#'     \item{mcse_target, event, notes}{The requested precision target, event
#'       definition, and interpretation limits.}
#'     \item{metadata_establishes_eligibility}{Always \code{FALSE}. This report
#'       does not change the fit's sampling or derived-scale diagnostic gates.}
#'   }
#' @details Each event is evaluated on fitted log scales in map units, avoiding
#'   overflow from exponentiating distances. The \code{posterior} package
#'   computes mean ESS, MCSE and rank-normalized Rhat of the binary indicators.
#'   At least two chains with four retained iterations each are required for
#'   precision diagnostics. Constant indicators have nonestimable precision,
#'   even when their empirical probability is zero or one. Nonfinite derived
#'   log scales remain explicit failures; iterations are never dropped.
#'
#'   MCSE describes computational uncertainty in a probability estimate, not
#'   posterior uncertainty about the scale or observation sampling variation.
#'   \code{precision_target_met} only compares that MCSE with the requested
#'   target. It does not override Rhat, ESS, raw sampler diagnostics, or model
#'   checks. A small estimated MCSE cannot rule out an unvisited tail or mode.
#'   Read \code{\link{diagnose_bayes_scales}} alongside this report and compare
#'   independent runs when consequential probabilities remain uncertain.
#'
#'   The function does not fit a model, change priors, simulate future groups,
#'   or use stored importance candidates as MCMC chains. Private Gaussian
#'   projection blocks require a dedicated diagnostic adapter. Reading a full
#'   retained sampler can require substantial memory; the returned report keeps
#'   only summaries and retains neither a sampler nor chain arrays.
#' @seealso \code{\link{diagnose_bayes_scales}},
#'   \code{\link{posterior_scales}}, \code{\link{predict_scales}}
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' report <- diagnose_bayes_probabilities(example$fits$pooled,
#'   thresholds = list(habitat = c(100, 250)), mcse_target = 0.005)
#' report$available  # FALSE: this compact example retains no original sampler
#' report$reason
#' @export
diagnose_bayes_probabilities <- function(object, thresholds,
    scale = c("distance", "log"), mcse_target = NULL) {
  if (!inherits(object, "multiScaleR_bayes"))
    stop("object must be a multiScaleR_bayes fit.", call. = FALSE)
  scale <- match.arg(scale)
  if (!is.list(thresholds) || !length(thresholds) || is.null(names(thresholds)) ||
      anyNA(names(thresholds)) || any(!nzchar(names(thresholds))) || anyDuplicated(names(thresholds)) ||
      !all(names(thresholds) %in% object$specification$covariates) ||
      any(vapply(thresholds, function(x) !is.numeric(x) || is.complex(x) || !is.null(dim(x)) || !length(x) ||
        any(!is.finite(x)) || anyDuplicated(x) > 0L || (scale == "distance" && any(x <= 0)), TRUE)))
    stop("thresholds must uniquely name fitted kernel predictors with distinct finite numeric thresholds; distance thresholds must be positive.", call. = FALSE)
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
      precision_target_met = logical(), reason = character()),
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
  if (is.null(object$sampler)) {
    out$reason <- "Full unpermuted NUTS chains are not retained; stored permuted or thinned draws cannot recover them."
    return(out)
  }
  if (.msr_has_gaussian_blocks(object$specification)) {
    out$reason <- "Private Gaussian projection blocks require a dedicated diagnostic adapter."
    return(out)
  }
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Install posterior to compute probability diagnostics.", call. = FALSE)
  a <- .msr_bayes_sampler_array(object$sampler)
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
        metrics <- c(mcse = NA_real_, ess_mean = NA_real_, rhat = NA_real_)
        finite <- all(is.finite(log_scale)); reason <- ""
        indicator <- matrix(NA_real_, out$iterations, out$chains)
        if (finite) indicator[] <- as.numeric(log_scale <= if (scale == "log") threshold else log(threshold))
        constant <- if (finite) length(unique(as.vector(indicator))) == 1L else NA
        if (!finite) reason <- "Nonfinite derived log scales; no iterations were dropped." else
          if (constant) reason <- "Constant indicator; probability precision is nonestimable." else
            if (out$chains < 2L || out$iterations < 4L)
              reason <- "At least two chains with four iterations each are required." else {
                metrics <- tryCatch(c(mcse = posterior::mcse_mean(indicator),
                  ess_mean = posterior::ess_mean(indicator), rhat = posterior::rhat(indicator)),
                  error = function(e) { reason <<- conditionMessage(e); metrics })
                if (any(!is.finite(metrics)) || metrics[["mcse"]] <= 0 || metrics[["ess_mean"]] <= 0)
                  if (!nzchar(reason)) reason <- "Indicator precision diagnostics are nonestimable; inspect chain mixing."
              }
        usable <- !nzchar(reason)
        summaries[[length(summaries) + 1L]] <- cbind(id,
          estimate = if (finite) mean(indicator) else NA_real_,
          mcse = unname(metrics[["mcse"]]), ess_mean = unname(metrics[["ess_mean"]]),
          rhat = unname(metrics[["rhat"]]), constant = constant,
          precision_available = usable, precision_target_met = if (usable && !is.null(mcse_target))
            metrics[["mcse"]] <= mcse_target else NA, reason = reason)
        chains[[length(chains) + 1L]] <- cbind(id[rep(1L, out$chains), , drop = FALSE],
          chain = seq_len(out$chains), estimate = if (finite) colMeans(indicator) else rep(NA_real_, out$chains),
          constant = if (finite) apply(indicator, 2L, function(x) length(unique(x)) == 1L) else
            rep(NA, out$chains))
      }
    }
  }
  out$probabilities <- do.call(rbind, summaries)
  out$chain_probabilities <- do.call(rbind, chains)
  rownames(out$probabilities) <- rownames(out$chain_probabilities) <- NULL
  out$available <- TRUE
  out
}

#' Check Sampling Precision of Fitted Bayesian Scales
#'
#' Check chain mixing and simulation precision for fitted log distances and
#' between-group SDs. These derived quantities can mix differently from the
#' individual coefficients checked by the original sampler diagnostics.
#'
#' @param object A \code{\link{multiScale_bayes}} fit. Computing new diagnostics
#'   requires its original NUTS sampler with unpermuted post-warmup chains.
#'   Prepared, approximate, and compact fits return an unavailable report.
#' @param probs Numeric vector of distinct probabilities strictly between zero
#'   and one. Defaults to \code{c(0.025, 0.5, 0.975)}, the endpoints and median
#'   of a central 95\% interval. Changing this vector changes the quantiles
#'   whose simulation error is reported, not any diagnostic threshold.
#' @return A list with the following components:
#'   \describe{
#'     \item{available, reason}{Whether full chain diagnostics could be computed
#'       and an explanation when they are unavailable.}
#'     \item{raw_diagnostics}{The original, unchanged sampler or approximation
#'       diagnostics stored in the fit.}
#'     \item{parameters}{A data frame identifying each derived \code{variable},
#'       \code{quantity}, landscape \code{covariate}, \code{group}, representative
#'       observation \code{row}, and \code{units}, followed by \code{rhat},
#'       \code{ess_bulk}, and \code{ess_tail}. Quantities are \code{log_distance},
#'       \code{scale_group_sd}, and \code{response_group_sd}.}
#'     \item{quantiles}{A data frame with the same identifying columns plus
#'       \code{probability}, \code{estimate}, and \code{mcse}, the estimated Monte
#'       Carlo standard error of that quantile. Both estimate and MCSE use the
#'       stated units. Undefined diagnostics remain \code{NA}.}
#'     \item{chains, iterations}{Actual numbers of retained chains and
#'       post-warmup iterations per chain, or \code{NA} when unavailable.}
#'     \item{derived_passed}{Whether at least two chains, Rhat no greater than
#'       1.01, and bulk and tail effective sample sizes of at least 100 per
#'       chain hold for every derived quantity. \code{NA} when unavailable.}
#'     \item{passed}{Whether both the original sampler checks and the derived
#'       checks passed. \code{NA} when derived diagnostics are unavailable.
#'       This does not replace or change \code{object$diagnostics$passed}.}
#'     \item{criteria, notes}{The fixed derived checks and their limitations.}
#'   }
#' @details Landscape distances are diagnosed on the natural log axis in map
#'   units, avoiding overflow and the potentially infinite arithmetic distance
#'   moments permitted by the hierarchical priors. Exponentiate an estimated
#'   log-distance quantile to express it in map units. Its MCSE is on the log
#'   axis; it is computational uncertainty, not a posterior interval or a
#'   biological SD. Scale-group SDs measure population variation in log
#'   distance. Response-group SDs use the response linear-predictor axis.
#'
#'   The function uses \code{posterior} for rank-normalized Rhat, effective
#'   sample sizes, and quantile MCSE. A passed screen does not guarantee a
#'   precise tail endpoint or discovery of every mode. Compare MCSE with the
#'   precision your scientific decision needs, inspect independent runs when
#'   tails matter, and assess the likelihood and priors separately.
#'
#'   Permuted or thinned stored draws do not reconstruct the original chain
#'   order or full-chain precision. Compact examples therefore retain their
#'   original diagnostic evidence but cannot provide new chain diagnostics.
#'   MCMC diagnostics are not defined for Laplace draws or weighted importance
#'   candidates; inspect those fits' original numerical checks instead.
#' @seealso \code{\link{multiScale_bayes}}, \code{\link{posterior_scales}}
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'                               package = "multiScaleR"))
#' report <- diagnose_bayes_scales(example$fits$pooled)
#' report$available  # FALSE: the compact example contains no full sampler
#' report$reason
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
  if (is.null(object$sampler)) {
    out$reason <- paste("Full unpermuted NUTS chains are not retained.",
      "Stored permuted or thinned draws cannot recover original chain diagnostics.")
    return(out)
  }
  if (!requireNamespace("posterior", quietly = TRUE))
    stop("Install posterior to compute derived sampling diagnostics.", call. = FALSE)
  a <- .msr_bayes_sampler_array(object$sampler)
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
  out$passed <- isTRUE(object$diagnostics$passed) && out$derived_passed
  out$available <- TRUE
  out
}

.msr_bayes_sampler_array <- function(sampler) {
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("Install rstan to read the retained NUTS sampler.", call. = FALSE)
  rstan::extract(sampler, permuted = FALSE, inc_warmup = FALSE)
}

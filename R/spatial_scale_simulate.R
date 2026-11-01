# Complete stored parameter blocks are shared by both private adapters.
.msr_spatial_response_joint_ids <- function(object, draws = NULL) {
  if (!identical(object$target, "posterior"))
    stop("Stored joint predictions require the posterior target.", call. = FALSE)
  ids <- .msr_bayes_draw_ids(object, draws)
  n <- nrow(object$draws$beta)
  if (!length(ids) || !all(vapply(names(object$specification$parameters), function(nm) {
    x <- object$draws[[nm]]
    is.matrix(x) && nrow(x) == n && ncol(x) == length(object$specification$parameters[[nm]]) &&
      all(is.finite(x))
  }, TRUE))) stop("Retain complete finite joint posterior parameter blocks.", call. = FALSE)
  ids
}

#' Simulate Responses from a Private Spatial Scale Prototype
#'
#' Generate response replicates from complete stored joint vectors and the
#' frozen spatial scale field, including family-specific observation noise.
#'
#' @inheritParams .msr_spatial_response_predict
#' @param nsim Positive integer number of replicate datasets. The default is one.
#' @param trials Optional positive integer binomial totals, one scalar or one
#'   per observation. NULL retains fitted totals at training rows and uses one
#'   trial at new rows from a binary fit. New grouped-binomial rows require
#'   explicit totals. This argument applies only to binomial models.
#' @return A data frame with observation rows and one column per replicate.
#'   Gaussian responses are in response units; binomial responses are success
#'   counts out of their trial totals; Poisson and negative-binomial responses
#'   are counts.
#'   Attributes retain point IDs, unseen-group policies and the prediction
#'   contract. Replicates include observation noise, unlike expected means.
#' @details Each replicate selects one complete stored joint vector. Ordinary
#'   unseen-group effects are shared within that replicate, and observation
#'   noise is generated separately for each row. The field uses the frozen
#'   finite basis without an independent new-location innovation. Zeroing
#'   scale effects gives a zero-deviation reference, not field integration.
#'   No year column is needed when the prepared model has no temporal terms.
#'   Original sampler/approximation labels and failed diagnostics are retained.
#'   These outputs do not establish calibration or diagnostic eligibility.
#' @seealso \code{\link{.msr_spatial_response_predict}}
#' @keywords internal
.msr_spatial_response_simulate <- function(object, nsim = 1L, seed = NULL,
    newdata = NULL, kernel_inputs = NULL, row_id = NULL, re.form = NULL,
    scale_re.form = NULL, new_scale_levels = c("error", "sample", "zero"),
    new_response_levels = c("error", "sample", "zero"), trials = NULL,
    include_fixed = FALSE) {
  validate_scalar_numeric(nsim, "nsim", integerish = TRUE, lower = 1)
  ids <- .msr_spatial_response_joint_ids(object)
  include_response <- .msr_bayes_prediction_re_form(re.form, "re.form")
  include_scale <- .msr_bayes_prediction_re_form(scale_re.form, "scale_re.form")
  context <- .msr_bayes_prediction_seed(seed,
    .msr_spatial_response_context(object, newdata, kernel_inputs, row_id,
      match.arg(new_scale_levels), match.arg(new_response_levels),
      include_scale, include_response, include_fixed))
  spec <- context$spec
  if (!is.null(trials)) {
    if (spec$family != "binomial") stop("trials applies only to binomial models.", call. = FALSE)
    if (!is.numeric(trials) || !is.null(dim(trials)) ||
        !length(trials) %in% c(1L, length(spec$y)) ||
        any(!is.finite(trials) | trials < 1 | trials != floor(trials)))
      stop("trials must be positive integer totals, one or one per prediction row.", call. = FALSE)
    spec$trials <- rep(trials, length.out = length(spec$y))
  } else if (!is.null(newdata) && spec$family == "binomial" &&
      any(object$specification$trials != 1))
    stop("Supply trials explicitly for new grouped-binomial observations.", call. = FALSE)
  ans <- .msr_bayes_prediction_seed(seed, {
    values <- matrix(NA_real_, length(spec$y), nsim)
    selected <- sample(ids, nsim, replace = TRUE)
    for (j in seq_len(nsim)) {
      p <- .msr_bayes_prediction_parameters(.msr_bayes_draw(object, selected[j]), context)
      eta <- .msr_bayes_components(spec, p, include_response)$eta
      values[, j] <- switch(spec$family,
        gaussian = stats::rnorm(length(eta), eta, exp(p$log_dispersion)),
        poisson = stats::rpois(length(eta), exp(eta)),
        binomial = stats::rbinom(length(eta), spec$trials, stats::plogis(eta)),
        nbinom2 = stats::rnbinom(length(eta), mu = exp(eta), size = exp(p$log_dispersion)))
      if (any(!is.finite(values[, j])))
        stop("Nonfinite spatial response replicate; inspect inputs and joint vectors.", call. = FALSE)
    }
    values
  })
  colnames(ans) <- paste0("sim_", seq_len(nsim)); ans <- as.data.frame(ans)
  attr(ans, "point_ids") <- spec$observation_ids
  attr(ans, "new_groups") <- context$new_groups
  contract <- .msr_bayes_prediction_contract(object, context, "replicate_response", "response",
    include_scale, include_response)
  contract$effects$spatial_scale <- context$spatial_field
  attr(ans, "prediction_contract") <- contract
  ans
}

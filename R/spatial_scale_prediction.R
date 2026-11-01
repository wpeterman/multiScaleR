# Prepare the ordinary response design before restoring the one frozen field.
# Other Gaussian blocks still reach the shared validator and remain rejected.
.msr_spatial_response_context <- function(object, newdata, kernel_inputs, row_id,
    scale_policy, response_policy, include_scale, include_response, include_fixed) {
  if (!inherits(object, "msr_spatial_scale_prototype"))
    stop("Use a prepared private spatial scale prototype.", call. = FALSE)
  # Validate the target/prior identity without evaluating or refitting it.
  .msr_spatial_scale_objective(object)
  nm <- object$prototype$covariate
  term <- object$specification$scales[[nm]]$random[[".spatial_field"]]
  if (is.null(term$gaussian_block$basis))
    stop("The spatial prototype needs its original frozen field basis.", call. = FALSE)
  ordinary <- object
  ordinary$specification$scales[[nm]]$random[[".spatial_field"]] <- NULL
  context <- .msr_bayes_prediction_prepare(ordinary, newdata, kernel_inputs, row_id,
    scale_policy, response_policy, include_scale, include_response, include_fixed)
  if (include_scale) {
    if (!is.null(newdata)) {
      xy <- .msr_spatial_coordinates(newdata, term$gaussian_block$coordinates)
      term$gaussian_block$design <- .msr_spatial_basis_predict(
        term$gaussian_block$basis, xy, context$spec$scales[[nm]]$X)
    }
    context$spec$scales[[nm]]$random[[".spatial_field"]] <- term
  }
  context$spatial_field <- list(covariate = nm, included = include_scale,
    coordinates = term$gaussian_block$coordinates,
    basis = "frozen_training_projection_and_normalization", range = "fixed",
    new_location_innovation = FALSE)
  context
}

#' Predict Responses from a Private Spatial Scale Prototype
#'
#' Evaluate a response mean or linear predictor with a frozen spatial scale
#' field. This experimental adapter uses the shared Bayesian input validation
#' without admitting spatial fields to public prediction or fitting methods.
#'
#' @param object A private \code{msr_spatial_scale_prototype} object. Stored
#'   posterior draws are required unless \code{parameters} is supplied.
#' @param newdata Optional observation data. NULL uses the training rows. New
#'   rows need projected coordinate columns when the spatial field is included.
#' @param kernel_inputs Matching prepared binned kernel means or raw weighted-SQ inputs for
#'   \code{newdata}. NULL is appropriate only for training predictions.
#' @param row_id Optional point-ID column in \code{newdata}. NULL uses the
#'   training ID column. Repeated IDs retain every observation row.
#' @param parameters Optional named list in the original parameter layout.
#'   This returns a point evaluation, with no parameter uncertainty. NULL uses
#'   complete stored joint posterior vectors, including stored joint resamples.
#' @param type Output scale. \code{"response"} returns the expected Gaussian
#'   response, binomial probability, or Poisson/negative-binomial mean.
#'   \code{"link"} returns the linear predictor before the inverse link.
#' @param scale_re.form,re.form NULL includes fitted scale and response effects,
#'   respectively. NA sets their deviations to zero. Zeroing the scale effects
#'   also omits the spatial field; it does not integrate over that field.
#' @param new_scale_levels,new_response_levels Policy for unseen ordinary
#'   groups: \code{"error"}, \code{"sample"}, or \code{"zero"}. Sampling
#'   shares one effect per group within each complete parameter vector.
#' @param draws Optional positive number of evenly spaced stored vectors to use.
#'   NULL uses all vectors. Do not combine this with \code{parameters}.
#' @param seed Optional positive integer seed for sampled new groups. Supplying
#'   it restores the caller's random-number state, including after an error.
#' @param include_fixed Logical. FALSE uses supplied fixed model columns. TRUE
#'   reevaluates fixed landscape features using their fitted reference scales.
#' @return A numeric vector for explicit parameters, or a draws-by-observations
#'   matrix for stored joint vectors, on the requested output scale. The
#'   \code{new_groups} attribute describes unseen ordinary groups. The
#'   \code{prediction_contract} attribute records the target, effects, fixed
#'   spatial basis, omitted observation noise and original diagnostic status.
#'   It never establishes fit eligibility or interval calibration.
#' @details Predictions keep the original training center, spread, spatial
#'   projection, range and prior normalization. New coordinates use the same
#'   finite basis without an independent field innovation. No year column is
#'   needed when the prepared model has no temporal terms. Expected responses
#'   omit observation noise; these are not posterior predictive replicates.
#'   Stored Gaussian approximation vectors remain approximations and stored
#'   sampler vectors retain their original diagnostics. Neither this adapter
#'   nor its metadata qualifies failed diagnostics or extreme-tail precision.
#' @seealso \code{\link{.msr_spatial_scale_prepare}},
#'   \code{\link{prediction_contract}}
#' @keywords internal
.msr_spatial_response_predict <- function(object, newdata = NULL,
    kernel_inputs = NULL, row_id = NULL, parameters = NULL,
    type = c("response", "link"), scale_re.form = NULL, re.form = NULL,
    new_scale_levels = c("error", "sample", "zero"),
    new_response_levels = c("error", "sample", "zero"), draws = NULL,
    seed = NULL, include_fixed = FALSE) {
  type <- match.arg(type)
  include_scale <- .msr_bayes_prediction_re_form(scale_re.form, "scale_re.form")
  include_response <- .msr_bayes_prediction_re_form(re.form, "re.form")
  if (!is.null(parameters) && !is.null(draws))
    stop("Use parameters or stored draws, not both.", call. = FALSE)
  context <- .msr_bayes_prediction_seed(seed,
    .msr_spatial_response_context(object, newdata, kernel_inputs, row_id,
      match.arg(new_scale_levels), match.arg(new_response_levels),
      include_scale, include_response, include_fixed))
  spec <- context$spec
  ids <- if (is.null(parameters)) .msr_spatial_response_joint_ids(object, draws) else 1L
  ans <- .msr_bayes_prediction_seed(seed, {
    values <- matrix(NA_real_, length(ids), length(spec$y))
    for (i in seq_along(ids)) {
      p <- if (is.null(parameters)) .msr_bayes_draw(object, ids[i]) else parameters
      valid <- is.list(p) && all(vapply(names(spec$parameters), function(nm) {
        x <- p[[nm]]
        is.numeric(x) && length(x) == length(spec$parameters[[nm]]) &&
          all(is.finite(x) | (nm %in% c("log_scale_sd", "log_response_sd") & x == -Inf))
      }, TRUE))
      if (!isTRUE(valid)) stop("Supply the original numeric parameter blocks.", call. = FALSE)
      p <- .msr_bayes_prediction_parameters(p, context)
      eta <- .msr_bayes_components(spec, p, include_response)$eta
      values[i, ] <- if (type == "link") eta else switch(spec$family,
        gaussian = eta, binomial = stats::plogis(eta), exp(eta))
      if (any(!is.finite(values[i, ])))
        stop("Nonfinite spatial response prediction; inspect inputs and parameter vectors.", call. = FALSE)
    }
    values
  })
  quantity <- if (type == "link") "linear_predictor" else "expected_response"
  if (is.null(parameters)) {
    contract <- .msr_bayes_prediction_contract(object, context, quantity, type,
      include_scale, include_response)
    # A raw matrix provides joint vectors, not an interval calculation.
    contract$uncertainty$intervals <- "not_provided"
  } else {
    contract <- .msr_prediction_contract(
      if (object$target == "posterior") "posterior" else "likelihood",
      "supplied_parameter_evaluation", quantity, type,
      population = "fixed_parameters", kernel = "fixed_parameters",
      conditioning = c("supplied_parameter_vector", "frozen_training_design_and_normalization"),
      effects = list(scale = list(included = include_scale, unseen = context$scale_policy),
        response = list(included = include_response, unseen = context$response_policy)),
      omitted = c("parameter_uncertainty", "observation_noise"),
      coordinates = "original_parameter_blocks", joint_map = TRUE,
      diagnostics = object$diagnostics)
    ans <- as.numeric(ans[1L, ])
  }
  contract$effects$spatial_scale <- context$spatial_field
  attr(ans, "new_groups") <- context$new_groups
  attr(ans, "prediction_contract") <- contract
  ans
}

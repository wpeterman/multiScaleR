#' Compare Private Spatial Scales at Paired Locations
#'
#' Evaluate ratios of ecological kernel scales on the same complete joint
#' vector at both locations, preserving their shared spatial dependence.
#'
#' @inheritParams .msr_spatial_scale_summary
#' @param pairs A nonempty numeric matrix with two columns of positive integer
#'   location-row indices. Column one is the denominator (from); column two is
#'   the numerator (to). Rows refer to newdata, or training data when it is NULL.
#'   Input order, repeated pairs and optional row labels are retained.
#' @param scale \code{"log_ratio"} (default) returns log(scale_to/scale_from).
#'   \code{"ratio"} returns scale_to/scale_from. Both are dimensionless kernel
#'   scale comparisons, not ratios of automatically transformed weight radii.
#' @return With \code{summary = TRUE}, a data frame with one row per pair:
#'   \code{pair} is input pair order, \code{pair_name} retains its supplied row
#'   label or numeric position, \code{covariate} names the selected predictor,
#'   \code{from_row}/\code{to_row} are original location-row indices, and
#'   \code{from_x}/\code{from_y}/\code{to_x}/\code{to_y} retain projected
#'   coordinates. \code{mean}, \code{median}, \code{sd}, \code{lower} and
#'   \code{upper} describe the joint contrast on the requested output scale.
#'   SD is posterior variation of that contrast, not the field-SD parameter.
#'   Otherwise, a joint-vectors-by-pairs matrix with pair labels as column names.
#'   Attributes retain the original \code{location_pairs} matrix, the identified
#'   \code{contrast_pairs} table, \code{coordinate_columns}, and evaluated
#'   \code{prediction_locations} with their original row indices.
#'   The \code{prediction_contract} identifies the scale-contrast target and
#'   its uncertainty. Summaries additionally retain the interval \code{level}.
#' @details All supplied coordinates are validated; only locations referenced
#'   by pairs require derived-vector storage. The original frozen basis, range,
#'   training projection and prior normalization remain unchanged. Ordinary
#'   scale-group effects remain guarded. Response parameters are validated as
#'   part of complete joint vectors but do not enter ecological scales.
#'
#'   Each pair uses the same stored vector at both locations. Equal locations
#'   therefore give log ratio zero or ratio one, even when each scale separately
#'   has substantial posterior variation. This algebraic result is not a test
#'   of field constancy or a Monte Carlo precision estimate. Summaries are
#'   computed from contrast vectors on the requested scale; finite-sample
#'   interpolated summaries need not equal exponentiated log-scale summaries.
#'
#'   Ratios that underflow to zero or overflow cause an error. Use log ratios
#'   instead; no vectors or pairs are clipped or dropped. Raw matrices have no
#'   computed intervals. Pointwise contrast intervals are not obtained by
#'   subtracting marginal interval endpoints. Observation noise is omitted.
#'   Original method and failed/unavailable diagnostic labels are preserved.
#'   No new-location field innovation, MCSE, calibration, fitting admission or
#'   model eligibility is introduced.
#' @seealso \code{\link{.msr_spatial_scale_summary}},
#'   \code{\link{prediction_contract}}
#' @keywords internal
.msr_spatial_scale_contrasts <- function(object, pairs, newdata = NULL,
    scale = c("log_ratio", "ratio"), level = .95, draws = NULL,
    summary = TRUE, seed = NULL) {
  if (!inherits(object, "msr_spatial_scale_prototype"))
    stop("Use a prepared private spatial scale prototype.", call. = FALSE)
  .msr_spatial_scale_objective(object)
  scale <- match.arg(scale)
  validate_scalar_logical(summary, "summary")
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
    inclusive_lower = FALSE, inclusive_upper = FALSE)
  data <- if (is.null(newdata)) object$specification$data else newdata
  if (!is.data.frame(data) || nrow(data) < 1L || anyDuplicated(names(data)))
    stop("newdata must be a nonempty data frame with distinct column names.", call. = FALSE)
  if (!is.matrix(pairs) || !is.numeric(pairs) || is.complex(pairs) ||
      ncol(pairs) != 2L || nrow(pairs) < 1L || any(!is.finite(pairs)) ||
      any(pairs != floor(pairs) | pairs < 1 | pairs > nrow(data)))
    stop("pairs must have two columns of positive integer indices within the supplied location rows.", call. = FALSE)
  nm <- object$prototype$covariate
  field <- object$specification$scales[[nm]]$random[[".spatial_field"]]$gaussian_block
  .msr_spatial_coordinates(data, field$coordinates)
  # Evaluate each referenced location once, without a full draws-by-map array.
  rows <- unique(as.integer(pairs))
  from <- match(pairs[, 1L], rows); to <- match(pairs[, 2L], rows)
  logs <- .msr_spatial_scale_summary(object, data[rows, , drop = FALSE], level,
    draws, scale = "log", summary = FALSE, seed = seed)
  values <- logs[, to, drop = FALSE] - logs[, from, drop = FALSE]
  if (any(!is.finite(values)))
    stop("Nonfinite spatial log contrast; no vectors or pairs were dropped.", call. = FALSE)
  if (scale == "ratio") {
    values <- exp(values)
    if (any(!is.finite(values) | values <= 0))
      stop("Spatial ratios exceed numeric range; use scale = 'log_ratio'.", call. = FALSE)
  }
  labels <- if (is.null(rownames(pairs))) as.character(seq_len(nrow(pairs))) else rownames(pairs)
  colnames(values) <- labels
  locations <- attr(logs, "prediction_locations"); locations$row <- rows
  identified <- data.frame(pair = seq_len(nrow(pairs)), pair_name = labels, covariate = nm,
    from_row = as.integer(pairs[, 1L]), to_row = as.integer(pairs[, 2L]),
    from_x = locations$x[from], from_y = locations$y[from],
    to_x = locations$x[to], to_y = locations$y[to])
  ans <- if (summary) cbind(identified, .msr_bayes_intervals(values, level)) else values
  context <- list(n_scale = 0L, n_response = 0L, scale_policy = "error", response_policy = "error")
  contract <- .msr_bayes_prediction_contract(object, context, "kernel_scale_contrast", scale)
  contract$target_population <- "specified_location_pairs"
  contract$effects$spatial_scale <- attr(logs, "prediction_contract")$effects$spatial_scale
  contract$conditioning <- c(contract$conditioning, "same_joint_vector_at_both_locations")
  contract$component <- list(name = "paired_location_kernel_scales", direction = "to_over_from")
  contract$uncertainty$intervals <- if (summary) "conditional_posterior_contrast_quantiles" else "not_provided"
  attr(ans, "location_pairs") <- pairs
  attr(ans, "contrast_pairs") <- identified
  attr(ans, "coordinate_columns") <- field$coordinates
  attr(ans, "prediction_locations") <- locations
  attr(ans, "prediction_contract") <- contract
  if (summary) attr(ans, "level") <- level
  ans
}

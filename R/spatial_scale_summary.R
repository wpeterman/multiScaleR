#' Summarize a Private Spatial-Only Scale Surface
#'
#' Report the selected ecological kernel scale at each supplied location using
#' complete stored joint vectors and the original frozen spatial basis.
#'
#' @param object A private \code{msr_spatial_scale_prototype} with complete
#'   posterior vectors. Ordinary scale-group effects are not supported here.
#' @param newdata Optional data frame of projected coordinates and any other
#'   required scale covariates. NULL uses every original training row.
#' @param level Central pointwise interval probability, strictly between zero
#'   and one. The default is 0.95. These are marginal location intervals.
#' @param draws Optional positive number of evenly spaced stored vectors.
#'   NULL uses all vectors. This does not repair failed diagnostics.
#' @param scale \code{"distance"} returns physical map units; \code{"log"}
#'   returns their natural logarithm. Training normalization stays fixed.
#' @param summary Logical. TRUE returns pointwise summaries. FALSE returns the
#'   complete draws-by-locations matrix for joint calculations.
#' @param seed Optional positive integer seed. Supplying it restores the
#'   caller's random-number state, including after an error.
#' @return With \code{summary = TRUE}, a data frame retaining every location row:
#'   \code{row} is input order, \code{covariate} names the selected kernel,
#'   \code{x} and \code{y} are the supplied projected coordinates, and
#'   \code{mean}, \code{median}, \code{sd}, \code{lower}, \code{upper} describe
#'   the ecological scale on the requested output scale. The SD is posterior
#'   variation in ecological scale at that location, not the field-SD parameter.
#'   \code{p_below_min_D} and \code{p_above_max_D} are fractions of retained
#'   vectors outside the original reference distances, not precision estimates.
#'   Otherwise, a joint draws-by-locations matrix. Attributes
#'   \code{coordinate_columns}, \code{prediction_locations} and
#'   \code{prediction_contract} retain the original coordinate names, row
#'   coordinates and uncertainty description, respectively.
#'   Summaries additionally retain their interval \code{level}.
#' @details This private adapter admits only the spatial field and fixed scale
#'   designs, with no ordinary annual or other scale-group effects anywhere in
#'   the model. Response effects do not enter the ecological-scale quantity.
#'   New locations use the original finite basis without independent field
#'   innovation. Repeated locations retain their joint dependence and rows.
#'   Range, prior normalization, original method and diagnostics are unchanged.
#'   Observation noise is omitted. Marginal intervals cannot be added to obtain
#'   an interval for a joint contrast. Use the complete matrix for that purpose.
#'   This output does not establish calibration, Monte Carlo precision or fit
#'   eligibility, including when stored vectors use a Gaussian approximation.
#' @seealso \code{\link{.msr_spatial_scale_prepare}},
#'   \code{\link{.msr_spatial_scale_plot}}, \code{\link{prediction_contract}}
#' @keywords internal
.msr_spatial_scale_summary <- function(object, newdata = NULL, level = .95,
    draws = NULL, scale = c("distance", "log"), summary = TRUE, seed = NULL) {
  if (!inherits(object, "msr_spatial_scale_prototype"))
    stop("Use a prepared private spatial scale prototype.", call. = FALSE)
  .msr_spatial_scale_objective(object)
  scale <- match.arg(scale)
  validate_scalar_logical(summary, "summary")
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
    inclusive_lower = FALSE, inclusive_upper = FALSE)
  nm <- object$prototype$covariate; spec <- object$specification
  field <- spec$scales[[nm]]$random[[".spatial_field"]]$gaussian_block
  allowed <- vapply(names(spec$scales), function(covariate) {
    effects <- spec$scales[[covariate]]$random
    if (identical(covariate, nm))
      identical(names(effects), ".spatial_field") && length(effects) == 1L else !length(effects)
  }, TRUE)
  if (is.null(field$basis) || !all(allowed))
    stop("Spatial-only summaries require the frozen field and no ordinary scale-group effects.", call. = FALSE)
  ids <- .msr_spatial_response_joint_ids(object, draws)
  if (is.null(newdata)) newdata <- spec$data
  if (!is.data.frame(newdata) || nrow(newdata) < 1L || anyDuplicated(names(newdata)))
    stop("newdata must be a nonempty data frame with distinct column names.", call. = FALSE)
  xy <- .msr_spatial_coordinates(newdata, field$coordinates)
  values <- .msr_bayes_prediction_seed(seed,
    .msr_spatial_scale_predict(object, newdata, scale = scale, draws = draws))
  if (!identical(dim(values), c(length(ids), nrow(newdata))) || any(!is.finite(values)))
    stop("Nonfinite spatial scale prediction; inspect vectors or use the log scale.", call. = FALSE)
  locations <- data.frame(row = seq_len(nrow(newdata)), x = xy[, 1], y = xy[, 2])
  ans <- values
  if (summary) {
    bounds <- vapply(c("min_D", "max_D"), function(nm) as.numeric(spec$support[[nm]]), 0)
    if (scale == "log") bounds <- log(bounds)
    ans <- cbind(locations, covariate = nm, .msr_bayes_intervals(values, level),
      p_below_min_D = colMeans(values < bounds[["min_D"]]),
      p_above_max_D = colMeans(values > bounds[["max_D"]]))
    attr(ans, "level") <- level
  }
  context <- list(n_scale = 0L, n_response = 0L, scale_policy = "error", response_policy = "error")
  contract <- .msr_bayes_prediction_contract(object, context, "kernel_scale", scale)
  if (!summary) contract$uncertainty$intervals <- "not_provided"
  contract$effects$spatial_scale <- list(covariate = nm, included = TRUE,
    coordinates = field$coordinates, basis = "frozen_training_projection_and_normalization",
    range = "fixed", new_location_innovation = FALSE)
  attr(ans, "coordinate_columns") <- field$coordinates
  attr(ans, "prediction_locations") <- locations
  attr(ans, "prediction_contract") <- contract
  ans
}

#' Plot a Private Spatial-Only Scale Summary
#'
#' Display location points colored by a posterior scale summary without
#' interpolating marginal summaries into an additional spatial surface.
#'
#' @inheritParams .msr_spatial_scale_summary
#' @param statistic Summary used for point colors: median, SD, lower or upper
#'   interval limit. The default is median. SD describes ecological-scale
#'   uncertainty at each location, not the prior field-SD parameter.
#' @return Invisibly, the summary rows that were plotted, including their
#'   prediction contract. The current graphics device receives a point plot.
#' @details Coordinate axes use the original column names. The color legend
#'   states the ecological scale or log scale. The title retains the original
#'   computation label. The subtitle retains the original
#'   diagnostic status and states that interval calibration is not established.
#'   All input rows are retained, including repeated coordinates. No new field
#'   innovation, spatial interpolation, covariance repair or fit qualification
#'   is introduced. Nonfinite colors cause an error rather than clipping.
#' @keywords internal
.msr_spatial_scale_plot <- function(object, newdata = NULL,
    statistic = c("median", "sd", "lower", "upper"), scale = c("distance", "log"),
    level = .95, draws = NULL, seed = NULL) {
  statistic <- match.arg(statistic); scale <- match.arg(scale)
  table <- .msr_spatial_scale_summary(object, newdata, level, draws, scale, seed = seed)
  values <- table[[statistic]]
  if (any(!is.finite(values)))
    stop("Nonfinite spatial summary colors; retain diagnostics or use another output scale.", call. = FALSE)
  limits <- range(values); span <- diff(limits)
  if (!is.finite(span)) stop("The summary range is too wide to plot.", call. = FALSE)
  palette <- grDevices::hcl.colors(32L, "Viridis")
  color_index <- if (span == 0) rep(16L, length(values)) else
    1L + floor(31 * (values - limits[1]) / span)
  contract <- attr(table, "prediction_contract")
  coordinates <- attr(table, "coordinate_columns")
  units <- if (scale == "distance") "map units" else "log map units"
  margin <- graphics::par("mar")
  on.exit(graphics::par(mar = margin), add = TRUE)
  outside_margin <- margin; outside_margin[4] <- max(margin[4], 8)
  outside_margin[3] <- max(margin[3], 6.1)
  graphics::par(mar = outside_margin)
  graphics::plot(table$x, table$y, asp = 1, pch = 21, bg = palette[color_index],
    col = "grey25", xlab = coordinates[1], ylab = coordinates[2],
    main = paste("Private spatial scale:", statistic, "\n", object$inference),
    sub = paste("Fit diagnostics:", contract$diagnostics_status, "| Calibration not established"))
  legend_index <- if (span == 0) 16L else c(1L, 32L)
  legend_values <- if (span == 0) limits[1] else limits
  region <- graphics::par("usr")
  graphics::legend(region[2] + .04 * diff(region[1:2]), region[4],
    legend = format(legend_values, digits = 3), xjust = 0, yjust = 1, xpd = NA,
    pch = 21, pt.bg = palette[legend_index], bty = "n", cex = .8,
    title = paste0(statistic, " (", units, ")"))
  invisible(table)
}

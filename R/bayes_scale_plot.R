#' Plot Bayesian Landscape Scales and Group Variation
#'
#' Display posterior medians and central credible intervals for fitted scales
#' or the standard deviations that describe variation among groups. Plotting
#' does not refit the model or change its posterior draws.
#'
#' @param x A fitted \code{\link{multiScale_bayes}} object with posterior draws.
#' @param type \code{"scales"} (default) plots one estimate per distinct
#'   scale-formula group. \code{"group_sd"} plots between-group SDs, with
#'   scale and response effects in separate panels.
#' @param covariate Optional character vector selecting kernel predictors for
#'   \code{type = "scales"}. \code{NULL} (default) includes all predictors.
#' @param level Central posterior interval probability; default 0.95.
#' @param scale \code{"distance"} (default) plots scales in projected map
#'   units; \code{"log"} uses natural log map units. Used only for
#'   \code{type = "scales"}. These are kernel parameters, not effective radii.
#' @param estimator \code{"auto"} (default) uses the weighted production
#'   candidates for importance-corrected fits and stored draws otherwise.
#'   \code{"draws"} always uses stored draws; \code{"weighted"} requires
#'   retained importance candidates. Weighted plotting avoids resampling noise
#'   but does not improve effective sample size or correct a missed mode.
#' @param ... Reserved for extensions; unsupported arguments cause an error.
#' @return A ggplot object. Its \code{data} contains the plotted summaries,
#'   including \code{median}, \code{lower}, and \code{upper}. Scale panels also
#'   contain \code{covariate} and \code{group}; SD panels contain
#'   \code{component} and \code{term}. Print the object inside a function or
#'   loop to display it. A failed or unavailable stored diagnostic flag is
#'   shown in the subtitle. This does not compute new chain diagnostics.
#' @details Annual intervals describe uncertainty in each fitted year's scale.
#'   The between-year SD describes population variation on the log-distance
#'   axis; it is not an annual scale's standard error. Response SDs use the
#'   response link scale. Shared-scale models without group effects have no
#'   group SD to plot. No plot establishes calibration or model adequacy.
#'   Continuous half-normal SD priors have no point mass at zero. A positive
#'   lower SD interval limit does not establish between-group heterogeneity;
#'   assess meaningful variation and compare a shared-scale model instead.
#'   For NUTS fits the stored flag checks raw parameters and sampler behavior;
#'   use \code{\link{diagnose_bayes_scales}} to check derived annual quantities.
#' @seealso \code{\link{posterior_scales}}, \code{\link{multiScale_bayes}},
#'   \code{\link{diagnose_bayes_scales}};
#'   \code{vignette("bayesian-scales", package = "multiScaleR")}.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' plot(example$fits$pooled)                    # fitted annual distances
#' plot(example$fits$pooled, type = "group_sd") # between-group variation
#' @method plot multiScaleR_bayes
#' @export
plot.multiScaleR_bayes <- function(x, type = c("scales", "group_sd"), covariate = NULL,
                                  level = 0.95, scale = c("distance", "log"),
                                  estimator = c("auto", "draws", "weighted"), ...) {
  if (length(list(...))) stop("Unused plotting arguments.", call. = FALSE)
  type <- match.arg(type); scale <- match.arg(scale); estimator <- match.arg(estimator)
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
    inclusive_lower = FALSE, inclusive_upper = FALSE)
  .msr_bayes_draw_ids(x, NULL)
  if (estimator == "auto") estimator <- if (!is.null(x$importance)) "weighted" else "draws"
  subtitle <- paste0("Central ", format(100 * level, trim = TRUE), "% credible intervals")
  if (isFALSE(x$diagnostics$passed)) subtitle <- paste(subtitle, "| Diagnostics need attention")
  else if (!isTRUE(x$diagnostics$passed))
    subtitle <- paste(subtitle, "| Stored numerical checks unavailable")
  if (type == "scales") {
    table <- posterior_scales(x, level = level, scale = scale, estimator = estimator)
    if (!is.null(covariate)) {
      if (!is.character(covariate) || !length(covariate) || anyNA(covariate) ||
          anyDuplicated(covariate) || any(!covariate %in% table$covariate))
        stop("covariate must select distinct fitted kernel predictor names.", call. = FALSE)
      table <- table[table$covariate %in% covariate, , drop = FALSE]
    }
    if (any(!is.finite(as.matrix(table[c("median", "lower", "upper")]))))
      stop("Nonfinite distance quantiles cannot be plotted; try scale = 'log'.", call. = FALSE)
    table$group <- factor(table$group, levels = unique(table$group))
    .msr_bayes_interval_plot(table, "group", "covariate", "free_x",
      "Scale-formula group", if (scale == "log") "Log kernel scale (log map units)" else
        "Kernel scale (map units)", subtitle)
  } else {
    if (!is.null(covariate)) stop("covariate selection applies only to type = 'scales'.", call. = FALSE)
    table <- .msr_bayes_plot_sd(x, level, estimator)
    if (!nrow(table)) stop("This model has no between-group SDs to plot.", call. = FALSE)
    if (any(!is.finite(as.matrix(table[c("median", "lower", "upper")]))))
      stop("Nonfinite group-SD quantiles cannot be plotted.", call. = FALSE)
    table$term <- factor(table$term, levels = unique(table$term))
    .msr_bayes_interval_plot(table, "term", "component", "free",
      "Grouping term", "Group standard deviation", subtitle)
  }
}

# Keep the plotting environment separate from the fitted object. A returned
# plot should not retain a stanfit or a large importance-candidate population.
.msr_bayes_interval_plot <- function(table, xvar, facet, facet_scales, xlabel, ylabel, subtitle) {
  ggplot2::ggplot(table, ggplot2::aes(x = .data[[xvar]], y = .data$median,
    ymin = .data$lower, ymax = .data$upper)) + ggplot2::geom_pointrange() +
    ggplot2::facet_wrap(stats::reformulate(facet), scales = facet_scales) +
    ggplot2::labs(x = xlabel, y = ylabel, subtitle = subtitle) + ggplot2::theme_bw()
}

.msr_bayes_plot_sd <- function(object, level, estimator) {
  blocks <- object$draws; weights <- NULL
  if (estimator == "weighted") {
    if (is.null(object$importance$candidates) || is.null(object$importance$weights))
      stop("Weighted plots require retained importance candidates.", call. = FALSE)
    weights <- object$importance$weights
    if (!is.numeric(weights) || length(weights) != nrow(object$importance$candidates) ||
        any(!is.finite(weights) | weights < 0) || !is.finite(sum(weights)) || sum(weights) <= 0)
      stop("Invalid retained importance weights.", call. = FALSE)
    keep <- weights > 0; weights <- weights[keep] / sum(weights[keep])
    blocks <- .msr_bayes_parameter_blocks(object$importance$candidates[keep, , drop = FALSE],
      object$parameters)
  }
  rows <- list()
  for (name in c("log_scale_sd", "log_response_sd")) {
    if (!ncol(blocks[[name]])) next
    values <- exp(blocks[[name]])
    table <- if (is.null(weights)) .msr_bayes_intervals(values, level) else
      .msr_bayes_weighted_intervals(values, weights, level)
    table$component <- if (name == "log_scale_sd") "Log-scale group SD" else "Response group SD"
    table$term <- if (name == "log_scale_sd") unlist(lapply(
      names(object$specification$scales), function(nm) {
        terms <- names(object$specification$scales[[nm]]$random)
        if (length(terms)) paste(nm, terms, sep = ":") else character()
      })) else
      names(object$specification$response$random)
    rows[[name]] <- table
  }
  if (!length(rows)) return(data.frame())
  table <- do.call(rbind, rows);rownames(table) <- NULL
  table
}

utils::globalVariables(".data")

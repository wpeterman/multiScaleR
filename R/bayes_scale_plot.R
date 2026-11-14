#' Plot Bayesian Landscape Scales and Group Variation
#'
#' Draw posterior medians and central credible intervals for the fitted
#' scales, or for the standard deviations that describe variation among
#' groups. Plotting does not refit the model or change its draws.
#'
#' @param x A fitted [multiScale_bayes()] object with posterior draws.
#' @param type `"scales"` (default) plots one interval for each distinct
#'   combination of scale covariates and groups, such as each year.
#'   `"group_sd"` plots the between-group SDs, with scale and response
#'   effects in separate panels.
#' @param covariate Optional character vector of landscape predictors to plot
#'   with `type = "scales"`. The default `NULL` plots every predictor.
#' @param level Probability of the central credible interval; default 0.95.
#' @param scale `"distance"` (default) plots scales in projected map units.
#'   `"log"` uses natural log map units, which helps when the upper limits are
#'   very large. Used only with `type = "scales"`. These are kernel scale
#'   parameters, not effective radii.
#' @param estimator `"auto"` (default) uses the weighted production
#'   candidates of an importance-corrected fit and the stored draws
#'   otherwise. `"draws"` always uses the stored draws, and `"weighted"`
#'   requires retained importance candidates. Weighted plotting avoids
#'   resampling noise but does not raise the effective sample size or recover
#'   a missed mode.
#' @param ... Not used. Supplying any further argument is an error.
#'
#' @return A ggplot object. Its `data` holds the plotted summaries, including
#'   `median`, `lower`, and `upper`. Scale panels also hold `covariate` and
#'   `group`, and SD panels hold `component` and `term`. Print the object to
#'   display it inside a function or loop. When the stored diagnostic flag of
#'   the fit failed or is unavailable, the subtitle says so. Plotting computes
#'   no new chain diagnostics.
#'
#' @details
#' # Read the plots
#'
#' With `type = "scales"`, each point is the posterior median of one group's
#' scale and each bar its central credible interval. Annual intervals describe
#' uncertainty in each fitted year's scale. They are separate intervals, not a
#' joint band. With `type = "group_sd"`, the scale panel shows the SD of log
#' scale among groups (tau), and `exp(tau)` is the multiplier on the scale
#' for a group one SD above typical. It is not the standard error of any
#' year's scale. The response panel shows random-intercept SDs on the link
#' scale. A shared-scale model without group effects has no group SD to plot.
#'
#' The half-normal prior on an SD has no point mass at zero, so a positive
#' lower limit does not establish that groups differ. Judge whether the
#' variation is meaningful, and compare a shared-scale model. No plot
#' establishes calibration or model adequacy.
#'
#' For NUTS and slice fits, the stored flag describes checks on the original
#' parameters, and for NUTS also on sampler behavior. Use
#' [diagnose_bayes_scales()] to check the mixing and simulation precision of
#' the derived log scales and group SDs.
#'
#' # Draw chain plots
#'
#' For chain plots, pass `fit$chains` for a slice fit, or
#' `rstan::extract(fit$sampler, permuted = FALSE, inc_warmup = FALSE)` for a
#' NUTS fit, to bayesplot. Use `mcmc_trace()` for trajectories,
#' `mcmc_rank_overlay()` for rank distributions by chain, and
#' `mcmc_dens_overlay()` for posterior densities by chain. Select
#' `regex_pars = "^scale_beta"` to see the scale-formula coefficients on the
#' internal log scale. Chain plots need the original post-warmup draws in
#' iteration and chain order. Approximate fits and compact examples support
#' only interval plots and posterior predictive checks.
#'
#' @seealso [posterior_scales()], [multiScale_bayes()],
#'   [diagnose_bayes_scales()], [simulate.multiScaleR_bayes()],
#'   [multiScaleR_bayes-class];
#'   `vignette("v13-bayesian-scales", package = "multiScaleR")`.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' plot(example$fits$pooled)                    # sigma for each year, meters
#' plot(example$fits$pooled, type = "group_sd") # tau and the site SD
#'
#' # The plotted summaries are kept in the ggplot object
#' p <- plot(example$fits$pooled, scale = "log")
#' p$data[c("group", "median", "lower", "upper")]   # log meters
#'
#' # Trace plots from a stored slice fit that keeps its chains
#' shared <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
#'   package = "multiScaleR"))$fit
#' if (requireNamespace("bayesplot", quietly = TRUE)) {
#'   print(bayesplot::mcmc_trace(shared$chains, regex_pars = "^scale_beta"))
#' }
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

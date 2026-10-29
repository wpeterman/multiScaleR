.msr_bayes_draw_ids <- function(object, draws) {
  if (is.null(object$draws))
    stop("This object has no posterior draws. Fit with method = 'sample', 'laplace', or 'marginal'.", call. = FALSE)
  n <- nrow(object$draws$beta)
  if (is.null(draws)) return(seq_len(n))
  validate_scalar_numeric(draws, "draws", integerish = TRUE, lower = 1, upper = n)
  unique(as.integer(round(seq(1, n, length.out = draws))))
}

.msr_bayes_draw <- function(object, i) {
  lapply(object$draws, function(z) unname(z[i, ]))
}

.msr_bayes_intervals <- function(x, level = 0.95) {
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
                           inclusive_lower = FALSE, inclusive_upper = FALSE)
  limits <- apply(x, 2L, stats::quantile, probs = c((1 - level) / 2, (1 + level) / 2),
                   names = FALSE)
  data.frame(mean = colMeans(x), median = apply(x, 2L, stats::median),
             sd = apply(x, 2L, stats::sd), lower = limits[1L, ],
             upper = limits[2L, ], row.names = NULL)
}

# Direct candidate summaries avoid adding multinomial resampling error. These
# are moments/quantiles of a finite weighted population, not exact posterior
# moments or additional effective independent draws.
.msr_bayes_weighted_intervals <- function(x, weights, level) {
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)
  limits <- apply(x, 2L, .msr_bayes_weighted_quantile, weights = weights,
                  probs = c((1 - level) / 2, .5, (1 + level) / 2))
  means <- colSums(x * weights)
  deviations <- sweep(x, 2L, means, "-")
  data.frame(mean = means, median = limits[2L, ],
    sd = sqrt(colSums(deviations^2 * weights)), lower = limits[1L, ],
    upper = limits[3L, ], row.names = NULL)
}

#' Summarize Posterior Landscape Scales
#'
#' Report shared or group-specific landscape scales in projected map units,
#' including their posterior uncertainty and mass outside the reference range.
#'
#' @param object A \code{\link{multiScale_bayes}} fit with posterior draws
#'   from sampling or either approximation. Prepared objects have no draws.
#' @param level Central posterior interval probability; default 0.95.
#' @param draws Optional positive integer limiting the number of stored
#'   posterior draws used. Draws are selected at evenly spaced positions.
#'   The default \code{NULL} uses all stored posterior draws.
#' @param scale Summary scale: \code{"distance"} (default) gives projected
#'   map units; \code{"log"} gives natural log distances in those map units.
#'   Log summaries avoid exponentiating very large distances and provide
#'   finite-moment summaries under the documented group-SD priors.
#' @param estimator \code{"draws"} (default) summarizes the stored posterior
#'   draws. For importance-corrected fits, \code{"weighted"} instead summarizes
#'   all retained production candidates using their normalized importance
#'   weights. This avoids additional resampling error. It requires
#'   \code{draws = NULL} and is unavailable for NUTS or uncorrected Laplace fits.
#' @return A data frame with one row per predictor and distinct scale-formula
#'   covariate combination. \code{covariate} identifies the landscape predictor;
#'   \code{group} labels the scale predictors; \code{row} is a representative
#'   observation row. \code{mean}, \code{median}, \code{sd}, \code{lower}, and
#'   \code{upper} summarize the stored draws in map units, or log map units
#'   when \code{scale = "log"}. The \code{mean} and \code{sd} are empirical
#'   draw summaries; a finite sample does not establish that the corresponding
#'   posterior moment exists. Exponentiate log-scale medians and interval
#'   limits to express them as distances. Exponentiating a log mean gives a
#'   geometric summary, not an arithmetic posterior mean.
#'   Weighted summaries use empirical inverse-CDF quantiles and the population
#'   SD of the finite weighted candidate distribution. They do not increase
#'   importance effective sample size, correct an unseen mode, or establish
#'   that an arithmetic posterior moment exists.
#'   \code{p_below_min_D} and \code{p_above_max_D} give posterior mass outside
#'   the input reference range. These are truncation and resolution diagnostics,
#'   not automatic evidence of a biological boundary.
#' @details For \code{~1 + (1 | year)}, each row describes one fitted year's
#'   scale after partial pooling. A median is not a separately optimized annual
#'   estimate. The shared log-scale intercept describes a population median;
#'   it differs from the arithmetic mean of annual scales. Intervals describe
#'   uncertainty conditional on the likelihood, priors, maps, and binning.
#'
#'   A random scale effect is \eqn{\tau z}, with standard-normal \eqn{z}
#'   and half-normal \eqn{\tau} whose prior scale is \eqn{s}. Its prior
#'   exponential moment of order \eqn{r} exists only when \eqn{rs < 1}.
#'   Consequently, \code{scale_group_sd = 0.5} places annual distance variance
#'   at the infinite-moment boundary; \code{scale_group_sd >= 1} also permits
#'   an infinite arithmetic distance mean. A binned likelihood can retain
#'   these tails because very large distances approach uniform weighting
#'   within the extracted neighborhood. Prefer medians and interval quantiles
#'   for distance inference, and \code{scale = "log"} for means and SDs.
#'   Log-scale SD here measures uncertainty in a fitted year's log distance;
#'   it differs from the between-year population SD in \code{summary(object)}.
#' @seealso \code{\link{multiScale_bayes}}, \code{\link{plot.multiScaleR_bayes}},
#'   \code{\link{diagnose_bayes_scales}};
#'   \code{vignette("bayesian-scales", package = "multiScaleR")} and
#'   \code{vignette("bayesian-computation", package = "multiScaleR")}.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' fit <- example$fits$pooled  # compact draws from an actual annual-scale fit
#' annual <- posterior_scales(fit)
#' annual[c("group", "median", "lower", "upper")]
#' # Distances and central 95% credible intervals are in meters in this example.
#' posterior_scales(fit, scale = "log")[c("group", "mean", "sd")]
#' # Log-distance SD measures uncertainty in one year's scale, not variation among years.
#' @export
posterior_scales <- function(object, level = 0.95, draws = NULL,
                             scale = c("distance", "log"),
                             estimator = c("draws", "weighted")) {
  if (!inherits(object, "multiScaleR_bayes"))
    stop("object must be a multiScaleR_bayes fit.", call. = FALSE)
  scale <- match.arg(scale)
  estimator <- match.arg(estimator)
  weights <- NULL
  if (estimator == "weighted") {
    if (is.null(object$importance$candidates) || is.null(object$importance$weights))
      stop("Weighted summaries require an importance-corrected fit with retained candidates.", call. = FALSE)
    if (!is.null(draws)) stop("Weighted summaries require draws = NULL.", call. = FALSE)
    weights <- object$importance$weights
    if (!is.numeric(weights) || length(weights) != nrow(object$importance$candidates) ||
        any(!is.finite(weights) | weights < 0) || !is.finite(sum(weights)) || sum(weights) <= 0)
      stop("Invalid retained importance weights.", call. = FALSE)
    # Zero target mass contributes nothing and may accompany extreme candidates
    # whose exponentiation would otherwise create 0 * Inf in a weighted moment.
    keep <- weights > 0
    weights <- weights[keep] / sum(weights[keep])
    object$draws <- .msr_bayes_parameter_blocks(
      object$importance$candidates[keep, , drop = FALSE], object$parameters)
  }
  ids <- .msr_bayes_draw_ids(object, draws)
  spec <- object$specification
  tables <- lapply(spec$covariates, function(nm) {
    sm <- spec$scales[[nm]]
    rows <- sm$representatives
    x <- matrix(NA_real_, length(ids), length(rows))
    if (!is.null(weights)) {
      log_s <- object$draws$scale_beta[, sm$beta_index, drop = FALSE] %*%
        t(sm$X[rows, , drop = FALSE])
      for (term in sm$random)
        log_s <- log_s + object$draws$scale_z[, term$z_index[term$index[rows]], drop = FALSE] *
          exp(object$draws$log_scale_sd[, term$sd_index])
      x <- if (scale == "log") log_s + log(spec$unit_conv) else exp(log_s) * spec$unit_conv
    } else for (i in seq_along(ids)) {
      p <- .msr_bayes_draw(object, ids[i])
      log_s <- drop(sm$X[rows, , drop = FALSE] %*% p$scale_beta[sm$beta_index])
      for (term in sm$random)
        log_s <- log_s + exp(p$log_scale_sd[term$sd_index]) *
          p$scale_z[term$z_index][term$index[rows]]
      x[i, ] <- if (scale == "log") log_s + log(spec$unit_conv) else
        exp(log_s) * spec$unit_conv
    }
    below <- x < if (scale == "log") log(spec$support[["min_D"]]) else spec$support[["min_D"]]
    above <- x > if (scale == "log") log(spec$support[["max_D"]]) else spec$support[["max_D"]]
    cbind(data.frame(covariate = nm, group = sm$labels, row = rows),
      if (is.null(weights)) .msr_bayes_intervals(x, level) else
        .msr_bayes_weighted_intervals(x, weights, level),
      p_below_min_D = if (is.null(weights)) colMeans(below) else colSums(below * weights),
      p_above_max_D = if (is.null(weights)) colMeans(above) else colSums(above * weights))
  })
  ans <- do.call(rbind, tables)
  rownames(ans) <- NULL
  ans
}

#' @rdname multiScale_bayes
#' @param object A \code{multiScaleR_bayes} object.
#' @param ... Reserved for method extensions; unsupported arguments cause an error.
#' @export
summary.multiScaleR_bayes <- function(object, ...) {
  if (length(list(...))) stop("Unused summary arguments.", call. = FALSE)
  if (is.null(object$draws))
    return(list(inference = object$inference, family = object$specification$family,
      observations = length(object$specification$y),
      scale_formulas = lapply(object$specification$scales, function(z) z[["formula"]])))
  b <- .msr_bayes_intervals(object$draws$beta)
  b$coefficient <- object$specification$beta_names
  scale_coefficients <- lapply(object$specification$scales, function(sm) {
    x <- object$draws$scale_beta[, sm$beta_index, drop = FALSE]
    # Report intercepts in log map units, not the internal normalized distance.
    location <- if (is.null(sm$distance_location))
      c(TRUE, rep(FALSE, ncol(x) - 1L)) else sm$distance_location
    x[, location] <- x[, location, drop = FALSE] + log(object$specification$unit_conv)
    table <- .msr_bayes_intervals(x)
    table$coefficient <- colnames(sm$X)
    table
  })
  group_sd <- lapply(c("log_scale_sd", "log_response_sd"), function(nm) {
    if (!ncol(object$draws[[nm]])) return(NULL)
    .msr_bayes_intervals(exp(object$draws[[nm]]))
  })
  names(group_sd) <- c("log_scale_group_sd", "response_group_sd")
  if (!is.null(group_sd$log_scale_group_sd))
    group_sd$log_scale_group_sd$term <- unlist(lapply(
      names(object$specification$scales), function(nm) {
        groups <- names(object$specification$scales[[nm]]$random)
        if (!length(groups)) return(character())
        paste(nm, groups, sep = ":")
      }))
  if (!is.null(group_sd$response_group_sd))
    group_sd$response_group_sd$term <- names(object$specification$response$random)
  list(inference = object$inference, family = object$specification$family,
       observations = length(object$specification$y), coefficients = b,
       scales = posterior_scales(object), scale_coefficients = scale_coefficients,
       group_sd = group_sd,
       dispersion = if (ncol(object$draws$log_dispersion))
         .msr_bayes_intervals(exp(object$draws$log_dispersion)) else NULL,
       diagnostics = object$diagnostics)
}

#' @rdname multiScale_bayes
#' @param x A \code{multiScaleR_bayes} object to print.
#' @export
print.multiScaleR_bayes <- function(x, ...) {
  if (length(list(...))) stop("Unused print arguments.", call. = FALSE)
  cat("Experimental Bayesian landscape-scale model\n")
  cat("Inference:", x$inference, "\n")
  cat("Family:", x$specification$family,
      "| observations:", length(x$specification$y), "\n")
  cat("Scale formulas:\n")
  for (nm in x$specification$covariates)
    cat(" ", nm, ":", paste(deparse(x$specification$scales[[nm]]$formula),
                             collapse = " "), "\n")
  if (!is.null(x$diagnostics) || !is.null(x$draws)) {
    # Printing reports stored evidence only. Undefined flags must not raise an
    # error or suggest that derived annual quantities have been diagnosed.
    status <- if (isTRUE(x$diagnostics$passed)) "passed" else
      if (isFALSE(x$diagnostics$passed)) "need attention" else "unavailable"
    cat(if (!is.null(x$importance)) "Marginal importance diagnostics:" else
        if (!is.null(x$marginal)) "Marginal Laplace numerical checks:" else
        if (is.null(x$laplace)) "Raw-parameter sampling checks:" else "Laplace numerical checks:",
        status, "\n")
    if (is.null(x$importance) && is.null(x$marginal) && is.null(x$laplace))
      cat("Check derived annual scales with diagnose_bayes_scales().\n")
  }
  invisible(x)
}

#' Predict Means from a Bayesian Scale Model
#'
#' Predict fitted or new-observation means using joint scale, coefficient,
#' and random-effect uncertainty. New groups require an explicit decision
#' about sampling their effects or fixing their deviations at zero.
#'
#' @param object A \code{\link{multiScale_bayes}} fit with posterior draws
#'   from sampling or either approximation.
#' @param newdata \code{NULL} (default) uses the fitted observations. Otherwise
#'   a nonempty data frame containing fixed covariates, offsets and included
#'   grouping variables. Response values are neither needed nor used.
#' @param type \code{"response"} (default) returns expected response values,
#'   or \code{"link"} returns the linear predictor. Binomial response values
#'   are success probabilities, including for grouped-binomial data.
#' @param re.form \code{NULL} (default) includes fitted response random effects.
#'   \code{NA} sets response random effects to zero while retaining fitted
#'   group-specific landscape scales. This does not integrate over random effects.
#' @param summary Logical. Return posterior summaries (default \code{TRUE})
#'   or a draws-by-observations matrix (\code{FALSE}).
#' @param level Central posterior interval probability; default 0.95.
#' @param draws Optional number of evenly spaced posterior draws to use;
#'   \code{NULL} uses all draws.
#' @param kernel_inputs Matching binned inputs for \code{newdata}, from
#'   \code{\link{kernel_prep}} or \code{\link{kernel_prep_by_group}}.
#'   Required for new observations and otherwise \code{NULL}. Use the fitted
#'   kernel, predictor definitions and physical map units. Numerical
#'   \code{unit_conv} normalizers may differ and are converted internally;
#'   this does not convert coordinate systems or physical units.
#' @param row_id Name of a \code{newdata} column matching unique prepared
#'   \code{binned$point_ids}. Observation keys may repeat. \code{NULL} uses
#'   the fitted row-ID column when retained, or new-data row names otherwise.
#'   Supply it explicitly for older saved objects. Used only with new data.
#' @param scale_re.form \code{NULL} (default) includes fitted or explicitly
#'   sampled scale random effects. \code{NA} fixes all scale deviations at zero.
#'   This is separate from \code{re.form}, which controls response effects.
#' @param new_scale_levels,new_response_levels Separate policies for unseen
#'   included scale and response group levels: \code{"error"} (default),
#'   \code{"sample"}, or \code{"zero"}. Sampling draws one new standard-normal
#'   effect per group and complete posterior vector, multiplied by that vector's
#'   group SD. Rows in the same group share it. Zeroing new effects gives a
#'   conditional reference, not a population average. Known effects retain
#'   their fitted draws. These policies do not create new fixed-factor coefficients.
#' @param seed Optional positive integer seed for new effects. When supplied,
#'   the caller's random-number state is restored, including on failure.
#' @param ... Reserved for future extensions; unused arguments cause an error.
#' @return With \code{summary = TRUE}, a data frame in the supplied observation
#'   order containing \code{row}, \code{point_id}, \code{mean}, \code{median},
#'   \code{sd}, \code{lower}, and \code{upper}. Otherwise a numeric matrix.
#'   New-data summaries also include \code{new_scale_level} and
#'   \code{new_response_level}, indicating unseen groups in included components.
#'   New-data outputs have a \code{new_groups} attribute with \code{component},
#'   \code{covariate}, \code{group}, \code{level} and \code{action}, recording
#'   which new levels were sampled or fixed at zero.
#'   Intervals describe uncertainty in expected responses or linear predictors;
#'   they do not include a new observation's sampling variation. These are
#'   marginal intervals for individual observations, not a simultaneous band.
#'   When new effects are sampled, they also include variation among new groups.
#'   The mean and SD are empirical draw summaries. A finite Monte Carlo result
#'   does not establish that the corresponding posterior or new-group moment
#'   exists; medians and quantiles are preferable for heavy-tailed predictions.
#' @details Fitting standardization, fixed-effect bases and contrasts are reused.
#'   Prediction does not refit the model, learn a new scale reference or condition
#'   on response values in \code{newdata}. An unseen fixed year factor remains an
#'   error, even when random-effect sampling is enabled. A hierarchical response
#'   year term or a scientifically justified numeric trend is needed for future
#'   response years. Exchangeable annual effects alone impose no temporal trend
#'   or autocorrelation. \code{\link{predict_scales}} can predict annual scales
#'   separately when the response formula still contains fixed year effects.
#'
#'   Use one batch for joint predictions across rows. Within that batch, all
#'   new observations in one grouping term and level share their new effect
#'   for each posterior vector. Distinct terms and groups follow the model's
#'   conditional independence, while joint posterior parameter dependence is
#'   retained. Separate calls generate separate new effects; do not combine
#'   independently generated batches to infer joint spatial or annual uncertainty.
#'
#'   Setting effects to zero is not integration over group variation. With a
#'   nonlinear link or kernel, the two predictions generally differ. Grouping
#'   columns used only by omitted components may be absent. New fixed-factor
#'   levels are always rejected. Formula transforms must remain available in
#'   base/stats or through explicit package namespaces; precompute caller-defined
#'   transforms as covariate columns. Older fits without design metadata must
#'   reproduce their stored training matrix before new predictions are allowed.
#'   Prediction-aware bases such as \code{scale()} and \code{poly()} retain
#'   their training definitions. Arbitrary expressions such as
#'   \code{I(x - mean(x))} evaluate their mean on the new batch by R formula
#'   rules; precompute with training constants when that is not intended.
#'
#'   New kernels use the neighborhoods actually prepared. Larger predicted
#'   scales cannot recover missing distant landscape cells. The method uses
#'   stored joint draws, including resamples from corrected fits, and does not
#'   validate failed fit diagnostics. Nonfinite new predictions cause an error
#'   without clipping or replacing draws. Large draws-by-observations matrices
#'   require memory; \code{draws} reduces size at a cost in simulation precision.
#' @seealso \code{\link{multiScale_bayes}}, \code{\link{posterior_scales}},
#'   \code{\link{predict_scales}}, \code{\link{simulate.multiScaleR_bayes}};
#'   \code{vignette("bayesian-prediction", package = "multiScaleR")}.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' fit <- example$fits$pooled
#' fitted_means <- predict(fit, draws = 100)
#' head(fitted_means[c("point_id", "median", "lower", "upper")])
#' # Expected counts at the fitted rows, including fitted site effects.
#' # The draw limit keeps this demonstration short; use all draws for final summaries.
#' # simulate() additionally includes observation variation.
#' new_surveys <- example$data[1:3, ]
#' new_surveys$count <- NULL
#' predict(fit, new_surveys, kernel_inputs = example$inputs,
#'   row_id = "site_year", draws = 100)
#' # These rows reuse known annual maps and group effects for new surveys.
#' @exportS3Method stats::predict
predict.multiScaleR_bayes <- function(object, newdata = NULL,
                                      type = c("response", "link"), re.form = NULL,
                                      summary = TRUE, level = 0.95, draws = NULL,
                                      kernel_inputs = NULL, row_id = NULL,
                                      scale_re.form = NULL,
                                      new_scale_levels = c("error", "sample", "zero"),
                                      new_response_levels = c("error", "sample", "zero"),
                                      seed = NULL, ...) {
  if (length(list(...))) stop("Unused prediction arguments.", call. = FALSE)
  include_response <- .msr_bayes_prediction_re_form(re.form, "re.form")
  include_scale <- .msr_bayes_prediction_re_form(scale_re.form, "scale_re.form")
  validate_scalar_logical(summary, "summary")
  type <- match.arg(type)
  ids <- .msr_bayes_draw_ids(object, draws)
  # Formula evaluation can fail after touching RNG state. Protect preparation
  # as well as new-effect sampling when the caller supplies a seed.
  context <- .msr_bayes_prediction_seed(seed,
    .msr_bayes_prediction_prepare(object, newdata, kernel_inputs, row_id,
      match.arg(new_scale_levels), match.arg(new_response_levels), include_scale, include_response))
  spec <- context$spec
  ans <- .msr_bayes_prediction_seed(seed, {
    values <- matrix(NA_real_, length(ids), length(spec$y))
    for (i in seq_along(ids)) {
      p <- .msr_bayes_prediction_parameters(.msr_bayes_draw(object, ids[i]), context)
      eta <- .msr_bayes_components(spec, p, include_response = include_response)$eta
      values[i, ] <- if (type == "link") eta else switch(spec$family,
        gaussian = eta, binomial = stats::plogis(eta), exp(eta))
      if (!is.null(newdata) && any(!is.finite(values[i, ])))
        stop("Nonfinite new-observation prediction; inspect the design, scale support and parameter draws.", call. = FALSE)
    }
    values
  })
  if (summary) {
    ans <- cbind(data.frame(row = seq_along(spec$y), point_id = spec$observation_ids),
      .msr_bayes_intervals(ans, level))
    if (!is.null(newdata)) {
      ans$new_scale_level <- Reduce(`|`, context$scale_flags)
      ans$new_response_level <- context$response_flags
    }
  }
  if (!is.null(newdata)) attr(ans, "new_groups") <- context$new_groups
  ans
}

#' Predict Landscape Scales for New Covariates and Groups
#'
#' Propagate joint posterior uncertainty into scales at new scale-formula
#' covariates, including new years when the model contains annual random effects.
#'
#' @param object A fitted \code{\link{multiScale_bayes}} object with stored
#'   joint posterior draws from sampling or either approximation.
#' @param newdata A nonempty data frame containing the fixed scale-formula
#'   covariates and any included scale grouping variables. Responses and
#'   response-only covariates are not needed. Prepared landscape bins are not
#'   required because this function predicts the scale distribution itself.
#' @param new_levels Policy for unseen scale random-effect levels:
#'   \code{"error"} (default) requires an explicit decision;
#'   \code{"sample"} draws a new standardized group effect for each joint
#'   posterior vector; \code{"zero"} fixes new deviations at zero.
#'   Known levels retain their fitted effects. This does not create coefficients
#'   for an unseen fixed-factor level.
#' @param re.form \code{NULL} (default) includes scale random effects;
#'   \code{NA} sets all scale deviations to zero. This is a reference prediction,
#'   not integration over group variation. Response effects are not used here.
#' @param level Central interval probability; default 0.95.
#' @param draws Optional positive number of evenly spaced stored posterior
#'   vectors; \code{NULL} uses all stored vectors.
#' @param summary Logical; default \code{TRUE} returns interval summaries.
#'   Repeated active scale combinations are computed once, avoiding a large
#'   draws-by-all-observations matrix when many rows share a year's scale.
#'   \code{FALSE} returns one draws-by-new-observations matrix per predictor.
#' @param scale \code{"distance"} (default) returns projected map units;
#'   \code{"log"} returns natural log distances in those map units.
#' @param seed Optional positive integer seed for new effects. Supplying it
#'   restores the caller's random-number state, including on failure.
#' @return With \code{summary = TRUE}, a data frame with one row per predictor
#'   and distinct active scale-formula combination. \code{covariate} names the
#'   predictor, \code{group} labels the covariates and groups, and \code{row}
#'   indexes a representative row of \code{newdata}. \code{new_level} identifies
#'   rows containing an unseen included scale group. \code{distribution} is
#'   \code{"conditional_posterior"}, \code{"new_group_predictive"}, or
#'   \code{"zero_deviation_reference"}. \code{mean}, \code{median}, \code{sd},
#'   \code{lower}, and \code{upper} summarize the finite draw distribution on
#'   the requested scale. Prefer distance medians and quantiles when arithmetic
#'   moments are poorly defined; see \code{\link{posterior_scales}}.
#'   \code{p_below_min_D} and \code{p_above_max_D} use the fitted input reference
#'   range, not an unsupplied future map's coverage.
#'
#'   With \code{summary = FALSE}, a named list of numeric matrices retains every
#'   new observation column in its supplied order. Repeated rows within a group
#'   share its effect. Both forms have a small \code{new_groups} attribute with
#'   \code{component}, \code{covariate}, \code{group}, \code{level} and
#'   \code{action}; it identifies the explicitly sampled or zeroed new levels.
#' @details A new year's sampled effect is shared by that year's rows within
#'   one posterior vector. Different new years are independent conditional on
#'   the population parameters, whose posterior dependence is retained.
#'   Effects for separate model components follow their specified independence.
#'   Separate calls sample independently unless their complete random-number
#'   streams coincide; use a single batch to retain within-group dependence.
#'
#'   New-group intervals include variation among future groups as well as
#'   uncertainty about population parameters. They are marginal prediction
#'   intervals, not simultaneous bands or estimates of a realized future year's
#'   scale after observing its responses. Exchangeable annual effects do not
#'   impose temporal autocorrelation. A fitted numeric trend extrapolates under
#'   that trend assumption; an unseen fixed year factor is an error.
#'   Zeroing effects does not compute a population-averaged distance.
#'
#'   The method uses stored joint draws, including stored resamples for an
#'   importance-corrected fit. It does not add effective posterior information
#'   or validate a failed fit. Nonfinite distance calculations cause an error;
#'   use log summaries for extreme scales. Formula transforms must remain
#'   available in base/stats or through an explicit package namespace. Fitted
#'   contrasts and prediction-aware bases such as \code{scale()} and
#'   \code{poly()} are reused. An expression such as
#'   \code{I(x - mean(x))} evaluates its mean on the supplied new batch by
#'   ordinary R formula rules. Use a prediction-aware basis or precomputed
#'   training constants when the original centering must be retained.
#' @seealso \code{\link{posterior_scales}},
#'   \code{\link{predict.multiScaleR_bayes}},
#'   \code{\link{simulate.multiScaleR_bayes}},
#'   \code{vignette("bayesian-prediction", package = "multiScaleR")}.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' future <- data.frame(year = c("2026", "2027"))
#' predict_scales(example$fits$pooled, future,
#'   new_levels = "sample", draws = 100, seed = 12)
#' # These intervals include new-year variation and pooling-SD uncertainty.
#' # A fixed response factor(year) still prevents predicting future counts;
#' # forecasting the scale distribution alone does not invent that coefficient.
#' @export
predict_scales <- function(object, newdata,
    new_levels = c("error", "sample", "zero"), re.form = NULL,
    level = 0.95, draws = NULL, summary = TRUE,
    scale = c("distance", "log"), seed = NULL) {
  ids <- .msr_bayes_draw_ids(object, draws)
  if (!is.data.frame(newdata) || nrow(newdata) < 1L || anyDuplicated(names(newdata)))
    stop("newdata must be a nonempty data frame with distinct column names.", call. = FALSE)
  validate_scalar_logical(summary, "summary")
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
    inclusive_lower = FALSE, inclusive_upper = FALSE)
  include <- .msr_bayes_prediction_re_form(re.form, "re.form")
  policy <- match.arg(new_levels); scale <- match.arg(scale)
  context <- .msr_bayes_prediction_seed(seed,
    .msr_bayes_prediction_scales(object$specification, newdata, policy, include))
  context$n_response <- 0L; context$response_policy <- "error"
  spec <- context$spec
  prediction_rows <- stats::setNames(lapply(spec$covariates, function(nm)
    if (summary) spec$scales[[nm]]$representatives else seq_len(nrow(newdata))), spec$covariates)
  # Build the full prediction basis before collapsing repeated combinations.
  # This preserves formula evaluation and avoids copying designs per draw.
  designs <- stats::setNames(lapply(spec$covariates, function(nm) {
    sm <- spec$scales[[nm]]
    if (summary) {
      rows <- prediction_rows[[nm]]
      sm$X <- sm$X[rows, , drop = FALSE]
      for (k in seq_along(sm$random))
        sm$random[[k]]$index <- sm$random[[k]]$index[rows]
    }
    sm
  }), spec$covariates)
  values <- .msr_bayes_prediction_seed(seed, {
    out <- stats::setNames(lapply(spec$covariates, function(nm)
      matrix(NA_real_, length(ids), length(prediction_rows[[nm]]),
        dimnames = list(NULL, rownames(newdata)[prediction_rows[[nm]]]))), spec$covariates)
    for (i in seq_along(ids)) {
      p <- .msr_bayes_prediction_parameters(.msr_bayes_draw(object, ids[i]), context)
      for (nm in spec$covariates) {
        sm <- designs[[nm]]
        x <- drop(sm$X %*% p$scale_beta[sm$beta_index])
        for (term in sm$random)
          x <- x + exp(p$log_scale_sd[term$sd_index]) * p$scale_z[term$z_index][term$index]
        out[[nm]][i, ] <- if (scale == "log") x + log(spec$unit_conv) else exp(x) * spec$unit_conv
        if (any(!is.finite(out[[nm]][i, ])))
          stop("Nonfinite predicted scales; inspect the draws and try scale = 'log'.", call. = FALSE)
      }
    }
    out
  })
  if (!summary) {
    attr(values, "new_groups") <- context$new_groups
    return(values)
  }
  tables <- lapply(spec$covariates, function(nm) {
    sm <- spec$scales[[nm]]; rows <- sm$representatives
    x <- values[[nm]]
    kind <- if (!include) rep("zero_deviation_reference", length(rows)) else
      ifelse(context$scale_flags[[nm]][rows], if (policy == "sample")
        "new_group_predictive" else "zero_deviation_reference", "conditional_posterior")
    cbind(data.frame(covariate = nm, group = sm$labels, row = rows,
      new_level = context$scale_flags[[nm]][rows], distribution = kind),
      .msr_bayes_intervals(x, level),
      p_below_min_D = colMeans(x < if (scale == "log") log(spec$support[["min_D"]]) else spec$support[["min_D"]]),
      p_above_max_D = colMeans(x > if (scale == "log") log(spec$support[["max_D"]]) else spec$support[["max_D"]]))
  })
  ans <- do.call(rbind, tables); rownames(ans) <- NULL
  attr(ans, "new_groups") <- context$new_groups
  ans
}

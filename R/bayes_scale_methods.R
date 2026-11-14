.msr_bayes_draw_ids <- function(object, draws) {
  if (is.null(object$draws))
    stop("This object has no posterior draws. Fit with method = 'sample', 'slice', 'laplace', or 'marginal'.", call. = FALSE)
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
#' Report the scale of each landscape predictor in projected map units, for
#' every distinct combination of scale-formula covariates and groups, with its
#' posterior uncertainty and the posterior mass outside the range of the
#' prepared inputs.
#'
#' @param object A fitted [multiScale_bayes()] object with posterior draws,
#'   from sampling or either approximation. A prepared object has no draws and
#'   gives an error.
#' @param level Probability of the central posterior interval; default 0.95.
#' @param draws Optional positive integer number of stored draws to use, taken
#'   at evenly spaced positions. The default `NULL` uses every stored draw.
#'   Fewer draws run faster and add simulation noise.
#' @param scale `"distance"` (default) reports scales in projected map units.
#'   `"log"` reports natural log distances in those units. Log summaries avoid
#'   exponentiating very large distances, and their means and SDs stay finite
#'   under the default group-SD prior (see "Read means of long-tailed
#'   distances").
#' @param estimator `"draws"` (default) summarizes the stored draws, for every
#'   method. `"weighted"` summarizes all retained production candidates of an
#'   importance-corrected fit with their normalized importance weights, which
#'   avoids the extra noise of resampling. It requires `draws = NULL` and a
#'   fit made with `control = list(correction = "importance")`.
#'
#' @return A data frame with one row per landscape predictor and distinct
#'   combination of scale-formula covariates and groups.
#'   * `covariate`: the landscape predictor.
#'   * `group`: a label for the combination, such as `"year=2021"`, or
#'     `"shared"` for a predictor with a single scale.
#'   * `row`: a representative observation row with that combination.
#'   * `mean`, `median`, `sd`, `lower`, `upper`: the posterior mean, median,
#'     and SD, and the limits of the central `level` credible interval. Values
#'     are in map units, or in log map units when `scale = "log"`. For log
#'     summaries, `exp()` of the median and the limits gives distances, but
#'     `exp()` of the mean gives a geometric mean, not the arithmetic mean
#'     distance.
#'   * `p_below_min_D`, `p_above_max_D`: the posterior probability that the
#'     scale lies below `min_D` or above `max_D` of the prepared inputs. These
#'     flag truncation and resolution limits. They are not evidence of a
#'     biological boundary.
#'
#'   `mean` and `sd` summarize a finite set of draws. A finite value does not
#'   show that the corresponding posterior moment exists. With
#'   `estimator = "weighted"`, quantiles are empirical inverse-CDF quantiles of
#'   the weighted candidates and `sd` is the SD of that finite weighted
#'   population.
#'
#' @details
#' # Read annual scales after partial pooling
#'
#' For a scale formula such as `~1 + (1 | year)`, each row describes one
#' fitted year's scale after partial pooling. Its median is not a separately
#' estimated annual value. The shared log-scale intercept describes a
#' population median, which differs from the arithmetic mean of the annual
#' scales. Intervals are conditional on the likelihood, the priors, the maps,
#' and the binning. With `scale = "log"`, `sd` measures uncertainty in one
#' year's log distance. It differs from the between-year SD reported in
#' `summary(object)$group_sd`.
#'
#' # Read means of long-tailed distances
#'
#' A random scale effect is \eqn{\tau z}{tau z}, with standard-normal z and a
#' half-normal tau whose prior scale is s (`scale_group_sd`). Its prior
#' exponential moment of order r exists only when \eqn{r s < 1}{r s < 1}. The
#' default `scale_group_sd = 0.5` therefore places the variance of annual
#' distances at the boundary where it becomes infinite, and
#' `scale_group_sd >= 1` also allows an infinite arithmetic mean distance. The
#' binned likelihood can keep these tails, because very large distances
#' approach uniform weighting within the extracted neighborhood. Report medians
#' and interval limits for distances, and use `scale = "log"` when you need
#' means and SDs.
#'
#' # Summarize importance-corrected fits directly
#'
#' For a fit made with `correction = "importance"`, `estimator = "weighted"`
#' summarizes the retained candidate population with its normalized weights
#' instead of the resampled draws. Both target the same posterior. The weighted
#' version avoids multinomial resampling noise, but it neither raises the
#' importance effective sample size nor recovers a posterior region the
#' candidates missed. Candidates with zero weight are dropped. The cautions
#' about infinite distance moments still apply.
#'
#' [plot.multiScaleR_bayes()] uses the weighted estimator by default for
#' corrected fits, while this function and `summary()` default to the draws.
#' Request `estimator = "weighted"` in both to keep a table and a plot
#' consistent.
#'
#' Weighted summaries are only as precise as the candidate population.
#' `importance_max` in the `control` list of [multiScale_bayes()] is a cap, and
#' a population that passes its checks at a smaller size stops there. To
#' require more candidates for precise quantiles, raise `importance_min` as
#' well, as in `control = list(correction = "importance",
#' importance_min = 65536, importance_max = 131072)`. Larger populations cost
#' more and do not guarantee a chosen quantile error. Check the retained count
#' in `fit$diagnostics$candidates`, and compare independently seeded fits.
#'
#' @seealso [multiScale_bayes()], [multiScaleR_bayes-class],
#'   [plot.multiScaleR_bayes()], [predict_scales()],
#'   [diagnose_bayes_scales()], [kernel_dist()];
#'   `vignette("v13-bayesian-scales", package = "multiScaleR")` and
#'   `vignette("v15-bayesian-computation", package = "multiScaleR")`.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' fit <- example$fits$pooled  # compact draws from an annual-scale NUTS fit
#' annual <- posterior_scales(fit)
#' annual[c("group", "median", "lower", "upper", "p_below_min_D")]
#' # Sigma for each year in meters, with central 95% credible intervals.
#' # Part of each posterior lies below min_D (100 m, one raster cell).
#'
#' posterior_scales(fit, scale = "log")[c("group", "mean", "sd")]
#' # The log-distance SD measures uncertainty in one year's scale, not the
#' # variation among years.
#'
#' # Convert the median sigma to the radius holding 90% of the kernel weight
#' kernel_dist(sigma = annual$median[1], kernel = "gaussian")
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
#' @param object A `multiScaleR_bayes` object from `multiScale_bayes()`,
#'   prepared or fitted.
#' @param ... Not used. Supplying any further argument is an error.
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
#' @param x A `multiScaleR_bayes` object to print, prepared or fitted.
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
      cat("Check derived scales with diagnose_bayes_scales().\n")
  }
  invisible(x)
}

#' Predict Means from a Bayesian Scale Model
#'
#' Predict expected responses, or linear predictors, at fitted or new
#' observations. Each prediction carries the joint posterior uncertainty in
#' the scales, the coefficients, and the random effects. A site or year that
#' was not in the fitted data needs an explicit choice, either drawing its
#' effect from the estimated spread among groups or setting it to zero.
#'
#' @param object A fitted [multiScale_bayes()] object with posterior draws,
#'   from sampling or either approximation.
#' @param newdata `NULL` (default) uses the fitted observations.
#'   Otherwise a nonempty data frame of new observations that holds the fixed
#'   covariates, the offset variables, and the grouping variables of every
#'   component you include. Response values are neither needed nor used.
#' @param type `"response"` (default) returns expected responses, which are
#'   expected counts for Poisson and NB2 models, success probabilities for
#'   binomial models (also for grouped binomial data), and means for Gaussian
#'   models. `"link"` returns the linear predictor.
#' @param re.form `NULL` (default) includes the response random effects. `NA`
#'   sets them all to zero while keeping group-specific scales. This is a
#'   reference prediction, not an average over the random effects.
#' @param summary Logical. `TRUE` (default) returns posterior summaries.
#'   `FALSE` returns the matrix of draws by observations.
#' @param level Probability of the central interval; default 0.95.
#' @param draws Optional positive integer number of evenly spaced posterior
#'   draws to use. The default `NULL` uses every draw. Fewer draws save time
#'   and memory at a cost in simulation precision.
#' @param kernel_inputs Prepared landscape inputs for `newdata`, from
#'   [kernel_prep()] or [kernel_prep_by_group()]. Required with new data and
#'   `NULL` otherwise. Use the fitted kernel, predictor definitions, and
#'   physical map units. See "Match new rows to their landscapes" in
#'   [predict.multiScaleR_bayes()].
#' @param row_id Name of the column in `newdata` that matches each row to an
#'   entry of `kernel_inputs$binned$point_ids`, or to a row name of
#'   `kernel_inputs$kernel_dat` for surface-only inputs. IDs may repeat. The
#'   default `NULL` uses the fitted row-ID column when the fit kept it, and
#'   the row names of `newdata` otherwise. Supply it for fits saved by older
#'   versions of the package. Used only with new data.
#' @param scale_re.form `NULL` (default) includes the fitted scale random
#'   effects and any sampled new ones. `NA` sets every scale deviation to
#'   zero. It is separate from `re.form`, which controls response effects.
#' @param new_scale_levels,new_response_levels What to do with a group level
#'   that was not in the fitted data, set separately for scale and response
#'   random effects. `"error"` (default) stops, `"sample"` draws a new effect,
#'   and `"zero"` sets the new effect to zero. A sampled effect is one
#'   standard-normal value per new group and posterior draw, multiplied by
#'   that draw's group SD, and shared by every row in the group. Zeroing gives
#'   a typical-group reference, not a population average. Known groups keep
#'   their fitted effects. Neither choice creates a coefficient for a new
#'   level of a fixed factor.
#' @param seed Optional positive integer seed for sampled new effects. Your
#'   own random-number state is restored afterward, even when prediction
#'   fails.
#' @param include_fixed Logical. With new data, `TRUE` recomputes declared
#'   fixed-radius landscape predictors from `kernel_inputs`, using their
#'   original preparation means and SDs. The default `FALSE` requires those
#'   columns in `newdata`. Source, metric, radius, kernel, and resolution must
#'   match the fitted preparation, including its `max_D` buffer extent. This
#'   adds no scale parameters or uncertainty.
#' @param ... Not used. Supplying any further argument is an error.
#'
#' @return With `summary = TRUE`, a data frame with one row per observation,
#'   in the supplied order.
#'   * `row`: the position of the observation.
#'   * `point_id`: its prepared landscape ID.
#'   * `mean`, `median`, `sd`, `lower`, `upper`: summaries of the predicted
#'     quantity on the scale set by `type`. `lower` and `upper` bound the
#'     central `level` credible interval.
#'   * `new_scale_level`, `new_response_level` (new data only): whether the
#'     row belongs to an unseen group of an included scale or response
#'     component.
#'
#'   With `summary = FALSE`, a numeric matrix with one row per posterior draw
#'   and one column per observation.
#'
#'   Results for new data carry a `new_groups` attribute, a data frame with
#'   `component`, `covariate`, `group`, `level`, and `action` that lists each
#'   unseen level and whether it was sampled or set to zero. Summaries also
#'   carry a `prediction_contract` attribute that describes the target, the
#'   included effects, the uncertainty, and the stored parameter map (see
#'   [prediction_contract]). Its diagnostic status repeats the fit's own
#'   status. Prediction does not make a failed fit eligible.
#'
#' @details
#' # What the intervals include
#'
#' Intervals describe uncertainty in the expected response or the linear
#' predictor. They include uncertainty in the scales, the coefficients, and
#' the fitted group effects, and, when new effects are sampled, the variation
#' among new groups. They leave out the scatter of a single new observation
#' around its mean. Use [simulate.multiScaleR_bayes()] for that. Each interval
#' is marginal for one observation, not a simultaneous band. `mean` and `sd`
#' summarize a finite set of draws and do not show that the corresponding
#' moment exists, so prefer the median and the interval limits for
#' long-tailed predictions.
#'
#' # Match new rows to their landscapes
#'
#' `kernel_inputs` supplies the prepared bins for the prediction locations and
#' maps. Each prepared location has one unique entry in `binned$point_ids`,
#' and `row_id` matches every row of `newdata` to one of them. Rows can come
#' in any order and can repeat, as for repeated visits to one location. Build
#' the inputs with [kernel_prep()] or [kernel_prep_by_group()], as in
#' `vignette("v06-year-specific-rasters", package = "multiScaleR")`. Keep the
#' fitted kernel, predictor definitions, and projected map units. The
#' numerical distance normalizer `unit_conv` may differ and is converted
#' internally, but that does not convert coordinate systems or physical
#' units. Kernel means need bins. Weighted SQ needs matching source and
#' metric specifications and the raw cell buffers.
#'
#' Prediction reuses the fitted reference center and spread, the contrasts,
#' and the formula bases. It never standardizes the new rows on their own,
#' never learns a new scale reference, and never refits the model or uses
#' responses in `newdata`. A new kernel uses only the neighborhood you
#' prepared, so a predicted scale larger than `max_D` cannot draw on cells
#' beyond it.
#'
#' # Choose how to treat new groups
#'
#' `vignette("v14-bayesian-prediction", package = "multiScaleR")` compares
#' the choices in a table. A few rules matter for every call. An unseen level
#' of a fixed year factor, such as `factor(year)`, is always an error, under
#' every policy. To predict response years you have not surveyed, use a
#' random intercept for year or a scientifically justified numeric trend.
#' Exchangeable annual effects impose no trend or autocorrelation.
#' [predict_scales()] can still predict a new year's scale when the response
#' formula has fixed year effects. Setting new effects to zero is not
#' integration over group variation, and with a nonlinear link or kernel the
#' two predictions generally differ. Grouping columns used only by a
#' component you leave out with `re.form = NA` or `scale_re.form = NA` may be
#' absent from `newdata`.
#'
#' # Predict jointly across rows
#'
#' Put every row whose joint prediction matters into one call. Within a call,
#' all new rows in the same grouping term and level share one new effect for
#' each posterior draw. Distinct terms and groups follow the model's
#' conditional independence, and the joint posterior dependence among
#' parameters is kept. Separate calls draw separate new effects, so do not
#' combine independently generated batches to infer joint spatial or annual
#' uncertainty or a total. For a total of expected values, set
#' `summary = FALSE` and sum each row of the matrix with `rowSums()`.
#'
#' The full matrix of draws by observations is built even when
#' `summary = TRUE`, at 8 bytes per value. When you need only an interval for
#' each location of a large map, predict smaller batches and keep their
#' summaries, or lower `draws` at a cost in simulation precision.
#'
#' # Use formula transformations safely
#'
#' `scale()`, `poly()`, `splines::ns()`, and `splines::bs()` keep their
#' training definitions inside ordinary expressions such as `I(scale(temp))`
#' and in transformed offsets, so single rows, batches, and effect curves use
#' the same constants. An expression such as `I(x - mean(x))` instead takes
#' its mean from the new batch, by ordinary R formula rules. Compute such a
#' covariate with training constants before fitting when you need the
#' original centering. Transforms must come from base R or stats, or carry an
#' explicit package namespace. Custom functions defined only in the fitting
#' workspace are not saved with the model, so compute them as covariate
#' columns instead. New levels of a fixed factor are always rejected. For
#' fits saved without design metadata, the rebuilt formula must reproduce the
#' stored training matrix from the original data, or prediction stops. Give
#' a namespace prefix to any custom function that reuses a learner name such
#' as `scale`.
#'
#' # Know what prediction cannot fix
#'
#' Prediction passes the fitted posterior forward. It does not repair failed
#' sampling, an inaccurate approximation, or a misspecified observation
#' model, and it does not check the fit's diagnostics, so inspect
#' `summary(object)$diagnostics` first. Every method uses its stored joint
#' draws. For an importance-corrected fit these are resampled candidates, so
#' prediction has the finite precision of those resamples rather than a fresh
#' weighted integration.
#'
#' Prediction stops on unmatched IDs, new fixed-factor levels, unavailable
#' formula transforms, incompatible kernel inputs, and nonfinite results. It
#' never clips scales to the input limits or replaces a problematic draw.
#'
#' To test predictions, hold out data at the level you intend to predict. A
#' held-out visit at a known site and year is easier to predict than a whole
#' new site or year, and repeated sites or years are not independent
#' validation replicates. [cross_validate_bayes()] refits the model without
#' the held-out rows and scores them.
#'
#' @seealso [multiScale_bayes()], [multiScaleR_bayes-class],
#'   [posterior_scales()], [prediction_contract], [predict_scales()],
#'   [simulate.multiScaleR_bayes()], [cross_validate_bayes()],
#'   [multiScaleR-prediction] for prediction across the package;
#'   `vignette("v14-bayesian-prediction", package = "multiScaleR")`.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' fit <- example$fits$pooled
#' fitted_means <- predict(fit, draws = 100)
#' head(fitted_means[c("point_id", "median", "lower", "upper")])
#' # Expected counts at the fitted rows, including fitted site effects.
#' # The draw limit keeps this demonstration short. Use all draws for
#' # final summaries.
#'
#' # New surveys at known site-years reuse their annual maps and group effects
#' new_surveys <- example$data[1:3, c("site", "year", "site_year")]
#' predict(fit, new_surveys, kernel_inputs = example$inputs,
#'   row_id = "site_year", draws = 100)
#' # simulate() would add the scatter of individual counts.
#'
#' # A fit with random year intercepts can predict new sites and years
#' forecast <- readRDS(system.file("extdata", "bayesian-prediction-example.rds",
#'   package = "multiScaleR"))
#' future <- predict(forecast$fit, forecast$newdata,
#'   kernel_inputs = forecast$kernel_inputs, row_id = "site_year",
#'   new_scale_levels = "sample", new_response_levels = "sample",
#'   draws = 200, seed = 802)
#' future[c("point_id", "median", "lower", "upper", "new_response_level")]
#' attr(future, "new_groups")   # the unseen levels and what was done to them
#' @exportS3Method stats::predict
predict.multiScaleR_bayes <- function(object, newdata = NULL,
                                      type = c("response", "link"), re.form = NULL,
                                      summary = TRUE, level = 0.95, draws = NULL,
                                      kernel_inputs = NULL, row_id = NULL,
                                      scale_re.form = NULL,
                                      new_scale_levels = c("error", "sample", "zero"),
                                      new_response_levels = c("error", "sample", "zero"),
                                      seed = NULL, include_fixed = FALSE, ...) {
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
      match.arg(new_scale_levels), match.arg(new_response_levels), include_scale, include_response,
      include_fixed))
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
  if (summary) attr(ans, "prediction_contract") <- .msr_bayes_prediction_contract(
    object, context, if (type == "link") "linear_predictor" else "expected_response",
    type, include_scale, include_response)
  ans
}

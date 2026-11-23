#' Plot marginal effects from a fitted model
#'
#' @description
#' Draw one panel per covariate that shows the model's expected response as
#' that covariate sweeps across its observed range, with every other input
#' held at a reference value, and a 95% confidence band. Read each curve as
#' the direction and size of one effect in the covariate's own units, at the
#' selected scales of effect.
#'
#' `vignette("v02-user-guide", package = "multiScaleR")` reads these plots for
#' a count model and an N-mixture model. [multiScaleR-prediction] explains
#' what the bands include.
#'
#' @param x A fitted `multiScaleR` object from [multiScale_optim()]. It must
#'   hold `opt_mod` (the model refit at the selected scales) and `scl_params`
#'   (the mean and standard deviation of each optimized covariate, used to put
#'   the x-axis back in original units). Also accepts a `multiScaleR_bayes`
#'   fit from [multiScale_bayes()], as described in "Bayesian fits".
#' @param ylab Character scalar. Y-axis label for every panel. Default
#'   `"Estimated response"`. When you leave it unset, Cox fits use
#'   `"Log relative risk"` (or `"Relative risk"` with `link = TRUE`), and GAMs
#'   with several linear predictors name the native output column.
#' @param length.out Positive integer. Number of equally spaced values at
#'   which each numeric curve is evaluated. Default `100`. Larger values give
#'   smoother curves at a small cost in time.
#' @param type Character. Component to plot for `unmarked` models, ignored
#'   for other models. Default `"state"`.
#'   * `occu()` and `pcount()` fits: `"state"` (occupancy probability or
#'     expected abundance), `"det"` (detection probability), or `"lambda"`,
#'     an alias for `"state"`.
#'   * `gdistremoval()` fits: an active `"lambda"`, `"phi"`, `"dist"`, or
#'     `"rem"` component, with `"state"` as an alias for `"lambda"`.
#'   * Other `unmarked` classes: `"state"` or `"lambda"`, as the native model
#'     names its state component.
#' @param link Logical. If `TRUE`, predictions are passed through the
#'   model's inverse link, so the curve is on the response scale. Default
#'   `FALSE`. The setting matters only for some models:
#'   * `glm` fits, which include the model in `opt_mod` from the analytic
#'     backends, and single-predictor GAMs and BAMs always plot the response
#'     scale. spaMM `HLfit` curves are always transformed too.
#'   * Cox fits plot the log relative risk with `FALSE` and the relative risk
#'     with `TRUE`.
#'   * GAMs with several linear predictors plot the native link-scale output
#'     with `FALSE` and the native response-scale output with `TRUE`.
#'     `component` chooses the column.
#'   * `unmarked`, glmmTMB, lme4, and Bayesian curves use their native
#'     response scale for both settings.
#'   * For other models, set `TRUE` if the curve appears to be on the link
#'     scale.
#' @param draws For Bayesian fits, either a positive number of stored draws
#'   to select evenly across the retained sample, or explicit draw indices.
#'   The default `NULL` uses at most 1000 evenly spaced draws. Ignored for
#'   other fits.
#' @param reference_data Optional data frame of raw predictor values for the
#'   rows used to fit the model. Default `NULL`, which recovers the data from
#'   the fitted model. Supply it when the fitted model frame keeps only
#'   transformed columns, such as `poly(temperature, 2)`. It needs unique row
#'   names that match the native training frame (rows may be reordered), and
#'   its values must reproduce the fitted features, design, factor coding, and
#'   offsets. Available for native LM, GLM, GAM, and Cox fits. See "Supply raw
#'   reference rows".
#' @param component Optional positive integer selecting a native output
#'   column of a GAM with several linear predictors. Required for those fits
#'   and otherwise must be `NULL`. With `link = FALSE` it selects a linear
#'   predictor. With `link = TRUE` it selects a native response column, whose
#'   meaning and number depend on the family. For `mgcv::ziplss()` response
#'   curves, select `component = 1`, the overall expected count, which
#'   combines both linear predictors.
#' @param add_label Logical. If `TRUE`, a caption under each panel states what
#'   the plotted limits include, such as conditioning on the selected scales
#'   or the treatment of random effects. Bayesian captions also report the
#'   number of draws and the stored status of the fit checks. Default `FALSE`,
#'   which draws the panels without a caption, as suits a figure for
#'   publication.
#'
#' @return Invisibly, a named list of `ggplot2` objects, one per plotted
#'   covariate and named after it. Each plot is also printed. A panel shows:
#'   * the predicted response as a line for a numeric covariate, or as points
#'     for the levels of a factor. The x-axis is in the covariate's original
#'     units, with optimized covariates back-transformed by `x$scl_params`.
#'   * a shaded band, or error bars for a factor, for the 95% limits when the
#'     model's prediction method supplies standard errors. Bayesian panels show
#'     the posterior median and central 95% credible interval.
#'   * a subtitle naming any `I()` term (such as `I(x^2)`) or interaction that
#'     involves the covariate.
#'
#'   For a frequentist fit, each plot's `data` element holds the plotted
#'   values: `x`, `fit`, `lwr`, `upr`, and `variable`. The list can carry
#'   attributes: a
#'   [prediction_contract] for `unmarked`, glmmTMB, lme4, multi-component GAM,
#'   and Bayesian fits, `reference_data_contract` when you supply
#'   `reference_data`, and `draw_ids` for Bayesian fits.
#'
#' @details
#' # How each panel is built
#'
#' The focal covariate sweeps from its minimum to its maximum observed value.
#' A factor uses its observed levels, with one point and interval per level.
#' Every other numeric input is held at its sample mean. For a standardized
#' kernel-weighted covariate the mean is 0, the average landscape value at the
#' survey points. Other non-numeric inputs are held at their first observed
#' value, which keeps an actual fitted level. Variables that enter only
#' through an offset are held at their sample mean and get no panel. Numeric
#' predictors that [multiScale_optim()] did not optimize keep their own units.
#'
#' Panels come from the active model terms. For ordinary `lm`, `glm`,
#' `MASS::glm.nb`, and time-fixed, unpenalized Cox and `clogit` fits, the
#' fitted terms expand dot shortcuts, keep active offset covariates, and
#' annotate interactions. Removed terms add no panel or prediction input.
#' Quoted predictor names keep their original spelling. Factor contrasts and
#' fixed transformations are evaluated by the fitted model, and the plotting
#' grid never relearns a learned transformation. When `x$opt_mod` wraps
#' another fitted model, the nested model is used to find the predictors,
#' recover the data, and predict.
#'
#' With an interaction, each panel is one slice of the interaction surface,
#' at the mean of the interacting variable. A message says so, and the panel
#' subtitle names the variables held at their mean.
#'
#' # Scale of the curves and their limits
#'
#' * **GLMs.** Predictions and their standard errors are computed on the link
#'   scale. The limits are the estimate plus and minus 1.96 standard errors,
#'   transformed by the inverse link (for example `exp()` for a Poisson log
#'   link) and ordered for increasing or decreasing links. Limits with
#'   nonfinite standard errors, invalid endpoints for the link or family, or
#'   an inverse-link pole inside the interval are not shown, with a warning.
#'   Fitted values stay native, and the method does not clip or repair
#'   invalid limits.
#' * **Models without a link**, such as `lm`: the estimate plus and minus 1.96
#'   standard errors.
#' * **glmmTMB and lme4:** native response predictions, including
#'   zero-inflation and truncated-family means, with native response-scale
#'   Wald limits when available.
#' * **spaMM:** the linear predictor and its native fixed-effect variance
#'   (`fixefVar`), transformed through the inverse link. The response variance
#'   is not used as an interval for the mean.
#' * **pscl `zeroinfl`:** class-specific predictions with their own limits.
#' * **`mgcv::ziplss()`:** response-scale Wald limits that propagate the full
#'   native coefficient covariance through both linear predictors. They can
#'   extend below zero.
#'
#' Every frequentist band is a confidence band for the expected response. It
#' holds the selected scales of effect at their estimates, even for the
#' analytic backends, and it adds no observation noise, so it is not an
#' interval for a future count.
#'
#' # Random effects
#'
#' Mixed-model and `unmarked` curves set the random effects to zero
#' (`re.form = NA`). These are predictions at the stated covariate values for
#' a group with zero deviation, not averages over sites, years, or the
#' random-effect distribution. To show a curve for each year, build your own
#' prediction grid conditioned on that year with the fitted factor coding.
#' Zero-effect curves need only the fixed design and offsets. They do not
#' evaluate unused learned random slopes or need new grouping variables.
#'
#' # Unmarked models
#'
#' Panels use the variables of the selected component's formula, recovered
#' by column name. Factor contrasts such as `Year2` are coefficient labels,
#' not data columns, so the factor itself gets one point per observed level.
#' Site keys that are not fixed predictors, and detection-only predictors,
#' get no state panels. Other factors are held at their first observed
#' level. Curves stay on the native response scale for both `link` settings,
#' and their limits are the native `lower` and `upper` confidence limits.
#'
#' * **`occu()` and `pcount()`.** State curves use the original site rows,
#'   and detection curves use the original observation rows. Detection is
#'   conditional on occupancy for `occu()`. For `pcount()` it is the
#'   probability of detecting each individual on one survey. Curves for
#'   Poisson, negative binomial, and ZIP `pcount()` fits use the same checked
#'   native target as point, snapshot, and raster predictions. ZIP state
#'   curves keep the overall expected abundance, including structural zeros.
#'   Each component prepares only its own training basis, and fixed scaling,
#'   polynomial, and spline calls and transformed offsets use the original
#'   raw covariates before site removal.
#' * **`gdistremoval()`.** Curves use the component's own site,
#'   primary-period, or removal-interval rows and its fitted design, with
#'   spline knots, polynomial recipes, nested scaling, and factor
#'   interactions fixed at training. `"lambda"` is the native abundance or
#'   density parameter, `"phi"` is availability, `"dist"` is the native
#'   distance-detection parameter, and `"rem"` is the removal parameter.
#'   Inactive components have no curve. These are native component values,
#'   not survey averages or random-effect marginal means. GDR models with
#'   active formula offsets or TMB zero inflation are rejected, as
#'   [aic_tab()] explains.
#' * **Other classes** that would need a different basis or a derived state
#'   prediction stop rather than plot one.
#'
#' Factor coding uses retained or unambiguously recoverable original
#' metadata. A logical pipe inside `I()`, such as `I(a | b)`, is a fixed
#' operation, not a random-effect declaration, but the native `unmarked`
#' parser can reject it, or a column name containing a pipe, before fitting.
#' Then precompute the logical predictor as a named column, or rename the
#' column. For complete logical inputs, `I(as.integer(a) + as.integer(b) > 0)`
#' expresses the same OR without a pipe. With missing values, precompute
#' `a | b` to keep R's missing-value logic.
#'
#' # Supply raw reference rows
#'
#' A fitted model frame can keep `poly(temperature, 2)`, `I(elevation^2)`, or
#' a compound offset without the raw columns behind them. A square cannot
#' recover the sign of its input, and an orthogonal polynomial needs its
#' fitted prediction metadata. If you kept the raw predictors, pass them as
#' `reference_data`, with row names that identify the fitted rows:
#'
#' ```
#' fit_rows <- row.names(fit$opt_mod$model)
#' reference <- training_data[match(fit_rows, row.names(training_data)), ]
#' row.names(reference) <- fit_rows
#' plot_marginal_effects(fit, reference_data = reference)
#' ```
#'
#' The rows are matched by ID and checked against the evaluated training
#' features, the full native coefficient design or smooth basis, the factor
#' contrasts, and the formula offsets. Missing rows or inputs, nonfinite
#' values, new levels, and mismatches raise an error. Transformations and
#' contrasts are taken from the fit, not relearned from the supplied rows or
#' the plotting grid. Supported nested calls use a single standard learner's
#' original retained attributes, and the default refitter can also recover
#' supported Cox recipes while the original training rows are available.
#' Precompute a predictor whose recipe is ambiguous or missing. For native
#' conditional logistic fits with learned predictors, keep the model frame
#' with `model = TRUE` and fit complete choice sets with
#' `na.action = stats::na.fail`, because some survival versions drop training
#' attributes when they omit rows.
#'
#' For native refits, the reference must also hold the final fitted spatial
#' columns in their fitted units, standardized by default. The analytic
#' backends instead rebuild the selected kernel values from their verified
#' preparation. Passing the check establishes that the supplied rows
#' reproduce the fitted features. It cannot prove a unique raw value behind a
#' square, polynomial, or other many-to-one transformation, and it does not
#' refit the model. The `reference_data_contract` attribute records this.
#' Offsets supplied outside the formula, unsupported native structures,
#' matrix-valued raw inputs, and nested learned transformations without
#' recoverable metadata are rejected. `unmarked` and Bayesian curves keep
#' their own data requirements.
#'
#' # GAMs with several linear predictors
#'
#' Select `component` explicitly. Link curves use the selected predictor's
#' active variables. Response curves consider every component's inputs,
#' because one response column can depend on several predictors. Native
#' response meanings are kept. For example, `mgcv::gaulss()` returns the mean
#' and the reciprocal standard deviation, and multinomial response columns
#' are category probabilities. Limits use native standard errors on the
#' selected output scale and condition on the fitted smoothing parameters and
#' scales of effect. Native GAM curves keep the saved factor coding of
#' interaction-only factors when session contrast options change. The
#' returned [prediction_contract] identifies the native output column.
#'
#' # Cox and step-selection models
#'
#' Cox panels, including native `survival::clogit()` fits and the conditional
#' logistic backend, default to the log relative risk. `link = TRUE`
#' exponentiates it to a relative risk or relative selection score. Neither
#' is a probability of occurrence or survival, or a probability within a
#' choice set. Strata and cluster-only inputs are held at their first
#' observed value, which keeps an actual fitted group.
#'
#' # Bayesian fits
#'
#' Bayesian panels use the stored marginal response-coefficient draws and the
#' original pooled reference normalization. Each landscape covariate sweeps a
#' hypothetical grid of values while the others stay at their reference
#' values. The coefficient draws keep their posterior variation, but the
#' landscape is not reevaluated at each draw's kernel scales. Response random
#' effects are zero, and observation noise is excluded. Use
#' [predict.multiScaleR_bayes()] with matched spatial inputs for
#' location-specific posterior expectations. A credible band does not
#' establish convergence or calibration, so check the fit diagnostics. Their
#' stored status appears in the caption when `add_label = TRUE`.
#'
#' @seealso [multiScaleR-prediction], [predict.multiScaleR()],
#'   [prediction_contract], [multiScale_optim()], and
#'   `vignette("v02-user-guide", package = "multiScaleR")`.
#'
#' @examples
#' ## Running example: counts at 100 survey points and two landscape layers
#' data("landscape_counts")
#' data("surv_pts")
#' pts <- terra::vect(surv_pts)
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))[[c("land1", "land2")]]
#' inputs <- kernel_prep(pts, land, max_D = 1700, kernel = "gaussian",
#'                       bin = TRUE, verbose = FALSE)
#' dat <- data.frame(landscape_counts, inputs$kernel_dat)
#' start <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
#' fit <- multiScale_optim(start, inputs, backend = "profile_glm",
#'                         n_cores = 1, verbose = FALSE)
#'
#' ## Expected count across each covariate, others at their means. The bands
#' ## are 95% limits that hold the scales of effect at their estimates.
#' curves <- plot_marginal_effects(fit, ylab = "Expected count")
#' names(curves)
#'
#' ## The values behind one panel. land1 is back on its original scale, the
#' ## kernel-weighted proportion of habitat around a point.
#' head(curves$land1$data)
#'
#' \donttest{
#' ## A caption stating what the limits include, useful while exploring
#' plot_marginal_effects(fit, ylab = "Expected count", add_label = TRUE)
#' }
#'
#' ## Expected abundance from a stored N-mixture fit, with native 95%
#' ## confidence limits
#' nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
#'                             package = "multiScaleR"))
#' plot_marginal_effects(nmix, type = "state", ylab = "Expected abundance")
#'
#' @export
#' @importFrom unmarked coef predict
#' @importFrom insight get_parameters link_inverse find_predictors get_predicted
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line xlab ylab labs
#' @importFrom cowplot theme_cowplot
#' @importFrom utils globalVariables
plot_marginal_effects <- function(x,
                                  ylab = "Estimated response",
                                  length.out = 100,
                                  type = "state",
                                  link = FALSE, draws = NULL,
                                  reference_data = NULL, component = NULL,
                                  add_label = FALSE) {
  validate_scalar_logical(link, "link")
  validate_scalar_logical(add_label, "add_label")
  if (inherits(x, "multiScaleR_bayes")) {
    if (!is.null(reference_data) || !is.null(component))
      stop("reference_data and component do not redefine Bayesian effect curves.", call. = FALSE)
    return(.msr_bayes_marginal_effects(x, ylab, length.out, draws, add_label))
  }
  if (!is.null(draws)) stop("draws applies only to Bayesian effect curves.", call. = FALSE)
  if(!inherits(x, "multiScaleR")){
    stop("`x` must be a fitted `multiScaleR` object.")
  }
  validate_character_scalar(ylab, "ylab")
  validate_scalar_numeric(length.out, "length.out", integerish = TRUE, positive = TRUE)

  if(is.null(x$opt_mod)){
    stop("`x` must contain an `opt_mod` fitted model.")
  }
  if(!is.list(x$scl_params) ||
     !all(c("mean", "sd") %in% names(x$scl_params))){
    stop("`x` must contain `scl_params` with `mean` and `sd` elements.")
  }

  mod <- x$opt_mod
  .msr_gdr_native_guard(.analysis_model(mod))
  analysis_mod <- .msr_native_active_metadata(
    .msr_native_fixed_metadata(.analysis_model(mod), nested = TRUE))
  scl <- x$scl_params
  native_gam_components <- inherits(analysis_mod, "gam") &&
    is.list(analysis_mod$pterms) && !inherits(analysis_mod$pterms, "terms")
  if (native_gam_components) {
    if (is.null(component))
      stop("Multi-component GAM effect curves require an explicit native output selection with component.", call. = FALSE)
    validate_scalar_numeric(component, "component", integerish = TRUE, positive = TRUE)
    if (!link && component > length(analysis_mod$pterms))
      stop("component must select an available native GAM linear predictor.", call. = FALSE)
  } else if (!is.null(component))
    stop("component applies only to multi-component GAM effect curves.", call. = FALSE)
  reference <- if (is.null(reference_data)) NULL else
    .msr_curve_reference_context(x, analysis_mod, reference_data)
  if (!is.null(reference)) analysis_mod <- reference$model
  namespace(analysis_mod)
  native_mixed_response <- inherits(analysis_mod, c("glmmTMB", "merMod"))
  if (missing(ylab) && inherits(analysis_mod, "coxph"))
    ylab <- if (link) "Relative risk" else "Log relative risk"

  if(inherits(analysis_mod, "glm") && !native_gam_components){
    link <- TRUE
  }
  if (missing(ylab) && native_gam_components)
    ylab <- paste("Native GAM", if (link) "response column" else "linear predictor", component)

  # These are populated for non-unmarked models; used when building plot subtitles
  poly_terms  <- character(0)
  poly_vars   <- character(0)
  inter_terms <- character(0)
  inter_vars  <- character(0)

  ordinary_unmarked_curve <- methods::is(mod, "unmarkedFit") &&
    class(mod)[1L] %in% c("unmarkedFitOccu", "unmarkedFitPCount")
  gdr_unmarked_curve <- class(mod)[1L] %in% c("unmarkedFitGDR", "unmarkedFitGDRField")
  trained_unmarked_curve <- ordinary_unmarked_curve || gdr_unmarked_curve
  # Extract variables from fitted model
  if (any(grepl("^unmarked", class(mod)))) {
    allowed <- if (gdr_unmarked_curve) c("state", names(.msr_gdr_active_formulas(mod))) else
      if (ordinary_unmarked_curve) c("lambda", "state", "det") else c("lambda", "state")
    if (is.null(type) || !type %in% allowed) {
      if (gdr_unmarked_curve)
        stop("`type` must name an active native GDR component: ", paste(setdiff(allowed, "state"), collapse = ", "), ".")
      if (ordinary_unmarked_curve)
        stop("`type` must be specified as 'state', 'lambda', or 'det'.")
      stop("`type` must be specified as either 'lambda' or 'state'.")
    }
    # Keep the documented plotting alias; ordinary providers name it state.
    if (ordinary_unmarked_curve && identical(type, "lambda")) type <- "state"
    if (gdr_unmarked_curve && identical(type, "state")) type <- "lambda"
    component <- .msr_unmarked_marginal_formula(mod, type)
    dat <- if (trained_unmarked_curve)
      get("get_covariates", asNamespace("unmarked"))(mod, type) else
        unmarked::siteCovs(unmarked::getData(mod))
    active_terms <- stats::terms(component, data = dat)
    vars <- .msr_marginal_variables(active_terms, offsets = TRUE)
    c_vars <- .msr_marginal_variables(active_terms, offsets = FALSE)
    # Legacy class paths may evaluate inactive frame context. The checked
    # ordinary zero-effect predictor uses only its active fixed design/offsets.
    if (!trained_unmarked_curve) {
      frame_variables <- as.list(attr(active_terms, "variables"))[-1L]
      vars <- unique(c(vars, unlist(lapply(frame_variables, all.vars), use.names = FALSE)))
    }
    missing_vars <- setdiff(vars, names(dat))
    if (length(missing_vars)) {
      level <- if (trained_unmarked_curve && type %in% c("det", "rem")) "observation" else
        if (gdr_unmarked_curve && type %in% c("phi", "dist")) "primary-period" else "site"
      stop("Could not recover ", level,
           " covariates for marginal effects: ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }

    # Use column names, not coefficient names or column positions. Repeated
    # site-year keys may coexist with factors, offsets, and unscaled predictors.
    pred_defaults <- lapply(dat[vars], function(values) {
      if (is.numeric(values)) return(mean(values, na.rm = TRUE))
      values[which(!is.na(values))[1]]
    })
    factor_vars <- vars[vapply(dat[vars], is.factor, logical(1))]
    custom_contrasts <- any(vapply(dat[factor_vars], function(values) {
      fitted <- attr(values, "contrasts")
      if (is.null(fitted)) return(FALSE)
      attr(values, "contrasts") <- NULL
      !isTRUE(all.equal(fitted, stats::contrasts(values)))
    }, logical(1)))
    # Ordinary native curves use the same checked target as points and maps.
    # Zero effects need no new random slopes, groups or random-design mapping.
    native_curve <- if (trained_unmarked_curve)
      .msr_unmarked_prediction_function(mod, re.form = NA, component = type) else NULL
    basis <- if (is.null(native_curve)) .msr_unmarked_curve_basis(active_terms, dat) else NULL
    plot_list <- lapply(c_vars, function(v) {
      values <- dat[[v]]
      if (all(is.na(values))) {
        stop("Covariate '", v, "' has no observed values.", call. = FALSE)
      }
      if (is.numeric(values)) {
        x_seq <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE),
                     length.out = length.out)
      } else if (is.factor(values)) {
        # Retain factor levels, ordering, and any fitted contrast attributes.
        x_seq <- values[match(levels(droplevels(values)), as.character(values))]
      } else {
        x_seq <- unique(values[!is.na(values)])
      }
      x_unscaled <- x_seq
      if (v %in% names(scl$mean)) {
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      newdata <- as.data.frame(lapply(pred_defaults, function(value) {
        rep(value, length(x_seq))
      }), stringsAsFactors = FALSE, check.names = FALSE)
      newdata[[v]] <- x_seq

      # Keep original factor attributes for the legacy class fallback. The
      # ordinary predictor applies its independently checked training coding.
      for (factor_var in factor_vars) {
        attr(newdata[[factor_var]], "contrasts") <-
          if (custom_contrasts) attr(dat[[factor_var]], "contrasts") else NULL
      }
      preds <- tryCatch(
        if (!is.null(native_curve)) {
          native_curve(mod, type = type, newdata = newdata, re.form = NA)
        } else if (custom_contrasts || basis$learned) {
          .msr_unmarked_contrast_prediction(mod, basis$terms, newdata, type)
        } else {
          unmarked::predict(mod, newdata = newdata, type = type, re.form = NA)
        },
        error = function(e) {
          stop(
            sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
            call. = FALSE
          )
        }
      )

      fit <- preds$Predicted
      se <- preds$SE
      lwr <- preds$lower
      upr <- preds$upper

      data.frame(x = x_unscaled, fit = fit, lwr = lwr, upr = upr, variable = v)
    })

  } else {
    # A conditional Cox score uses exp(lp), even when generic metadata treats
    # the clogit subclass as binomial. It is not an occurrence probability.
    curve_inverse <- if (inherits(analysis_mod, "coxph")) exp else
      if (!native_mixed_response && !native_gam_components &&
          !inherits(analysis_mod, "zeroinfl")) link_inverse(analysis_mod) else NULL
    curve_predictor <- if (inherits(analysis_mod, "gam")) .msr_predict_gam else predict
    native_fixed_terms <- if (class(analysis_mod)[1L] %in%
        c("lm", "glm", "negbin", "coxph", "clogit")) stats::terms(analysis_mod) else NULL
    cox_groups <- if (inherits(analysis_mod, "coxph"))
      .msr_cox_curve_conditioning(stats::terms(analysis_mod)) else
        list(inputs = character(), fixed_panel_inputs = character())
    if (identical(x$backend, "conditional_logit")) {
      vars <- .msr_active_formula_variables(x$profile_fit$contract$terms)
      c_vars <- vars
    } else if (native_gam_components) {
      forms <- stats::formula(analysis_mod)
      vars <- unique(unlist(lapply(forms, .msr_marginal_variables, offsets = TRUE),
        use.names = FALSE))
      panels <- lapply(forms, .msr_marginal_variables, offsets = FALSE)
      c_vars <- if (link) unique(unlist(panels, use.names = FALSE)) else panels[[component]]
      context <- unique(unlist(lapply(analysis_mod$pterms, function(tt) {
        expressions <- as.list(attr(tt, "variables"))[-1L]
        if (attr(tt, "response") > 0L) expressions <- expressions[-attr(tt, "response")]
        unlist(lapply(expressions, all.vars), use.names = FALSE)
      }), use.names = FALSE))
      vars <- union(vars, context)
    } else if (!is.null(native_fixed_terms) || inherits(analysis_mod, c("glm", "lm"))) {
      # Native terms retain active variables and quoted names without a
      # provider's optional name-sanitizing interpretation.
      vars <- .msr_active_formula_variables(stats::terms(analysis_mod))
      c_vars <- vars
    } else {
      c_vars <- find_predictors(analysis_mod)[[1]]
      vars <- unlist(find_predictors(analysis_mod))
    }
    vars <- union(vars, cox_groups$inputs)
    c_vars <- setdiff(c_vars, setdiff(cox_groups$inputs, cox_groups$fixed_panel_inputs))
    # Offset-only variables describe effort, not covariate effects. Keep them
    # in vars so prediction holds them at their mean, but omit their panels.
    # A variable that also enters as a model term remains a plotted covariate.
    offset_vars <- tryCatch(insight::find_offset(analysis_mod), error = function(e) NULL)
    if (!native_gam_components && length(offset_vars) > 0) {
      term_vars <- tryCatch({
        labels <- attr(if (is.null(native_fixed_terms)) stats::terms(formula(analysis_mod)) else
          native_fixed_terms, "term.labels")
        if (length(labels)) all.vars(stats::reformulate(labels)) else character(0)
      }, error = function(e) character(0))
      c_vars <- setdiff(c_vars, setdiff(offset_vars, term_vars))
    }
    if (!is.null(reference)) vars <- union(vars, reference$required_inputs)
    dat_all <- if (!is.null(reference)) reference$data else
      if (identical(x$backend, "joint_mixed") &&
        is.data.frame(x$mixed_fit$prediction_data)) x$mixed_fit$prediction_data else
          .model_data(analysis_mod)
    if(is.null(dat_all)){
      dat_all <- extract_model_data(analysis_mod)
    }
    if(is.null(dat_all)){
      stop("Could not recover the original model data needed to plot marginal effects.")
    }

    design <- x$profile_fit$contract$objective$response_design
    if (is.null(reference) && !is.null(design$expression_programs)) {
      # Evaluated powers cannot identify the underlying signed kernel values.
      # Rebuild those values from the selected, verified spatial preparation.
      e <- x$profile_fit$contract$objective$engine
      if (inherits(e, "msr_surface_engine") && isTRUE(e$rebuild_from_inputs))
        e <- .msr_prepare_surface_engine(x$kernel_inputs, e$covariates, e$rows,
                                         e$standardize, e$implementation)
      value <- .msr_eval_kernel_engine(e,
        x$profile_fit$par[seq_along(e$parameter_names)])$value
      if (nrow(dat_all) != nrow(value)) stop("Effect-curve training rows do not match the fitted kernels.")
      for (nm in e$covariates) dat_all[[nm]] <- value[, nm]
      missing <- setdiff(vars, names(dat_all))
      if (length(missing)) stop("Raw fixed-transform inputs were not retained for effect curves: ",
        paste(missing, collapse = ", "), ". Refit with recoverable raw training data.", call. = FALSE)
    }

    absent <- setdiff(vars, names(dat_all))
    if (length(absent)) stop("Raw inputs were not retained for effect curves: ",
      paste(absent, collapse = ", "), ". Supply verified reference_data.", call. = FALSE)

    # Parse formula term labels to detect I() polynomial and interaction terms.
    # predict() handles I(x^2) automatically when newdata contains the base column,
    # so no special prediction logic is needed. Detection is used only for plot
    # annotations and the interaction warning below.
    curve_terms <- if (native_gam_components) lapply(forms, stats::terms) else
      list(tryCatch(if (is.null(native_fixed_terms)) stats::terms(formula(analysis_mod)) else
        native_fixed_terms, error = function(e) NULL))
    fterms <- unlist(lapply(curve_terms, attr, "term.labels"), use.names = FALSE)
    term_order <- unlist(lapply(curve_terms, attr, "order"), use.names = FALSE)
    unique_terms <- !duplicated(fterms)
    fterms <- fterms[unique_terms]; term_order <- term_order[unique_terms]
    poly_terms  <- fterms[grepl("^I\\(", fterms)]
    inter_terms <- fterms[term_order > 1L]
    poly_variables <- lapply(poly_terms, function(term) all.vars(str2lang(term)))
    inter_variables <- lapply(inter_terms, function(term) all.vars(str2lang(term)))

    # Which base variables appear inside an I() term?
    poly_vars <- c_vars[sapply(c_vars, function(v) {
      any(vapply(poly_variables, function(names) v %in% names, logical(1)))
    })]

    # Which base variables appear in at least one interaction term?
    inter_vars <- if(length(inter_terms) > 0) {
      c_vars[sapply(c_vars, function(v) {
        any(vapply(inter_variables, function(names) v %in% names, logical(1)))
      })]
    } else character(0)

    if(length(inter_terms) > 0) {
      message(
        "Interaction term(s) detected: ", paste(inter_terms, collapse = ", "), ".\n",
        "Each panel shows the effect of that variable holding all others at their ",
        "sample mean. This represents one conditional slice of the interaction ",
        "surface and may not capture the full effect of each variable."
      )
    }

    # Hold non-focal predictors at representative observed values. Numeric
    # columns use their sample mean; non-numeric columns keep the first
    # non-missing observed value so special terms like strata() can still be
    # evaluated during prediction.
    pred_defaults <- stats::setNames(vector("list", length(vars)), vars)
    # A bare strata(id) input is held at the first fitted ID, taken from the
    # fitted strata labels. Recovered training data can store numeric IDs as
    # a factor (retained model frames), which survival reads as a new level.
    strata_ids <- .msr_strata_ids(analysis_mod)
    for (var in vars) {
      if (!is.null(strata_ids) && var == strata_ids$name &&
          var %in% setdiff(cox_groups$inputs, c_vars)) {
        pred_defaults[[var]] <- strata_ids$ids[1L]
        next
      }
      if (!var %in% names(dat_all)) {
        pred_defaults[[var]] <- 0
        next
      }

      values <- dat_all[[var]]
      if (var %in% setdiff(cox_groups$inputs, c_vars)) {
        pred_defaults[[var]] <- values[1L]
      } else if (is.numeric(values)) {
        pred_defaults[[var]] <- mean(values, na.rm = TRUE)
      } else {
        keep <- which(!is.na(values))[1]
        pred_defaults[[var]] <- if (is.na(keep)) values[1] else values[keep]
      }
    }

    plot_list <- lapply(c_vars, function(v) {
      values <- dat_all[[v]]
      if (all(is.na(values)))
        stop("Covariate '", v, "' has no observed values.", call. = FALSE)
      if (is.numeric(values)) {
        x_seq <- seq(min(values, na.rm = TRUE), max(values, na.rm = TRUE),
                     length.out = length.out)
      } else if (is.factor(values)) {
        # Preserve fitted levels, ordering, and contrasts in native prediction.
        x_seq <- values[match(levels(droplevels(values)), as.character(values))]
      } else {
        x_seq <- unique(values[!is.na(values)])
      }
      x_unscaled <- x_seq
      if (v %in% names(scl$mean)) {
        x_unscaled <- (x_seq * scl$sd[v]) + scl$mean[v]
      }

      # Construct newdata: focal variable sweeps its range; all others at sample mean.
      # R's predict() evaluates I(x^2), log(x), etc. from the base column, so
      # newdata only needs to contain the base variable names.
      newdata <- as.data.frame(
        lapply(pred_defaults, function(value) rep(value, length(x_seq))),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      newdata[[v]] <- x_seq

      if (native_gam_components) {
        preds <- .msr_gam_curve_component(analysis_mod, newdata, component, link)
      } else if(inherits(analysis_mod, "HLfit")){
        preds_ <- tryCatch(
          predict(analysis_mod, newdata = newdata, type = "link",
            variances = list(fixefVar = TRUE), re.form = NA),
          error = function(e) {
            stop(
              sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
              call. = FALSE
            )
          }
        )
        preds <- list(preds = as.vector(preds_),
                      se = sqrt(attr(preds_, "fixefVar")))
      } else if(inherits(analysis_mod, "zeroinfl")) {
        preds <- tryCatch(
          as.data.frame(get_predicted(analysis_mod, data = newdata)),
          error = function(e) {
            stop(
              sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
              call. = FALSE
            )
          }
        )
      } else if (native_mixed_response) {
        # Native response prediction also handles zero-inflation and truncated
        # families. Transforming only the conditional count link would not.
        preds <- tryCatch(
          predict(analysis_mod, newdata = newdata, type = "response",
                  re.form = NA, se.fit = TRUE),
          error = function(e) stop(sprintf(
            "Failed to compute marginal effects for covariate '%s': %s", v,
            e$message), call. = FALSE)
        )
      } else if (inherits(analysis_mod, "coxph")) {
        preds <- predict(analysis_mod, newdata = newdata, type = "lp",
                         reference = "zero", se.fit = TRUE)
      } else {
        preds <- tryCatch(
          curve_predictor(analysis_mod, newdata = newdata, se.fit = TRUE),
          error = function(e) {
            stop(
              sprintf("Failed to compute marginal effects for covariate '%s': %s", v, e$message),
              call. = FALSE
            )
          }
        )
      }

      if (native_gam_components) {
        if (is.list(preds)) {
          fit_ <- preds$fit
          se <- preds$se.fit
          invalid <- !is.finite(se) | se < 0
          if (any(invalid)) {
            warning("Native GAM confidence limits unavailable for nonfinite or negative standard errors.", call. = FALSE)
            se[invalid] <- NA_real_
          }
          lwr <- fit_ + qnorm(.025) * se
          upr <- fit_ + qnorm(.975) * se
        } else {
          fit_ <- preds; lwr <- upr <- NA_real_
        }
      } else if(inherits(analysis_mod, "zeroinfl")) {
        fit_ <- preds$Predicted
        lwr  <- preds$CI_low
        upr  <- preds$CI_high
        if(is.null(lwr) && is.null(upr)) lwr <- upr <- NA

      } else if(!native_mixed_response && !is.null(curve_inverse) && (link || inherits(analysis_mod, "HLfit"))){
        if(!inherits(preds, 'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- curve_inverse(as.data.frame(preds)[, 1])
          lwr  <- upr <- NA
        } else {
          fit_ <- curve_inverse(preds[[1]])
          bounds <- .msr_native_link_ci(analysis_mod, preds[[1]], preds[[2]],
            curve_inverse)
          lwr <- bounds$lower
          upr <- bounds$upper
        }

      } else {
        if(!inherits(preds, 'list') || (is.list(preds) && length(preds) == 1)){
          fit_ <- as.data.frame(preds)[[1]]
          lwr  <- upr <- NA
        } else {
          fit_ <- preds[[1]]
          lwr  <- preds[[1]] + qnorm(0.025) * preds[[2]]
          upr  <- preds[[1]] + qnorm(0.975) * preds[[2]]
        }
      }

      data.frame(x = x_unscaled, fit = fit_, lwr = lwr, upr = upr, variable = v)
    })
  }

  names(plot_list) <- c_vars

  # Build ggplots, adding subtitles for polynomial and interaction terms
  gg_list <- lapply(c_vars, function(v) {
    df <- plot_list[[v]]

    subtitle_parts <- character(0)

    # Note any I() terms involving this variable
    if(v %in% poly_vars) {
      rel_poly <- poly_terms[vapply(poly_variables, function(names) v %in% names, logical(1))]
      subtitle_parts <- c(subtitle_parts,
                          paste0("Includes: ", paste(rel_poly, collapse = ", ")))
    }

    # Note the mean-conditioning for interaction terms
    if(v %in% inter_vars) {
      partners <- inter_variables[vapply(inter_variables,
        function(names) v %in% names, logical(1))]
      others <- setdiff(unique(unlist(partners, use.names = FALSE)), v)
      subtitle_parts <- c(subtitle_parts,
                          paste0("At mean of: ", paste(others, collapse = ", ")))
    }

    subtitle <- if(length(subtitle_parts) > 0) paste(subtitle_parts, collapse = "; ") else NULL

    p <- ggplot(df, aes(x = x, y = fit)) +
      xlab(df$variable[1]) +
      ylab(ylab) +
      theme_cowplot()
    if (all(is.na(df$lwr)) && all(is.na(df$upr)))
      p <- p + labs(caption = paste(c(
        if (trained_unmarked_curve || native_mixed_response || inherits(analysis_mod, "HLfit"))
          "Random effects are set to zero, not integrated.",
        "Confidence limits unavailable; point predictions condition on spatial scales."), collapse = " "))
    else if (trained_unmarked_curve)
      p <- p + labs(caption = paste(
        "Random effects are set to zero, not integrated.",
        "Native component confidence limits condition on the selected spatial scales."))
    else if (native_mixed_response)
      p <- p + labs(caption = paste(
        "Random effects are set to zero, not integrated.",
        "Intervals use native response-scale standard errors and condition on spatial scales."))
    else if (inherits(analysis_mod, "HLfit"))
      p <- p + labs(caption = paste(
        "Random effects are set to zero, not integrated.",
        "Limits transform native link-scale fixed-effect uncertainty and condition on spatial scales."))
    else if (!all(is.na(df$lwr)))
      p <- p + labs(caption = "Confidence limits condition on the selected spatial scales.")

    if (native_gam_components)
      p <- p + labs(caption = paste("Native GAM", if (link) "response column" else "linear predictor",
        component, if (link && .msr_gam_ziplss_model(analysis_mod))
          ". Limits propagate full coefficient covariance through both predictors;" else
          ". Limits use native standard errors on this output scale;",
        "kernel scales and smoothing parameters are fixed at their fitted values."))

    discrete <- !is.numeric(df$x)
    if (discrete) {
      p <- p + ggplot2::geom_point()
    } else {
      p <- p + geom_line(linewidth = 1)
    }

    if(!is.null(subtitle))
      p <- p + labs(subtitle = subtitle)

    if(!all(is.na(df$lwr)) && !all(is.na(df$upr))) {
      if (discrete) {
        p <- p + ggplot2::geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15)
      } else {
        p <- p + geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25)
      }
    }

    p
  })

  names(gg_list) <- c_vars

  # Explanatory captions are opt-in; the subtitle for interaction terms stays.
  if (!add_label)
    gg_list <- lapply(gg_list, function(p) p + labs(caption = NULL))

  lapply(gg_list, print)
  if (!is.null(reference)) attr(gg_list, "reference_data_contract") <- reference$contract
  if (native_gam_components) {
    contract <- .msr_native_prediction_contract(x, "gam", if (link) "response" else "link",
      list(), se.fit = TRUE)
    contract$native_output <- list(type = if (link) "response" else "link", column = component)
    attr(gg_list, "prediction_contract") <- contract
  }
  if (trained_unmarked_curve)
    attr(gg_list, "prediction_contract") <- .msr_native_prediction_contract(x,
      "unmarked", type, list(re.form = NA), se.fit = TRUE)
  if (native_mixed_response)
    attr(gg_list, "prediction_contract") <- .msr_native_prediction_contract(x,
      if (inherits(analysis_mod, "glmmTMB")) "glmmTMB" else "merMod",
      "response", list(re.form = NA), se.fit = TRUE)
  invisible(gg_list)
}


# Recover the state/abundance formula without depending on coefficient names
# or unmarked's private formula helpers. Older fits retain a double formula;
# newer classes may retain named component formulas in a formlist slot.
.msr_unmarked_marginal_formula <- function(mod, type) {
  slots <- methods::slotNames(mod)
  if ("formlist" %in% slots && inherits(mod@formlist[[type]], "formula")) {
    return(mod@formlist[[type]])
  }
  if ("formula" %in% slots && inherits(mod@formula, "formula")) {
    f <- mod@formula
    if (length(f) == 3L) {
      return(stats::as.formula(call("~", f[[3]]), env = environment(f)))
    }
    if (length(f) == 2L) return(f)
  }
  stop("Could not recover the '", type, "' formula for this unmarked model.",
       call. = FALSE)
}


# Formula variables are data columns; factor contrasts are coefficient names.
# Grouping variables describe random effects, and offset-only variables need
# prediction defaults but do not describe a separate effect to plot.
# A quoted name may contain a pipe. Logical pipes inside I() or another
# fixed function are not formula random-effect terms.
.msr_marginal_random_call <- function(expr) {
  while (is.call(expr) && identical(expr[[1L]], as.name("(")) && length(expr) == 2L)
    expr <- expr[[2L]]
  is.call(expr) && is.symbol(expr[[1L]]) &&
    as.character(expr[[1L]]) %in% c("|", "||")
}

.msr_marginal_variables <- function(expr, offsets) {
  if (inherits(expr, c("terms", "formula"))) {
    tt <- stats::terms(expr)
    expressions <- lapply(attr(tt, "term.labels"), str2lang)
    if (offsets) {
      variables <- as.list(attr(tt, "variables"))[-1L]
      expressions <- c(expressions, variables[attr(tt, "offset")])
    }
    return(as.character(unique(unlist(lapply(expressions,
      .msr_marginal_variables, offsets = offsets), use.names = FALSE))))
  }
  if (!is.call(expr)) return(all.vars(expr))
  fun <- if (is.symbol(expr[[1L]])) as.character(expr[[1L]]) else ""
  if (.msr_marginal_random_call(expr) || (!offsets && fun == "offset"))
    return(character())
  # Only formula operators can introduce a random-effect term. Preserve every
  # raw input in a fixed function, including logical |/|| inside I().
  if (!fun %in% c("~", "+", "-", "*", "/", ":", "^", "%in%", "("))
    return(all.vars(expr))
  as.character(unique(unlist(lapply(as.list(expr)[-1L], .msr_marginal_variables,
    offsets = offsets), use.names = FALSE)))
}


# Reconstruct the native fixed training frame before site removal. Store only
# its learned prediction terms; never relearn a basis on a hypothetical grid.
.msr_unmarked_curve_basis <- function(component, data) {
  tt <- stats::terms(component, data = data)
  labels <- attr(tt, "term.labels")
  keep <- !vapply(lapply(labels, str2lang), .msr_marginal_random_call, logical(1))
  variables <- as.list(attr(tt, "variables"))[-1L]
  offsets <- vapply(variables[attr(tt, "offset")], deparse1, character(1), backtick = TRUE)
  labels <- c(labels[keep], offsets)
  intercept <- attr(tt, "intercept") == 1L
  form <- if (length(labels)) stats::reformulate(labels, intercept = intercept,
    env = attr(tt, ".Environment")) else stats::as.formula(if (intercept) "~1" else "~0",
      env = attr(tt, ".Environment"))
  frame <- stats::model.frame(form, data, na.action = stats::na.pass,
    drop.unused.levels = FALSE)
  fixed <- .msr_fixed_prediction_terms(attr(frame, "terms"), frame)
  fixed <- .msr_cv_prediction_terms(fixed, data)
  list(terms = fixed,
    learned = !identical(attr(fixed, "predvars"), attr(fixed, "variables")))
}


# Preserve explicit contrasts when native data-frame prediction drops them.
# Restrict the fallback to models whose state is a single linked estimate;
# derived-state predictions (for example ZIP abundance) need their own method.
.msr_unmarked_contrast_prediction <- function(mod, component, newdata, type) {
  supported <- identical(class(mod)[1], "unmarkedFitOccu") ||
    (identical(class(mod)[1], "unmarkedFitPCount") && mod@mixture %in% c("P", "NB"))
  if (!supported) {
    stop("Explicit nondefault factor contrasts cannot safely be predicted for ",
         "this unmarked class by plot_marginal_effects(). Use its fitted design ",
         "matrix and class-specific prediction method.", call. = FALSE)
  }
  terms <- stats::terms(component)
  random <- which(vapply(lapply(attr(terms, "term.labels"), str2lang),
    .msr_marginal_random_call, logical(1)))
  if (length(random)) terms <- stats::drop.terms(terms, random, keep.response = FALSE)
  frame <- stats::model.frame(terms, newdata, na.action = stats::na.pass)
  X <- stats::model.matrix(terms, frame)
  offset <- stats::model.offset(frame)
  coefficient_names <- names(unmarked::coef(mod, type = type, altNames = FALSE))
  if (!identical(sub("^\\(Intercept\\)$", "Int", colnames(X)), coefficient_names)) {
    stop("The prediction design does not match the fitted unmarked coefficients.",
         call. = FALSE)
  }
  # One row at a time avoids constructing a dense covariance matrix over the
  # whole plotting grid. Confidence bounds retain the fitted link function.
  rows <- lapply(seq_len(nrow(X)), function(i) {
    lc <- unmarked::linearComb(mod, X[i, , drop = FALSE], type = type,
      offset = if (is.null(offset)) NULL else offset[i], re.form = NA)
    prediction <- unmarked::backTransform(lc)
    ci <- unmarked::confint(prediction, level = 0.95)
    data.frame(Predicted = as.numeric(unmarked::coef(prediction)),
      SE = as.numeric(unmarked::SE(prediction)), lower = ci[1], upper = ci[2])
  })
  do.call(rbind, rows)
}


safe_predict <- function(mod, newdata) {
  pkg <- extract_namespace(mod)

  if(is.null(pkg)){
    mod_classes <- class(mod)

    # Try S3 generic method (e.g., predict.lm)
    pkg <- NULL
    for (cls in mod_classes) {
      mthd <- tryCatch(
        getS3method("predict", cls, optional = TRUE),
        error = function(e) NULL
      )
      if (!is.null(mthd)) {
        env <- environment(mthd)
        if (!is.null(env) && !is.null(pkgname <- environmentName(env)) &&
            pkgname %in% loadedNamespaces()) {
          pkg <- pkgname
          break
        }
      }
    }

    # Try S4 class definition
    if (is.null(pkg) && isS4(mod)) {
      for (cls in mod_classes) {
        s4class <- tryCatch(getClassDef(cls), error = function(e) NULL)
        if (!is.null(s4class) && !is.null(pkgname <- s4class@package)) {
          pkg <- pkgname
          break
        }
      }
    }

    if (is.null(pkg)) {
      stop("Could not determine the package for the model object.")
    }
  }

  # Load namespace if needed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  if (!paste0("package:", pkg) %in% search()) {
    .msr_attach_package(pkg)
  }

  # Dispatch predict
  predict(mod, newdata, se.fit = TRUE)
}

namespace <- function(x) {

  pkg <- extract_namespace(x)

  if (is.null(pkg)) {
    mod_classes <- class(x)

    # Try S3 generic method (e.g., predict.lm)
    pkg <- NULL
    for (cls in mod_classes) {
      mthd <- tryCatch(
        getS3method("predict", cls, optional = TRUE),
        error = function(e) NULL
      )
      if (!is.null(mthd)) {
        env <- environment(mthd)
        if (!is.null(env) && !is.null(pkgname <- environmentName(env)) &&
            pkgname %in% loadedNamespaces()) {
          pkg <- pkgname
          break
        }
      }
    }

    # Try S4 class definition
    if (is.null(pkg) && isS4(x)) {
      for (cls in mod_classes) {
        s4class <- tryCatch(getClassDef(cls), error = function(e) NULL)
        if (!is.null(s4class) && !is.null(pkgname <- s4class@package)) {
          pkg <- pkgname
          break
        }
      }
    }

    if (is.null(pkg)) {
      stop("Could not determine the package for the model object.")
    }
  }

  # Check if namespace is already loaded
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  if (!paste0("package:", pkg) %in% search()) {
    .msr_attach_package(pkg)
  }

  invisible(pkg)
}

.msr_attach_package <- function(pkg) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  )
}

utils::globalVariables(c("fit", "lwr", "upr"))


# Transform native mean limits only on a valid branch. Sorting handles a
# decreasing inverse link; it cannot repair a pole, invalid domain or bad SE.
.msr_native_link_ci <- function(model, fit, se, inverse) {
  eta_lower <- fit + qnorm(.025) * se
  eta_upper <- fit + qnorm(.975) * se
  valid <- is.finite(fit) & is.finite(se) & se >= 0 &
    is.finite(eta_lower) & is.finite(eta_upper)
  family <- if (inherits(model, "coxph")) NULL else
    tryCatch(stats::family(model), error = function(e) NULL)
  if (is.function(family$valideta)) for (i in which(valid))
    valid[i] <- isTRUE(family$valideta(c(eta_lower[i], eta_upper[i])))
  if (isTRUE(family$link %in% c("inverse", "1/mu^2")))
    valid <- valid & !(eta_lower <= 0 & eta_upper >= 0)
  if (identical(family$link, "sqrt")) valid <- valid & eta_lower >= 0
  lower <- upper <- rep(NA_real_, length(fit))
  index <- which(valid)
  if (length(index)) {
    a <- as.numeric(inverse(eta_lower[index]))
    b <- as.numeric(inverse(eta_upper[index]))
    acceptable <- is.finite(a) & is.finite(b)
    if (is.function(family$validmu)) for (j in which(acceptable))
      acceptable[j] <- isTRUE(family$validmu(c(a[j], b[j])))
    valid[index[!acceptable]] <- FALSE
    lower[index[acceptable]] <- pmin(a[acceptable], b[acceptable])
    upper[index[acceptable]] <- pmax(a[acceptable], b[acceptable])
  }
  if (any(!valid)) warning(sprintf(
    "Native mean confidence limits were not shown for %d row(s): nonfinite SEs or limits outside the native link/family domain.",
    sum(!valid)), call. = FALSE)
  list(lower = lower, upper = upper)
}

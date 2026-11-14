# Record outcomes only for matching future comparisons, after fitting/scoring.
# This metadata never enters training, normalization, bins or predictions.
.msr_cv_targets <- function(formula, data, family, weights, folds) {
  used <- unique(unlist(lapply(folds, function(s) c(s$analysis, s$assessment)), use.names = FALSE))
  rows <- match(sort(used), rownames(data))
  tryCatch({
    y <- .msr_cv_outcome(formula, data[rows, , drop = FALSE], family)
    data.frame(row_id = rownames(data)[rows], observed = y$observed,
      trials = y$trials, weight = as.numeric(weights[rows]), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}

.msr_cv_finish <- function(output, folds, ids, family, metrics) {
  for (item in c("predictions", "metrics", "fold_results")) {
    output[[item]] <- do.call(rbind, output[[item]])
    rownames(output[[item]]) <- NULL
  }
  output$complete <- all(output$fold_results$status == "scored")
  output$summary <- do.call(rbind, lapply(metrics, function(metric)
    .msr_cv_metric(output$predictions, metric, family, output$complete, nrow(output$predictions))))
  output$folds <- folds
  output$assessment_counts <- stats::setNames(tabulate(match(output$predictions$row_id, ids), length(ids)), ids)
  class(output) <- if (identical(family, "conditional_choice")) c("multiScaleR_choice_cv", "multiScaleR_cv") else "multiScaleR_cv"
  output
}

#' Cross-Validate an Ordinary GLM Reference Model
#'
#' Fit a GLM to the analysis rows of each fold and score its predictions for
#' the held-out rows. No scales of effect are estimated. Use this function to
#' build null and nonspatial reference models on the same folds you give
#' [cross_validate()], so that a scale-of-effect model has a baseline to beat.
#'
#' @inheritParams cross_validate
#' @param formula Two-sided model formula, written as for [stats::glm()]. The
#'   response must be a column of `data`, optionally wrapped in `I()`, or
#'   `cbind(successes, failures)` for grouped binomial data. Predictors can
#'   include factors, interactions, offsets, and the data-dependent terms
#'   `scale()`, `poly()`, `splines::ns()`, and `splines::bs()`. Use `y ~ 1`
#'   for a null model. Every variable must be a column of `data`.
#' @param data Data frame with one row per observation, holding the raw
#'   responses and every predictor. Its row names are the row IDs that the
#'   folds refer to, and they must be unique. Do not supply predictors
#'   standardized or transformed with all rows (see "What each fold fits").
#' @param glm_control Named list of [stats::glm.control()] settings for the
#'   GLM fitted in each fold. The default keeps the settings of `glm()`.
#' @param theta Optional NB2 size parameter fixed in advance, one positive
#'   number. The default `NULL` estimates theta in each fold with
#'   [MASS::glm.nb()]. A fixed theta is fitted with
#'   `glm(family = MASS::negative.binomial(theta))` and adds no parameter to
#'   K. Fix it only from information independent of the assessment outcomes.
#'   Other families reject it.
#' @return A list of class `multiScaleR_cv`, with class
#'   `multiScaleR_distribution_cv` added for `"nbinom2"` and `"Gamma"`. It has
#'   the elements of a [cross_validate()] result, with these differences.
#'   `settings$backend` is `"glm"`. `training` holds only the coefficients of
#'   each fold, plus the frozen nuisance parameter for NB2 and Gamma.
#'   `diagnostics` holds the convergence record of each fold's GLM, and `fits`
#'   holds ordinary `glm` objects when `keep_fits = TRUE`. See
#'   [multiScaleR_cv-class] for every element and how to read each score.
#' @details
#' # What each fold fits
#'
#' Each fold fits `glm()` to its analysis rows with the declared family and
#' frequency weights, then predicts the expected response of the assessment
#' rows. The families, response forms, weights, and scores are those of
#' [cross_validate()]. Factors, interactions, and offsets use the native GLM
#' machinery. The constants of `scale()`, `poly()`, `splines::ns()`, and
#' `splines::bs()` come from the analysis rows, even inside expressions and
#' offsets. Any other preprocessing that learns from the data must also
#' happen within the analysis rows of each fold, so do not supply predictors
#' standardized with all rows. The section "Formula terms that learn from the
#' data" in [cross_validate()] lists the supported terms and their limits.
#' Missing values are never dropped silently, and native warnings and failed
#' folds are kept.
#'
#' # Binomial responses
#'
#' A binomial response can be numeric or logical 0/1 values, a two-level
#' factor (the second level counts as a success), or
#' `cbind(successes, failures)` for several trials per row. Keep trial counts
#' in the response and case frequencies in `weights`, because the two are
#' scored differently.
#'
#' * The negative log score of a row is the negative log binomial probability
#'   of its observed number of successes, given its trials. Rows are averaged
#'   with their frequency weights.
#' * The Brier score of a row is the squared error of the predicted success
#'   probability, averaged over the row's trials. Rows are averaged with
#'   weights equal to frequency weight x trials, so the pooled Brier score is
#'   a mean per trial. A row with many trials counts more in the Brier score
#'   than in the log score.
#' * A row with zero trials or zero weight stays in `predictions` but
#'   contributes to neither score.
#' * An observed outcome given zero probability produces an infinite log
#'   score. It is kept as is, not clipped.
#'
#' # Negative binomial and Gamma responses
#'
#' `"nbinom2"` fits [MASS::glm.nb()] when theta is estimated and
#' `glm(family = MASS::negative.binomial(theta))` when it is fixed, so it
#' needs the MASS package. `"Gamma"` uses a log link. In each fold the
#' nuisance parameter (theta, or the Gamma dispersion) is estimated from the
#' analysis rows, frozen, and used to score the assessment rows. Coefficients
#' and their native standard errors are those of the native fit. A fixed
#' Gamma `dispersion` is built into a normalized Gamma family before fitting.
#' The section "Negative binomial and Gamma responses" in [cross_validate()]
#' describes how these parameters enter K and the rules they follow.
#'
#' ## How to read the Gamma dispersion
#'
#' `distribution[[fold]]$dispersion` is the training deviance divided by the
#' training sum of case frequencies, the value native Gamma AIC uses. It is
#' held fixed when the held-out Gamma densities are scored. Native
#' coefficient standard errors use the Pearson dispersion instead, so they
#' can rest on a different number. The density dispersion does not replace
#' that coefficient uncertainty. `convention` records the rule that produced
#' the value, and `estimated` is `FALSE` when you fixed it with `dispersion`.
#' Supply `theta` or `dispersion` only when the value is known or was fixed
#' independently of the assessment outcomes.
#'
#' # What the scores do and do not include
#'
#' Predictions use each fold's fitted parameters as if they were known. They
#' do not integrate parameter uncertainty, provide prediction intervals, or
#' establish calibration. Use [cross_validate()] to estimate scales of effect
#' within the same folds.
#' @examples
#' data("landscape_counts")   # counts and a site-level covariate
#' data("surv_pts")           # survey locations
#' # Three spatial folds. The simulated coordinates are in meters.
#' folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
#'                           crs = 32617)
#'
#' # Null and site-only reference models on the same folds
#' null_cv <- cross_validate_glm(counts ~ 1, landscape_counts, folds,
#'                               family = "poisson")
#' site_cv <- cross_validate_glm(counts ~ site, landscape_counts, folds,
#'                               family = "poisson")
#' site_cv   # pooled held-out scores, lower is better
#' cv_tab(list(null = null_cv, site = site_cv), metric = "neg_log_score")
#'
#' # Negative binomial: theta is estimated in each fold and frozen for scoring
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   nb_cv <- cross_validate_glm(counts ~ site, landscape_counts, folds,
#'                               family = "nbinom2")
#'   sapply(nb_cv$distribution, function(d) d$theta)   # small theta: overdispersed
#' }
#'
#' # Gamma model for a strictly positive response
#' i <- 1:96
#' positive <- data.frame(
#'   y = exp(0.2 + 0.35 * sin(i * 0.73)) *
#'     rep(c(0.3, 0.6, 1.1, 2.4, 0.8, 1.8, 0.4, 1.6), 12),
#'   x = sin(i * 0.73))
#' gamma_folds <- grouped_cv_folds(rep(letters[1:12], each = 8), k = 3)
#' gamma_cv <- cross_validate_glm(y ~ x, positive, gamma_folds, family = "Gamma")
#' gamma_cv$summary
#' # Training deviance / sum of case frequencies, frozen for scoring each fold
#' sapply(gamma_cv$distribution, function(d) d$dispersion)
#' @seealso [cross_validate()], [cv_tab()], [select_cv()],
#'   [multiScaleR_cv-class], and
#'   `vignette("v08-cross-validation", package = "multiScaleR")`.
#' @export
cross_validate_glm <- function(formula, data, folds,
    family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"), weights = NULL,
    glm_control = stats::glm.control(), metrics = NULL, keep_fits = FALSE,
    keep_selection = FALSE, theta = NULL, dispersion = NULL) {
  if (identical(family, "gamma")) family <- "Gamma"
  family <- match.arg(family)
  if (.msr_cv_distribution_family(family)) {
    out <- .msr_cv_distribution(formula, data, folds, family, weights, glm_control, metrics,
      keep_fits, keep_selection, theta, dispersion)
    out$call <- match.call(); return(out)
  }
  if (!is.null(theta) || !is.null(dispersion)) stop("theta/dispersion controls apply only to NB2/Gamma CV.", call. = FALSE)
  validate_scalar_logical(keep_fits, "keep_fits")
  validate_scalar_logical(keep_selection, "keep_selection")
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      length(setdiff(all.vars(formula), c(names(data), "."))))
    stop("Supply a two-sided GLM formula and nonempty data with every formula variable and unique columns.", call. = FALSE)
  ids <- .msr_fold_ids(rownames(data), nrow(data))
  folds <- .msr_cv_folds(folds, ids)
  defaults <- list(gaussian = c("rmse", "mae"),
    poisson = c("rmse", "mae", "neg_log_score", "deviance"), binomial = c("neg_log_score", "brier"))
  if (is.null(metrics)) metrics <- defaults[[family]]
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) ||
      anyDuplicated(metrics) || any(!metrics %in% defaults[[family]]))
    stop("metrics must uniquely name scores supported for this family.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (!is.numeric(weights) || !is.null(dim(weights)) || length(weights) != nrow(data) ||
      any(!is.finite(weights) | weights < 0 | weights != floor(weights)) ||
      (!is.null(names(weights)) && !identical(names(weights), ids)) ||
      (family == "gaussian" && any(weights != 1)))
    stop("weights must align with data and be nonnegative integer frequencies; Gaussian weights must equal one.", call. = FALSE)
  if (!is.list(glm_control) || is.null(names(glm_control)) || anyDuplicated(names(glm_control)) ||
      any(!names(glm_control) %in% names(formals(stats::glm.control))))
    stop("glm_control must be a named glm.control() list.", call. = FALSE)
  glm_control <- do.call(stats::glm.control, glm_control)
  fam <- switch(family, gaussian = stats::gaussian(), poisson = stats::poisson(), binomial = stats::binomial())
  output <- list(predictions = list(), metrics = list(), fold_results = list(), training = list(),
    diagnostics = list(), warnings = list(), errors = list(), prediction_contracts = list(),
    fits = if (keep_fits) list() else NULL,
    selection_records = if (keep_selection) list() else NULL)
  for (name in names(folds)) {
    split <- folds[[name]]; ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    pred <- data.frame(fold = name, row_id = split$assessment, observed = NA_real_, trials = NA_real_,
      prediction = NA_real_, weight = unname(weights[ti]), eligible = FALSE, stringsAsFactors = FALSE)
    for (metric in metrics) pred[[metric]] <- NA_real_
    fit <- NULL; messages <- character(); stage <- "preparation"
    error <- tryCatch(withCallingHandlers({
      outcome <- .msr_cv_outcome(formula, data[ti, , drop = FALSE], family)
      pred$observed <- outcome$observed; pred$trials <- outcome$trials
      .msr_cv_outcome(formula, data[ai, , drop = FALSE], family)
      stage <- "training_glm"
      fit <- do.call(stats::glm, list(formula = formula, data = data[ai, , drop = FALSE],
        family = fam, weights = unname(weights[ai]), control = glm_control,
        na.action = stats::na.fail, model = TRUE))
      output$training[[name]] <- list(coefficients = stats::coef(fit))
      if (!isTRUE(fit$converged)) stop("Training GLM did not converge.", call. = FALSE)
      if (is.factor(stats::model.response(fit$model)) && length(unique(stats::model.response(fit$model))) != 2L)
        stop("A training binomial factor must contain both declared response levels.", call. = FALSE)
      stage <- "prediction"
      fit <- .msr_cv_prediction_metadata(fit, data[ai, , drop = FALSE])
      pred$prediction <- as.numeric(stats::predict(fit, newdata = data[ti, , drop = FALSE], type = "response"))
      stage <- "scoring"
      pred <- .msr_cv_score(pred, family, metrics)
      NULL
    }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(stage = stage, message = conditionMessage(e)))
    complete <- is.null(error)
    if (!complete) {
      pred$prediction <- NA_real_; pred$eligible <- FALSE
      for (metric in metrics) pred[[metric]] <- NA_real_
    }
    output$predictions[[name]] <- pred
    output$metrics[[name]] <- do.call(rbind, lapply(metrics, function(metric)
      cbind(fold = name, .msr_cv_metric(pred, metric, family, complete, nrow(pred)))))
    output$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai), n_assessment = length(ti),
      status = if (complete) "scored" else "failed", stringsAsFactors = FALSE)
    output$warnings[name] <- list(messages); output$errors[name] <- list(error)
    output$diagnostics[name] <- list(list(fitted_glm = if (is.null(fit)) NULL else fit[c("converged", "boundary", "iter")]))
    output$prediction_contracts[name] <- list(.msr_prediction_contract(
      "likelihood", "glm", "expected_response", "response",
      population = "fixed_at_estimate", kernel = "not_applicable",
      conditioning = "training_fit_and_fixed_prediction_rows",
      omitted = "parameter_uncertainty", fixed_inputs = TRUE))
    if (keep_selection) output$selection_records[name] <- list(.msr_cv_selection_record(fit,
      output$diagnostics[[name]], .msr_cv_selection_controls(list(backend = "glm", glm_control = glm_control))))
    if (keep_fits) output$fits[name] <- list(fit)
  }
  output <- .msr_cv_finish(output, folds, ids, family, metrics)
  output$targets <- .msr_cv_targets(formula, data, family, weights, folds)
  output$settings <- list(family = family, link = fam$link, backend = "glm", metrics = metrics,
    glm_control = glm_control, normalization = "training_only", prediction_target = "held_out_GLM_expected_response",
    parameter_uncertainty = "not_integrated", observation_intervals = "not_provided", calibration = "not_established", workers = 1L)
  output$call <- match.call()
  output
}

# Comparisons bind observation definitions and analysis/assessment membership.
# Old objects are usable only if their prediction records identify every used
# target, including analysis-only rows. Hidden training responses are not guessed.
.msr_cv_comparison <- function(mod_list, metric = NULL) {
  if (!is.list(mod_list) || !length(mod_list) || is.null(names(mod_list)) ||
      anyNA(names(mod_list)) || any(!nzchar(names(mod_list))) || anyDuplicated(names(mod_list)))
    stop("mod_list must be a nonempty, uniquely named list of CV results.", call. = FALSE)
  if (any(vapply(mod_list, inherits, logical(1), "multiScaleR_choice_cv")))
    return(.msr_choice_comparison(mod_list, metric))
  if (any(vapply(mod_list, inherits, logical(1), "multiScaleR_distribution_cv")))
    return(.msr_cv_distribution_comparison(mod_list, metric))
  signatures <- lapply(mod_list, function(x) {
    if (!inherits(x, "multiScaleR_cv") || !isTRUE(x$complete) ||
        !is.data.frame(x$fold_results) || any(x$fold_results$status != "scored"))
      stop("Every candidate must be a complete multiScaleR_cv result; failed folds cannot be dropped.", call. = FALSE)
    s <- x$settings
    if (!s$family %in% c("gaussian", "poisson", "binomial") ||
        !identical(s$prediction_target, "held_out_GLM_expected_response") ||
        !identical(s$parameter_uncertainty, "not_integrated") ||
        !identical(s$observation_intervals, "not_provided"))
      stop("Candidate prediction targets or uncertainty contracts are not qualified for comparison.", call. = FALSE)
    if (!is.null(metric) && !metric %in% s$metrics)
      stop("Every candidate must record the requested metric.", call. = FALSE)
    used <- sort(unique(unlist(lapply(x$folds, function(f) c(f$analysis, f$assessment)), use.names = FALSE)))
    .msr_cv_folds(x$folds, used)
    folds <- lapply(x$folds[sort(names(x$folds))], function(f) lapply(f[c("analysis", "assessment")], sort))
    attributes(folds) <- list(names = names(folds))
    p <- x$predictions
    if (!is.data.frame(p) || !all(c("fold", "row_id", "observed", "trials", "weight", "prediction", "eligible") %in% names(p)))
      stop("Missing recorded assessment targets.", call. = FALSE)
    for (name in names(folds)) {
      rows <- p[p$fold == name, , drop = FALSE]
      if (anyNA(rows$row_id) || anyDuplicated(rows$row_id) || !setequal(rows$row_id, folds[[name]]$assessment))
        stop("Recorded assessment rows do not match the split plan.", call. = FALSE)
    }
    if (anyNA(p$fold) || any(!p$fold %in% names(folds))) stop("Unknown recorded folds.", call. = FALSE)
    targets <- x$targets
    if (is.null(targets)) {
      targets <- unique(p[c("row_id", "observed", "trials", "weight")])
      if (anyDuplicated(targets$row_id) || !setequal(targets$row_id, used))
        stop("Legacy CV targets do not identify every training/assessment row; rerun cross-validation for comparison.", call. = FALSE)
    }
    columns <- c("row_id", "observed", "trials", "weight")
    if (!is.data.frame(targets) || !all(columns %in% names(targets)) ||
        anyDuplicated(targets$row_id) || !setequal(targets$row_id, used))
      stop("Invalid recorded target identities.", call. = FALSE)
    targets <- targets[match(used, targets$row_id), columns, drop = FALSE]; rownames(targets) <- NULL
    for (nm in columns[-1L]) {
      if (!is.numeric(targets[[nm]]) || any(!is.finite(targets[[nm]])))
        stop("Targets must record finite numeric responses, trials and frequency weights.", call. = FALSE)
      targets[[nm]] <- as.numeric(targets[[nm]])
      if (!identical(as.numeric(p[[nm]]), targets[[nm]][match(p$row_id, targets$row_id)]))
        stop("Observed responses, trial counts and frequency weights must match recorded targets.", call. = FALSE)
    }
    if (any(targets$weight < 0 | targets$weight != floor(targets$weight)) ||
        (s$family == "gaussian" && any(targets$weight != 1)) ||
        any(targets$trials < 0 | targets$trials != floor(targets$trials)) ||
        (s$family != "binomial" && any(targets$trials != 1)) ||
        (s$family == "binomial" && any(targets$observed < 0 | targets$observed > targets$trials | targets$observed != floor(targets$observed))) ||
        (s$family == "poisson" && any(targets$observed < 0 | targets$observed != floor(targets$observed))))
      stop("Targets must preserve valid outcomes, trials and integer frequency weights.", call. = FALSE)
    check <- .msr_cv_score(p, s$family, if (is.null(metric)) s$metrics else metric)
    if (!identical(check$eligible, p$eligible)) stop("Recorded score eligibility is inconsistent.", call. = FALSE)
    list(folds = folds, targets = targets, family = s$family, link = s$link)
  })
  first <- signatures[[1L]]
  for (s in signatures) {
    if (!identical(first$folds, s$folds) || !identical(first$targets, s$targets))
      stop("CV comparison requires identical splits, observed responses, trial counts and weights.", call. = FALSE)
    if ((is.null(metric) || !metric %in% c("rmse", "mae")) &&
        (!identical(first$family, s$family) || !identical(first$link, s$link)))
      stop("This comparison requires the same response family and link.", call. = FALSE)
  }
  invisible(first)
}

#' Compare Cross-Validated Prediction Scores
#'
#' Rank candidate models by one held-out score, after checking that every
#' candidate was scored on the same folds and the same observations. Lower
#' scores are better.
#'
#' @param mod_list Named list of complete cross-validation results to compare,
#'   for example `list(landscape = cv, null = null_cv)`. The names label the
#'   rows of the table and must be unique. Every candidate must use identical
#'   folds, responses, trials, and frequency weights, and the candidates must
#'   be of compatible kinds (see "What can be compared").
#' @param metric One score name, chosen explicitly: `"rmse"`, `"mae"`,
#'   `"neg_log_score"`, `"deviance"`, or `"brier"`. Every candidate must have
#'   recorded it. There is no default, because the scores answer different
#'   questions.
#' @param by_fold Logical. The default `FALSE` ranks the pooled scores. `TRUE`
#'   ranks the candidates separately within each fold, which shows whether a
#'   ranking is consistent across folds.
#' @return A data frame with one row per candidate (per candidate and fold
#'   when `by_fold = TRUE`), sorted from best to worst. Its columns are:
#'   * `Modnames`: the candidate's name in `mod_list`.
#'   * `Rank`: 1 for the best score. Tied scores share the lowest rank.
#'   * `Score`: the pooled score, or the fold score when `by_fold = TRUE`.
#'   * `Delta`: `Score` minus the best score, in the units of the score. It is
#'     `NA` when every score is infinite.
#'   * `Metric`, `Unit`: the score and what one score refers to, for example
#'     `"observation"`, `"trial"`, or `"choice_set"`. See
#'     [multiScaleR_cv-class].
#'   * `N_assessment`, `N_scored`, `Weight_sum`: the planned assessment
#'     occurrences, the number that entered the score, and their total
#'     weight.
#'   * `Warnings`: the number of warnings recorded while fitting the
#'     candidate (in that fold when `by_fold = TRUE`). Read them in the
#'     candidate's `warnings` element.
#'   * `Fold`: the fold, when `by_fold = TRUE`.
#'
#'   The candidate results themselves are not modified.
#' @details
#' # How to read the table
#'
#' The best candidate is listed first. `Delta` is each candidate's distance
#' from it, in the units of the score. For RMSE and MAE that is a difference
#' in response units, such as animals per survey. For the negative log score,
#' a `Delta` of D means the best model gave the observed outcomes exp(D)
#' times the probability, as a geometric mean per observation (per choice set
#' for conditional choice). The per-fold table from `by_fold = TRUE` is
#' descriptive. Use it to see whether one fold drives a ranking, never to pick
#' a model for that fold from its own assessment score.
#'
#' # What can be compared
#'
#' Every candidate must be complete, with no failed fold. An incomplete
#' candidate stops the comparison, because a model scored only on its easier
#' folds would look better than it is. No failed model or fold is removed
#' silently. Results of different kinds cannot be mixed.
#'
#' | Results | Comparisons allowed |
#' |:---|:---|
#' | Gaussian, Poisson, and binomial from [cross_validate()] and [cross_validate_glm()] | RMSE and MAE across families on identical observations. The log score, deviance, and Brier score need the same family and link. |
#' | NB2 or Gamma, alone or with Gaussian or Poisson | RMSE and MAE across these families on identical observations. Log scores between Poisson and NB2 on the same counts. A Gamma log score compares only with another Gamma. |
#' | [cross_validate_clogit()] | Other conditional-choice results with the same complete sets and partial-likelihood `method`, using `neg_log_score` only. |
#' | [cross_validate_field()] | Other field results with the same responses, folds, and projected coordinates. Poisson and NB2 share a count target. |
#' | [cross_validate_mixed()] | Other mixed-model results with the same groups, prediction target, and folds. |
#' | [cross_validate_bayes()] | Other Bayesian results with the same targets, folds, and blocks, using `neg_log_score` only. |
#'
#' A binomial log score refers to the probability of the observed count of
#' successes, while a Brier score is a mean over individual trials. Field,
#' mixed-model, and Bayesian log scores are joint scores for blocks of
#' held-out observations, so they are kept apart from the per-observation
#' scores of ordinary GLMs.
#'
#' # What the table does not do
#'
#' The scores carry no penalty for model complexity, and `Delta` is not an
#' AIC difference, an Akaike weight, or a model probability. The spread of
#' fold scores is not a standard error, because folds share training data
#' and are few. If you use the table to choose a model or a tuning rule, the
#' winning score is optimistic. Evaluate that choice with an independent test
#' set or an outer validation layer. To score a predeclared rule that picks a
#' model by training AIC, AICc, or BIC, use [select_cv()].
#' @examples
#' data("landscape_counts")
#' data("surv_pts")
#' # Three spatial folds. The simulated coordinates are in meters.
#' folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
#'                           crs = 32617)
#' candidates <- list(
#'   null = cross_validate_glm(counts ~ 1, landscape_counts, folds,
#'                             family = "poisson"),
#'   site = cross_validate_glm(counts ~ site, landscape_counts, folds,
#'                             family = "poisson"))
#'
#' # Pooled ranking. Delta is the gap from the best model in animals per survey.
#' cv_tab(candidates, metric = "rmse")
#'
#' # exp(Delta) is how many times more probability the best model gave the
#' # observed counts, as a geometric mean per survey
#' cv_tab(candidates, metric = "neg_log_score")
#'
#' # Is the ranking the same in every fold?
#' cv_tab(candidates, metric = "rmse", by_fold = TRUE)
#' @seealso [cross_validate()], [cross_validate_glm()], [select_cv()],
#'   [multiScaleR_cv-class], [aic_tab()], and
#'   `vignette("v08-cross-validation", package = "multiScaleR")`.
#' @export
cv_tab <- function(mod_list, metric, by_fold = FALSE) {
  validate_scalar_logical(by_fold, "by_fold")
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      !metric %in% c("rmse", "mae", "neg_log_score", "deviance", "brier"))
    stop("Select one supported metric explicitly.", call. = FALSE)
  if (is.list(mod_list) && any(vapply(mod_list, inherits, logical(1), "multiScaleR_bayes_cv")))
    return(.msr_bayes_cv_tab(mod_list, metric, by_fold))
  if (is.list(mod_list) && any(vapply(mod_list, inherits, logical(1), "multiScaleR_field_cv")))
    return(.msr_field_cv_tab(mod_list, metric, by_fold))
  if (is.list(mod_list) && any(vapply(mod_list, inherits, logical(1), "multiScaleR_mixed_cv")))
    return(.msr_mixed_cv_tab(mod_list, metric, by_fold))
  .msr_cv_comparison(mod_list, metric)
  groups <- if (by_fold) names(mod_list[[1L]]$folds) else "pooled"
  tables <- lapply(groups, function(fold) {
    rows <- lapply(names(mod_list), function(name) {
      x <- mod_list[[name]]; p <- x$predictions
      if (by_fold) p <- p[p$fold == fold, , drop = FALSE]
      p <- .msr_cv_score(p, x$settings$family, metric)
      score <- .msr_cv_metric(p, metric, x$settings$family, TRUE, nrow(p))
      if (is.na(score$estimate)) stop("Every candidate must have a defined score with positive scoring weight.", call. = FALSE)
      messages <- if (by_fold) x$warnings[fold] else x$warnings
      data.frame(Modnames = name, Score = score$estimate, Metric = metric, Unit = score$unit,
        N_assessment = score$n_assessment, N_scored = score$n_scored,
        Weight_sum = score$weight_sum, Warnings = length(unlist(messages, use.names = FALSE)))
    })
    tab <- do.call(rbind, rows)
    tab$Rank <- rank(tab$Score, ties.method = "min")
    best <- min(tab$Score); tab$Delta <- if (is.finite(best)) tab$Score - best else NA_real_
    if (by_fold) tab$Fold <- fold
    tab <- tab[order(tab$Rank), c("Modnames", "Rank", "Score", "Delta", "Metric", "Unit",
      "N_assessment", "N_scored", "Weight_sum", "Warnings", if (by_fold) "Fold"), drop = FALSE]
    tab
  })
  out <- do.call(rbind, tables); rownames(out) <- NULL
  out
}

.msr_cv_training_identity <- function(fit, split, targets, family, method = NULL) {
  if (identical(family, "conditional_choice")) return(.msr_choice_identity(fit, split, targets, method))
  if (.msr_cv_distribution_family(family)) return(.msr_cv_distribution_training_identity(fit, split, targets, family))
  model <- if (inherits(fit, "multiScaleR")) fit$opt_mod else fit
  if (!inherits(model, "glm") || !is.data.frame(model$model) ||
      !identical(rownames(model$model), split$analysis))
    stop("Retained GLMs must match their analysis row IDs exactly; full-data fits cannot replace training fits.", call. = FALSE)
  link <- switch(family, gaussian = "identity", poisson = "log", binomial = "logit")
  if (!identical(model$family$family, family) || !identical(model$family$link, link))
    stop("Retained training family and link must match the CV scoring target.", call. = FALSE)
  data <- data.frame(row.names = split$analysis)
  data$outcome <- stats::model.response(model$model)
  y <- .msr_cv_outcome(outcome ~ 1, data, family)
  w <- stats::model.weights(model$model)
  if (is.null(w)) w <- rep(1, nrow(data))
  expected <- targets[match(split$analysis, targets$row_id), , drop = FALSE]
  if (!identical(y$observed, expected$observed) || !identical(y$trials, expected$trials) ||
      !identical(as.numeric(w), expected$weight))
    stop("Retained training responses, trials or weights disagree with recorded targets.", call. = FALSE)
  invisible(TRUE)
}

#' Validate Training-Only Model Selection Within CV Folds
#'
#' Score a model-selection rule rather than a single model. In each fold the
#' rule picks the candidate with the lowest training AIC, AICc, or BIC, and
#' the result keeps the held-out predictions that candidate already made. No
#' models are refitted.
#'
#' @param mod_list Named list of complete cross-validation results, all scored
#'   on the same folds and observations. Each must keep its training records
#'   (`keep_selection = TRUE`) or full fits (`keep_fits = TRUE`). Results from
#'   [cross_validate_field()] and [cross_validate_mixed()] always keep their
#'   records.
#' @param criterion Training criterion, declared in advance: `"AICc"`
#'   (default), `"AIC"`, or `"BIC"`. It is computed from each candidate's
#'   training log-likelihood and parameter count K with the same rules as
#'   [aic_tab()] and [bic_tab()]. AICc needs more training observations than
#'   K + 1.
#' @param ties What to do when two candidates have exactly the same training
#'   criterion. The default `"error"` marks the fold as failed. `"first"`
#'   picks the tied candidate listed first in `mod_list`. No tolerance is
#'   applied, so only exact ties count.
#' @param keep_fits Logical. The default `FALSE` returns compact output.
#'   `TRUE` also returns the selected training fit of each fold in `fits`.
#'   This requires every candidate to have been run with `keep_fits = TRUE`,
#'   because records cannot rebuild a fit.
#' @return A cross-validation result of the same class as the candidates
#'   that scores the selection procedure. Each fold's predictions, training
#'   values, and diagnostics come from the candidate chosen in that fold, and
#'   `settings$backend` is `"training_selection"`. The result adds
#'   `selection`, which holds `chosen` (the candidate picked in each fold,
#'   `NA` where selection failed), `tables` (each fold's training table, with
#'   every candidate's K, log-likelihood, and criterion), and the `criterion`,
#'   `ties`, and `candidates`. Results from [cross_validate()],
#'   [cross_validate_glm()], and [cross_validate_clogit()] also record the
#'   `metrics` scored and `metrics_dropped`, and each candidate's
#'   `candidate_settings`. See [multiScaleR_cv-class] for the other elements.
#' @details
#' # How the procedure is scored
#'
#' Suppose your plan is to fit several candidate models and keep the one with
#' the lowest AICc. Scoring only the winner on held-out data would flatter it,
#' because it was chosen after seeing all of the data. `select_cv()` scores
#' the whole plan instead. In each fold it computes the training criterion of
#' every candidate from its fit to that fold's analysis rows, picks the
#' lowest, and keeps the held-out predictions that candidate made. Assessment
#' outcomes and CV scores never choose the winner. The pooled score therefore
#' evaluates the declared candidate set and rule. If you then tune the
#' candidate set or the rule using these results, evaluate the final choice
#' with an independent test set or an outer validation layer.
#'
#' # What the candidates must keep
#'
#' Run every candidate with `keep_selection = TRUE` (compact records) or
#' `keep_fits = TRUE` (full fits). A record holds the native training
#' log-likelihood, the native and ecological parameter counts, the rows that
#' entered the likelihood, their responses, trials, and frequency weights,
#' the formula as text, the controls, and the original diagnostics. NB2 and
#' Gamma records also hold the estimated or fixed nuisance parameter. A
#' record holds no fitted model, design matrix, raster, or formula
#' environment, so it is much smaller than a fit, although its row metadata
#' still grows with the number of observations and folds. Keeping records
#' adds no fitting and does not change which folds pass their diagnostics.
#'
#' Before selecting, the function checks each record or fit against the
#' fold's analysis rows and the recorded responses and weights. If a
#' candidate keeps both records and fits, they must agree exactly. A record
#' that could not be built makes selection fail in that fold, without
#' changing the candidate's own result. Results that kept fits work without
#' records. A result that kept neither must be run again with
#' `keep_selection = TRUE`.
#'
#' # Which candidates can be combined
#'
#' The candidates must pass the same checks as in [cv_tab()], with identical
#' folds and observations and comparable likelihoods. NB2 and Poisson
#' candidates can share a count target. Gamma and Gaussian candidates can
#' share a continuous target, scored with RMSE or MAE. An information
#' criterion cannot compare a count likelihood with a continuous density.
#' Conditional-choice candidates must also share their partial-likelihood
#' `method`.
#'
#' For results from [cross_validate()], [cross_validate_glm()], and
#' [cross_validate_clogit()], the procedure is scored on the metrics that
#' every candidate recorded. A default Poisson candidate records `"deviance"`
#' and a default NB2 candidate does not, so their selection keeps `"rmse"`,
#' `"mae"`, and `"neg_log_score"`. `selection$metrics` lists the metrics kept
#' and `selection$metrics_dropped` those left out. Field and mixed-model
#' candidates must record the same metrics. Bayesian results have no
#' training likelihood and are rejected.
#'
#' Penalties follow the native parameter counts, and scales and shapes
#' estimated by [multiScale_optim()] count toward K.
#'
#' * Conditional choice: the sample size in AICc and BIC is the number of
#'   choice sets (one event each), not the number of alternatives.
#' * Spatial fields: the coefficients, field SD, estimated range, and any
#'   observation dispersion count toward K. The integrated field values add
#'   nothing. The held-out joint block score is that of the model chosen in
#'   the fold.
#' * Mixed models: the sample size is the number of analysis observations.
#'   The covariance parameters count toward K, and the integrated group
#'   effects add nothing.
#'
#' # Read the result
#'
#' Read `complete` and `selection$chosen` first. A fold fails, and the pooled
#' scores become `NA`, when a criterion cannot be computed, is not finite, or
#' ties under `ties = "error"`. `selection$tables` shows the training table
#' behind each choice. For GLM and conditional-choice results, the shared
#' `settings` are copied from the first candidate, named in
#' `settings$controls_from`, and `selection$candidate_settings` holds every
#' candidate's controls. NB2 and Gamma selections also keep the settings
#' (`settings$fold_settings`) and `distribution` of the model chosen in each
#' fold, including when a Poisson or Gaussian reference wins. GLM,
#' conditional-choice, and mixed-model selections also keep every
#' candidate's warnings and diagnostics in `selection`, and grouped-map
#' labels (`raster_groups`) are kept for the chosen model. The score does not
#' integrate parameter or selection uncertainty, add intervals, or establish
#' calibration.
#' @examples
#' data("landscape_counts")
#' data("surv_pts")
#' # Three spatial folds. The simulated coordinates are in meters.
#' folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
#'                           crs = 32617)
#'
#' # keep_selection = TRUE stores the training likelihoods that select_cv() reads
#' candidates <- list(
#'   null = cross_validate_glm(counts ~ 1, landscape_counts, folds,
#'                             family = "poisson", keep_selection = TRUE),
#'   site = cross_validate_glm(counts ~ site, landscape_counts, folds,
#'                             family = "poisson", keep_selection = TRUE))
#' selected <- select_cv(candidates, criterion = "AICc")
#' selected$selection$chosen         # model picked by training AICc in each fold
#' selected$selection$tables$fold1   # the training table behind the first pick
#' selected$summary                  # held-out scores of the whole procedure
#' @seealso [cv_tab()], [cross_validate()], [cross_validate_glm()],
#'   [multiScaleR_cv-class], [aic_tab()], [bic_tab()], and
#'   `vignette("v08-cross-validation", package = "multiScaleR")`.
#' @export
select_cv <- function(mod_list, criterion = c("AICc", "AIC", "BIC"),
    ties = c("error", "first"), keep_fits = FALSE) {
  criterion <- match.arg(criterion); ties <- match.arg(ties)
  validate_scalar_logical(keep_fits, "keep_fits")
  if (is.list(mod_list) && any(vapply(mod_list, inherits, logical(1), "multiScaleR_bayes_cv")))
    stop("Bayesian CV has no training ML selection criterion. Compare predictive scores with cv_tab().", call. = FALSE)
  if (is.list(mod_list) && any(vapply(mod_list, inherits, logical(1), "multiScaleR_field_cv")))
    return(.msr_field_select_cv(mod_list, criterion, ties, keep_fits))
  if (is.list(mod_list) && any(vapply(mod_list, inherits, logical(1), "multiScaleR_mixed_cv")))
    return(.msr_mixed_select_cv(mod_list, criterion, ties, keep_fits))
  signature <- .msr_cv_comparison(mod_list)
  first <- mod_list[[1L]]; folds <- first$folds
  # Score the selection procedure on the metrics every candidate recorded.
  # Default metric sets differ by family (Poisson adds `deviance`; NB2 and
  # Gamma do not), so requiring identical sets rejected the default
  # Poisson-versus-NB2 comparison. Candidates with identical sets are unaffected.
  candidate_metrics <- lapply(mod_list, function(x) x$settings$metrics)
  metrics <- Reduce(intersect, candidate_metrics)
  if (!length(metrics))
    stop("Candidate results share no scoring metric; rerun cross-validation with a common `metrics` set.", call. = FALSE)
  dropped <- setdiff(unique(unlist(candidate_metrics, use.names = FALSE)), metrics)
  choice <- inherits(first, "multiScaleR_choice_cv")
  distribution <- any(vapply(mod_list, inherits, logical(1), "multiScaleR_distribution_cv"))
  all_fits <- all(vapply(mod_list, function(x) !is.null(x$fits), logical(1)))
  records <- lapply(mod_list, function(x) {
    has_fits <- !is.null(x$fits) && setequal(names(x$fits), names(folds))
    has_records <- !is.null(x$selection_records) && setequal(names(x$selection_records), names(folds))
    if (!is.null(x$fits) && !has_fits)
      stop("Retain every training fit when supplying fits.", call. = FALSE)
    if (!is.null(x$selection_records) && !has_records)
      stop("Retain every training selection record when supplying records.", call. = FALSE)
    if (!has_fits && !has_records)
      stop("Candidate results must retain every training fit or selection record: use keep_fits = TRUE or keep_selection = TRUE.", call. = FALSE)
    if (keep_fits && !has_fits)
      stop("Selected full fits require keep_fits = TRUE in every candidate; selection records cannot reconstruct fits.", call. = FALSE)
    out <- stats::setNames(vector("list", length(folds)), names(folds))
    for (name in names(folds)) {
      fold_settings <- .msr_cv_fold_settings(x, name)
      if (has_fits) .msr_cv_training_identity(x$fits[[name]], x$folds[[name]], signature$targets, fold_settings$family, fold_settings$method)
      native <- if (has_fits && (has_records || !all_fits)) {
        if (choice) .msr_choice_record(x$fits[[name]], x$diagnostics[[name]], .msr_cv_selection_controls(fold_settings),
          x$folds[[name]], signature$targets) else
          .msr_cv_selection_record(x$fits[[name]], x$diagnostics[[name]], .msr_cv_selection_controls(fold_settings))
      } else NULL
      if (has_records) {
        .msr_cv_record_identity(x$selection_records[[name]], x, name, signature$targets)
        if (has_fits && !identical(x$selection_records[[name]], native))
          stop("Selection records disagree with their retained native training fits; rerun CV.", call. = FALSE)
      }
      out[name] <- list(if (has_records) x$selection_records[[name]] else native)
    }
    out
  })
  output <- list(predictions = list(), metrics = list(), fold_results = list(), training = list(),
    diagnostics = list(), warnings = list(), errors = list(), prediction_contracts = list(),
    fits = if (keep_fits) list() else NULL)
  if (distribution) output$distribution <- list()
  selected_settings <- list()
  selection <- list(criterion = criterion, ties = ties, candidates = names(mod_list),
    metrics = metrics, metrics_dropped = dropped,
    candidate_settings = lapply(mod_list, `[[`, "settings"),
    chosen = stats::setNames(rep(NA_character_, length(folds)), names(folds)), tables = list(),
    candidate_warnings = lapply(mod_list, `[[`, "warnings"), candidate_diagnostics = lapply(mod_list, `[[`, "diagnostics"))
  for (name in names(folds)) {
    messages <- character(); winner <- NA_character_
    error <- tryCatch(withCallingHandlers({
      if (all_fits) {
        models <- lapply(mod_list, function(x) x$fits[[name]])
        tab <- if (criterion == "BIC") bic_tab(models, mod_names = names(mod_list), verbose = FALSE) else
          aic_tab(models, AICc = criterion == "AICc", mod_names = names(mod_list), verbose = FALSE)
      } else {
        tab <- .msr_cv_record_table(lapply(records, `[[`, name), criterion, names(mod_list))
      }
      selection$tables[[name]] <- tab
      values <- stats::setNames(tab[[criterion]], tab$Modnames)[names(mod_list)]
      if (length(values) != length(mod_list) || any(!is.finite(values)))
        stop("Training information criteria must all be finite.", call. = FALSE)
      best <- names(values)[values == min(values)]
      if (length(best) > 1L && ties == "error") stop("Exact training-criterion tie; predeclare ties = 'first' to resolve it.", call. = FALSE)
      winner <- best[1L]
      NULL
    }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(stage = "training_selection", message = conditionMessage(e)))
    selection$chosen[name] <- winner
    x <- if (is.null(error)) mod_list[[winner]] else first
    pred <- if (distribution) .msr_cv_distribution_rows(x, name) else x$predictions[x$predictions$fold == name, , drop = FALSE]
    # Remove recorded losses outside the shared set so folds won by different
    # candidates keep identical prediction columns.
    pred[intersect(dropped, names(pred))] <- NULL
    if (is.null(error)) pred <- .msr_cv_score(pred, .msr_cv_fold_settings(x, name)$family, metrics)
    if (!is.null(error)) {
      pred$prediction <- NA_real_; pred$eligible <- FALSE
      if (choice) pred$linear_predictor <- pred$log_probability <- NA_real_
      for (metric in metrics) pred[[metric]] <- NA_real_
    }
    output$predictions[[name]] <- pred
    output$metrics[[name]] <- do.call(rbind, lapply(metrics, function(metric)
      cbind(fold = name, .msr_cv_metric(pred, metric, x$settings$family, is.null(error), nrow(pred)))))
    output$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(folds[[name]]$analysis),
      n_assessment = length(folds[[name]]$assessment), status = if (is.null(error)) "scored" else "failed")
    output$errors[name] <- list(error)
    output$warnings[name] <- list(c(messages, if (is.null(error)) x$warnings[[name]]))
    if (is.null(error)) {
      for (item in c("training", "diagnostics", "prediction_contracts")) output[[item]][name] <- x[[item]][name]
      if (!is.null(x$raster_groups)) output$raster_groups[name] <- x$raster_groups[name]
      if (keep_fits) output$fits[name] <- x$fits[name]
      if (distribution) {
        output$distribution[name] <- list(.msr_cv_distribution_metadata(x, name))
        selected_settings[name] <- list(.msr_cv_fold_settings(x, name))
      }
    }
  }
  output <- .msr_cv_finish(output, folds, names(first$assessment_counts), first$settings$family, metrics)
  output$targets <- signature$targets; output$selection <- selection
  output$settings <- first$settings
  output$settings$metrics <- metrics
  output$settings$backend <- "training_selection"
  output$settings$controls_from <- names(mod_list)[1L]
  output$settings$selection_uncertainty <- "not_integrated"
  if (distribution) {
    output$settings$fold_settings <- selected_settings
    class(output) <- c("multiScaleR_distribution_cv", "multiScaleR_cv")
  }
  output$call <- match.call()
  output
}

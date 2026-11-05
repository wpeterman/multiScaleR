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
#' Fit a GLM using only each fold's analysis rows, then score its held-out
#' response predictions. This provides null and nonspatial reference models for
#' the same folds used to validate estimated ecological scales.
#'
#' @inheritParams cross_validate
#' @param data Data frame of raw responses and predictor variables, with unique
#'   row IDs matching the fold plan. Supply every formula variable explicitly.
#' @param glm_control Named native \code{stats::glm.control()} controls.
#' @return A \code{multiScaleR_cv} object as described in
#'   \code{\link{cross_validate}}, with backend \code{"glm"}, no estimated
#'   spatial scales, and ordinary GLM fits when \code{keep_fits = TRUE}.
#'   Compact \code{targets} metadata records used row IDs, observed outcomes,
#'   trials and frequency weights for guarded model comparison. Optional
#'   \code{selection_records} have the same contract as \code{cross_validate()}.
#' @details This adapter supports the same qualified GLM families, outcome and
#'   frequency-weight definitions as \code{cross_validate()}. Formulas, factors,
#'   interactions, offsets and learned \code{poly()} terms use native GLM
#'   machinery. All learned preprocessing must occur within the training fit;
#'   do not supply predictors standardized with assessment data.
#'
#' NB2 uses \code{MASS::glm.nb()} when theta is estimated and
#'   \code{glm(family = MASS::negative.binomial(theta))} when it is fixed.
#'   Gamma retains its native deviance-based AIC scale when dispersion is
#'   estimated. Known Gamma dispersion uses the normalized fixed-scale family
#'   constructed before fitting. Coefficients and native coefficient SEs are
#'   retained; density dispersion is not a replacement for the SE convention.
#'   Nuisance parameters and supported training transformations are frozen
#'   before assessment prediction. See \code{cross_validate()} for their counts,
#'   support, compact metadata and uncertainty distinctions.
#'
#' No ecological scale is estimated here. Use \code{cross_validate()} for
#'   training-only spatial scale refitting. Predictions use fitted parameters;
#'   parameter uncertainty is not integrated, observation intervals are not
#'   provided and calibration is not established. Fold failures and native
#'   warnings are retained. No missing rows are silently dropped.
#' @examples
#' dat <- data.frame(y = c(1, 2, 0, 3, 2, 1, 4, 2), x = 1:8)
#' folds <- grouped_cv_folds(rep(letters[1:4], each = 2), k = 2)
#' cv <- cross_validate_glm(y ~ x, dat, folds, family = "poisson")
#' cv
#' @seealso \code{\link{cross_validate}}, \code{\link{cv_tab}},
#'   \code{\link{select_cv}}
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
#' Rank candidate predictions using identical recorded targets and splits.
#' Lower scores are better. Scores are pooled from assessment occurrences,
#' preserving frequency weights and binomial trial definitions.
#'
#' @param mod_list Uniquely named list of complete \code{multiScaleR_cv} results.
#'   Use identical explicit folds and response/weight definitions for every model.
#'   Complete \code{\link{cross_validate_clogit}} results compare only with the
#'   same conditional-choice target, complete sets and native partial method.
#' @param metric One recorded score, selected explicitly: \code{"rmse"},
#'   \code{"mae"}, \code{"neg_log_score"}, \code{"deviance"} or \code{"brier"}.
#' @param by_fold Logical. Default \code{FALSE} compares pooled scores.
#'   \code{TRUE} gives descriptive ranks within each assessment fold.
#' @return A data frame with \code{Modnames}, \code{Rank}, \code{Score},
#'   \code{Delta}, \code{Metric}, \code{Unit}, \code{N_assessment},
#'   \code{N_scored}, \code{Weight_sum}, \code{Warnings} and, for per-fold
#'   comparisons, \code{Fold}. Ties share their minimum rank. If every score
#'   is infinite, deltas are \code{NA}. Original candidate fits and diagnostics
#'   are not modified.
#' @details Incomplete candidates are rejected as a group; no failed model or
#'   fold is silently removed. Comparisons match analysis/assessment IDs and all
#'   used responses, trials and frequency weights. RMSE/MAE can compare qualified
#'   conditional-mean predictions across qualified Gaussian, Poisson, NB2 and
#'   Gamma targets on the same observations. NB2 and Poisson log scores can
#'   compare their count distributions. Continuous/discrete density scores
#'   cannot be mixed, and Gaussian density scores remain unqualified. Gamma
#'   density comparisons preserve positive-response and training-dispersion
#'   guards. Other canonical scores retain their family/link requirements. Conditional-choice comparisons use
#'   \code{neg_log_score} only, per held-out set, preserving every alternative
#'   and its membership. They cannot be mixed with Bernoulli/GLM targets.
#'   Binomial negative log scores refer to
#'   count probabilities; Brier scores average over individual trials.
#'
#' There is no added complexity penalty, Akaike weight or model probability.
#'   Fold variation is not an independence-based standard error or calibration
#'   certificate. Per-fold assessment ranks must not select that fold's model.
#'   Use \code{\link{select_cv}} for training-only information-criterion selection.
#'   Choosing a model or tuning rule from pooled CV results requires independent
#'   assessment or an outer validation layer to evaluate that choice.
#' @examples
#' dat <- data.frame(y = c(1, 2, 0, 3, 2, 1, 4, 2), x = 1:8)
#' folds <- grouped_cv_folds(rep(letters[1:4], each = 2), k = 2)
#' candidates <- list(null = cross_validate_glm(y ~ 1, dat, folds, "poisson"),
#'                    trend = cross_validate_glm(y ~ x, dat, folds, "poisson"))
#' cv_tab(candidates, metric = "neg_log_score")
#' @seealso \code{\link{cross_validate}}, \code{\link{cross_validate_glm}},
#'   \code{\link{select_cv}}, \code{\link{aic_tab}}
#' @export
cv_tab <- function(mod_list, metric, by_fold = FALSE) {
  validate_scalar_logical(by_fold, "by_fold")
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      !metric %in% c("rmse", "mae", "neg_log_score", "deviance", "brier"))
    stop("Select one supported metric explicitly.", call. = FALSE)
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
#' Select a candidate using each fold's training likelihood, then evaluate its
#' already recorded held-out predictions. No additional model fits are run.
#'
#' @inheritParams cv_tab
#' @param criterion Predeclared training criterion: \code{"AICc"} (default),
#'   \code{"AIC"} or \code{"BIC"}, calculated with existing
#'   \code{\link{aic_tab}} or \code{\link{bic_tab}} parameter counts and guards.
#' @param ties How to handle an exact training-criterion tie. Default
#'   \code{"error"} retains the fold as failed. \code{"first"} selects the
#'   first tied candidate in the declared list. No tolerance is added.
#' @param keep_fits Logical. Default \code{FALSE} keeps compact selection
#'   records. \code{TRUE} retains the selected training fit per fold in memory.
#' @return A \code{multiScaleR_cv} result for the selection procedure, with
#'   \code{selection} metadata containing the criterion, tie policy, candidate
#'   names, chosen model and native training table per fold, including original
#'   parameter penalties, and \code{candidate_settings} with each candidate's
#'   controls. Shared output settings inherit the first candidate, identified
#'   by \code{settings$controls_from}; use \code{candidate_settings} for the
#'   complete procedure. NB2/Gamma procedures also retain the actual selected
#'   \code{settings$fold_settings} and \code{distribution} metadata per fold,
#'   including when a Poisson or Gaussian reference wins. Selected
#'   \code{raster_groups} are retained when present.
#'   Candidate warnings/diagnostics and selection failures
#'   are retained. Failed selection makes pooled scores \code{NA}.
#' @details Inputs must retain native training-selection records via
#'   \code{keep_selection = TRUE} or full fold fits via \code{keep_fits = TRUE},
#'   share compatible likelihood targets and splits, and record the same scoring
#'   metrics. NB2/Poisson candidates may share a count target. Gamma/Gaussian
#'   candidates may share a continuous target using RMSE/MAE. Density scores
#'   cannot compare different observation measures. Other canonical family/link
#'   and conditional-choice guards are unchanged.
#'   Analysis-row identities and response/weight definitions are checked against
#'   each retained GLM or record. Compact records store native log likelihood,
#'   native and ecological parameter counts, effective row IDs, outcomes, trials,
#'   frequencies, formula text, controls and diagnostics. NB2/Gamma records also
#'   preserve estimated/fixed nuisance metadata and native nuisance counts.
#'   They retain no model
#'   environment, design matrix or raw raster. If both representations are
#'   retained, they must agree exactly. Record-only selection uses the same
#'   likelihood-signature guards and information-criterion table constructors.
#'   Selection records cannot reconstruct fits requested by \code{keep_fits}.
#'   Metadata failures remain selection failures without discarding scored folds.
#'   Conditional-choice records retain event-based nobs and native
#'   partial-method/stratum signatures. Alternative counts do not substitute
#'   for native events in AICc/BIC. Native likelihood comparability, parameter
#'   counts and numerical guards remain unchanged. Invalid inputs are rejected; criterion
#'   failures, nonfinite criteria and unresolved ties are retained per fold.
#'
#' Assessment outcomes and CV scores never choose the winner inside a fold.
#'   The resulting held-out score evaluates the predeclared candidate set and
#'   selection rule. If those are tuned using these results, use an independent
#'   test set or outer validation layer. Selection does not integrate parameter
#'   or selection uncertainty, add intervals or establish calibration.
#' @examples
#' dat <- data.frame(y = c(1, 2, 0, 3, 2, 1, 4, 2), x = 1:8)
#' folds <- grouped_cv_folds(rep(letters[1:4], each = 2), k = 2)
#' candidates <- list(null = cross_validate_glm(y ~ 1, dat, folds, "poisson", keep_selection = TRUE),
#'                    trend = cross_validate_glm(y ~ x, dat, folds, "poisson", keep_selection = TRUE))
#' selected <- select_cv(candidates, criterion = "AIC", ties = "first")
#' selected$selection$chosen
#' @seealso \code{\link{cv_tab}}, \code{\link{cross_validate}},
#'   \code{\link{cross_validate_glm}}, \code{\link{aic_tab}}, \code{\link{bic_tab}}
#' @export
select_cv <- function(mod_list, criterion = c("AICc", "AIC", "BIC"),
    ties = c("error", "first"), keep_fits = FALSE) {
  criterion <- match.arg(criterion); ties <- match.arg(ties)
  validate_scalar_logical(keep_fits, "keep_fits")
  signature <- .msr_cv_comparison(mod_list)
  first <- mod_list[[1L]]; metrics <- first$settings$metrics; folds <- first$folds
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
    if (!setequal(metrics, x$settings$metrics) || (!has_fits && !has_records))
      stop("Candidate results must share metrics and retain every training fit or selection record: use keep_fits = TRUE or keep_selection = TRUE.", call. = FALSE)
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

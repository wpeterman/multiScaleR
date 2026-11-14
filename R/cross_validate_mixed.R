.msr_mixed_cv_formula <- function(formula, data, family) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("Supply a two-sided mixed-model formula.", call. = FALSE)
  if ("." %in% all.vars(formula))
    stop("Use explicit predictors in mixed CV formulas.", call. = FALSE)
  terms <- stats::terms(formula, data = data)
  labels <- attr(terms, "term.labels")
  expressions <- lapply(labels, str2lang)
  bars <- which(vapply(expressions, .msr_marginal_random_call, logical(1)))
  if (length(bars) != 1L)
    stop("Use one random-effect block with a bare grouping column, such as (1 | site).", call. = FALSE)
  bar <- expressions[[bars]]
  while (identical(bar[[1L]], as.name("("))) bar <- bar[[2L]]
  if (!identical(bar[[1L]], as.name("|")) || !is.symbol(bar[[3L]]))
    stop("Use a single standard random-effect block with a bare grouping column.", call. = FALSE)
  group <- as.character(bar[[3L]])
  if (!group %in% names(data) || !is.atomic(data[[group]]) || !is.null(dim(data[[group]])) ||
      anyNA(data[[group]]) || any(!nzchar(as.character(data[[group]]))))
    stop("The grouping column must contain a nonmissing group ID for every row.", call. = FALSE)
  random <- stats::as.formula(call("~", bar[[2L]]), env = environment(formula))
  random_vars <- all.vars(random)
  if (family != "gaussian" && (!is.numeric(bar[[2L]]) || !identical(as.numeric(bar[[2L]]), 1)))
    stop("Non-Gaussian mixed CV uses one random intercept. Gaussian CV also supports numeric random slopes.", call. = FALSE)
  variables <- as.list(attr(terms, "variables"))[-1L]
  offsets <- vapply(variables[attr(terms, "offset")], deparse1, character(1), backtick = TRUE)
  fixed <- stats::reformulate(c(labels[-bars], offsets), response = formula[[2L]],
    intercept = attr(terms, "intercept") == 1L, env = environment(formula))
  if (length(setdiff(c(all.vars(fixed), random_vars), names(data))) ||
      any(all.vars(fixed[[2L]]) %in% c(group, all.vars(fixed[[3L]]), random_vars)) ||
      group %in% random_vars)
    stop("Supply explicit raw formula variables, with a separate group and no response in the predictors.", call. = FALSE)
  if (any(!vapply(data[random_vars], function(x) is.numeric(x) && is.null(dim(x)), logical(1))))
    stop("Gaussian random slopes must use numeric predictor columns.", call. = FALSE)
  list(fixed = fixed, random = random, group = group)
}

.msr_mixed_cv_basis <- function(formula, training, contrasts = NULL) {
  tt <- stats::delete.response(stats::terms(formula, data = training))
  tt <- .msr_cv_prediction_terms(tt, training, qualify = TRUE)
  frame <- stats::model.frame(tt, training, na.action = stats::na.fail, drop.unused.levels = TRUE)
  X <- stats::model.matrix(tt, frame, contrasts.arg = contrasts)
  list(terms = tt, xlevels = lapply(frame[vapply(frame, is.factor, logical(1))], levels),
    contrasts = attr(X, "contrasts"), columns = colnames(X), X = X)
}

.msr_mixed_cv_design <- function(basis, data) {
  result <- .msr_unmarked_fixed_design(basis, data)
  if (nrow(result$X) != nrow(data) || any(!is.finite(result$X)) ||
      any(!is.finite(result$offset)))
    stop("Mixed-model predictions require a complete, finite design and offset.", call. = FALSE)
  if (is.null(result$offset)) result$offset <- numeric(nrow(data))
  result
}

.msr_mixed_cv_native <- function(model, spec, training, assessment, family) {
  frame <- stats::model.frame(model)
  if (!identical(rownames(frame), rownames(training)) ||
      !isTRUE(model$fit$convergence == 0L) || !isTRUE(model$sdr$pdHess) ||
      !is.finite(.msr_model_loglik(model)))
    stop("The native mixed fit must retain all analysis rows, converge and have a positive-definite Hessian.", call. = FALSE)
  expected <- eval(spec$fixed[[2L]], training, environment(spec$fixed))
  actual <- stats::model.response(frame)
  if (!identical(as.vector(expected), as.vector(actual)))
    stop("The native mixed fit changed the recorded analysis responses.", call. = FALSE)
  reference <- as.matrix(glmmTMB::getME(model, "X"))
  fixed <- .msr_mixed_cv_basis(spec$fixed, training, model$modelInfo$contrasts)
  random <- .msr_mixed_cv_basis(spec$random, training)
  tr <- .msr_mixed_cv_design(fixed, training)
  te <- .msr_mixed_cv_design(fixed, assessment)
  Ztrain <- .msr_mixed_cv_design(random, training)$X
  Ztest <- .msr_mixed_cv_design(random, assessment)$X
  beta <- glmmTMB::fixef(model)$cond
  covariance <- glmmTMB::VarCorr(model)$cond[[spec$group]]
  effects <- glmmTMB::ranef(model)$cond[[spec$group]]
  groups <- match(as.character(training[[spec$group]]), rownames(effects))
  close <- function(x, y) isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-9, check.attributes = FALSE))
  if (is.null(covariance) || any(!is.finite(covariance)) || anyNA(groups) ||
      !identical(colnames(Ztrain), rownames(covariance)) ||
      !identical(colnames(Ztrain), names(effects)) ||
      !identical(colnames(tr$X), names(beta)) || !identical(colnames(tr$X), colnames(reference)) ||
      !close(tr$X, reference) || any(!is.finite(beta)))
    stop("The training-only recipes did not reproduce the native mixed-model design.", call. = FALSE)
  eta_train <- drop(tr$X %*% beta) + tr$offset
  eta_test <- drop(te$X %*% beta) + te$offset
  native <- stats::predict(model, type = "link")
  fitted_random <- rowSums(Ztrain * as.matrix(effects)[groups, , drop = FALSE])
  if (!close(eta_train + fitted_random, native) || any(!is.finite(c(eta_train, eta_test))))
    stop("Training fixed effects, random effects and offsets do not reproduce the native link predictions.", call. = FALSE)
  dispersion <- if (family %in% c("gaussian", "nbinom2", "Gamma")) stats::sigma(model) else 1
  if (!is.finite(dispersion) || dispersion <= 0)
    stop("The native dispersion must be positive and finite.", call. = FALSE)
  recipe <- function(x) list(predvars = deparse1(attr(x$terms, "predvars")),
    xlevels = x$xlevels, contrasts = x$contrasts, columns = x$columns)
  list(eta_train = eta_train, eta_test = eta_test, Ztrain = Ztrain, Ztest = Ztest,
    covariance = covariance, dispersion = dispersion, beta = beta,
    recipe_sha256 = .msr_gdr_hash(list(recipe(fixed), recipe(random))))
}

.msr_mixed_cv_prepare <- function(inputs, ids, covariates, transform) {
  train <- .msr_cv_inputs(inputs, ids, covariates)
  raw <- .msr_cv_initial_values(train)
  normalization <- .msr_cv_raw_identity(colMeans(raw), apply(raw, 2L, stats::sd), transform,
    train$scale_vars)
  if (any(!is.finite(normalization$sd) | normalization$sd <= 0))
    stop("Training spatial predictors need positive finite spread.", call. = FALSE)
  train$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, normalization$mean), 2L,
    normalization$sd, "/"))
  train$scl_params <- normalization
  train
}

#' Cross-Validate Mixed-Effects Models
#'
#' Refit a mixed model in every analysis split and score the observations
#' held out from it. Choose whether prediction is for groups represented in
#' training or for wholly new groups. Integrate over the shared group effect
#' when calculating response means and joint group log scores.
#'
#' This function is experimental. Its arguments and results may change in a
#' later release without a deprecation period.
#'
#' @param formula Two-sided native mixed-model formula with one grouping
#'   column, for example \code{count ~ habitat + (1 | site)}. Gaussian models
#'   also support one block of numeric random slopes such as
#'   \code{(1 + habitat | site)}. Use explicit predictors and formula offsets.
#'   The response is a data column, or \code{cbind(successes, failures)} for
#'   grouped binomial data with positive trial counts.
#' @param data Raw data frame with unique row names. The names identify rows
#'   in the split plan and, if supplied, in the prepared raw-cell buffers.
#' @param folds Named list of explicit analysis and assessment row-ID splits.
#'   Use \code{\link{grouped_cv_folds}} or \code{\link{spatial_cv_folds}} to hold
#'   out entire groups. Use explicit splits within groups for existing-group
#'   prediction.
#' @param family Character response family. Choose \code{"gaussian"} (identity,
#'   default), \code{"poisson"} (log), \code{"binomial"} (logit),
#'   \code{"nbinom2"} (log), or \code{"Gamma"} (log).
#' @param prediction Character prediction target. \code{"new_group"} (default)
#'   requires each assessment group to be absent from analysis in that split.
#'   \code{"existing_group"} requires every assessment group to have analysis
#'   observations in that split.
#' @param kernel_inputs Optional raw-cell \code{multiScaleR_data} object for
#'   training-only ecological scale fitting with native refits. The default
#'   \code{NULL} fits the supplied formula without scale optimization.
#' @param optim_args Named list of extra serial \code{\link{multiScale_optim}}
#'   controls. Use only with \code{kernel_inputs}.
#' @param metrics Character vector of \code{"neg_log_score"}, \code{"rmse"},
#'   \code{"mae"}, or \code{"brier"} (binomial only). By default, report the
#'   joint group log score and RMSE, or the joint group log score and Brier
#'   score for binomial responses. The log score is summed across assessment
#'   groups and divided by the number of assessment rows. Binomial RMSE and
#'   MAE compare the expected success probability with the observed fraction
#'   of successes. Brier score averages squared probability error over trials.
#' @param control A \code{glmmTMB::glmmTMBControl()} list for the native
#'   analysis fits. The default uses native controls.
#' @param integration_tol Relative tolerance for the one-dimensional adaptive
#'   integration in non-Gaussian groups. The default is \code{1e-8}.
#' @param keep_fits Logical. Retain the fitted model for each successfully
#'   scored fold. Defaults to \code{FALSE}.
#'
#' @return An object of class \code{multiScaleR_mixed_cv} and
#'   \code{multiScaleR_cv}. \code{predictions} contains the fold and row ID,
#'   group ID, observed response, trials, integrated response mean, and an
#'   eligibility flag. \code{joint_scores} has one table per scored fold with
#'   the number of held-out rows and the negative log predictive density
#'   for each group. \code{summary} and \code{metrics} give pooled and
#'   fold-specific scores. \code{training} retains analysis IDs, native
#'   coefficients, covariance, dispersion, and any fitted ecological scales.
#'   \code{selection_records} stores the native analysis log likelihood,
#'   its row count and parameter counts. \code{settings}, \code{folds},
#'   \code{targets}, and \code{prediction_contracts} identify the comparison.
#'   \code{fold_results}, \code{errors}, \code{diagnostics}, and \code{warnings}
#'   retain failed folds and native convergence information.
#'   \code{fits} is optional. A failed fold makes pooled scores unavailable.
#'
#' @details Mixed CV uses maximum-likelihood fits with constant dispersion
#'   and a Gaussian distribution for the group effects. Given the fitted
#'   parameters, new groups use the estimated population distribution of
#'   effects. Existing groups condition that distribution on their analysis
#'   responses. Gaussian response models use exact multivariate Normal
#'   conditioning. For other families the function integrates over one group
#'   intercept with adaptive quadrature. Each assessment group contributes one
#'   joint predictive density. Parameter, ecological-scale and model-selection
#'   uncertainty are outside these plug-in predictions.
#'
#'   Use \code{\link{cv_tab}} to compare complete results with the same
#'   observation support, group assignments, prediction target and splits.
#'   Use \code{\link{select_cv}} to evaluate a training information-criterion
#'   selection rule. Its sample size for the criterion is the native number
#'   of analysis observations. Estimated covariance parameters count toward
#'   the penalty. Integrated group effects do not count individually.
#'   Native fitting, preprocessing and selection use analysis observations.
#'
#' @seealso \code{\link{cross_validate}}, \code{\link{cross_validate_field}},
#'   \code{\link{cv_tab}}, \code{\link{select_cv}}
#' @examples
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   surveys <- glmmTMB::Salamanders
#'   folds <- grouped_cv_folds(surveys$site, k = 3, ids = rownames(surveys))
#'   cv <- cross_validate_mixed(count ~ spp + (1 | site), surveys,
#'     folds[1], family = "poisson")
#'   cv$summary
#' }
#' @export
cross_validate_mixed <- function(formula, data, folds,
    family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"),
    prediction = c("new_group", "existing_group"), kernel_inputs = NULL,
    optim_args = list(), metrics = NULL, control = NULL, integration_tol = 1e-8,
    keep_fits = FALSE) {
  family <- match.arg(family)
  prediction <- match.arg(prediction)
  validate_scalar_logical(keep_fits, "keep_fits")
  if (!requireNamespace("glmmTMB", quietly = TRUE))
    stop("Mixed CV requires the glmmTMB package.", call. = FALSE)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("Supply nonempty raw data with unique column names.", call. = FALSE)
  if (is.null(control)) control <- glmmTMB::glmmTMBControl()
  if (!is.list(control)) stop("Supply native glmmTMBControl settings.", call. = FALSE)
  validate_scalar_numeric(integration_tol, "integration_tol", lower = 50 * .Machine$double.eps, upper = 0.01)
  if (is.null(metrics)) metrics <- c("neg_log_score", if (family == "binomial") "brier" else "rmse")
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) || anyDuplicated(metrics) ||
      any(!metrics %in% c("neg_log_score", "rmse", "mae", if (family == "binomial") "brier")))
    stop("Mixed CV metrics are neg_log_score, rmse, mae and, for binomial models, brier.", call. = FALSE)
  ids <- .msr_fold_ids(rownames(data), nrow(data))
  folds <- .msr_cv_folds(folds, ids)
  inputs <- kernel_inputs
  spatial <- character()
  if (!is.null(inputs)) {
    if (!inherits(inputs, "multiScaleR_data") || !is.data.frame(inputs$kernel_dat) ||
        !is.list(inputs$raw_cov) || !is.list(inputs$d_list))
      stop("Mixed scale CV requires raw-cell kernel inputs.", call. = FALSE)
    prepared_ids <- rownames(inputs$kernel_dat)
    n <- nrow(inputs$kernel_dat)
    aligned <- function(values) length(values) == n &&
      (is.null(names(values)) || identical(names(values), prepared_ids))
    if (!.msr_valid_observation_ids(prepared_ids, n) || !setequal(ids, prepared_ids) ||
        !aligned(inputs$raw_cov) || !aligned(inputs$d_list) ||
        (!is.null(inputs$raster_group) && !aligned(inputs$raster_group)))
      stop("Mixed scale CV requires matching raw-cell buffers and row IDs.", call. = FALSE)
    spatial <- inputs$scale_vars$covariate
    for (nm in setdiff(spatial, names(data))) data[[nm]] <- 0
  }
  spec <- .msr_mixed_cv_formula(formula, data, family)
  outcomes <- .msr_cv_outcome(spec$fixed, data, family)
  if (any(outcomes$trials <= 0)) stop("Mixed binomial CV requires positive trials on every row.", call. = FALSE)
  group <- as.character(data[[spec$group]])
  covariates <- intersect(spatial, union(all.vars(spec$fixed[[3L]]), all.vars(spec$random)))
  if (any(spatial %in% c(spec$group, all.vars(formula[[2L]]))))
    stop("Spatial inputs must not contain a model response or grouping column.", call. = FALSE)
  if (!is.null(inputs) && (!length(covariates) || !any(inputs$scale_vars$optimize[match(covariates, spatial)])))
    stop("Mixed scale CV needs an optimized spatial predictor in its formula.", call. = FALSE)
  allowed <- setdiff(names(formals(multiScale_optim)), c("fitted_mod", "kernel_inputs",
    "join_by", "n_cores", "PSOCK", "refit_fn", "backend", "verbose", "profile_control"))
  if (!is.list(optim_args) || (length(optim_args) && (is.null(inputs) ||
      is.null(names(optim_args)) || anyNA(names(optim_args)) || anyDuplicated(names(optim_args)) ||
      any(!names(optim_args) %in% allowed))))
    stop("Use supported serial refit optim_args only with kernel_inputs.", call. = FALSE)
  transform <- .msr_cv_covariate_transform(optim_args, "refit", !is.null(inputs))
  for (name in names(folds)) {
    split <- folds[[name]]
    new <- !group[match(split$assessment, ids)] %in% group[match(split$analysis, ids)]
    if ((prediction == "new_group" && !all(new)) || (prediction == "existing_group" && any(new)))
      stop("Fold ", name, " does not match prediction = '", prediction, "'.", call. = FALSE)
  }
  family_object <- switch(family, gaussian = stats::gaussian(), poisson = stats::poisson(),
    binomial = stats::binomial(), nbinom2 = glmmTMB::nbinom2(), Gamma = stats::Gamma("log"))
  settings <- list(family = family, link = family_object$link, prediction = prediction,
    metrics = metrics, provider = "glmmTMB", backend = if (is.null(inputs)) "native" else "refit",
    normalization = "training_only", integration_tol = integration_tol, control = control,
    optim_args = optim_args, parameter_uncertainty = "not_integrated", observation_intervals = "not_provided")
  out <- .msr_mixed_cv_output(folds, data.frame(row_id = ids, group = group,
    observed = outcomes$observed, trials = outcomes$trials), settings, keep_fits)
  out$call <- match.call()
  for (name in names(folds)) {
    split <- folds[[name]]
    ai <- match(split$analysis, ids)
    ti <- match(split$assessment, ids)
    training <- data[ai, , drop = FALSE]
    assessment <- data[ti, , drop = FALSE]
    pred <- cbind(fold = name, out$targets[ti, , drop = FALSE],
      prediction = NA_real_, eligible = FALSE)
    train <- model <- fit <- NULL
    messages <- character()
    stage <- "preparation"
    error <- tryCatch(withCallingHandlers({
      if (!is.null(inputs)) {
        train <- .msr_mixed_cv_prepare(inputs, split$analysis, covariates, transform)
        training[covariates] <- train$kernel_dat[covariates]
      }
      stage <- "native_fit"
      # Bind the native call's local inputs where update() can find them
      # during a scale refit. Only the analysis data enter that environment.
      fit_env <- new.env(parent = environment(formula))
      .msr_mixed_formula <- formula
      environment(.msr_mixed_formula) <- fit_env
      .msr_mixed_family <- family_object
      .msr_mixed_control <- control
      fit_env$.msr_mixed_formula <- .msr_mixed_formula
      fit_env$.msr_mixed_family <- .msr_mixed_family
      fit_env$.msr_mixed_control <- .msr_mixed_control
      fit_env$training <- training
      model <- glmmTMB::glmmTMB(.msr_mixed_formula, data = training,
        family = .msr_mixed_family, REML = FALSE, na.action = stats::na.fail,
        control = .msr_mixed_control)
      fit <- model
      if (!is.null(train)) {
        stage <- "scale_fit"
        fit <- do.call(multiScale_optim, c(list(fitted_mod = model, kernel_inputs = train,
          backend = "refit", n_cores = 1L, PSOCK = FALSE, verbose = FALSE), optim_args))
        if (!identical(as.integer(fit$optim_results$convergence), 0L))
          stop("The training ecological scale optimizer did not converge.", call. = FALSE)
        model <- fit$opt_mod
        training <- .msr_refit_prediction_data(fit, training, train, include_fixed = TRUE)
        assess <- .msr_cv_inputs(inputs, split$assessment, covariates,
          if (is.null(train$binned)) NULL else train$binned$cv_reference)
        assessment <- .msr_refit_prediction_data(fit, assessment, assess, include_fixed = TRUE)
      }
      stage <- "group_prediction"
      design <- .msr_mixed_cv_native(model, spec, training, assessment, family)
      scored <- .msr_mixed_cv_predict(design, out$targets[ai, , drop = FALSE],
        out$targets[ti, , drop = FALSE], family, integration_tol)
      pred$prediction <- scored$prediction
      pred$eligible <- TRUE
      out$joint_scores[[name]] <- scored$joint
      out$selection_records[[name]] <- .msr_mixed_cv_selection_record(fit, split$analysis)
      out$training[[name]] <- list(row_ids = split$analysis, coefficients = design$beta,
        covariance = design$covariance, dispersion = design$dispersion, recipe_sha256 = design$recipe_sha256,
        normalization = if (!is.null(train)) fit$scl_params,
        scale_est = if (!is.null(train)) fit$scale_est,
        bin_reference = if (!is.null(train)) train$binned$cv_reference)
      out$prediction_contracts[[name]] <- .msr_prediction_contract("native", "glmmTMB_mixed_CV",
        "native_component", "response", population = "fixed_at_training_ML_estimates",
        kernel = "fixed_at_training_fold_estimates", effects = list(group = prediction),
        conditioning = c("analysis_responses_only", "fitted_training_parameters", "frozen_training_recipes"),
        integration = list(group = if (family == "gaussian") "exact_Gaussian" else "adaptive_quadrature",
          assessment = "joint_within_group"), omitted = c("parameter_uncertainty", "selection_uncertainty"))
      if (keep_fits) out$fits[[name]] <- fit
      NULL
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(stage = stage, message = conditionMessage(e)))
    if (!is.null(error)) {
      pred$prediction <- NA_real_
      pred$eligible <- FALSE
      for (item in c("joint_scores", "training", "selection_records", "prediction_contracts"))
        out[[item]][[name]] <- NULL
    }
    out$diagnostics[[name]] <- list(convergence = model$fit$convergence,
      positive_hessian = model$sdr$pdHess, native_message = model$fit$message,
      scale_convergence = if (inherits(fit, "multiScaleR")) fit$optim_results$convergence)
    out$predictions[[name]] <- pred
    out$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai),
      n_assessment = length(ti), status = if (is.null(error)) "scored" else "failed")
    out$warnings[name] <- list(messages)
    out$errors[name] <- list(error)
  }
  .msr_mixed_cv_finish(out)
}

.msr_mixed_cv_predict <- function(design, train, test, family, tolerance) {
  result <- numeric(nrow(test))
  groups <- split(seq_len(nrow(test)), test$group, drop = TRUE)
  rows <- lapply(names(groups), function(group) {
    ai <- which(train$group == group)
    ti <- groups[[group]]
    analysis <- list(eta = design$eta_train[ai], Z = design$Ztrain[ai, , drop = FALSE],
      y = train$observed[ai], n = train$trials[ai])
    assessment <- list(eta = design$eta_test[ti], Z = design$Ztest[ti, , drop = FALSE],
      y = test$observed[ti], n = test$trials[ti])
    group_fit <- if (family == "gaussian")
      .msr_mixed_cv_gaussian_group(analysis, assessment, design$covariance, design$dispersion) else
        .msr_mixed_cv_glmm_group(analysis, assessment, family, design$dispersion,
          sqrt(design$covariance[1L, 1L]), tolerance)
    result[ti] <<- group_fit$prediction
    data.frame(group = group, n = length(ti), neg_log_score = group_fit$neg_log_score)
  })
  joint <- do.call(rbind, rows)
  if (any(!is.finite(result)) || any(!is.finite(joint$neg_log_score)))
    stop("Mixed-model integration did not yield finite predictions and group log scores.", call. = FALSE)
  list(prediction = result, joint = joint)
}

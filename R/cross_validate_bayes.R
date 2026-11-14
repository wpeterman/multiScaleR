# Keep one copy of a prepared buffer per point, while retaining every observation.
.msr_bayes_cv_inputs <- function(inputs, point_ids, keys, bin_reference, reference = NULL) {
  selected <- unique(keys)
  rows <- match(selected, point_ids)
  out <- inputs
  out$kernel_dat <- if (is.null(inputs$kernel_dat)) data.frame(row.names = selected) else
    inputs$kernel_dat[rows, , drop = FALSE]
  out$raw_cov <- inputs$raw_cov[rows]
  out$d_list <- inputs$d_list[rows]
  if (!is.null(inputs$raster_group)) out$raster_group <- inputs$raster_group[selected]
  out$scl_params <- list(mean = numeric(), sd = numeric())
  if (!is.null(inputs$binned)) {
    if (bin_reference == "training") {
      out$binned <- .msr_cv_bins(out, reference)
    } else {
      for (nm in c("vsum", "csum"))
        out$binned[[nm]] <- lapply(inputs$binned[[nm]], function(z) z[rows, , drop = FALSE])
      out$binned$point_ids <- selected
      out$binned$n_points <- length(selected)
    }
  }
  out
}

.msr_bayes_cv_fixed <- function(inputs, data, keys, formula, scales) {
  specs <- inputs$scale_vars
  if (is.null(specs)) return(list(inputs = inputs, data = data))
  used <- unique(c(all.vars(formula[[3L]]),
    unlist(lapply(scales, all.vars), use.names = FALSE)))
  if ("." %in% used) used <- union(used, setdiff(names(data), all.vars(formula[[2L]])))
  fixed <- specs$covariate[!specs$optimize & specs$covariate %in% used]
  if (!length(fixed)) return(list(inputs = inputs, data = data))
  if (is.null(inputs$raw_cov))
    stop("Fixed spatial predictors need raw cells for training-only CV standardization.", call. = FALSE)
  context <- list(covs = character(), n_covs = 0L, scale_vars = specs,
    unit_conv = inputs$unit_conv, resolution = inputs$resolution, n_cols = inputs$n_cols,
    binned = inputs$binned)
  raw <- .msr_kernel_cov_w(numeric(), inputs$d_list, inputs$raw_cov,
    inputs$kernel, context, covariates = fixed)
  rows <- match(keys, rownames(inputs$kernel_dat))
  center <- colMeans(raw[rows, , drop = FALSE])
  spread <- apply(raw[rows, , drop = FALSE], 2L, stats::sd)
  if (any(!is.finite(center) | !is.finite(spread) | spread <= 0))
    stop("Fixed spatial predictors must vary and be finite on the training observations.", call. = FALSE)
  inputs$scl_params <- list(mean = center, sd = spread)
  inputs$kernel_dat[fixed] <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"))
  data[fixed] <- inputs$kernel_dat[rows, fixed, drop = FALSE]
  list(inputs = inputs, data = data)
}

# Weighted log-sum-exp, including a component with exactly zero density.
.msr_bayes_cv_logadd <- function(a, b) {
  out <- pmax(a, b)
  finite <- is.finite(out)
  out[finite] <- out[finite] + log1p(exp(-abs(a[finite] - b[finite])))
  out
}

.msr_bayes_cv_scores <- function(fit, context, target, new_group_draws, seed) {
  if (is.null(fit$importance)) {
    draws <- fit$draws
    weight <- rep(1 / nrow(draws$beta), nrow(draws$beta))
  } else {
    weight <- fit$importance$weights
    draws <- .msr_bayes_parameter_blocks(fit$importance$candidates, fit$parameters)
  }
  if (!length(weight) || length(weight) != nrow(draws$beta) ||
      any(!is.finite(weight) | weight < 0) || !any(weight > 0))
    stop("Invalid posterior draws or importance weights.", call. = FALSE)
  indices <- which(weight > 0)
  weight <- weight / sum(weight)
  # Each augmented vector supplies one value for each new group in the whole
  # assessment set. Joint blocks sum conditional log densities before averaging.
  spec <- context$spec
  n <- nrow(target)
  group <- match(target$assessment_group, unique(target$assessment_group))
  n_group <- max(group)
  integrate_new <- (context$scale_policy == "sample" && context$n_scale > 0L) ||
    (context$response_policy == "sample" && context$n_response > 0L)
  repeats <- if (integrate_new) new_group_draws else 1L
  count <- spec$family %in% c("poisson", "nbinom2")
  .msr_bayes_prediction_seed(seed, {
    marginal <- rep(-Inf, n)
    joint <- rep(-Inf, n_group)
    mean_value <- if (count) rep(-Inf, n) else numeric(n)
    for (i in indices) {
      parameters <- lapply(draws, function(z) unname(z[i, ]))
      w <- weight[i] / repeats
      logw <- log(weight[i]) - log(repeats)
      for (j in seq_len(repeats)) {
        p <- .msr_bayes_prediction_parameters(parameters, context)
        eta <- .msr_bayes_components(spec, p)$eta
        if (any(!is.finite(eta))) stop("Nonfinite assessment linear predictors.", call. = FALSE)
        ll <- .msr_scale_response_loglik(spec$family, target$observed, eta,
          p$log_dispersion, target$trials)
        if (anyNA(ll) || any(ll == Inf))
          stop("Invalid assessment response log densities.", call. = FALSE)
        marginal <- .msr_bayes_cv_logadd(marginal, ll + logw)
        block <- drop(rowsum(ll, group, reorder = FALSE))
        joint <- .msr_bayes_cv_logadd(joint, block + logw)
        mean_value <- if (count) .msr_bayes_cv_logadd(mean_value, eta + logw) else
          mean_value + w * if (spec$family == "binomial") stats::plogis(eta) else eta
      }
    }
    list(prediction = unname(if (count) exp(mean_value) else mean_value),
      marginal_neg_log_score = unname(-marginal), joint = unname(-joint),
      block_size = tabulate(group),
      posterior_vectors = length(indices), new_group_draws = repeats)
  })
}

.msr_bayes_cv_diagnostics <- function(fit, method) {
  raw <- fit$diagnostics
  derived <- if (method %in% c("sample", "slice") && !is.na(raw$passed))
    diagnose_bayes_scales(fit) else NULL
  passed <- if (is.null(derived)) raw$passed else derived$passed
  list(passed = passed, raw = raw, derived = derived,
    inference = fit$inference, importance_history = fit$importance$history,
    target_evaluation = fit$importance$target_evaluation)
}

#' Cross-Validate Bayesian Scale Models
#'
#' Refit a Bayesian scale model in each explicit analysis split. Score the
#' posterior predictive density of the held-out observations, with coherent
#' random effects shared by observations in the same predicted group.
#'
#' This function is experimental. Its arguments and results may change in a
#' later release without a deprecation period.
#'
#' @param formula Two-sided Bayesian response formula. Use a data column as the
#'   response, or \code{cbind(successes, failures)} for grouped binomial data.
#'   An \code{I()} wrapper is also accepted. Predictor and scale formulas use
#'   \code{\link{multiScale_bayes}} syntax.
#' @param data Data frame of raw observations and formula variables. Unique row
#'   names identify observations in \code{folds}, even if several observations
#'   use the same prepared point.
#' @param kernel_inputs Spatial inputs from \code{\link{kernel_prep}} or
#'   \code{\link{kernel_prep_by_group}}. For training-only rebinning, prepare
#'   kernel means with \code{bin = TRUE, store_cell_data = TRUE}. Weighted SQ and
#'   fixed spatial metrics use raw cells.
#' @param folds Uniquely named list of analysis and assessment splits. Each
#'   split contains nonempty, disjoint character vectors \code{analysis} and
#'   \code{assessment} of observation row names in \code{data}.
#' @param scales Named scale formulas as in \code{\link{multiScale_bayes}}.
#'   Default \code{NULL} estimates one shared scale per spatial predictor.
#' @param family Response family. Use \code{"nbinom2"} (default) or
#'   \code{"poisson"} for counts, \code{"binomial"} for successes and trials,
#'   or \code{"gaussian"} for continuous outcomes.
#' @param row_id Optional name of the column holding prepared point IDs. The
#'   default uses the observation row names as point IDs. Repeated point IDs
#'   are allowed and do not replace the unique observation row names.
#' @param priors Named prior settings passed to \code{\link{multiScale_bayes}}.
#'   Specify them before looking at assessment outcomes. Default \code{list()}
#'   uses the model's documented defaults.
#' @param method Posterior method. The default \code{"slice"} runs without
#'   optional sampling packages. \code{"sample"} runs NUTS. \code{"laplace"}
#'   and \code{"marginal"} use the respective posterior approximations,
#'   including direct importance-weighted scoring when correction is requested.
#' @param fit_args Named list of additional fitting arguments, chosen from
#'   \code{chains}, \code{iter}, \code{warmup}, \code{control},
#'   \code{laplace_draws} and \code{cores}. Defaults follow
#'   \code{\link{multiScale_bayes}}. Folds run sequentially.
#' @param assessment_group Optional name of a data column defining joint
#'   scoring blocks within each assessment split. Default \code{NULL} scores
#'   the whole assessment split jointly. Choose the blocks to represent the
#'   prediction question before observing assessment outcomes.
#' @param new_scale_levels,new_response_levels Separate unseen-level policies
#'   for scale and response random effects. Default \code{"sample"} integrates
#'   new group variation by simulation. \code{"error"} refuses an unseen level
#'   and \code{"zero"} predicts at its zero deviation. Existing groups retain
#'   their fitted posterior effects.
#' @param new_group_draws Positive integer number of new-effect simulations per
#'   posterior vector when sampling unseen levels. Default one. Increase it to
#'   reduce conditional Monte Carlo error from new-effect simulation. For rare
#'   joint responses or concentrated posterior weights, check score stability
#'   at larger budgets with the same seed and fitting settings.
#' @param bin_reference \code{"training"} (default) computes representative
#'   distances from unique analysis points and freezes them for assessment.
#'   Use \code{"fixed"} only when supplied bins and representative distances
#'   were fixed independently of all analysis and assessment data.
#' @param seed Optional positive integer seed for distinct fold fitting and
#'   scoring streams. A supplied seed preserves the caller's R random state.
#' @param keep_fits Logical. Default \code{FALSE} stores compact fold records.
#'   Use \code{TRUE} to retain complete Bayesian fits for investigation.
#' @return A \code{multiScaleR_bayes_cv} object inheriting
#'   \code{multiScaleR_cv}. \code{predictions} contains the posterior mean
#'   response and marginal negative log score per observation. For binomial
#'   models, \code{prediction} is a per-trial probability and log scores use
#'   the observed count and trial total. \code{joint_scores} contains each
#'   fold's block sizes and joint negative log scores. \code{summary} and
#'   \code{metrics} report pooled and per-fold joint negative log score divided
#'   by the number of assessment observations. \code{fold_results} records
#'   which folds scored, failed fitting or scoring, or failed diagnostics.
#'   \code{training} records analysis and point IDs, frozen normalization,
#'   bin references, priors and random seeds. \code{diagnostics},
#'   \code{prediction_contracts}, \code{warnings}, \code{errors},
#'   \code{targets}, \code{folds}, \code{settings}, \code{assessment_counts},
#'   \code{complete}, and optional \code{fits} support audit and comparison.
#' @details The joint score is the negative log of the posterior average of
#'   the conditional product of observation densities in each assessment
#'   block. For grouped binomial data, the block product uses one normalized
#'   count probability per observation. New effects are drawn once per
#'   group in each joint parameter vector and shared by all corresponding
#'   assessment observations. The sum of marginal row scores is a different
#'   target. Pooled scores weight each joint block by its observation count
#'   after expressing its score per observation.
#'
#'   Each fit learns normalization, formula bases, coefficients and ecological
#'   scales only from analysis rows. Repeated visits can share one point
#'   buffer. Use whole-site or whole-year folds if prediction concerns an
#'   unseen site or year. Where several groups cross, the assessment block
#'   should cover the observations whose joint response you want to predict.
#'
#'   Failed native or derived diagnostics keep the fold incomplete, with any
#'   computed scores available for investigation. Slice sampling without the
#'   optional posterior package can score folds with diagnostics recorded as
#'   unavailable. Install posterior to check mixing and simulation precision.
#'   Laplace and marginal scores inherit their fitted approximation and
#'   diagnostic limits. Importance correction uses normalized candidate
#'   weights directly. Numerical diagnostic checks do not establish
#'   predictive calibration.
#'
#'   Use \code{\link{cv_tab}} with \code{metric = "neg_log_score"} for complete
#'   Bayesian results on the same targets, folds, blocks and new-level policies.
#'   Poisson and NB2 can compare predictive mass on the same count observations.
#'   Training AIC, AICc and BIC selection does not apply to this posterior
#'   workflow. Use an outer split to evaluate a model-selection rule.
#' @seealso \code{\link{multiScale_bayes}}, \code{\link{grouped_cv_folds}},
#'   \code{\link{diagnose_bayes_scales}}, \code{\link{cv_tab}}
#' @export
cross_validate_bayes <- function(formula, data, kernel_inputs, folds, scales = NULL,
    family = c("nbinom2", "poisson", "binomial", "gaussian"), row_id = NULL,
    priors = list(), method = c("slice", "sample", "laplace", "marginal"),
    fit_args = list(), assessment_group = NULL, new_scale_levels = c("sample", "error", "zero"),
    new_response_levels = c("sample", "error", "zero"), new_group_draws = 1L,
    bin_reference = c("training", "fixed"), seed = NULL, keep_fits = FALSE) {
  family <- match.arg(family); method <- match.arg(method)
  bin_reference <- match.arg(bin_reference)
  new_scale_levels <- match.arg(new_scale_levels); new_response_levels <- match.arg(new_response_levels)
  validate_scalar_numeric(new_group_draws, "new_group_draws", integerish = TRUE, lower = 1)
  validate_scalar_logical(keep_fits, "keep_fits")
  if (!is.data.frame(data) || nrow(data) < 3L || anyDuplicated(names(data)))
    stop("Supply a data frame of at least three observations with distinct columns.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !all(all.vars(formula[[2L]]) %in% names(data)))
    stop("Supply a two-sided formula with explicit response columns in data.", call. = FALSE)
  response <- all.vars(formula[[2L]])
  ids <- rownames(data)
  folds <- .msr_cv_folds(folds, ids)
  if (any(vapply(folds, function(f) length(f$analysis) < 2L, logical(1))))
    stop("Bayesian CV needs at least two analysis observations in every fold.", call. = FALSE)
  if (!is.list(fit_args) || (length(fit_args) && (is.null(names(fit_args)) ||
      anyDuplicated(names(fit_args)) ||
      !all(names(fit_args) %in% c("chains", "iter", "warmup", "control", "laplace_draws", "cores")))))
    stop("fit_args must be a uniquely named list of Bayesian fitting controls.", call. = FALSE)
  if (!is.null(scales) && !is.list(scales))
    stop("scales must be a named list of scale formulas.", call. = FALSE)
  point_ids <- .msr_scale_point_ids(kernel_inputs)
  if (!is.null(row_id) && (!is.character(row_id) || length(row_id) != 1L ||
      !row_id %in% names(data) || row_id %in% response))
    stop("row_id must name a prepared point-ID column in data.", call. = FALSE)
  point_keys <- if (is.null(row_id)) ids else data[[row_id]]
  if (!is.atomic(point_keys) || !is.null(dim(point_keys)))
    stop("Prepared point keys must be an atomic vector.", call. = FALSE)
  point_keys <- as.character(point_keys)
  used <- sort(unique(unlist(folds, use.names = FALSE)))
  used_rows <- match(used, ids)
  if (!all(point_keys[used_rows] %in% point_ids))
    stop("Every used observation must match a prepared point ID.", call. = FALSE)
  if (!is.null(kernel_inputs$kernel_dat) &&
      !identical(rownames(kernel_inputs$kernel_dat), point_ids))
    stop("Prepared point IDs must match the spatial covariate rows in order.", call. = FALSE)
  raw <- !is.null(kernel_inputs$raw_cov) && !is.null(kernel_inputs$d_list)
  if (raw && (!identical(names(kernel_inputs$raw_cov), point_ids) ||
      !identical(names(kernel_inputs$d_list), point_ids)))
    stop("Raw spatial buffers must be named in prepared point-ID order.", call. = FALSE)
  if (!is.null(kernel_inputs$raster_group) &&
      !identical(names(kernel_inputs$raster_group), point_ids))
    stop("Raster groups must be named in prepared point-ID order.", call. = FALSE)
  if (bin_reference == "training" && !is.null(kernel_inputs$binned) &&
      (!raw || is.null(kernel_inputs$scale_vars)))
    stop("Training binning needs declared spatial predictors and raw cells. Supply store_cell_data = TRUE.", call. = FALSE)
  if (!is.null(assessment_group) && (!is.character(assessment_group) ||
      length(assessment_group) != 1L || !assessment_group %in% names(data) ||
      assessment_group %in% response))
    stop("assessment_group must name a column defining assessment blocks.", call. = FALSE)
  group <- if (is.null(assessment_group)) rep("whole_fold", length(used)) else
    data[[assessment_group]][used_rows]
  if (!is.atomic(group) || !is.null(dim(group)) || anyNA(group) || any(!nzchar(as.character(group))))
    stop("Assessment block labels must be nonmissing and nonempty.", call. = FALSE)
  observed <- .msr_cv_outcome(formula, data[used_rows, , drop = FALSE], family)
  if (any(observed$trials < 1L))
    stop("Binomial CV requires a positive number of trials per observation.", call. = FALSE)
  targets <- data.frame(row_id = used, point_id = point_keys[used_rows],
    observed = observed$observed, trials = observed$trials, assessment_group = as.character(group))
  seeds <- .msr_bayes_prediction_seed(seed, matrix(sample.int(.Machine$integer.max,
    2L * length(folds)), ncol = 2L, dimnames = list(names(folds), c("fit", "score"))))
  out <- list(predictions = list(), joint_scores = list(), fold_results = list(), training = list(),
    diagnostics = list(), prediction_contracts = list(), warnings = list(), errors = list(),
    fits = if (keep_fits) list() else NULL)
  for (name in names(folds)) {
    f <- folds[[name]]
    a <- match(f$analysis, ids); b <- match(f$assessment, ids)
    target <- targets[match(f$assessment, used), , drop = FALSE]
    train_target <- targets[match(f$analysis, used), , drop = FALSE]
    pred <- cbind(data.frame(fold = name), target, prediction = NA_real_,
      marginal_neg_log_score = NA_real_, new_scale_level = NA, new_response_level = NA,
      eligible = FALSE)
    blocks <- unique(target$assessment_group)
    joint <- data.frame(fold = name, assessment_group = blocks,
      n = tabulate(match(target$assessment_group, blocks)), neg_log_score = NA_real_)
    fit <- NULL; messages <- character(); stage <- "preparation"; status <- "failed"
    out$diagnostics[name] <- out$prediction_contracts[name] <- list(NULL)
    out$training[[name]] <- list(analysis = f$analysis, point_ids = point_keys[a],
      fit_seed = seeds[name, "fit"], score_seed = seeds[name, "score"])
    error <- tryCatch(withCallingHandlers({
      train <- .msr_bayes_cv_inputs(kernel_inputs, point_ids, point_keys[a], bin_reference)
      analysis <- data[a, , drop = FALSE]
      # Preserve the declared success level even when only one outcome occurs
      # in an analysis split and native model.frame drops unused factor levels.
      if (family == "binomial" && length(response) == 1L && is.factor(analysis[[response]]))
        analysis[[response]] <- train_target$observed
      prepared <- .msr_bayes_cv_fixed(train, analysis, point_keys[a], formula, scales)
      train <- prepared$inputs
      stage <- "fit"
      fit <- .msr_bayes_prediction_seed(seeds[name, "fit"],
        do.call(multiScale_bayes, c(list(formula = formula, inputs = train,
          data = prepared$data, scales = scales, family = family, row_id = row_id,
          priors = priors, method = method, seed = seeds[name, "fit"]), fit_args)))
      spec <- fit$specification
      if (!identical(as.numeric(spec$y), train_target$observed) ||
          !identical(as.numeric(spec$trials), train_target$trials) ||
          !identical(spec$observation_ids, point_keys[a]))
        stop("The fit did not preserve the analysis response, trials and point mapping.", call. = FALSE)
      out$training[[name]] <- c(out$training[[name]], list(center = spec$center, spread = spec$spread,
        fixed_center = spec$fixed_reference$center, fixed_spread = spec$fixed_reference$spread,
        bin_distances = train$binned$bin_rep, bin_reference = train$binned$cv_reference,
        priors = spec$priors))
      stage <- "diagnostics"
      out$diagnostics[[name]] <- .msr_bayes_cv_diagnostics(fit, method)
      stage <- "prediction"
      assess <- .msr_bayes_cv_inputs(kernel_inputs, point_ids, point_keys[b], bin_reference,
        train$binned$cv_reference)
      # Response columns cannot enter the assessment design. Grouped trial
      # counts are supplied separately to the normalized response density.
      future <- data[b, setdiff(names(data), response), drop = FALSE]
      context <- .msr_bayes_prediction_prepare(fit, future, assess, row_id,
        new_scale_levels, new_response_levels, TRUE, TRUE, include_fixed = TRUE)
      stage <- "score"
      score <- .msr_bayes_cv_scores(fit, context, target, new_group_draws, seeds[name, "score"])
      pred$prediction <- score$prediction
      pred$marginal_neg_log_score <- score$marginal_neg_log_score
      pred$new_scale_level <- Reduce(`|`, context$scale_flags)
      pred$new_response_level <- context$response_flags
      joint$neg_log_score <- score$joint
      stopifnot(identical(joint$n, score$block_size))
      out$prediction_contracts[[name]] <- list(new_groups = context$new_groups,
        posterior_vectors = score$posterior_vectors, new_group_draws = score$new_group_draws,
        importance_weighted = !is.null(fit$importance))
      passed <- out$diagnostics[[name]]$passed
      if (isFALSE(passed)) {
        status <- "diagnostic_failed"
        list(stage = "diagnostics", message = "The fold did not pass its posterior diagnostics.")
      } else {
        status <- "scored"
        pred$eligible <- TRUE
        NULL
      }
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(stage = stage, message = conditionMessage(e)))
    out$predictions[[name]] <- pred
    out$joint_scores[[name]] <- joint
    out$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(a),
      n_assessment = length(b), status = status)
    out$warnings[name] <- list(unique(messages))
    out$errors[name] <- list(error)
    if (keep_fits) out$fits[name] <- list(fit)
  }
  out$targets <- targets
  out$folds <- folds
  out$settings <- list(family = family, link = switch(family, gaussian = "identity",
    binomial = "logit", "log"), method = method, fit_args = fit_args, seed = seed,
    metrics = "neg_log_score", prediction_target = "posterior_predictive_joint_blocks",
    new_scale_levels = new_scale_levels, new_response_levels = new_response_levels,
    assessment_group = assessment_group, new_group_draws = new_group_draws,
    normalization = "training_only", bin_reference = bin_reference)
  out$call <- match.call()
  .msr_bayes_cv_finish(out, ids)
}

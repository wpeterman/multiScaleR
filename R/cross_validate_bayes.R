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
#' Refit a [multiScale_bayes()] model in each analysis split, then score the
#' held-out observations with the posterior predictive density. Observations
#' that share a predicted group, such as a held-out site or year, also share
#' its random effect when they are scored.
#'
#' This function is experimental. Its arguments and results may change in a
#' later release without a deprecation period.
#'
#' @param formula Two-sided response formula in [multiScale_bayes()] syntax.
#'   Use a data column as the response, or `cbind(successes, failures)` for
#'   grouped binomial data. An `I()` wrapper is also accepted.
#' @param data Data frame of raw observations and formula variables. Its
#'   unique row names identify observations in `folds`, even when several
#'   observations share one prepared point.
#' @param kernel_inputs Prepared landscape inputs from [kernel_prep()] or
#'   [kernel_prep_by_group()]. For the default training-only rebinning,
#'   prepare kernel means with `bin = TRUE` and `store_cell_data = TRUE` (both
#'   defaults of [kernel_prep()]). Weighted SQ and fixed-radius metrics also
#'   use the raw cells.
#' @param folds Uniquely named list of splits. Each split holds nonempty,
#'   disjoint character vectors `analysis` and `assessment` of row names of
#'   `data`, as returned by [grouped_cv_folds()] or [spatial_cv_folds()].
#'   Every analysis set needs at least two observations.
#' @param scales Named list of scale formulas, as in [multiScale_bayes()]. The
#'   default `NULL` estimates one shared scale per landscape predictor.
#' @param family Response family, one of `"nbinom2"` (default) or `"poisson"`
#'   for counts, `"binomial"` for successes and trials, or `"gaussian"` for
#'   continuous outcomes.
#' @param row_id Optional name of the column of `data` that holds prepared
#'   point IDs. The default `NULL` uses the row names of `data` as point IDs.
#'   Point IDs may repeat. They do not replace the unique row names that
#'   identify observations.
#' @param priors Named list of prior settings passed to [multiScale_bayes()].
#'   Set them before you look at any assessment outcome. The default `list()`
#'   uses the documented defaults.
#' @param method Computation for every fold. The default `"slice"` needs no
#'   optional sampling packages. `"sample"` runs NUTS. `"laplace"` and
#'   `"marginal"` use the approximations, and an importance-corrected
#'   marginal fit is scored directly with its candidate weights.
#' @param fit_args Named list of further fitting arguments for
#'   [multiScale_bayes()], chosen from `chains`, `iter`, `warmup`, `control`,
#'   `laplace_draws`, and `cores`. Settings you leave out take the defaults of
#'   [multiScale_bayes()]. Folds always run one after another.
#' @param assessment_group Optional name of a column of `data` that defines
#'   joint scoring blocks within each assessment split. The default `NULL`
#'   scores each whole assessment split as one block. Choose the blocks to
#'   match your prediction question before you see any assessment outcome.
#' @param new_scale_levels,new_response_levels What to do with a group level
#'   that is absent from an analysis split, set separately for scale and
#'   response random effects. The default `"sample"` averages over new-group
#'   variation by simulation. `"error"` refuses an unseen level, and `"zero"`
#'   predicts at a zero deviation. Groups present in the analysis split keep
#'   their fitted effects.
#' @param new_group_draws Positive integer number of new-effect simulations
#'   per posterior vector when unseen levels are sampled; default 1. Raise it
#'   to reduce the simulation error of the average over new groups. See
#'   "Compare models".
#' @param bin_reference `"training"` (default) computes the representative
#'   bin distances from the unique analysis points of each fold and freezes
#'   them for its assessment points. It needs the raw cells, so it refuses
#'   inputs prepared with `store_cell_data = FALSE`, such as lean inputs from
#'   [kernel_prep_by_group()]. `"fixed"` uses the prepared bin distances as
#'   given, so lean inputs of kernel means work with it. Use `"fixed"` only
#'   when the supplied bins were fixed independently of all analysis and
#'   assessment data.
#' @param seed Optional positive integer seed. It generates distinct fitting
#'   and scoring seeds for every fold, and your own random-number state is
#'   restored afterward.
#' @param keep_fits Logical. The default `FALSE` stores compact fold records
#'   only. `TRUE` also keeps every fold's complete fit in `fits`, for
#'   investigation, at a cost in memory.
#'
#' @return An object of class `multiScaleR_bayes_cv`, which inherits from
#'   `multiScaleR_cv` (see [multiScaleR_cv-class]). It is a list with these
#'   elements.
#'   * `summary`: a one-row data frame with the pooled score. `estimate` is the
#'     sum of joint negative log scores divided by the number of assessment
#'     observations (`unit` is `"joint_block_density_per_observation"`), so
#'     lower is better. `n_assessment`, `n_scored`, `weight_sum`, and
#'     `complete` describe its basis. `estimate` is `NA` unless every fold
#'     scored.
#'   * `metrics`: the same score for each `fold`.
#'   * `predictions`: one row per assessed observation, with `fold`,
#'     `row_id`, `point_id`, `observed`, `trials`, `assessment_group`, the
#'     posterior mean response `prediction` (a per-trial probability for
#'     binomial models), `marginal_neg_log_score` for that observation alone,
#'     `new_scale_level` and `new_response_level` flags, and `eligible`, which
#'     is `TRUE` only in a scored fold.
#'   * `joint_scores`: one row per fold and assessment block, with its size
#'     `n` and joint `neg_log_score`.
#'   * `fold_results`: one row per fold, with `n_analysis`, `n_assessment`,
#'     and `status`, which is `"scored"`, `"diagnostic_failed"`, or
#'     `"failed"` (an error at any stage, recorded in `errors`).
#'   * `diagnostics`: for each fold, the fit's stored checks (`raw`), the
#'     [diagnose_bayes_scales()] report for sampled fits (`derived`), the
#'     overall `passed` flag, and the importance records of corrected fits.
#'   * `errors`, `warnings`: for each fold, the failing `stage` and `message`
#'     (or `NULL`), and the warnings captured during the fold.
#'   * `training`: for each fold, the analysis rows and point IDs, the fit and
#'     score seeds, the frozen normalization (`center`, `spread`, and those of
#'     fixed-radius metrics), the bin distances and their reference, and the
#'     priors.
#'   * `prediction_contracts`: for each fold, the unseen groups, the number of
#'     posterior vectors, `new_group_draws`, and whether importance weights
#'     were used.
#'   * `targets`, `folds`, `settings`, `assessment_counts`, `complete`, `call`,
#'     `integrity_sha256`: the scored observations, the validated splits, the
#'     settings, how often each row was assessed, whether every fold scored,
#'     the call, and a checksum that [cv_tab()] uses to detect altered
#'     records.
#'   * `fits`: the complete fold fits when `keep_fits = TRUE`, and `NULL`
#'     otherwise.
#'
#' @details
#' # What the score measures
#'
#' The joint score of a block is the negative log of the posterior average of
#' the product of the observation densities in that block. Lower is better.
#' New effects are drawn once per group in each joint parameter vector and
#' shared by every assessment observation in that group. For grouped binomial
#' data, each observation contributes one normalized probability of its count
#' of successes. The sum of the marginal row scores is a different target.
#' The pooled score weights each block by its number of observations. It is
#' the only cross-validation score in the package that averages over
#' uncertainty in the coefficients, the scales, and the group effects.
#'
#' # Choose folds for your prediction question
#'
#' Each fit learns the normalization, formula bases, coefficients, and
#' scales from its analysis rows only. Repeated visits can share one point
#' buffer. To predict new sites or years, hold out whole sites or years, for
#' example with [grouped_cv_folds()] on the site labels. Where groups cross,
#' make each assessment block cover the observations whose joint response you
#' want to predict.
#' `vignette("v08-cross-validation", package = "multiScaleR")` discusses fold
#' design.
#'
#' # Handle failed folds
#'
#' A fold that fails its native or derived diagnostics is marked
#' `"diagnostic_failed"` and leaves the result incomplete, though its computed
#' scores stay available for investigation. Slice sampling without the
#' posterior package can score folds with diagnostics recorded as
#' unavailable, so install posterior to check mixing and simulation
#' precision. Laplace and marginal scores inherit the limits of their
#' approximation and its checks. Importance-corrected fits are scored with
#' their normalized candidate weights. Passing numerical checks does not
#' establish predictive calibration. Every fold runs the full fitting
#' method, so a run takes about as long as one fit per fold. When any fold
#' is not scored, the function returns the result and gives one warning that
#' counts the failed folds and quotes their first distinct reasons. Read
#' `errors` for the stage and message of every failed fold.
#'
#' # Compare models
#'
#' Use [cv_tab()] with `metric = "neg_log_score"` to rank complete Bayesian
#' results computed on the same targets, folds, blocks, and new-level
#' policies. Poisson and NB2 results can be compared, because both score
#' probability mass on the same counts. Training AIC, AICc, and BIC selection
#' does not apply to posterior fits, so use an outer split to evaluate a
#' model-selection rule. The average over unseen group effects uses a finite
#' number of simulated effects. When two candidates score very closely, rerun
#' both with a larger `new_group_draws` and the same `seed` and fitting
#' settings, and treat the comparison as unresolved if the ranking changes.
#'
#' @seealso [multiScale_bayes()], [grouped_cv_folds()],
#'   [spatial_cv_folds()], [diagnose_bayes_scales()], [cv_tab()],
#'   [multiScaleR_cv-class];
#'   `vignette("v08-cross-validation", package = "multiScaleR")`.
#' @examples
#' \donttest{
#' data("landscape_counts")
#' data("surv_pts")
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))
#' # Keep the raw cells (the default) so each fold rebins from its own
#' # analysis points
#' inputs <- kernel_prep(pts = terra::vect(surv_pts), raster_stack = land[["land1"]],
#'                       max_D = 1700, kernel = "gaussian", verbose = FALSE)
#'
#' # Two folds of alternating rows
#' folds <- grouped_cv_folds(rep(1:2, length.out = nrow(landscape_counts)),
#'                           k = 2, ids = rownames(landscape_counts))
#'
#' # Very short slice chains keep this example quick. With the posterior
#' # package installed, they usually fail the diagnostics, so the folds are
#' # marked "diagnostic_failed", the pooled score is NA, and a warning says
#' # so. Use the default chains and iter in an analysis.
#' cv <- cross_validate_bayes(counts ~ land1, landscape_counts, inputs, folds,
#'                            family = "poisson",
#'                            priors = list(scale_median = 300),
#'                            method = "slice",
#'                            fit_args = list(chains = 2, iter = 200), seed = 1)
#' cv$fold_results    # status of each fold
#' cv$joint_scores    # joint negative log score of each assessment block
#' cv$summary         # pooled score per observation, lower is better
#' }
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
    stop("Training binning needs declared spatial predictors and raw cells. Supply inputs prepared with store_cell_data = TRUE, or set bin_reference = 'fixed' to score with the prepared bins as given (appropriate only when those bins were fixed independently of the analysis and assessment data).", call. = FALSE)
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
  out <- .msr_bayes_cv_finish(out, ids)
  .msr_cv_warn_incomplete(out, "cross_validate_bayes")
  out
}

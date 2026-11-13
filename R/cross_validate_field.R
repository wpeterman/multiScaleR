# Reconstruct only protected training data through native kernel refits.
.msr_field_cv_refitter <- function(formula, original, xy, covariance, control,
    audit_dir, max_memory_mb, provider, spatial, max_refits, family = "gaussian") {
  state <- new.env(parent = emptyenv()); state$count <- 0L; state$exhausted <- FALSE
  callback <- function(model, data, context) {
    state$count <- state$count + 1L
    if (state$count > max_refits) {
      state$exhausted <- TRUE
      stop("The field CV fold exceeded max_refits.", call. = FALSE)
    }
    if (!is.data.frame(data) || !identical(rownames(data), rownames(original)) ||
        any(!spatial %in% names(data)))
      stop("Field refits must retain original training rows and declared ecological predictors.", call. = FALSE)
    protected <- intersect(setdiff(names(original), spatial), names(data))
    if (!identical(data[protected], original[protected]))
      stop("A field refit changed protected training responses or nonspatial inputs.", call. = FALSE)
    updated <- original
    updated[spatial] <- data[spatial]
    .msr_field_cv_fit(formula, updated, xy, covariance, control, audit_dir, max_memory_mb,
      provider, family)
  }
  scope <- list2env(list(formula = formula, original = original, xy = xy,
    covariance = covariance, control = control, audit_dir = audit_dir,
    max_memory_mb = max_memory_mb, provider = provider, spatial = spatial,
    max_refits = max_refits, family = family, state = state), parent = environment(.msr_field_cv_refitter))
  environment(callback) <- scope
  list(callback = callback, state = state)
}

#' Cross-Validate a Spatial Response Field
#'
#' Estimate a shared response field and its correlation range in each training
#' split, then score the joint distribution of the held-out block. Gaussian
#' responses use exact Gaussian conditioning at the fitted parameters. Count,
#' binary and positive responses use Laplace integration. Optional ecological
#' kernel scales are also refitted inside training.
#' @param formula Two-sided fixed-effects formula with a numeric response column.
#'   Products, interactions, factors, offsets and trained scale/poly/ns/bs terms
#'   retain their native meanings. Do not add random effects or penalized smooths;
#'   this adapter adds one shared response field to the fixed formula.
#' @param data Raw data frame with unique row names matching the fold IDs and
#'   coordinate order. Supply raw predictors, not full-data learned transforms.
#'   All observations have unit case weights.
#' @param coords Projected points or a numeric two-column table, one row per data
#'   row. Repeated locations identify shared field states.
#' @param folds Explicit named analysis/assessment ID splits, usually from
#'   [spatial_cv_folds()]. A location cannot occur in both parts of one split.
#'   Buffer rows excluded from that split enter neither its fit nor its score.
#' @param audit_dir Directory for compact write-once native fit/start records.
#'   Failures remain recorded. No models, tapes or dense covariance arrays are
#'   saved to these records.
#' @param crs Projected coordinate reference system for coordinate tables. Default
#'   `NULL` uses the existing CRS of supplied spatial points. Units are retained.
#' @param covariance Response correlation, `"exponential"` (default) or
#'   `"gaussian"`. Correlation is exp(-d/range) or exp(-(d/range)^2), respectively.
#'   Range and field SD are estimated in every fold, as is the observation
#'   dispersion parameter for Gaussian, NB2 and Gamma models.
#' @param kernel_inputs Optional original [kernel_prep()] inputs with raw cell
#'   data and matching IDs. Default `NULL` fits a spatial-field reference without
#'   estimating ecological scales. With inputs, use [msr_vars()] to identify the
#'   optimized predictors appearing in formula.
#' @param optim_args Named serial refit controls passed to [multiScale_optim()].
#'   Default an empty list. Backend, data, callback and worker overrides are
#'   reserved; this adapter uses the refit backend and one worker.
#' @param control Named native nlminb controls. Defaults are `iter.max = 1000`,
#'   `eval.max = 1500`, `rel.tol = 1e-10` and `x.tol = 1e-8`. Positive native controls
#'   are retained with the fit. Convergence and curvature gates remain mandatory.
#' @param metrics Recorded scores. Default `c("rmse", "mae", "neg_log_score")`.
#'   Negative log score uses one joint density or probability mass per
#'   assessment block. Pooled reporting sums the block losses and divides by
#'   the number of assessment observations.
#' @param max_refits Positive integer kernel-refit budget per fold, default 200.
#'   Each point checks two native starts. Exceeding the budget fails the fold.
#' @param max_memory_mb Positive finite conservative covariance-work budget in
#'   MiB, default 512. Node aggregation avoids a dense training-observation
#'   covariance for repeated surveys. This plan is not an operating-system cap.
#' @param keep_fits Logical, default `FALSE`. Retain complete training models
#'   when `TRUE`; compact training information-criterion records are always kept.
#' @param family Response family and link. Use `"gaussian"` (identity, default)
#'   for continuous observations, `"poisson"` or `"nbinom2"` (log) for
#'   nonnegative integer counts, `"binomial"` (logit) for numeric zero/one
#'   observations, or `"Gamma"` (log) for positive continuous observations.
#' @return A `multiScaleR_field_cv` result inheriting from `multiScaleR_cv`.
#'   `predictions` retains assessment IDs, means, observed values, squared and
#'   absolute losses, and separately labeled field and noise variances. The two
#'   field variance columns use the linear-predictor scale. Observation noise
#'   variance and observation predictive SD use the response scale. The noise
#'   column averages the conditional observation variance over the field.
#'   `joint_scores` retains each block's negative log density or mass. For
#'   Gaussian responses it includes the quadratic loss and log determinant.
#'   For the Laplace calculation it includes both marginal log likelihoods,
#'   their latent-mode diagnostics, and the analysis likelihood gap relative
#'   to the native fit. `training` retains the field SD, range in coordinate
#'   units, and the observation SD, NB2 size or Gamma shape when estimated.
#'   `selection_records`, `diagnostics`, `warnings` and `errors` retain training
#'   likelihoods, parameter counts, convergence information and failures.
#'   `summary` pools the recorded losses. `settings` and `prediction_contracts`
#'   specify the target and omitted uncertainty. Use [cv_tab()] for matching
#'   field-CV candidates and [select_cv()] for training-only AIC/AICc/BIC choices.
#' @details
#' The native fits use glmmTMB and TMB. Every fit checks two deterministic
#' training-only covariance starts and records package versions, convergence,
#' marginal score and fixed/joint/latent curvature. For Gaussian responses
#' both starts must pass and agree. For the other families, use the highest
#' likelihood among the starts that pass and retain the outcome of each start.
#'
#' For Gaussian responses, the score integrates the field and observation
#' noise by exact Gaussian conditioning. For the other families it subtracts
#' the Laplace log marginal likelihood for the joint analysis and assessment
#' responses from the analysis-only log marginal likelihood. Both evaluations
#' use parameters estimated solely in analysis, including field covariance,
#' fixed effects, observation dispersion and ecological scales. This ratio
#' retains dependence across held-out locations and repeated observations.
#' For discrete responses, a ratio above one beyond numerical roundoff is
#' recorded as a failed score.
#'
#' Bernoulli means and observation SDs integrate over a Gaussian approximation
#' to the field given the analysis responses. Log-link families instead use
#' Laplace ratios of exponentially tilted analysis likelihoods for the first
#' two response moments. Each distinct assessment location uses two additional
#' latent solves with the training covariance factor held fixed. This improves
#' response-moment accuracy in low-dimensional numerical comparisons, particularly
#' for sparse observations with strong fields. The joint block score uses its
#' own integration.
#'
#' Observation SDs include conditional observation noise and response-mean
#' variation. For a nonlinear link, obtain the total response variance by
#' squaring the observation predictive SD. The two field variance columns
#' describe the local Gaussian approximation on the linear-predictor scale.
#' Both routes condition on estimated parameters.
#' The field range is a response-correlation distance, separate from ecological
#' kernel bandwidth. Use the returned SDs for predictive variation rather than
#' interpreting them as native mean-prediction standard errors.
#'
#' A scored fold has an eligible native fit and a valid numerical score for the
#' specified calculation. Assess predictive reliability separately on the
#' held-out locations and retain failed folds in the availability denominator.
#' In simulations with buffered held-out locations, squared prediction errors
#' averaged 1.5 (Gaussian), 2.7 (Poisson), and 1.9 (NB2) times the predicted
#' observation variance, against 1 for a calibrated prediction. Treat the
#' returned predictive SDs as optimistic when parameters are estimated.
#'
#' Binning, normalization and supported learned design constants use training
#'   rows only. Custom cohort-dependent preprocessing must also be learned inside
#'   each training fold. Native formula evaluation cannot detect external data
#'   leaked through user-defined functions. No missing rows are silently dropped.
#'
#' The model contains one shared latent Gaussian field with exponential or
#' Gaussian correlation, constant observation dispersion when applicable, and
#' maximum-likelihood estimation.
#' \code{vignette("v08-cross-validation", package = "multiScaleR")} explains
#' fold design and how to read the scores.
#' @examples
#' # Start from raw data and explicit projected points, keeping all repeated
#' # surveys at a location in the same fold. Supply a compact audit directory.
#' # folds <- spatial_cv_folds(xy, block_size = 10000, k = 5, crs = 32617)
#' # cv <- cross_validate_field(y ~ habitat * group + offset(log(exposure)),
#' #   raw_data, xy, folds, audit_dir = "field-cv-audit", crs = 32617)
#' # count_cv <- cross_validate_field(count ~ habitat + offset(log(exposure)),
#' #   raw_data, xy, folds, audit_dir = "count-field-cv-audit", crs = 32617,
#' #   family = "poisson")
#' # cv_tab(list(count_field = count_cv), metric = "neg_log_score")
#' @seealso [cross_validate()], [cross_validate_glm()], [response_field_geometry()]
#' @md
#' @export
cross_validate_field <- function(formula, data, coords, folds, audit_dir, crs = NULL,
    covariance = c("exponential", "gaussian"), kernel_inputs = NULL, optim_args = list(),
    control = list(iter.max = 1000L, eval.max = 1500L, rel.tol = 1e-10, x.tol = 1e-8),
    metrics = c("rmse", "mae", "neg_log_score"), max_refits = 200L,
    max_memory_mb = 512, keep_fits = FALSE,
    family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma")) {
  covariance <- match.arg(covariance)
  family <- match.arg(family)
  provider <- .msr_field_cv_provider()
  validate_scalar_logical(keep_fits, "keep_fits")
  validate_scalar_numeric(max_refits, "max_refits", positive = TRUE, integerish = TRUE)
  validate_scalar_numeric(max_memory_mb, "max_memory_mb", positive = TRUE)
  if (!is.finite(max_memory_mb) || !is.finite(max_refits) || max_refits > .Machine$integer.max)
    stop("Use finite representable field CV budgets.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.symbol(formula[[2L]]) ||
      !is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)) ||
      any(startsWith(names(data), ".msr_field_")))
    stop("Supply a fixed-effects formula, numeric response column and uniquely named raw data without reserved .msr_field_ columns.", call. = FALSE)
  ids <- .msr_fold_ids(rownames(data), nrow(data)); response <- as.character(formula[[2L]])
  if (!is.numeric(data[[response]]) || !is.null(dim(data[[response]])) || any(!is.finite(data[[response]])))
    stop("Field CV responses must be finite numeric observations.", call. = FALSE)
  observed <- data[[response]]
  if (family %in% c("poisson", "nbinom2") && any(observed < 0 | observed != floor(observed)))
    stop("Poisson and NB2 field CV require nonnegative integer counts.", call. = FALSE)
  if (family == "binomial" && any(!observed %in% c(0, 1)))
    stop("Binomial field CV requires a numeric zero/one response.", call. = FALSE)
  if (family == "Gamma" && any(observed <= 0))
    stop("Gamma field CV requires strictly positive observations.", call. = FALSE)
  xy <- .msr_cv_coordinates(coords, crs, ids)
  .msr_response_field_crs(xy$crs)
  if (nrow(xy$xy) != nrow(data)) stop("Coordinate rows must match raw data exactly.", call. = FALSE)
  xy$xy <- unname(xy$xy)
  folds <- .msr_cv_folds(folds, ids)
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) || anyDuplicated(metrics) ||
      any(!metrics %in% c("rmse", "mae", "neg_log_score")))
    stop("Field CV metrics are rmse, mae and neg_log_score.", call. = FALSE)
  if (!is.list(control) || !length(control) || is.null(names(control)) || anyNA(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% c("iter.max", "eval.max", "rel.tol", "x.tol", "abs.tol")) ||
      any(!vapply(control, function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1))))
    stop("Use named positive native nlminb iteration/evaluation/tolerance controls.", call. = FALSE)
  if (any(vapply(control[intersect(names(control), c("iter.max", "eval.max"))],
      function(x) x != floor(x) || x > .Machine$integer.max, logical(1))))
    stop("Native iteration/evaluation budgets must be representable positive integers.", call. = FALSE)
  if (!is.character(audit_dir) || length(audit_dir) != 1L || is.na(audit_dir) || !nzchar(audit_dir))
    stop("Supply one explicit field CV audit_dir.", call. = FALSE)
  inputs <- kernel_inputs
  spatial <- if (is.null(inputs)) character() else inputs$scale_vars$covariate
  if (!is.null(inputs)) {
    if (!inherits(inputs, "multiScaleR_data") || !is.list(inputs$raw_cov) || !is.list(inputs$d_list) ||
        length(inputs$raw_cov) != nrow(inputs$kernel_dat) || length(inputs$d_list) != nrow(inputs$kernel_dat) ||
        !.msr_valid_observation_ids(rownames(inputs$kernel_dat), nrow(inputs$kernel_dat)) ||
        !setequal(ids, rownames(inputs$kernel_dat)))
      stop("Field scale CV requires matching raw-cell kernel inputs and row IDs.", call. = FALSE)
    for (nm in setdiff(spatial, names(data))) data[[nm]] <- 0
  }
  formula <- stats::formula(stats::terms(formula, data = data))
  has_random <- function(x) {
    if (!is.call(x)) return(FALSE)
    head <- x[[1L]]
    name <- if (is.symbol(head)) as.character(head) else
      if (is.call(head) && length(head) == 3L) as.character(head[[3L]]) else ""
    if (name %in% c("|", "||", "s", "te", "ti", "t2")) return(TRUE)
    parts <- as.list(x)
    any(vapply(seq_along(parts)[-1L], function(j)
      if (identical(unname(parts[j]), unname(alist(x = )))) FALSE else has_random(parts[[j]]), logical(1)))
  }
  if (has_random(formula[[3L]]) || response %in% all.vars(formula[[3L]]) ||
      length(setdiff(all.vars(formula), names(data))))
    stop("Use explicit raw fixed predictors without response leakage, additional random effects or penalized smooths.", call. = FALSE)
  covariates <- intersect(spatial, all.vars(formula[[3L]]))
  allowed <- setdiff(names(formals(multiScale_optim)), c("fitted_mod", "kernel_inputs", "join_by", "n_cores", "PSOCK", "refit_fn", "backend", "verbose"))
  if (!is.list(optim_args) || (length(optim_args) && (is.null(names(optim_args)) || anyNA(names(optim_args)) ||
      anyDuplicated(names(optim_args)) || any(!names(optim_args) %in% allowed))) ||
      (!is.null(optim_args$profile_control) && !is.list(optim_args$profile_control)) ||
      !is.null(optim_args$profile_control$worker_pool) || (is.null(inputs) && length(optim_args)))
    stop("Use supported serial refit optim_args only with kernel_inputs.", call. = FALSE)
  if (!is.null(inputs) && (!length(covariates) || !any(inputs$scale_vars$optimize[match(covariates, spatial)])))
    stop("Field scale CV needs an optimized predictor in its fixed formula.", call. = FALSE)
  transform <- .msr_cv_covariate_transform(optim_args, "refit", !is.null(inputs))
  position <- as.character(glmmTMB::numFactor(xy$xy))
  for (split in folds) {
    ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    if (length(intersect(position[ai], position[ti])))
      stop("Hold out complete locations; analysis and assessment cannot share a field location.", call. = FALSE)
    nodes <- if (family == "gaussian") length(unique(position[ai])) else
      length(unique(position[c(ai, ti)]))
    .msr_field_cv_budget(nodes, length(ti), length(all.vars(formula)) + 3L, max_memory_mb)
  }
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  audit_dir <- normalizePath(audit_dir, winslash = "/", mustWork = TRUE)
  out <- list(predictions = list(), joint_scores = list(), fold_results = list(), training = list(),
    selection_records = list(), diagnostics = list(), errors = list(), warnings = list(),
    prediction_contracts = list(), fits = if (keep_fits) list() else NULL)
  for (i in seq_along(folds)) {
    name <- names(folds)[i]; split <- folds[[i]]
    ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    directory <- file.path(audit_dir, sprintf("fold%03d", i)); dir.create(directory, showWarnings = FALSE)
    training <- data[ai, , drop = FALSE]; assessment <- data[ti, , drop = FALSE]
    train <- assess <- initial <- fit <- callback <- NULL; messages <- character(); stage <- "preparation"
    pred <- data.frame(fold = name, row_id = split$assessment, observed = data[[response]][ti],
      prediction = NA_real_, squared_error = NA_real_, absolute_error = NA_real_, eligible = FALSE)
    for (nm in c("fitted_training_field_variance", "new_location_field_variance",
        "observation_predictive_SE", "observation_noise_variance")) pred[[nm]] <- NA_real_
    error <- tryCatch(withCallingHandlers({
      if (!is.null(inputs)) {
        train <- .msr_cv_inputs(inputs, split$analysis, covariates)
        raw <- .msr_cv_initial_values(train)
        normalization <- .msr_cv_raw_identity(colMeans(raw), apply(raw, 2L, stats::sd), transform, train$scale_vars)
        if (any(!is.finite(normalization$sd) | normalization$sd <= 0)) stop("Training ecological predictors need positive finite spread.", call. = FALSE)
        train$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, normalization$mean), 2L, normalization$sd, "/"))
        train$scl_params <- normalization; training[covariates] <- train$kernel_dat[covariates]
      }
      stage <- "initial_field"
      initial <- .msr_field_cv_fit(formula, training, xy$xy[ai, , drop = FALSE], covariance,
        control, directory, max_memory_mb, provider, family)
      fit <- initial
      if (!is.null(train)) {
        stage <- "scale_fit"
        callback <- .msr_field_cv_refitter(formula, training, xy$xy[ai, , drop = FALSE], covariance,
          control, directory, max_memory_mb, provider, covariates, max_refits, family)
        fit <- do.call(multiScale_optim, c(list(fitted_mod = initial, kernel_inputs = train,
          backend = "refit", n_cores = 1L, PSOCK = FALSE, verbose = FALSE,
          refit_fn = callback$callback), optim_args))
        if (callback$state$exhausted || !identical(as.integer(fit$optim_results$convergence), 0L))
          stop("The training scale optimizer failed or exhausted max_refits.", call. = FALSE)
        training <- .msr_refit_prediction_data(fit, training, train, include_fixed = TRUE)
        assess <- .msr_cv_inputs(inputs, split$assessment, covariates,
          if (is.null(train$binned)) NULL else train$binned$cv_reference)
        assessment <- .msr_refit_prediction_data(fit, assessment, assess, include_fixed = TRUE)
      }
      model <- if (inherits(fit, "multiScaleR")) fit$opt_mod else fit
      stage <- "joint_prediction"
      result <- if (family == "gaussian") {
        .msr_field_cv_predict(model, formula, training, assessment,
          xy$xy[ai, , drop = FALSE], xy$xy[ti, , drop = FALSE], covariance, max_memory_mb)
      } else {
        .msr_field_cv_predict_laplace(model, formula, training, assessment,
          xy$xy[ai, , drop = FALSE], xy$xy[ti, , drop = FALSE], covariance, max_memory_mb, family)
      }
      # Native fit warnings fail inside .msr_field_cv_fit(). Keep any outer
      # scale-search diagnostics in the result without substituting new SEs.
      pred$prediction <- result$prediction
      pred$squared_error <- (pred$observed - pred$prediction)^2
      pred$absolute_error <- abs(pred$observed - pred$prediction); pred$eligible <- TRUE
      for (nm in c("fitted_training_field_variance", "new_location_field_variance", "observation_predictive_SE")) pred[[nm]] <- result[[nm]]
      pred$observation_noise_variance <- result$observation_noise_variance
      out$joint_scores[[name]] <- result$joint
      native_k <- .msr_parameter_count(model); total_k <- .msr_selection_parameter_count(fit, native_k)
      out$selection_records[[name]] <- list(logLik = .msr_model_loglik(model), K = total_k,
        native_K = native_k, ecological_K = total_k - native_k, nobs = length(ai),
        training_ids = split$analysis, training_response_sha256 = .msr_gdr_hash(data[[response]][ai]),
        training_coordinate_sha256 = .msr_gdr_hash(xy$xy[ai, , drop = FALSE]))
      out$training[[name]] <- list(row_ids = split$analysis, parameters = result$parameters,
        coefficients = glmmTMB::fixef(model)$cond, recipe_sha256 = result$recipe_sha256,
        normalization = if (inherits(fit, "multiScaleR")) fit$scl_params else NULL,
        scale_est = if (inherits(fit, "multiScaleR")) fit$scale_est else NULL,
        bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference,
        refits = if (is.null(callback)) 0L else callback$state$count)
      out$prediction_contracts[[name]] <- result$prediction_contract
      out$diagnostics[[name]] <- attr(model, "msr_field_cv_diagnostics")
      if (keep_fits) out$fits[[name]] <- fit
      NULL
    }, warning = function(w) {messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")}),
      error = function(e) list(stage = stage, message = conditionMessage(e)))
    out$predictions[[name]] <- pred
    out$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai), n_assessment = length(ti),
      status = if (is.null(error)) "scored" else "failed")
    out$warnings[name] <- list(messages); out$errors[name] <- list(error)
  }
  out$targets <- data.frame(row_id = ids, observed = data[[response]])
  out$coordinates <- data.frame(row_id = ids, x = xy$xy[, 1L], y = xy$xy[, 2L])
  out$folds <- folds; out$settings <- list(family = family,
    link = switch(family, gaussian = "identity", binomial = "logit", "log"), covariance = covariance,
    metrics = metrics, provider = provider, crs = xy$crs$wkt, backend = "native_field_refit",
    control = control, optim_args = optim_args, max_refits = max_refits,
    prediction_target = if (family == "gaussian") "held_out_Gaussian_response_joint_block" else
      "held_out_Laplace_response_joint_block", parameter_uncertainty = "not_integrated",
    ecological_scale_uncertainty = "not_integrated", observation_noise_included = TRUE,
    calibration = "not_established", workers = 1L)
  out$call <- match.call()
  .msr_field_cv_finish(out, ids)
}

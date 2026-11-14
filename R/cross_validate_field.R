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
#' Cross-validate a model that adds a spatial random field to its fixed
#' predictors. A *spatial random field* is a smooth surface of spatially
#' correlated variation that the predictors do not explain. In each fold the
#' function estimates the field, its correlation range, the regression
#' coefficients, and optionally the scales of effect from the analysis rows.
#' It then scores the joint distribution of each held-out block, so that
#' correlation among nearby held-out observations counts. Gaussian responses
#' use exact Gaussian conditioning at the fitted parameters. Count, binary,
#' and positive responses use Laplace integration.
#'
#' The count, binary, and positive response families are experimental. Their
#' arguments and results may change in a later release without a deprecation
#' period.
#'
#' @param formula Two-sided fixed-effects formula whose response is a numeric
#'   column of `data`. Products, interactions, factors, offsets, and the
#'   data-dependent terms `scale()`, `poly()`, `splines::ns()`, and
#'   `splines::bs()` keep their usual meanings, with their constants learned
#'   from the analysis rows. Do not add random effects or penalized smooths,
#'   because the function adds the one shared field itself.
#' @param data Data frame of raw responses and predictors, with unique row
#'   names that match the fold IDs, in the same row order as `coords`. Supply
#'   raw predictors, not transformations learned from all rows. Every row has
#'   a case weight of one.
#' @param coords Locations of the observations, one per row of `data`, as
#'   projected `sf` points, a `terra::SpatVector` of points, or a numeric
#'   two-column table. Rows with identical coordinates share one field value,
#'   so repeated surveys of a location are linked.
#' @param folds Named list of explicit splits, usually from
#'   [spatial_cv_folds()]. A location cannot appear in both the analysis and
#'   the assessment rows of one split, so keep repeated surveys of a location
#'   together. Rows excluded from a split, for example by a buffer, enter
#'   neither its fit nor its score.
#' @param audit_dir Directory where compact, write-once records of each
#'   native fit and its starting values are saved, one subfolder per fold.
#'   Records of failed fits are kept too. No models, tapes, or dense
#'   covariance matrices are saved. Give each candidate model its own
#'   directory to keep the records apart.
#' @param crs Projected coordinate reference system, in any form accepted by
#'   `sf::st_crs()`. It is required when `coords` is a coordinate table. The
#'   default `NULL` uses the CRS of spatial `coords`, which must be projected,
#'   and a supplied `crs` must match it. Coordinates are never transformed,
#'   so the correlation range is reported in the CRS units.
#' @param covariance Correlation function of the field, either
#'   `"exponential"` (default), with correlation exp(-d / range) at distance
#'   d, or `"gaussian"`, with correlation exp(-(d / range)^2). The range and
#'   the field SD are estimated in every fold, as is the observation
#'   dispersion for Gaussian, NB2, and Gamma models.
#' @param kernel_inputs Optional `multiScaleR_data` object from
#'   [kernel_prep()] with the cell-level data kept and row names matching
#'   `data`. The default `NULL` fits a field model without estimating scales.
#'   When it is supplied, each fold also estimates the scales of the
#'   optimized covariates that appear in `formula` (list them with
#'   [msr_vars()]) by refitting the field model.
#' @param optim_args Named list of extra [multiScale_optim()] arguments for the
#'   scale fit, used only with `kernel_inputs`. The default is an empty list.
#'   This function fixes the backend (`"refit"`), the data, the refit
#'   function, and the single worker, so those cannot be replaced.
#' @param control Named list of native `nlminb()` optimizer settings for each
#'   field fit. Accepted names are `iter.max`, `eval.max`, `rel.tol`, `x.tol`,
#'   and `abs.tol`, each a positive number. The defaults are
#'   `iter.max = 1000`, `eval.max = 1500`, `rel.tol = 1e-10`, and
#'   `x.tol = 1e-8`. The settings are stored with each fit, and the
#'   convergence and curvature checks apply whatever they are.
#' @param metrics Scores to record, any of `"rmse"`, `"mae"`, and
#'   `"neg_log_score"` (all three by default). The negative log score uses one
#'   joint density, or one joint probability for counts and 0/1 data, per
#'   assessment block. The pooled score sums the block losses and divides by
#'   the number of assessment observations.
#' @param max_refits Positive whole number, the most native refits allowed in
#'   one fold during the scale search (default 200). Each refit checks two
#'   starting values. A fold that exceeds the budget fails.
#' @param max_memory_mb Positive number, a conservative budget in MiB for the
#'   covariance calculations of one fold (default 512), checked before any
#'   fitting. Repeated surveys at one location share one field value, so they
#'   do not enlarge the training covariance matrix. This is a planning check,
#'   not an operating-system memory limit.
#' @param keep_fits Logical. `TRUE` keeps the complete training model of each
#'   fold in `fits`. The default `FALSE` saves memory. The compact training
#'   records that [select_cv()] reads are kept either way.
#' @param family Response family and link. Use `"gaussian"` (identity, the
#'   default) for continuous observations, `"poisson"` or `"nbinom2"` (log)
#'   for counts, `"binomial"` (logit) for numeric 0/1 observations, or
#'   `"Gamma"` (log) for strictly positive continuous observations. All
#'   families other than `"gaussian"` are experimental.
#' @return A list of class `multiScaleR_field_cv`, which inherits from
#'   `multiScaleR_cv`. `summary` holds the pooled scores, `predictions` the
#'   predicted mean, errors, and variance components of each held-out
#'   observation, `joint_scores` the joint loss of each held-out block, and
#'   `training` the field SD, correlation range (in coordinate units),
#'   dispersion, and coefficients of each fold. Compact `selection_records`
#'   are always kept. The spatial field section of [multiScaleR_cv-class]
#'   describes every element and its scale.
#' @details
#' # What each fold estimates
#'
#' The native fits use glmmTMB and TMB with maximum likelihood. The model has
#' one shared latent Gaussian field with exponential or Gaussian correlation,
#' and a constant observation dispersion where the family has one. From the
#' analysis rows each fold estimates the regression coefficients, the field
#' SD, the correlation range, the observation dispersion, and, with
#' `kernel_inputs`, the scales of effect. Distance bins, standardization, and
#' the constants of learned formula terms also come from the analysis rows
#' only. Custom preprocessing that depends on the data must be learned inside
#' each training fold too, because formula evaluation cannot detect outside
#' information passed in through your own functions. No missing rows are
#' dropped silently.
#'
#' The correlation range is the distance, in coordinate units, over which the
#' unexplained variation stays correlated. It is separate from the kernel
#' scale of effect, which describes how far the landscape around a point
#' affects it.
#'
#' Every fit tries two deterministic starting values chosen from the training
#' data and records the package versions, convergence, marginal score, and
#' curvature checks of each. For Gaussian responses both starts must pass and
#' agree, or the fold fails. For the other families the passing start with
#' the highest likelihood is used, and the outcome of each start is kept.
#' With few locations in a training fold, the field's standard deviation and
#' range are often weakly identified, and the two Gaussian starts can stop at
#' different points of a nearly flat likelihood. Those folds fail by design
#' ("starts disagree") rather than report one arbitrary solution. Use more
#' locations per fold, or fewer folds, and inspect `audit_dir` to see how far
#' apart the starts ended.
#'
#' # Build folds for a field model
#'
#' Hold out complete locations, so that no location has observations in both
#' parts of a split. Spatial blocks with a buffer test prediction at new
#' places. The field locations in `coords` can differ from the focal points
#' used for the landscape covariates, but each location must be given by its
#' actual, consistent coordinates, and `coords` is required even when
#' `kernel_inputs` is supplied. A template that compares two correlation
#' functions on the same folds:
#'
#' ```
#' folds <- spatial_cv_folds(field_xy, block_size = 10000, k = 5,
#'   buffer = 1000, crs = 32617, ids = rownames(raw_data))
#' covariances <- c(exponential = "exponential", gaussian = "gaussian")
#' candidates <- lapply(covariances, function(covariance) {
#'   cross_validate_field(y ~ habitat * cover + offset(log(exposure)),
#'     raw_data, field_xy, folds, crs = 32617, covariance = covariance,
#'     kernel_inputs = inputs,
#'     audit_dir = file.path("field-cv-audits", covariance))
#' })
#' cv_tab(candidates, metric = "neg_log_score")
#' selected <- select_cv(candidates, criterion = "AICc")
#' selected$selection$chosen   # chosen by the training criterion in each fold
#' ```
#'
#' # How the joint score works
#'
#' The negative log score treats the held-out observations of each fold as
#' one block. For Gaussian responses it is the negative log of one
#' multivariate normal density. Its covariance holds the field variation that
#' remains after conditioning on the analysis responses, plus the observation
#' noise. For the other families the score is the Laplace log marginal
#' likelihood of the analysis responses alone minus that of the analysis and
#' assessment responses together. Both calculations use parameters estimated
#' from the analysis rows only, including the field covariance, coefficients,
#' dispersion, and scales, so the score keeps the dependence among held-out
#' locations and repeated observations. For counts and 0/1 data, a block
#' probability above one beyond numerical roundoff is recorded as a failed
#' score.
#'
#' The pooled score is the total block loss divided by the number of
#' assessment observations, and smaller is better. RMSE and MAE summarize the
#' errors of the predicted response means. Multiplying separate marginal
#' densities would ignore the field dependence and define a different score,
#' so [cv_tab()] keeps these scores apart from the per-observation scores of
#' ordinary GLMs.
#'
#' # Predicted means and spreads
#'
#' For 0/1 responses, predicted means and observation SDs integrate over a
#' Gaussian approximation to the field given the analysis responses. For
#' log-link families, the first two response moments come from Laplace
#' ratios of exponentially tilted analysis likelihoods, with two additional
#' latent solves per distinct assessment location and the training
#' covariance factor held fixed. This improves the accuracy of the moments in
#' low-dimensional numerical comparisons, particularly for sparse
#' observations with strong fields. The joint block score uses its own
#' integration.
#'
#' `observation_predictive_SE` combines the conditional observation noise
#' with the variation of the response mean. For a nonlinear link, square it
#' to get the total response variance. The two field variance columns
#' describe the field on the linear-predictor scale. All of these condition
#' on the estimated parameters. Use them as predictive spreads for new
#' observations, not as standard errors of the predicted mean.
#'
#' # Check calibration separately
#'
#' A scored fold has an eligible native fit and a valid score. It does not
#' show that the predictive spreads are reliable. Assess that on the held-out
#' locations, and count failed folds in the denominator when you report how
#' often the method succeeds. In simulations with buffered held-out
#' locations, squared prediction errors averaged 1.5 (Gaussian), 2.7
#' (Poisson), and 1.9 (NB2) times the predicted observation variance, against
#' 1 for a calibrated prediction. Treat the returned spreads as optimistic
#' when parameters are estimated.
#'
#' # Compare and select field models
#'
#' Use [cv_tab()] to rank field results fitted to the same responses, folds,
#' and coordinates. Poisson and NB2 results share a count target. Use
#' [select_cv()] to score a training AIC, AICc, or BIC rule from the compact
#' records. The coefficients, field SD, estimated range, and any observation
#' dispersion count toward K, while the integrated field values add nothing.
#' Field scores do not integrate parameter or scale uncertainty and do not
#' establish calibration. `vignette("v08-cross-validation", package =
#' "multiScaleR")` explains fold design and how to read the scores, and
#' `vignette("v11-spatial-fields", package = "multiScaleR")` covers spatial
#' field models.
#' @examples
#' \donttest{
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   data("landscape_counts")
#'   data("surv_pts")
#'   xy <- crds(vect(surv_pts))   # simulated coordinates in meters
#'   folds <- spatial_cv_folds(xy, block_size = 1000, k = 3, crs = 32617)
#'
#'   # Poisson counts with a site covariate and an exponential field. The
#'   # Poisson family is experimental.
#'   field_cv <- cross_validate_field(counts ~ site, landscape_counts, xy,
#'     folds, audit_dir = file.path(tempdir(), "field-cv"), crs = 32617,
#'     family = "poisson")
#'   print(field_cv)   # the log score is a joint score per held-out survey
#'
#'   # Field SD (log-link scale) and correlation range (m) in each fold
#'   print(sapply(field_cv$training, function(fold)
#'     unlist(fold$parameters[c("field_SD", "range")])))
#' }
#' }
#' @seealso [cross_validate()], [cross_validate_glm()],
#'   [response_field_geometry()], [multiScaleR_cv-class], [cv_tab()], and
#'   [select_cv()].
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

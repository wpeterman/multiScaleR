# Validate the whole split plan before starting any fold fits. IDs, rather than
# changing row positions, bind raw data, spatial buffers and assessment output.
.msr_cv_folds <- function(folds, ids) {
  valid_ids <- function(x) is.character(x) && length(x) > 0L &&
    !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x) && all(x %in% ids)
  if (!is.list(folds) || !length(folds) || is.null(names(folds)) ||
      anyNA(names(folds)) || any(!nzchar(names(folds))) || anyDuplicated(names(folds)))
    stop("folds must be a nonempty, uniquely named list of explicit splits.", call. = FALSE)
  for (name in names(folds)) {
    split <- folds[[name]]
    if (!is.list(split) || !setequal(names(split), c("analysis", "assessment")) ||
        length(split) != 2L || !valid_ids(split$analysis) || !valid_ids(split$assessment))
      stop("Fold `", name, "` must contain analysis and assessment character row IDs, without missing, duplicate or unknown IDs.", call. = FALSE)
    if (length(intersect(split$analysis, split$assessment)))
      stop("Fold `", name, "` has overlapping analysis and assessment rows.", call. = FALSE)
  }
  folds
}

# Binning is itself learned preparation. Its physical domain is declared by
# max_D; representative distances use training buffers only and stay frozen.
.msr_cv_bins <- function(inputs, reference = NULL) {
  specs <- inputs$scale_vars[inputs$scale_vars$type == "kernel", , drop = FALSE]
  if (is.null(inputs$binned) || !nrow(specs)) return(NULL)
  n <- length(inputs$d_list)
  nb <- inputs$binned$nbins
  if (!is.numeric(nb) || length(nb) != 1L || !is.finite(nb) || nb < 2 || nb != floor(nb))
    stop("CV binning requires a valid prepared nbins value.", call. = FALSE)
  domain <- inputs$max_D / inputs$unit_conv
  if (is.null(reference)) {
    breaks <- seq(0, domain * (1 + 1e-9), length.out = nb + 1L)
  } else {
    breaks <- reference$breaks
    if (length(breaks) != nb + 1L || any(!is.finite(breaks)) || any(diff(breaks) <= 0) ||
        length(reference$bin_rep) != nb || any(!is.finite(reference$bin_rep)) ||
        !identical(reference$domain, domain))
      stop("Invalid frozen training bin reference.", call. = FALSE)
  }
  codes <- lapply(inputs$d_list, function(d) {
    if (any(!is.finite(d) | d < 0 | d > domain * (1 + 1e-9)))
      stop("Prepared cell distances must lie inside the declared physical buffer.", call. = FALSE)
    .bincode(d, breaks, include.lowest = TRUE)
  })
  if (is.null(reference)) {
    total <- count <- numeric(nb)
    for (i in seq_len(n)) {
      total <- total + .msr_bin_accumulate(inputs$d_list[[i]], codes[[i]], nb)
      count <- count + tabulate(codes[[i]], nb)
    }
    rep <- total / count
    empty <- !is.finite(rep)
    rep[empty] <- ((breaks[-1L] + breaks[-length(breaks)]) / 2)[empty]
    reference <- list(breaks = breaks, bin_rep = rep, domain = domain)
  }
  vsum <- csum <- stats::setNames(vector("list", nrow(specs)), specs$covariate)
  for (j in seq_len(nrow(specs))) {
    V <- C <- matrix(0, n, nb)
    for (i in seq_len(n)) {
      value <- as.numeric(inputs$raw_cov[[i]][, specs$source[j]])
      if (length(value) != length(codes[[i]]))
        stop("Prepared distances and source-cell values must align.", call. = FALSE)
      keep <- is.finite(value)
      V[i, ] <- .msr_bin_accumulate(value[keep], codes[[i]][keep], nb)
      C[i, ] <- tabulate(codes[[i]][keep], nb)
    }
    vsum[[j]] <- V; csum[[j]] <- C
  }
  list(bin_rep = reference$bin_rep, vsum = vsum, csum = csum,
    covariates = specs$covariate, sources = stats::setNames(specs$source, specs$covariate),
    nbins = nb, n_points = n, point_ids = rownames(inputs$kernel_dat),
    cv_reference = reference)
}

.msr_cv_inputs <- function(inputs, ids, covariates, reference = NULL) {
  rows <- match(ids, rownames(inputs$kernel_dat))
  out <- inputs
  out$kernel_dat <- inputs$kernel_dat[rows, covariates, drop = FALSE]
  out$d_list <- inputs$d_list[rows]
  out$raw_cov <- inputs$raw_cov[rows]
  names(out$d_list) <- names(out$raw_cov) <- ids
  specs <- inputs$scale_vars
  old_opt <- specs$covariate[specs$optimize]
  out$scale_vars <- specs[match(covariates, specs$covariate), , drop = FALSE]
  selected <- match(out$scale_vars$covariate[out$scale_vars$optimize], old_opt)
  out$sigma <- inputs$sigma[selected]
  out$shape <- if (is.null(inputs$shape)) NULL else inputs$shape[selected]
  out$n_covs <- length(selected)
  out$binned <- .msr_cv_bins(out, reference)
  out
}

.msr_cv_initial_values <- function(inputs) {
  specs <- inputs$scale_vars
  value <- matrix(vapply(seq_along(inputs$d_list), function(i)
    .msr_eval_scale_vars(inputs$d_list[[i]], inputs$raw_cov[[i]], specs,
      inputs$sigma, inputs$shape, inputs$kernel, inputs$unit_conv,
      inputs$resolution, inputs$n_cols), numeric(nrow(specs))),
    nrow = length(inputs$d_list), ncol = nrow(specs), byrow = TRUE)
  dimnames(value) <- list(rownames(inputs$kernel_dat), specs$covariate)
  if (any(!is.finite(value))) stop("Initial spatial covariates must be finite within each fold.", call. = FALSE)
  value
}

# Restrict outcomes to an explicit recorded observation definition. Predictor
# equations remain native GLM formulas, including learned poly/spline metadata.
.msr_cv_outcome <- function(formula, data, family) {
  expr <- formula[[2L]]
  while (is.call(expr) && identical(expr[[1L]], as.name("I"))) expr <- expr[[2L]]
  simple <- is.symbol(expr) || (is.call(expr) && identical(expr[[1L]], as.name("cbind")) &&
    length(expr) == 3L && all(vapply(as.list(expr)[-1L], is.symbol, logical(1))))
  if (!simple) stop("CV responses must be a data column or cbind(successes, failures); predefine other observation transformations explicitly.", call. = FALSE)
  y <- eval(expr, envir = data, enclos = environment(formula))
  if (is.matrix(y) && ncol(y) == 1L) y <- y[, 1L]
  trials <- rep(1, nrow(data))
  if (family == "binomial") {
    if (is.factor(y)) {
      if (nlevels(y) != 2L) stop("A binomial response factor must declare two levels.", call. = FALSE)
      y <- as.integer(y != levels(y)[1L])
    }
    if (is.matrix(y) && ncol(y) == 2L) {
      if ((!is.numeric(y) && !is.logical(y)) || any(!is.finite(y) | y < 0 | y != floor(y)))
        stop("Binomial successes and failures must be finite nonnegative integer counts.", call. = FALSE)
      trials <- rowSums(y)
      y <- y[, 1L]
    } else if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) || anyNA(y) || any(!y %in% c(0, 1))) {
      stop("Use binary binomial outcomes or cbind(successes, failures), not weighted proportions.", call. = FALSE)
    }
  } else {
    if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) || any(!is.finite(y)))
      stop("Gaussian and Poisson outcomes must be finite numeric vectors.", call. = FALSE)
    if (family == "poisson" && any(y < 0 | y != floor(y)))
      stop("Poisson outcomes must be nonnegative integer counts.", call. = FALSE)
  }
  if (length(y) != nrow(data)) stop("Response rows must match the explicit data rows.", call. = FALSE)
  list(observed = as.numeric(y), trials = as.numeric(trials))
}

.msr_cv_score <- function(predictions, family, metrics) {
  y <- predictions$observed; mu <- predictions$prediction
  n <- predictions$trials; w <- predictions$weight
  if (any(!is.finite(mu)) || (family == "binomial" && any(mu < 0 | mu > 1)) ||
      (family == "poisson" && any(mu < 0)))
    stop("Held-out predictions must be finite and inside the family's response domain.", call. = FALSE)
  eligible <- w > 0 & (family != "binomial" | n > 0)
  if (!any(eligible)) stop("The assessment set has no positive-weight observations or trials.", call. = FALSE)
  predictions$eligible <- eligible
  for (metric in metrics) {
    value <- switch(metric, rmse = (y - mu)^2, mae = abs(y - mu),
      neg_log_score = if (family == "poisson") -stats::dpois(y, mu, log = TRUE) else
        -stats::dbinom(y, n, mu, log = TRUE),
      deviance = stats::poisson()$dev.resids(y, mu, rep(1, length(y))),
      brier = (y * (1 - mu)^2 + (n - y) * mu^2) / pmax(n, 1))
    value[!eligible] <- NA_real_
    predictions[[metric]] <- value
  }
  predictions
}

.msr_cv_metric <- function(predictions, metric, family, complete, expected) {
  keep <- predictions$eligible & !is.na(predictions[[metric]])
  w <- predictions$weight[keep]
  if (metric == "brier") w <- w * predictions$trials[keep]
  value <- if (complete && length(w) && sum(w) > 0)
    sum(w * predictions[[metric]][keep]) / sum(w) else NA_real_
  if (metric == "rmse") value <- sqrt(value)
  data.frame(metric = metric, estimate = value,
    unit = if (metric == "brier") "trial" else "observation",
    n_assessment = expected, n_scored = sum(keep), weight_sum = sum(w), complete = complete,
    stringsAsFactors = FALSE)
}

#' Cross-Validate Estimated Scales with Explicit Folds
#'
#' Refit ecological scales and GLM coefficients using only each training set,
#' then score predictions for the explicitly held-out rows. A fold is one
#' training and assessment split. This first stage supports ordinary Gaussian,
#' Poisson and binomial GLMs and runs folds sequentially.
#'
#' @param formula Two-sided GLM formula. The response is a data column or
#'   \code{cbind(successes, failures)}. An \code{I()} wrapper is also accepted.
#'   Predictor expressions, interactions, factors, offsets and supported learned
#'   transformations use the existing backend's native formula machinery.
#'   Supply every variable in \code{data} or the prepared spatial inputs.
#' @param data Data frame containing raw responses and nonspatial formula
#'   variables. Unique row names must match the prepared spatial row IDs.
#'   Prepared spatial columns are rebuilt inside each fold, replacing any
#'   supplied values. No fitting rows are silently omitted for missing values.
#' @param kernel_inputs A \code{multiScaleR_data} object from
#'   \code{\link{kernel_prep}} with cell-level distances and source values
#'   retained. Use \code{store_cell_data = TRUE}. Its buffer, initial scales,
#'   metric definitions and other controls must be chosen independently of the
#'   assessment outcomes. Lean-only and grouped preparation are not supported.
#' @param folds Uniquely named list of explicit splits. Each split contains
#'   nonempty character vectors \code{analysis} and \code{assessment} of row
#'   IDs. They must be disjoint within a split. Splits can omit rows for spatial
#'   buffers and can repeat assessment rows across splits. Keep dependent sites,
#'   animals or other sampling groups together when defining the splits.
#' @param family Character response family: \code{"gaussian"} uses the identity
#'   link, \code{"poisson"} the log link, and \code{"binomial"} the logit link.
#'   The default is \code{"gaussian"}. Other families, links and quasi-likelihood
#'   targets require separately qualified scoring adapters.
#' @param backend Scale-fitting backend, \code{"profile_glm"} (default) or
#'   \code{"refit"}, passed to \code{\link{multiScale_optim}}. Its existing
#'   starts, convergence, curvature and kernel/metric guards remain in force.
#' @param weights Optional nonnegative integer frequency weights, one per data
#'   row in its supplied order. Names, when present, must match that order.
#'   Gaussian weights must all equal one. Binomial trials belong in the
#'   two-column response, separately from frequency weights; weighted proportions
#'   are not accepted. A backend may impose stricter native weight requirements.
#' @param optim_args Named list of additional \code{multiScale_optim()} controls,
#'   such as \code{par} or \code{profile_control}. Defaults preserve that backend's
#'   controls. Starts must be selected without assessment outcomes. Object/data
#'   arguments, custom refit functions, worker pools and parallel controls cannot
#'   be replaced. Each fold uses one fit worker.
#' @param glm_control Named list from \code{stats::glm.control()}, used for each
#'   initial training-only GLM. The default keeps native GLM controls.
#' @param metrics Optional character vector of score names. Defaults are
#'   \code{c("rmse", "mae")} for Gaussian, \code{c("rmse", "mae",
#'   "neg_log_score", "deviance")} for Poisson, and
#'   \code{c("neg_log_score", "brier")} for binomial. Lower values are better.
#' @param keep_fits Logical. Default \code{FALSE} keeps compact predictions,
#'   scales and diagnostics. Set \code{TRUE} to retain full fold models in memory
#'   for inspection; this can use substantially more memory.
#' @return An object of class \code{multiScaleR_cv}, a list containing:
#'   \describe{
#'     \item{\code{summary}}{Pooled scores with columns \code{metric},
#'       \code{estimate}, scoring \code{unit}, planned \code{n_assessment},
#'       \code{n_scored}, \code{weight_sum} and \code{complete}. Estimates are
#'       \code{NA} if any planned fold fails.}
#'     \item{\code{metrics}}{The same score columns per \code{fold}. RMSE is the
#'       square root of pooled squared errors, not the average of fold RMSEs.}
#'     \item{\code{predictions}}{Every planned assessment occurrence, including
#'       failed folds: \code{fold}, \code{row_id}, \code{observed} response/count,
#'       \code{trials}, response-scale \code{prediction}, frequency \code{weight},
#'       \code{eligible} scoring flag and requested score columns. The \code{rmse}
#'       column stores squared errors. Failed predictions/scores are \code{NA}.}
#'     \item{\code{fold_results}}{Per-fold \code{fold}, \code{n_analysis},
#'       \code{n_assessment} and \code{status} (\code{"scored"} or \code{"failed"}).}
#'     \item{\code{folds}, \code{assessment_counts}}{The explicit split plan and
#'       named counts of assessment occurrences for every input row, including
#'       unassessed rows. Repeated assessment occurrences each enter pooled scores.}
#'     \item{\code{training}}{Per-fold initial/fitted normalization,
#'       coefficients on the link scale, scales in map units, unitless shapes,
#'       frozen bin reference and optimizer/start
#'       records when available. Failed preparation can leave an element empty.}
#'     \item{\code{diagnostics}, \code{warnings}, \code{errors}}{Original scale-fit
#'       diagnostics, initial/final GLM convergence metadata and captured messages
#'       per fold. Failures are retained and never retried automatically.}
#'     \item{\code{prediction_contracts}}{Existing prediction contracts per fold,
#'       retaining the native conditioning and uncertainty distinctions.}
#'     \item{\code{fits}}{Full fold fits when \code{keep_fits = TRUE}; otherwise
#'       \code{NULL}. These fits use only their own training data.}
#'     \item{\code{complete}, \code{settings}, \code{call}}{Whether all planned
#'       folds were scored, family/backend/metric settings and the original call.
#'       Completion does not certify diagnostic eligibility or calibration.}
#'   }
#' @details
#' Raw spatial buffers are reused without changing them. Initial spatial
#' normalization and pooled representative bin distances are learned from the
#' training rows only. The bin domain remains the predeclared physical buffer.
#' Assessment values use the frozen training bins and normalization. Existing
#' coefficient, scale, shape and numerical-diagnostic calculations are unchanged.
#'
#' Negative log scores retain Poisson/binomial likelihood constants and average
#' over positive-weight observations using frequency weights. Binomial Brier
#' scores average squared prediction errors across trials, weighted by frequency
#' times trial count. Zero-trial rows contribute no score. Infinite log scores
#' correctly record an observed outcome assigned zero probability. Point scores
#' use fitted parameters; this function adds no parameter uncertainty, observation
#' intervals, calibration claim or new eligibility for guarded methods.
#'
#' Use identical folds and assessment targets when comparing candidate models.
#' If CV scores also choose a kernel, metric or model, assess the selected
#' procedure with outer validation or an independent test set. These scores
#' complement \code{\link{aic_tab}}; they are not AIC values or Akaike weights.
#'
#' Read \code{complete} first. Inspect failed folds and original diagnostics
#' before ranking models. Partial fold scores remain available for diagnosis,
#' while an incomplete pooled summary stays \code{NA}. This stage does not
#' generate random/spatial folds or score GAMs, mixed, unmarked or Bayesian fits.
#' @examples
#' # Explicit folds can be checked before any fitting is attempted.
#' ids <- paste0("site", 1:6)
#' folds <- list(first = list(analysis = ids[1:4], assessment = ids[5:6]),
#'               second = list(analysis = ids[3:6], assessment = ids[1:2]))
#' intersect(folds$first$analysis, folds$first$assessment) # no overlap
#' @seealso \code{\link{kernel_prep}}, \code{\link{multiScale_optim}},
#'   \code{\link{prepare_prediction}}, \code{\link{aic_tab}}
#' @export
cross_validate <- function(formula, data, kernel_inputs, folds,
    family = c("gaussian", "poisson", "binomial"),
    backend = c("profile_glm", "refit"), weights = NULL,
    optim_args = list(), glm_control = stats::glm.control(), metrics = NULL,
    keep_fits = FALSE) {
  family <- match.arg(family); backend <- match.arg(backend)
  validate_scalar_logical(keep_fits, "keep_fits")
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be a two-sided GLM formula.", call. = FALSE)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("data must be a nonempty data frame with unique columns.", call. = FALSE)
  ids <- rownames(data)
  if (!.msr_valid_observation_ids(ids, nrow(data)))
    stop("data requires unique, nonmissing row IDs.", call. = FALSE)
  inputs <- kernel_inputs
  if (!inherits(inputs, "multiScaleR_data") || is.null(inputs$scale_vars) ||
      !is.list(inputs$d_list) || !is.list(inputs$raw_cov) ||
      length(inputs$d_list) != nrow(inputs$kernel_dat) ||
      length(inputs$raw_cov) != nrow(inputs$kernel_dat) ||
      !.msr_valid_observation_ids(rownames(inputs$kernel_dat), nrow(inputs$kernel_dat)) ||
      !setequal(ids, rownames(inputs$kernel_dat)))
    stop("CV requires matching kernel_prep() row IDs and stored cell-level inputs; use store_cell_data = TRUE.", call. = FALSE)
  variables <- all.vars(formula)
  spatial <- inputs$scale_vars$covariate
  covariates <- if ("." %in% variables) spatial else intersect(spatial, variables)
  if (!any(inputs$scale_vars$optimize[match(covariates, spatial)]))
    stop("The formula must use an optimized spatial predictor.", call. = FALSE)
  if (length(setdiff(variables, c(names(data), spatial, "."))))
    stop("Supply every formula variable explicitly in data or spatial inputs.", call. = FALSE)
  if (length(intersect(all.vars(formula[[2L]]), spatial)))
    stop("A spatial predictor cannot be the CV response.", call. = FALSE)
  folds <- .msr_cv_folds(folds, ids)
  defaults <- list(gaussian = c("rmse", "mae"),
    poisson = c("rmse", "mae", "neg_log_score", "deviance"),
    binomial = c("neg_log_score", "brier"))
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
  allowed <- setdiff(names(formals(multiScale_optim)), c("fitted_mod", "kernel_inputs",
    "join_by", "n_cores", "PSOCK", "refit_fn", "backend", "verbose"))
  if (!is.list(optim_args) || (length(optim_args) &&
      (is.null(names(optim_args)) || anyNA(names(optim_args)) || anyDuplicated(names(optim_args)) ||
       any(!names(optim_args) %in% allowed))) ||
      (!is.null(optim_args$profile_control) && !is.list(optim_args$profile_control)) ||
      (is.list(optim_args$profile_control) && !is.null(optim_args$profile_control$worker_pool)))
    stop("optim_args must name supported serial scale controls; object/data/worker overrides are not allowed.", call. = FALSE)
  if (!is.list(glm_control) || is.null(names(glm_control)) || anyDuplicated(names(glm_control)) ||
      any(!names(glm_control) %in% names(formals(stats::glm.control))))
    stop("glm_control must be a named glm.control() list.", call. = FALSE)
  # Validate controls before any fit, without changing native defaults.
  glm_control <- do.call(stats::glm.control, glm_control)
  family_object <- switch(family, gaussian = stats::gaussian(), poisson = stats::poisson(), binomial = stats::binomial())
  output <- list(predictions = list(), metrics = list(), fold_results = list(),
    training = list(), diagnostics = list(), warnings = list(), errors = list(),
    prediction_contracts = list(), fits = if (keep_fits) list() else NULL)
  for (name in names(folds)) {
    split <- folds[[name]]
    ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    pred <- data.frame(fold = name, row_id = split$assessment, observed = NA_real_, trials = NA_real_,
      prediction = NA_real_, weight = unname(weights[ti]), eligible = FALSE, stringsAsFactors = FALSE)
    for (metric in metrics) pred[[metric]] <- NA_real_
    fit <- initial <- NULL; messages <- character(); stage <- "preparation"
    error <- tryCatch(withCallingHandlers({
      outcome <- .msr_cv_outcome(formula, data[ti, , drop = FALSE], family)
      pred$observed <- outcome$observed; pred$trials <- outcome$trials
      .msr_cv_outcome(formula, data[ai, , drop = FALSE], family)
      train <- .msr_cv_inputs(inputs, split$analysis, covariates)
      raw <- .msr_cv_initial_values(train)
      center <- colMeans(raw); spread <- apply(raw, 2L, stats::sd)
      if (any(!is.finite(spread) | spread <= 0))
        stop("Training spatial predictors must have positive finite standard deviations.", call. = FALSE)
      train$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"), optional = TRUE)
      train$scl_params <- list(mean = center, sd = spread)
      training_data <- data[ai, , drop = FALSE]
      for (nm in covariates) training_data[[nm]] <- train$kernel_dat[[nm]]
      stage <- "initial_glm"
      initial <- do.call(stats::glm, list(formula = formula, data = training_data,
        family = family_object, weights = unname(weights[ai]), control = glm_control,
        na.action = stats::na.fail, model = TRUE))
      output$training[[name]] <- list(initial_normalization = train$scl_params,
        initial_coefficients = stats::coef(initial),
        bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference)
      if (!isTRUE(initial$converged)) stop("Initial training GLM did not converge.", call. = FALSE)
      if (is.factor(stats::model.response(initial$model)) &&
          length(unique(stats::model.response(initial$model))) != 2L)
        stop("A training binomial factor must contain both declared response levels.", call. = FALSE)
      stage <- "scale_fit"
      fit <- do.call(multiScale_optim, c(list(fitted_mod = initial, kernel_inputs = train,
        n_cores = 1L, backend = backend, verbose = FALSE), optim_args))
      output$training[[name]] <- list(initial_normalization = train$scl_params,
        initial_coefficients = stats::coef(initial),
        fitted_normalization = fit$scl_params, coefficients = stats::coef(fit$opt_mod),
        scale_est = fit$scale_est, shape_est = fit$shape_est,
        bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference,
        optimizer = fit$optim_results[c("par", "value", "convergence", "message")],
        start_results = fit$profile_fit$start_results)
      if (!isTRUE(fit$opt_mod$converged) || !identical(as.integer(fit$optim_results$convergence), 0L))
        stop("Training GLM or scale optimizer did not converge.", call. = FALSE)
      stage <- "prediction"
      assess <- .msr_cv_inputs(inputs, split$assessment, covariates,
        if (is.null(train$binned)) NULL else train$binned$cv_reference)
      raw_assess <- .msr_cv_initial_values(assess)
      assess$kernel_dat <- as.data.frame(sweep(sweep(raw_assess, 2L, center), 2L, spread, "/"), optional = TRUE)
      assess$scl_params <- train$scl_params
      assessment_data <- data[ti, , drop = FALSE]
      for (nm in covariates) assessment_data[[nm]] <- assess$kernel_dat[[nm]]
      controls <- list(object = fit, newdata = assessment_data, kernel_inputs = assess,
        type = "response", se.fit = FALSE, uncertainty = "conditional", interval = "none")
      if (backend == "refit") controls$include_fixed <- TRUE
      mu <- do.call(stats::predict, controls)
      if (!is.numeric(mu) || length(mu) != nrow(pred)) stop("Prediction rows do not match the assessment set.", call. = FALSE)
      output$prediction_contracts[[name]] <- attr(mu, "prediction_contract", exact = TRUE)
      pred$prediction <- as.numeric(mu)
      stage <- "scoring"
      pred <- .msr_cv_score(pred, family, metrics)
      NULL
    }, warning = function(w) {
      messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) list(stage = stage, message = conditionMessage(e)))
    complete <- is.null(error)
    if (!complete) {
      pred$prediction <- NA_real_; pred$eligible <- FALSE
      for (metric in metrics) pred[[metric]] <- NA_real_
    }
    output$predictions[[name]] <- pred
    output$metrics[[name]] <- do.call(rbind, lapply(metrics, function(metric)
      cbind(fold = name, .msr_cv_metric(pred, metric, family, complete, nrow(pred)))))
    output$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai),
      n_assessment = length(ti), status = if (complete) "scored" else "failed", stringsAsFactors = FALSE)
    output$warnings[name] <- list(messages); output$errors[name] <- list(error)
    output$diagnostics[name] <- list(list(
      initial_glm = if (is.null(initial)) NULL else initial[c("converged", "boundary", "iter")],
      fitted_glm = if (is.null(fit)) NULL else fit$opt_mod[c("converged", "boundary", "iter")],
      scale = if (is.null(fit)) NULL else fit$diagnostics))
    if (keep_fits) output$fits[name] <- list(fit)
  }
  for (item in c("predictions", "metrics", "fold_results")) {
    output[[item]] <- do.call(rbind, output[[item]])
    rownames(output[[item]]) <- NULL
  }
  output$complete <- all(output$fold_results$status == "scored")
  output$summary <- do.call(rbind, lapply(metrics, function(metric)
    .msr_cv_metric(output$predictions, metric, family, output$complete, nrow(output$predictions))))
  output$folds <- folds
  output$assessment_counts <- stats::setNames(tabulate(match(output$predictions$row_id, ids), length(ids)), ids)
  output$settings <- list(family = family, link = family_object$link, backend = backend,
    metrics = metrics, optim_args = optim_args, glm_control = glm_control,
    normalization = "training_only", bin_reference = "training_only_declared_buffer_domain",
    prediction_target = "held_out_GLM_expected_response", parameter_uncertainty = "not_integrated",
    observation_intervals = "not_provided", calibration = "not_established", workers = 1L)
  output$call <- match.call()
  class(output) <- "multiScaleR_cv"
  output
}

#' @rdname cross_validate
#' @param x A \code{multiScaleR_cv} result.
#' @param ... Reserved for future print controls.
#' @export
print.multiScaleR_cv <- function(x, ...) {
  cat("Cross-validation:", sum(x$fold_results$status == "scored"), "of",
    nrow(x$fold_results), "folds scored;", sum(x$assessment_counts > 0), "of",
    length(x$assessment_counts), "input rows assessed.\n")
  if (!x$complete) cat("Incomplete: inspect errors and diagnostics before comparing models.\n")
  print(x$summary, row.names = FALSE)
  invisible(x)
}

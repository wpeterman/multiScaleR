# NB/Gamma CV keeps distribution parameters in the training fold, separately
# from regression/scale uncertainty and from native coefficient SE conventions.
.msr_cv_distribution_family <- function(x) is.character(x) && length(x) == 1L && x %in% c("nbinom2", "Gamma")

.msr_cv_fold_settings <- function(x, name) {
  if (!is.null(x$settings$fold_settings)) x$settings$fold_settings[[name]] else x$settings
}

.msr_cv_gamma_family <- function(dispersion = NULL) {
  family <- stats::Gamma(link = "log")
  if (!is.null(dispersion)) {
    # Native Gamma$aic estimates its own deviance scale even if the newer
    # dispersion field is fixed. Construct the declared known-scale density
    # before glm() fits, without changing any global method or existing fit.
    family$dispersion <- dispersion
    family$msr_cv_fixed_dispersion <- dispersion
    phi <- dispersion
    aic <- function(y, n, mu, wt, dev) {
      used <- wt > 0
      -2 * sum(wt[used] * stats::dgamma(y[used], shape = 1/phi,
        scale = mu[used] * phi, log = TRUE))
    }
    environment(aic) <- list2env(list(phi = dispersion), parent = baseenv())
    family$aic <- aic
  }
  family
}

.msr_cv_distribution_control <- function(family, theta, dispersion) {
  scalar <- function(x) is.numeric(x) && is.null(dim(x)) && length(x) == 1L && is.finite(x) && x > 0
  if (!.msr_cv_distribution_family(family) ||
      (!is.null(theta) && (!identical(family, "nbinom2") || !scalar(theta))) ||
      (!is.null(dispersion) && (!identical(family, "Gamma") || !scalar(dispersion))))
    stop("theta is a finite positive fixed NB2 size; dispersion is a finite positive fixed Gamma dispersion. Use NULL to estimate the applicable nuisance parameter.", call. = FALSE)
  list(family = family, theta = if (is.null(theta)) NULL else as.numeric(theta),
    dispersion = if (is.null(dispersion)) NULL else as.numeric(dispersion))
}

.msr_cv_distribution_outcome <- function(formula, data, family) {
  expr <- formula[[2L]]
  while (is.call(expr) && identical(expr[[1L]], as.name("I"))) expr <- expr[[2L]]
  if (!is.symbol(expr)) stop("NB/Gamma CV requires a direct response column, optionally wrapped in I().", call. = FALSE)
  y <- eval(expr, envir = data, enclos = environment(formula))
  if ((!is.numeric(y) && !is.logical(y)) || !is.null(dim(y)) || length(y) != nrow(data) || any(!is.finite(y)))
    stop("NB/Gamma outcomes must be finite numeric vectors matching data rows.", call. = FALSE)
  if ((family == "nbinom2" && any(y < 0 | y != floor(y))) || (family == "Gamma" && any(y <= 0)))
    stop("NB2 requires nonnegative integer counts; Gamma requires strictly positive outcomes, including zero-weight rows.", call. = FALSE)
  list(observed = as.numeric(y), trials = rep(1, nrow(data)))
}

.msr_cv_distribution_fit <- function(formula, data, weights, control, distribution) {
  common <- list(formula = formula, data = data, weights = unname(weights),
    control = control, na.action = stats::na.fail, model = TRUE)
  if (distribution$family == "nbinom2" && is.null(distribution$theta)) {
    if (!requireNamespace("MASS", quietly = TRUE)) stop("Estimated NB2 CV requires MASS.", call. = FALSE)
    do.call(MASS::glm.nb, common)
  } else {
    family <- if (distribution$family == "Gamma") .msr_cv_gamma_family(distribution$dispersion) else {
      if (!requireNamespace("MASS", quietly = TRUE)) stop("NB2 CV requires MASS.", call. = FALSE)
      MASS::negative.binomial(distribution$theta, link = "log")
    }
    do.call(stats::glm, c(common, list(family = family)))
  }
}

.msr_cv_distribution <- function(formula, data, folds, family, weights, glm_control, metrics,
    keep_fits, keep_selection, theta, dispersion, kernel_inputs = NULL,
    backend = "glm", optim_args = list()) {
  distribution <- .msr_cv_distribution_control(family, theta, dispersion)
  validate_scalar_logical(keep_fits, "keep_fits"); validate_scalar_logical(keep_selection, "keep_selection")
  if (!inherits(formula, "formula") || length(formula) != 3L || !is.data.frame(data) ||
      !nrow(data) || anyDuplicated(names(data))) stop("Supply a two-sided formula and nonempty data with unique columns.", call. = FALSE)
  ids <- .msr_fold_ids(rownames(data), nrow(data)); inputs <- kernel_inputs
  spatial <- if (is.null(inputs)) character() else inputs$scale_vars$covariate
  if (length(setdiff(all.vars(formula), c(names(data), spatial, "."))))
    stop("Supply every formula variable in data or prepared spatial inputs.", call. = FALSE)
  if (length(intersect(all.vars(formula[[2L]]), spatial))) stop("A spatial predictor cannot be the CV response.", call. = FALSE)
  folds <- .msr_cv_folds(folds, ids)
  if (is.null(metrics)) metrics <- c("rmse", "mae", "neg_log_score")
  if (!is.character(metrics) || !length(metrics) || anyNA(metrics) || anyDuplicated(metrics) ||
      any(!metrics %in% c("rmse", "mae", "neg_log_score")))
    stop("NB/Gamma metrics are rmse, mae and neg_log_score; nuisance-dependent deviances are not interchangeable comparison losses.", call. = FALSE)
  if (is.null(weights)) weights <- rep(1, nrow(data))
  if (!is.numeric(weights) || !is.null(dim(weights)) || length(weights) != nrow(data) ||
      any(!is.finite(weights) | weights < 0 | weights != floor(weights)) ||
      (!is.null(names(weights)) && !identical(names(weights), ids)))
    stop("weights must be nonnegative integer case frequencies aligned with data rows.", call. = FALSE)
  if (!is.list(glm_control) || is.null(names(glm_control)) || anyDuplicated(names(glm_control)) ||
      any(!names(glm_control) %in% names(formals(stats::glm.control)))) stop("Use named glm.control() settings.", call. = FALSE)
  glm_control <- do.call(stats::glm.control, glm_control)
  covariates <- if ("." %in% all.vars(formula)) spatial else intersect(spatial, all.vars(formula))
  if (!is.null(inputs)) {
    if (!inherits(inputs, "multiScaleR_data") || !is.list(inputs$raw_cov) || !is.list(inputs$d_list) ||
        length(inputs$raw_cov) != nrow(inputs$kernel_dat) || length(inputs$d_list) != nrow(inputs$kernel_dat) ||
        !.msr_valid_observation_ids(rownames(inputs$kernel_dat), nrow(inputs$kernel_dat)) ||
        !setequal(ids, rownames(inputs$kernel_dat)) || !length(covariates) ||
        !any(inputs$scale_vars$optimize[match(covariates, spatial)]))
      stop("Spatial CV requires matching raw-cell inputs and an optimized formula predictor.", call. = FALSE)
    if (!is.null(inputs$raster_group) && (!is.character(inputs$raster_group) ||
        !identical(names(inputs$raster_group), rownames(inputs$kernel_dat)) ||
        anyNA(inputs$raster_group) || any(!nzchar(inputs$raster_group))))
      stop("raster_group must align exactly with prepared row IDs.", call. = FALSE)
    if ((family == "Gamma" || !is.null(theta)) && backend != "refit")
      stop("Gamma and fixed-theta NB2 scale fits require backend = 'refit'; analytic NB backends estimate theta.", call. = FALSE)
    allowed <- setdiff(names(formals(multiScale_optim)), c("fitted_mod", "kernel_inputs", "join_by", "n_cores", "PSOCK", "refit_fn", "backend", "verbose"))
    if (!is.list(optim_args) || (length(optim_args) && (is.null(names(optim_args)) || anyNA(names(optim_args)) ||
        anyDuplicated(names(optim_args)) || any(!names(optim_args) %in% allowed))) ||
        (!is.null(optim_args$profile_control) && !is.list(optim_args$profile_control)) ||
        !is.null(optim_args$profile_control$worker_pool))
      stop("optim_args must name supported serial controls without data/worker overrides.", call. = FALSE)
  }
  settings <- list(family = family, link = "log", backend = backend, metrics = metrics,
    glm_control = glm_control, distribution_control = distribution, normalization = "training_only",
    prediction_target = "held_out_GLM_expected_response", density_target = "training_fitted_response_distribution",
    parameter_uncertainty = "not_integrated", observation_intervals = "not_provided", calibration = "not_established", workers = 1L)
  if (!is.null(inputs)) settings$optim_args <- optim_args
  output <- list(predictions = list(), metrics = list(), fold_results = list(), training = list(),
    diagnostics = list(), distribution = list(), warnings = list(), errors = list(), prediction_contracts = list(),
    fits = if (keep_fits) list() else NULL, selection_records = if (keep_selection) list() else NULL)
  if (!is.null(inputs$raster_group)) output$raster_groups <- list()
  targets <- .msr_cv_targets(formula, data, family, weights, folds)
  for (name in names(folds)) {
    split <- folds[[name]]; ai <- match(split$analysis, ids); ti <- match(split$assessment, ids)
    pred <- data.frame(fold = name, row_id = split$assessment, observed = NA_real_, trials = 1,
      prediction = NA_real_, weight = unname(weights[ti]), eligible = FALSE,
      distribution_family = family, theta = NA_real_, dispersion = NA_real_, poisson_limit = FALSE, stringsAsFactors = FALSE)
    for (metric in metrics) pred[[metric]] <- NA_real_
    fit <- initial <- NULL; nuisance <- NULL; messages <- character(); stage <- "preparation"
    error <- tryCatch(withCallingHandlers({
      observed <- .msr_cv_outcome(formula, data[ti, , drop = FALSE], family)
      pred$observed <- observed$observed
      .msr_cv_outcome(formula, data[ai, , drop = FALSE], family)
      training <- data[ai, , drop = FALSE]; assessment <- data[ti, , drop = FALSE]
      if (!is.null(inputs)) {
        train <- .msr_cv_inputs(inputs, split$analysis, covariates)
        raw <- .msr_cv_initial_values(train); center <- colMeans(raw); spread <- apply(raw, 2L, stats::sd)
        if (any(!is.finite(spread) | spread <= 0)) stop("Training spatial predictors require positive finite spread.")
        train$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"), optional = TRUE)
        train$scl_params <- list(mean = center, sd = spread)
        for (nm in covariates) training[[nm]] <- train$kernel_dat[[nm]]
        output$training[[name]] <- list(initial_normalization = train$scl_params,
          bin_reference = if (is.null(train$binned)) NULL else train$binned$cv_reference)
      }
      stage <- "initial_glm"
      initializer <- distribution
      # The existing joint-NB2 engine must assess its exact alpha=0 boundary
      # itself. A fixed-size native GLM initializes coefficients without first
      # requiring MASS's finite-theta alternation to converge at that boundary.
      if (!is.null(inputs) && backend == "joint_nb2" && family == "nbinom2") initializer$theta <- 1
      initial <- .msr_cv_distribution_fit(formula, training, weights[ai], glm_control, initializer)
      fit <- initial
      output$training[[name]]$initial_coefficients <- stats::coef(initial)
      output$training[[name]]$initial_distribution <- .msr_cv_distribution_parameters(initial, initializer)
      if (!is.null(inputs)) {
        stage <- "scale_fit"; fit <- NULL
        fit <- do.call(multiScale_optim, c(list(fitted_mod = initial, kernel_inputs = train, n_cores = 1L,
          backend = backend, verbose = FALSE), optim_args))
        if (!identical(as.integer(fit$optim_results$convergence), 0L)) stop("Training scale optimizer did not converge.")
        output$training[[name]] <- c(output$training[[name]], list(fitted_normalization = fit$scl_params,
          coefficients = stats::coef(fit$opt_mod), scale_est = fit$scale_est, shape_est = fit$shape_est,
          optimizer = fit$optim_results[c("par", "value", "convergence", "message")], start_results = fit$profile_fit$start_results))
      } else output$training[[name]]$coefficients <- stats::coef(fit)
      stage <- "training_identity"
      .msr_cv_distribution_training_identity(fit, split, targets, family)
      nuisance <- .msr_cv_distribution_parameters(fit, distribution)
      stage <- "prediction"
      if (is.null(inputs)) {
        pred$prediction <- as.numeric(stats::predict(.msr_native_fixed_metadata(fit), newdata = assessment, type = "response"))
      } else {
        assess <- .msr_cv_inputs(inputs, split$assessment, covariates, if (is.null(train$binned)) NULL else train$binned$cv_reference)
        raw <- .msr_cv_initial_values(assess)
        assess$kernel_dat <- as.data.frame(sweep(sweep(raw, 2L, center), 2L, spread, "/"), optional = TRUE)
        assess$scl_params <- train$scl_params
        for (nm in covariates) assessment[[nm]] <- assess$kernel_dat[[nm]]
        args <- list(object = fit, newdata = assessment, kernel_inputs = assess, type = "response",
          se.fit = FALSE, uncertainty = "conditional", interval = "none")
        if (backend == "refit") args$include_fixed <- TRUE
        mu <- do.call(stats::predict, args)
        output$prediction_contracts[[name]] <- attr(mu, "prediction_contract", exact = TRUE)
        pred$prediction <- as.numeric(mu)
      }
      if (length(pred$prediction) != nrow(pred)) stop("Prediction rows do not match assessment rows.")
      pred$theta <- nuisance$theta; pred$dispersion <- nuisance$dispersion; pred$poisson_limit <- nuisance$poisson_limit
      stage <- "scoring"; pred <- .msr_cv_distribution_score(pred, metrics); NULL
    }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(stage = stage, message = conditionMessage(e), start_results = e$start_results))
    complete <- is.null(error)
    if (!complete) {
      pred$prediction <- NA_real_; pred$eligible <- FALSE
      for (metric in metrics) pred[[metric]] <- NA_real_
    }
    output$predictions[[name]] <- pred
    output$metrics[[name]] <- do.call(rbind, lapply(metrics, function(metric) cbind(fold = name,
      .msr_cv_metric(pred, metric, family, complete, nrow(pred)))))
    output$fold_results[[name]] <- data.frame(fold = name, n_analysis = length(ai), n_assessment = length(ti),
      status = if (complete) "scored" else "failed", stringsAsFactors = FALSE)
    output$warnings[name] <- list(messages); output$errors[name] <- list(error); output$distribution[name] <- list(nuisance)
    native <- if (inherits(fit, "multiScaleR")) fit$opt_mod else fit
    output$diagnostics[name] <- list(list(initial_glm = if (is.null(initial)) NULL else initial[c("converged", "boundary", "iter", "th.warn", "SE.theta")],
      fitted_glm = if (is.null(native)) NULL else native[c("converged", "boundary", "iter", "th.warn", "SE.theta")],
      scale = if (inherits(fit, "multiScaleR")) fit$diagnostics else NULL))
    if (is.null(inputs)) output$prediction_contracts[[name]] <- .msr_prediction_contract("likelihood", "glm",
      "expected_response", "response", population = "fixed_at_estimate", kernel = "not_applicable",
      conditioning = "training_fit_and_fixed_prediction_rows", omitted = c("parameter_uncertainty", "selection_uncertainty"), fixed_inputs = TRUE)
    if (!is.null(inputs$raster_group)) output$raster_groups[[name]] <- lapply(split, function(rows) inputs$raster_group[rows])
    if (keep_fits) output$fits[name] <- list(fit)
    if (keep_selection) output$selection_records[name] <- list(.msr_cv_distribution_record(fit, output$diagnostics[[name]], .msr_cv_selection_controls(settings)))
  }
  output <- .msr_cv_finish(output, folds, ids, family, metrics)
  output$targets <- targets; output$settings <- settings; output$call <- match.call()
  class(output) <- c("multiScaleR_distribution_cv", "multiScaleR_cv")
  output
}

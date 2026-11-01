# A descriptor binds already computed rasters to their fitted transformation.
# It stores no cells, fit, formula environments, or prediction arrays.
.msr_projection_signature <- function(object) list(
  kernel = object$kernel_inputs$kernel, scales = object$scale_est$Mean,
  scale_names = rownames(object$scale_est), shapes = .msr_shape_values(object),
  optimized_center = object$scl_params$mean[rownames(object$scale_est)],
  optimized_spread = object$scl_params$sd[rownames(object$scale_est)],
  preparation_normalization = object$kernel_inputs$scl_params,
  specifications = object$kernel_inputs$scale_vars,
  support = c(object$min_D, object$max_D), resolution = object$kernel_inputs$resolution)

.msr_stamp_projection <- function(raster, object, scale_center, match_buffer,
                                   scale_vars, include_fixed, site_values, controls = NULL) {
  if (!inherits(object, "multiScaleR")) return(raster)
  fixed <- if (is.null(scale_vars)) character() else
    scale_vars$covariate[!scale_vars$optimize]
  attr(raster, "multiScaleR_projection") <- list(schema_version = 1L,
    signature = .msr_projection_signature(object), standardized = scale_center,
    match_buffer = match_buffer, specifications = scale_vars,
    unverified_fixed = if (include_fixed) character() else fixed,
    scenario_values = site_values, controls = controls)
  raster
}

.msr_native_raster_eval <- function(model, data) {
  for (nm in intersect(model$logical_names, names(data))) {
    if (any(!is.na(data[[nm]]) & !data[[nm]] %in% c(0, 1)))
      stop("A logical raster predictor must encode FALSE/TRUE as 0/1.", call. = FALSE)
    data[[nm]] <- as.logical(data[[nm]])
  }
  if (!is.null(model$analytic)) {
    contract <- model$analytic$contract
    terms <- stats::delete.response(contract$terms)
    absent <- setdiff(all.vars(terms), names(data))
    if (length(absent)) stop("Raster prediction is missing variables: ",
      paste(absent, collapse = ", "), call. = FALSE)
    frame <- stats::model.frame(terms, data = data, xlev = contract$xlevels,
                                na.action = stats::na.fail)
    X <- stats::model.matrix(terms, frame, contrasts.arg = contract$contrasts)
    X <- X[, contract$columns, drop = FALSE]
    offset <- stats::model.offset(frame)
    if (is.null(offset)) offset <- rep(0, nrow(X))
    if (any(!is.finite(X)) || any(!is.finite(offset)))
      stop("Raster prediction design and offsets must be finite.", call. = FALSE)
    prototype <- model$analytic$prototype
    prototype$design <- list(X = X, offset = offset, rows = rownames(frame))
    value <- predict.multiScaleR(prototype,
      type = if (is.null(model$type)) "response" else model$type,
      se.fit = model$se.fit, uncertainty = "conditional")
  } else value <- .msr_native_prediction_rows(model$object, data, model$type,
                                              model$se.fit, model$native_args)
  if (is.list(value) && !is.data.frame(value) && !is.null(value$fit))
    value <- data.frame(fit = value$fit, se.fit = value$se.fit)
  # Terra dispatches its writer on result class. Provider matrix subclasses
  # retain every numerical column, but need a plain matrix at this boundary.
  if (is.matrix(value) && is.object(value)) value <- unclass(value)
  if (NROW(value) != nrow(data))
    stop("The native component must return one prediction row per raster cell.", call. = FALSE)
  value
}

# Rasters are already smoothed and normalized. This adapter never replaces
# posterior scale integration with a single scale or creates joint scale SEs.
.msr_predict_raster <- function(object, raster, type, se.fit, uncertainty, interval, ...) {
  validate_scalar_logical(se.fit, "se.fit")
  if (interval != "none" || (se.fit && uncertainty != "conditional"))
    stop("Raster predictions support native conditional SEs only; use uncertainty = 'conditional' and interval = 'none'.", call. = FALSE)
  descriptor <- attr(raster, "multiScaleR_projection")
  if (is.null(descriptor) || !isTRUE(descriptor$standardized) ||
      !isTRUE(descriptor$match_buffer) ||
      !identical(descriptor$signature, .msr_projection_signature(object)))
    stop("Supply kernel_scale.raster() output from this spatial fit with scale_center = TRUE and match_buffer = TRUE.", call. = FALSE)
  if (length(descriptor$unverified_fixed))
    stop("Computed fixed-radius predictors need include_fixed = TRUE for verified training normalization.", call. = FALSE)
  expected <- object$kernel_inputs$scale_vars
  actual <- descriptor$specifications
  if (!is.null(actual)) {
    if (is.null(expected)) stop("The fit must retain matching spatial source specifications.", call. = FALSE)
    expected <- expected[match(actual$covariate, expected$covariate), , drop = FALSE]
    if (!isTRUE(all.equal(as.data.frame(actual), as.data.frame(expected), check.attributes = FALSE)))
      stop("Projected source and metric specifications must match the fitted model.", call. = FALSE)
    used <- if (methods::is(object$opt_mod, "unmarkedFit"))
      .unmarked_model_predictors(object$opt_mod) else .model_predictors(.analysis_model(object$opt_mod))
    required <- intersect(object$kernel_inputs$scale_vars$covariate, used)
    if (length(setdiff(required, c(actual$covariate, names(descriptor$scenario_values)))))
      stop("Every modeled spatial predictor needs a projected layer or an explicit site_values scenario.", call. = FALSE)
  }
  model <- .analysis_model(object$opt_mod)
  if (methods::is(model, "unmarkedFit") && is.null(type))
    stop("Specify the native unmarked prediction component explicitly.", call. = FALSE)
  if (identical(object$backend, "conditional_logit")) {
    if (is.null(type) || identical(type, "response"))
      stop("Choice-set probabilities need complete point choice sets; raster prediction requires type = 'link' or 'risk'.", call. = FALSE)
  }
  dots <- list(...)
  raster_args <- dots$raster_args
  dots$raster_args <- NULL
  if (is.null(raster_args)) raster_args <- list()
  if (!is.list(raster_args) || (length(raster_args) &&
      (is.null(names(raster_args)) || anyDuplicated(names(raster_args)) ||
       any(!names(raster_args) %in% c("filename", "overwrite", "wopt", "const", "cores", "na.rm")))))
    stop("raster_args must uniquely name supported terra prediction controls.", call. = FALSE)
  if (!is.null(raster_args$const) &&
      length(intersect(names(raster_args$const), names(raster))))
    stop("Raster const inputs must add variables; set existing predictor scenarios with site_values.",
      call. = FALSE)
  nested <- dots$native_args
  dots$native_args <- NULL
  if (!is.null(nested)) {
    if (!is.list(nested) || (length(nested) && (is.null(names(nested)) ||
        anyNA(names(nested)) || any(!nzchar(names(nested))))) || anyDuplicated(names(nested)) ||
        any(names(nested) %in% c(names(dots), "object", "newdata", "type", "se.fit")))
      stop("native_args must uniquely name native prediction controls.", call. = FALSE)
    dots <- c(dots, nested)
  }
  dots <- .msr_native_control_names(model, dots)
  if (.msr_native_lm_interval(model, dots, type) != "none" && uncertainty != "conditional")
    stop("Native intervals require uncertainty = 'conditional'; kernel uncertainty is not included.", call. = FALSE)
  if (isTRUE(as.logical(dots[["cov.fit"]])))
    stop("Cross-row native covariance requires point predictions; raster chunks cannot retain it.", call. = FALSE)
  if (inherits(model, "HLfit") && (length(dots[["intervals"]]) ||
      any(vapply(dots[["variances"]], function(x) !identical(x, FALSE), logical(1)))))
    stop("Native variance and interval attributes require point predictions; raster output cannot retain them.", call. = FALSE)
  if (any(names(dots) %in% c("join_by", "include_fixed")))
    stop("Point preparation controls do not apply to computed rasters.", call. = FALSE)
  if (identical(object$backend, "conditional_logit") && is.null(dots$reference))
    dots$reference <- "zero"
  frame <- if (methods::is(model, "unmarkedFit")) model@data@siteCovs else
    .msr_selection_frame(model)
  if (is.null(frame) && is.list(model) && is.data.frame(model$data)) frame <- model$data
  logical_names <- if (is.null(frame)) character() else
    names(frame)[vapply(frame, is.logical, logical(1))]
  thin <- list(object = list(opt_mod = model), type = type, se.fit = se.fit,
    native_args = dots, logical_names = logical_names)
  if (.msr_analytic_backend(object)) {
    if (length(dots)) stop("Native controls are unsupported for analytic raster prediction.", call. = FALSE)
    if (!is.null(type) && !type %in% c("response", "link"))
      stop("Analytic raster prediction type must be response or link.", call. = FALSE)
    if (se.fit) vcov.multiScaleR(object, type = "conditional")
    contract <- object$profile_fit$contract
    prototype <- structure(list(backend = object$backend, uncertainty = object$uncertainty,
      opt_mod = list(family = object$opt_mod$family),
      profile_fit = list(contract = list(covariates = contract$covariates)),
      kernel_inputs = list(kernel = object$kernel_inputs$kernel)), class = "multiScaleR_prediction")
    thin$analytic <- list(contract = list(terms = contract$terms, xlevels = contract$xlevels,
      contrasts = contract$contrasts,
      columns = if (is.null(contract$objective$response_design))
        c(contract$fixed_columns, contract$covariates) else
          contract$objective$response_design$columns),
      prototype = prototype)
    thin$object <- NULL
  }
  cores <- if (is.null(raster_args$cores)) 1L else raster_args$cores
  owned_cluster <- NULL
  if (is.numeric(cores)) {
    validate_scalar_numeric(cores, "raster_args$cores", integerish = TRUE, positive = TRUE)
    if (cores > 1L) {
      inputs <- object$kernel_inputs
      inputs$min_D <- object$min_D; inputs$max_D <- object$max_D
      .msr_stop_unsafe_parallel_ram(inputs, model, object$join_by, NULL, cores, TRUE)
      owned_cluster <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(owned_cluster), add = TRUE)
      cluster_prep(model, owned_cluster)
      raster_args$cores <- owned_cluster
    }
  }
  if (is.null(raster_args$na.rm)) raster_args$na.rm <- TRUE
  output <- do.call(terra::predict, c(list(object = raster, model = thin,
    fun = .msr_native_raster_eval), raster_args))
  attr(output, "prediction_contract") <- .msr_native_prediction_contract(object,
    if (methods::is(model, "unmarkedFit")) "unmarked" else class(model)[1L],
    if (is.null(type)) if (inherits(model, "coxph")) "lp" else "response" else
      if (inherits(model, "coxph") && type == "link") "lp" else type, dots)
  if (.msr_analytic_backend(object))
    attr(output, "prediction_contract") <- .msr_prediction_contract("likelihood",
      object$backend, if (identical(type, "link")) "linear_predictor" else "expected_response",
      output_scale = if (is.null(type)) "response" else type,
      population = "analytic_coefficient_information_conditional_on_kernel_estimate",
      kernel = "fixed_at_estimate", conditioning = "selected_kernel_parameters",
      intervals = "conditional_standard_errors_only",
      omitted = c("kernel_parameter_uncertainty", "observation_noise"), fixed_inputs = TRUE)
  attr(output, "prediction_contract") <- .msr_returned_prediction_uncertainty(
    attr(output, "prediction_contract"), se.fit)
  attr(output, "projection_contract") <- descriptor
  output
}

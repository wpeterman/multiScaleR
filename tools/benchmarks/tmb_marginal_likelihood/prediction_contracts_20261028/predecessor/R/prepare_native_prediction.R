# Freeze native numerical outputs. No native fit, formula environment or raw
# spatial inputs are retained, so an RDS snapshot remains independently usable.
.msr_prepare_native_prediction <- function(object, newdata, kernel_inputs,
    max_bytes, component, native_args, join_by = NULL, include_fixed = FALSE) {
  kind <- if (identical(object$backend, "joint_mixed")) "glmmTMB" else "unmarked"
  allowed <- if (kind == "glmmTMB") c("re.form", "allow.new.levels", "na.action") else
    c("backTransform", "na.rm", "appendData", "level", "re.form")
  if (!is.list(native_args) || (length(native_args) &&
      (is.null(names(native_args)) || anyNA(names(native_args)) ||
       any(!names(native_args) %in% allowed) || anyDuplicated(names(native_args)))))
    stop("native_args must uniquely name supported native prediction controls.", call. = FALSE)
  if (!is.null(component) && (!is.character(component) || length(component) != 1L ||
      is.na(component) || !nzchar(component)))
    stop("component must be one native prediction type.", call. = FALSE)
  if (kind == "unmarked") {
    if (is.null(component)) stop("Specify the native unmarked component explicitly.", call. = FALSE)
    model <- object$opt_mod
    nd <- if (is.null(kernel_inputs)) newdata else
      prepare_unmarked_data(object, newdata, kernel_inputs, join_by, include_fixed)
    if (is.null(kernel_inputs) && !is.null(join_by))
      stop("join_by requires explicit kernel_inputs and newdata.", call. = FALSE)
    if (is.null(kernel_inputs) && include_fixed)
      stop("include_fixed requires explicit kernel_inputs and newdata.", call. = FALSE)
    types <- component
  } else {
    model <- object$opt_mod
    nd <- .msr_mixed_prediction_data(object, newdata, kernel_inputs)
    types <- if (is.null(component)) c("response", "link") else component
  }
  capture <- function(args) {
    warnings <- character()
    # unmarked exports an S4 generic distinct from stats' S3 generic.
    predictor <- if (kind == "unmarked") unmarked::predict else stats::predict
    value <- withCallingHandlers(do.call(predictor, args), warning = function(w)
      warnings <<- c(warnings, conditionMessage(w)))
    list(value = value, warnings = warnings)
  }
  results <- stats::setNames(vector("list", length(types)), types)
  for (type in types) {
    args <- c(list(object = model, type = type), native_args)
    # unmarked distinguishes missing newdata from a supplied NULL.
    if (kind == "glmmTMB" || !is.null(nd)) args$newdata <- nd
    if (kind == "glmmTMB") {
      args$se.fit <- FALSE
      point <- capture(args)
      args$se.fit <- TRUE
      results[[type]] <- list(point = point, se = capture(args))
    } else results[[type]] <- list(table = capture(args))
  }
  first <- if (kind == "glmmTMB") results[[1L]]$point$value else results[[1L]]$table$value
  out <- structure(list(backend = if (kind == "unmarked") "unmarked" else object$backend,
    native_kind = kind, results = results, default_type = types[1L],
    n_rows = if (kind == "unmarked") nrow(first) else length(first),
    provenance = list(component = component, native_args = lapply(native_args, function(x)
      if (is.function(x) || is.language(x)) deparse(x) else x),
      scale_est = object$scale_est, native_class = class(model))),
    class = "multiScaleR_prediction")
  if (kind == "unmarked" && !is.null(kernel_inputs)) {
    # Retain row-to-site provenance, not buffers or a full scenario data frame.
    fixed <- if (include_fixed) {
      specs <- object$kernel_inputs$scale_vars
      specs$covariate[!specs$optimize &
        specs$covariate %in% .unmarked_model_predictors(model)]
    } else character()
    out$provenance$spatial <- list(kernel = kernel_inputs$kernel,
      unit_conv = kernel_inputs$unit_conv,
      include_fixed = include_fixed,
      covariates = c(names(object$scl_params$mean), fixed),
      fixed_covariates = fixed,
      row_names = rownames(newdata),
      site_keys = stats::setNames(lapply(names(join_by), function(nm) newdata[[nm]]),
                                 names(join_by)))
  }
  if (as.numeric(utils::object.size(out)) > max_bytes)
    stop("Prepared predictions exceed max_bytes; use fewer rows or an explicit larger budget.", call. = FALSE)
  out
}

.msr_predict_native_snapshot <- function(object, type, se.fit, uncertainty, interval) {
  if (!is.character(type) || length(type) != 1L || is.na(type) ||
      !type %in% names(object$results))
    stop("This native prediction type was not prepared; call prepare_prediction() again.", call. = FALSE)
  if (object$native_kind == "unmarked") item <- object$results[[type]]$table else {
    uncertainty <- match.arg(uncertainty, c("joint", "conditional"))
    interval <- match.arg(interval, c("none", "confidence"))
    if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit))
      stop("se.fit must be one logical value.", call. = FALSE)
    if (interval != "none" || (se.fit && uncertainty != "conditional"))
      stop("joint_mixed supports conditional prediction SEs only; use uncertainty = 'conditional' and interval = 'none'.", call. = FALSE)
    item <- object$results[[type]][[if (se.fit) "se" else "point"]]
  }
  for (message in item$warnings) warning(message, call. = FALSE)
  item$value
}

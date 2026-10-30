# Refit snapshots cache numerical outputs only. An unavailable native SE is
# recorded separately, so preparing point predictions never manufactures it.
.msr_prepare_refit_prediction <- function(object, newdata, kernel_inputs,
    max_bytes, component, native_args, join_by, include_fixed) {
  model <- .analysis_model(object$opt_mod)
  allowed <- c("re.form", "allow.new.levels", "na.action", "reference", "level")
  if (inherits(model, "gam")) allowed <- c(allowed,
    "exclude", "terms", "unconditional", "iterms.type", "block.size")
  if (!is.list(native_args) || (length(native_args) &&
      (is.null(names(native_args)) || anyNA(names(native_args)) ||
       any(!names(native_args) %in% allowed) || anyDuplicated(names(native_args)))))
    stop("native_args must uniquely name supported native prediction controls.", call. = FALSE)
  model <- .analysis_model(object$opt_mod)
  if (is.null(model)) stop("Prediction preparation requires a retained fitted native model.", call. = FALSE)
  if (is.null(component)) {
    if (inherits(model, "coxph")) types <- c("lp", "risk") else
    if (inherits(model, c("lm", "glm", "glmmTMB", "merMod", "gls", "lme")))
      types <- c("response", "link") else
      stop("Specify the native model's prediction component explicitly.", call. = FALSE)
  } else {
    validate_character_scalar(component, "component")
    types <- component
  }
  capture <- function(type, se) {
    warnings <- character()
    result <- tryCatch(withCallingHandlers(
      .msr_predict_refit(object, newdata, kernel_inputs, type, se,
        "conditional", "none", join_by = join_by, include_fixed = include_fixed,
        native_args = native_args),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(result, "error")) {
      if (!se) stop(conditionMessage(result), call. = FALSE)
      return(list(value = NULL, warnings = warnings, error = conditionMessage(result)))
    }
    list(value = result, warnings = warnings)
  }
  results <- stats::setNames(lapply(types, function(type)
    list(point = capture(type, FALSE), se = capture(type, TRUE))), types)
  first <- results[[1L]]$point$value
  out <- structure(list(backend = "refit", native_kind = "refit", results = results,
    default_type = types[1L], n_rows = if (is.matrix(first) || is.data.frame(first))
      nrow(first) else length(first), provenance = list(component = component,
      native_args = lapply(native_args, function(x)
        if (is.function(x) || is.language(x)) deparse(x) else x),
      scale_est = object$scale_est, native_class = class(model),
      include_fixed = include_fixed, explicit_site_join = !is.null(join_by),
      row_names = if (is.null(newdata)) names(first) else rownames(newdata),
      spatial = list(kernel = object$kernel_inputs$kernel,
        fitted_unit_conv = object$kernel_inputs$unit_conv,
        prediction_unit_conv = if (is.null(kernel_inputs)) NULL else kernel_inputs$unit_conv,
        source_specs = object$kernel_inputs$scale_vars,
        training_center = object$scl_params$mean,
        training_spread = object$scl_params$sd))),
    class = "multiScaleR_prediction")
  out$prediction_contract <- stats::setNames(lapply(types, function(type) {
    contract <- attr(results[[type]]$point$value, "prediction_contract")
    contract$uncertainty$conditional_se_available <- is.null(results[[type]]$se$error)
    contract
  }), types)
  if (as.numeric(utils::object.size(out)) > max_bytes)
    stop("Prepared predictions exceed max_bytes; use fewer rows or an explicit larger budget.", call. = FALSE)
  out
}

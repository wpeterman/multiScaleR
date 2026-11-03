# A supplied raw reference verifies observable fitted features. It cannot prove
# a unique inverse of a square, polynomial basis or another many-to-one map.
.msr_curve_reference_context <- function(object, model, reference_data) {
  if (!inherits(model, c("lm", "glm", "gam", "coxph")))
    stop("reference_data is not supported for this native model structure.", call. = FALSE)
  if (!is.data.frame(reference_data) || !nrow(reference_data) ||
      anyNA(names(reference_data)) || anyDuplicated(names(reference_data)) ||
      any(!nzchar(names(reference_data))))
    stop("reference_data must be a nonempty data frame with unique column names.", call. = FALSE)
  if (.row_names_info(reference_data, 1L) <= 0L ||
      anyNA(row.names(reference_data)) || anyDuplicated(row.names(reference_data)))
    stop("reference_data requires unique explicit fitted row IDs as row names.", call. = FALSE)
  model <- .msr_native_fixed_metadata(model)
  frame <- model$model
  original <- attr(frame, "terms")
  if (!is.data.frame(frame) || !inherits(original, "terms") || !nrow(frame))
    stop("reference_data requires the retained evaluated native training frame.", call. = FALSE)
  if (!is.null(model$call$offset))
    stop("reference_data does not support offsets supplied outside the model formula.", call. = FALSE)
  ids <- row.names(frame)
  if (length(ids) != nrow(reference_data) || !setequal(ids, row.names(reference_data)))
    stop("reference_data row IDs must match every retained fitted row exactly.", call. = FALSE)
  data <- reference_data[match(ids, row.names(reference_data)), , drop = FALSE]
  attr(data, "terms") <- attr(data, "na.action") <- NULL
  tt <- .msr_fixed_prediction_terms(stats::terms(model),
    .msr_response_frame_subset(original, frame, stats::terms(model)))
  if (inherits(model, "coxph") && (inherits(model, c("coxphms", "coxph.penal")) ||
      length(attr(tt, "specials")$tt)))
    stop("reference_data requires a checked time-fixed, unpenalized native Cox structure.", call. = FALSE)
  # Cox fits also use their retained prediction terms, without changing any
  # fitted coefficient, covariance, response, likelihood or diagnostic.
  if (inherits(model, "coxph")) {
    model$terms <- tt
    attr(model$model, "terms") <- tt
  }
  prediction_terms <- stats::delete.response(tt)
  calls <- as.list(attr(prediction_terms, "variables"))[-1L]
  required <- unique(unlist(lapply(calls, all.vars), use.names = FALSE))
  conditioning <- if (inherits(model, "coxph"))
    .msr_cox_curve_conditioning(tt)$inputs else character()
  kernel_source <- "caller_supplied_fitted_spatial_columns"
  if (.msr_analytic_backend(object) || identical(object$backend, "conditional_logit")) {
    engine <- object$profile_fit$contract$objective$engine
    if (inherits(engine, "msr_surface_engine") && isTRUE(engine$rebuild_from_inputs))
      engine <- .msr_prepare_surface_engine(object$kernel_inputs, engine$covariates,
        engine$rows, engine$standardize, engine$implementation)
    values <- .msr_eval_kernel_engine(engine,
      object$profile_fit$par[seq_along(engine$parameter_names)])$value
    if (nrow(values) != nrow(frame))
      stop("Fitted kernel rows do not match the retained effect-curve rows.", call. = FALSE)
    for (nm in engine$covariates) data[[nm]] <- values[, nm]
    kernel_source <- "verified_selected_preparation"
  }
  absent <- setdiff(required, names(data))
  if (length(absent)) stop("reference_data is missing raw inputs: ",
    paste(absent, collapse = ", "), call. = FALSE)
  for (nm in required) {
    value <- data[[nm]]
    if (!is.atomic(value) || !is.null(dim(value)) || is.complex(value) ||
        anyNA(value) || (is.numeric(value) && any(!is.finite(value))))
      stop("reference_data raw inputs must be finite scalar columns: ", nm, call. = FALSE)
  }
  .msr_curve_check_learned_calls(attr(prediction_terms, "predvars"))
  # Explicit supplied contrasts must agree with the fitted coding. Once
  # checked, use fitted levels and saved native contrasts in prediction grids.
  native_levels <- if (inherits(model, "gam")) .msr_gam_factor_context(model)$levels else
    model$xlevels
  for (nm in intersect(names(native_levels), required)) {
    value <- data[[nm]]
    if (is.factor(value) && length(setdiff(levels(value), native_levels[[nm]])))
      stop("reference_data contains unknown factor levels: ", nm, call. = FALSE)
    provided <- attr(value, "contrasts")
    saved <- model$contrasts[[nm]]
    if (!is.null(provided) && !is.null(saved)) {
      canonical <- factor(as.character(value), levels = native_levels[[nm]],
        ordered = is.ordered(value))
      expected <- canonical; stats::contrasts(expected) <- saved
      if (!isTRUE(all.equal(stats::contrasts(value), stats::contrasts(expected),
          tolerance = 1e-10)))
        stop("reference_data contrasts differ from the fitted coding: ", nm, call. = FALSE)
    }
    data[[nm]] <- factor(as.character(value), levels = native_levels[[nm]],
      ordered = is.ordered(value))
  }
  evaluated <- stats::model.frame(prediction_terms, data = data,
    xlev = model$xlevels, na.action = stats::na.fail)
  # Native model-frame reconstruction may drop factor attributes. Restore the
  # saved coding for native design methods that rely on those attributes.
  for (nm in intersect(names(model$contrasts), names(evaluated)))
    if (is.factor(evaluated[[nm]])) stats::contrasts(evaluated[[nm]]) <- model$contrasts[[nm]]
  retained <- .msr_response_frame_subset(original, frame, prediction_terms)
  if (!identical(names(evaluated), names(retained)) ||
      !identical(row.names(evaluated), ids))
    stop("reference_data evaluated features do not align to the fitted frame.", call. = FALSE)
  equivalent <- function(a, b) {
    if (is.list(a) || is.list(b))
      return(is.list(a) && is.list(b) && identical(names(a), names(b)) &&
        length(a) == length(b) && all(vapply(seq_along(a),
          function(i) equivalent(a[[i]], b[[i]]), logical(1))))
    if (is.numeric(a) || is.numeric(b))
      return(is.numeric(a) && is.numeric(b) && identical(dim(a), dim(b)) &&
        all(is.finite(a)) && all(is.finite(b)) &&
        isTRUE(all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE)))
    if (is.factor(a) || is.factor(b))
      return(is.factor(a) && is.factor(b) && identical(as.character(a), as.character(b)) &&
        identical(levels(a), levels(b)) && identical(is.ordered(a), is.ordered(b)))
    identical(a, b)
  }
  for (nm in names(evaluated)) if (!equivalent(evaluated[[nm]], retained[[nm]]))
    stop("reference_data does not reproduce the fitted feature: ", nm, call. = FALSE)
  offset <- function(mf) {
    value <- stats::model.offset(mf)
    if (is.null(value)) rep(0, nrow(mf)) else value
  }
  if (!equivalent(offset(evaluated), offset(frame)))
    stop("reference_data offsets differ from the fitted offsets.", call. = FALSE)
  if (inherits(model, "gam")) {
    # An explicit retained frame keeps only fitted rows under na.exclude.
    fitted_design <- .msr_predict_gam(model, newdata = frame, type = "lpmatrix")
    # Extra response columns must not activate a family's response-conditioned
    # prediction hook. Only the required predictor context is supplied.
    prediction_data <- if (length(required)) data[required] else frame
    reference_design <- .msr_predict_gam(model, newdata = prediction_data, type = "lpmatrix")
    if (!identical(attr(fitted_design, "lpi"), attr(reference_design, "lpi")))
      stop("reference_data GAM coefficient spans differ from the fitted basis.", call. = FALSE)
    if (!equivalent(attr(fitted_design, "model.offset"), attr(reference_design, "model.offset")) ||
        !equivalent(attr(fitted_design, "offset"), attr(reference_design, "offset")))
      stop("reference_data GAM offsets differ between fitted components.", call. = FALSE)
  } else if (inherits(model, "coxph")) {
    # The native method drops strata before constructing the coefficient design.
    fitted_design <- stats::model.matrix(model)
    reference_design <- stats::model.matrix(model, data = evaluated)
  } else {
    fitted_design <- stats::model.matrix(model)
    reference_design <- stats::model.matrix(prediction_terms, evaluated,
      contrasts.arg = model$contrasts)
  }
  if (!identical(dim(fitted_design), dim(reference_design)) ||
      !identical(colnames(fitted_design), colnames(reference_design)) ||
      !equivalent(fitted_design, reference_design))
    stop("reference_data does not reproduce the fitted native design or basis.", call. = FALSE)
  if (inherits(model, "gam")) {
    eta <- .msr_predict_gam(model, newdata = prediction_data, type = "link")
    saved_eta <- model$linear.predictors
    if (length(model$nsdf) == 1L) {
      # predict.gam returns a one-dimensional array while a fitted single
      # predictor can be stored as a plain vector. Keep row/value checks.
      single <- function(x) length(x) == nrow(frame) &&
        (is.null(dim(x)) || (length(dim(x)) <= 2L && dim(x)[1L] == nrow(frame) &&
          all(dim(x)[-1L] == 1L)))
      if (!single(eta) || !single(saved_eta))
        stop("Cannot align native GAM linear predictor rows.", call. = FALSE)
      eta <- as.numeric(eta); saved_eta <- as.numeric(saved_eta)
    }
    if (!equivalent(eta, saved_eta))
      stop("Cannot verify that native GAM prediction calls reproduce the fitted linear predictors.", call. = FALSE)
  }
  list(data = data, model = model, required_inputs = required,
    conditioning_variables = conditioning,
    contract = list(schema_version = 1L, fitted_rows = nrow(frame),
      alignment = "explicit_fitted_row_IDs", fitted_feature_equivalence = TRUE,
      native_design_and_offsets_verified = TRUE, equivalence_tolerance = 1e-10,
      original_raw_values_uniquely_established = FALSE, kernel_source = kernel_source,
      native_conditioning_variables = conditioning, metadata_establishes_eligibility = FALSE))
}

# Fixed metadata must already be in native prediction calls. Learned calls
# hidden inside other expressions cannot be relearned on a plotting grid.
.msr_curve_check_learned_calls <- function(expr) {
  walk <- function(ex) {
    if (!is.call(ex)) return(invisible(NULL))
    head <- ex[[1L]]
    name <- if (is.symbol(head)) as.character(head) else if (is.call(head) &&
        length(head) == 3L && as.character(head[[1L]]) %in% c("::", ":::"))
      as.character(head[[3L]]) else ""
    if (name %in% c("scale", "poly", "ns", "bs")) {
      if (name %in% c("ns", "bs") && !requireNamespace("splines", quietly = TRUE))
        stop("Saved spline prediction metadata requires the splines package.", call. = FALSE)
      definition <- switch(name, scale = base::scale, poly = stats::poly,
        ns = getExportedValue("splines", "ns"), bs = getExportedValue("splines", "bs"))
      call <- match.call(definition, ex, expand.dots = TRUE)
      fixed <- function(x) is.numeric(x) || identical(x, FALSE)
      known <- switch(name, scale = fixed(call$center) && fixed(call$scale),
        poly = isTRUE(call$raw) || is.list(call$coefs),
        ns = , bs = is.numeric(call$Boundary.knots) &&
          (is.numeric(call$knots) || (is.null(call$df) && is.null(call$knots))))
      if (!isTRUE(known)) stop("Effect curves require saved fixed prediction metadata; ",
        "nested learned transformations need fixed training calls: ", name, call. = FALSE)
    }
    for (child in as.list(ex)[-1L]) walk(child)
    invisible(NULL)
  }
  walk(expr)
}

# Keep a native family column on its stated scale. Response columns need not
# correspond one-to-one with linear predictors, for example multinomial fits.
.msr_gam_curve_component <- function(model, newdata, component, response) {
  value <- .msr_predict_gam(model, newdata = newdata,
    type = if (response) "response" else "link", se.fit = TRUE)
  fit <- if (is.list(value)) value$fit else value
  if (response && .msr_gam_ziplss_model(model) && !is.matrix(fit)) {
    n <- nrow(newdata)
    scalar <- function(x) is.numeric(x) && length(x) == n &&
      (is.null(dim(x)) || (length(dim(x)) == 1L && dim(x)[1L] == n))
    if (component != 1L || !scalar(fit))
      stop("ziplss has one overall expected-count response; select component = 1.", call. = FALSE)
    if (is.list(value) && !is.null(value$se.fit)) {
      if (!scalar(value$se.fit) || !identical(dim(fit), dim(value$se.fit)))
        stop("Native GAM standard errors do not match its scalar response.", call. = FALSE)
      return(list(fit = as.numeric(fit), se.fit = as.numeric(value$se.fit)))
    }
    return(as.numeric(fit))
  }
  if (!is.matrix(fit) || component > ncol(fit))
    stop("component must select an available native GAM output column.", call. = FALSE)
  if (is.list(value) && !is.null(value$se.fit)) {
    if (!identical(dim(fit), dim(value$se.fit)))
      stop("Native GAM standard errors do not match its output columns.", call. = FALSE)
    return(list(fit = fit[, component], se.fit = value$se.fit[, component]))
  }
  fit[, component]
}

# Native strata/cluster inputs condition a curve. Group-only inputs do not
# become effect panels, and their numeric labels must not be averaged.
.msr_cox_curve_conditioning <- function(terms) {
  special <- function(ex) {
    if (!is.call(ex)) return(FALSE)
    head <- ex[[1L]]
    name <- if (is.symbol(head)) as.character(head) else if (is.call(head) &&
        length(head) == 3L && as.character(head[[1L]]) %in% c("::", ":::"))
      as.character(head[[3L]]) else ""
    name %in% c("strata", "cluster")
  }
  calls <- as.list(attr(terms, "variables"))[-1L]
  selected <- calls[vapply(calls, special, logical(1))]
  inputs <- unique(unlist(lapply(selected, all.vars), use.names = FALSE))
  labels <- lapply(attr(terms, "term.labels"), str2lang)
  fixed <- labels[!vapply(labels, special, logical(1))]
  panels <- unique(unlist(lapply(fixed, .msr_marginal_variables,
    offsets = FALSE), use.names = FALSE))
  list(inputs = inputs, fixed_panel_inputs = panels)
}

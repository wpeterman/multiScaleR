# CV must retain learned transformations even when a wrapper such as I() or
# offset() prevents native makepredictcall() from seeing the learner itself.
# Only prediction calls change. Formula labels, fitted data and estimates do not.
.msr_cv_prediction_terms <- function(terms, data, qualify = FALSE, allow_learning = TRUE) {
  validate_scalar_logical(qualify, "qualify")
  validate_scalar_logical(allow_learning, "allow_learning")
  env <- environment(terms)
  known <- list(scale = "base", poly = "stats", ns = "splines", bs = "splines")
  pure <- .msr_prediction_pointwise_functions()
  # An identically named user function is not a base/stats preprocessing rule.
  trusted <- function(head, name, package) {
    if (is.symbol(head)) {
      if (!identical(head, as.name(name))) return(FALSE)
      return(.msr_standard_prediction_function(
        get0(name, envir = env, mode = "function", inherits = TRUE), name, package))
    }
    is.call(head) && length(head) == 3L && identical(head[[1L]], as.name("::")) &&
      identical(head[[2L]], as.name(package)) && identical(head[[3L]], as.name(name))
  }
  missing_arg <- function(parts, j) identical(unname(parts[j]), unname(alist(x = )))
  constant <- function(value) {
    if (is.object(value)) return(FALSE)
    if (is.null(value) || is.numeric(value) || is.logical(value)) return(TRUE)
    is.list(value) && all(vapply(value, constant, logical(1)))
  }
  normalize_scale <- function(ex) {
    parts <- as.list(ex)
    # Native makepredictcall can append exact names while retaining an original
    # abbreviation. Those aliases would otherwise become unused arguments.
    for (parameter in c("center", "scale")) {
      nms <- names(parts)
      if (parameter %in% nms) {
        aliases <- which(nzchar(nms) & startsWith(parameter, nms) & nms != parameter)
        if (length(aliases)) parts <- parts[-aliases]
      }
    }
    match.call(definition = base::scale, call = as.call(parts))
  }
  safe <- function(ex) {
    if (is.null(ex) || is.numeric(ex) || is.logical(ex) || is.character(ex)) return(TRUE)
    if (is.list(ex)) return(all(vapply(ex, safe, logical(1))))
    if (is.symbol(ex)) {
      name <- as.character(ex)
      value <- if (name %in% names(data)) data[[name]] else get0(name, envir = env, inherits = TRUE)
      return(constant(value))
    }
    if (!is.call(ex)) return(FALSE)
    name <- if (is.symbol(ex[[1L]])) as.character(ex[[1L]]) else
      if (is.call(ex[[1L]]) && length(ex[[1L]]) == 3L) as.character(ex[[1L]][[3L]]) else ""
    package <- if (name %in% names(known)) known[[name]] else "base"
    if (!name %in% c(pure, names(known)) || !trusted(ex[[1L]], name, package)) return(FALSE)
    parts <- as.list(ex)
    all(vapply(seq_along(parts)[-1L], function(j)
      missing_arg(parts, j) || safe(parts[[j]]), logical(1)))
  }
  walk <- function(ex) {
    if (!is.call(ex)) return(ex)
    parts <- as.list(ex)
    for (j in seq_along(parts)[-1L]) if (!missing_arg(parts, j)) parts[[j]] <- walk(parts[[j]])
    ex <- as.call(parts)
    if (qualify && is.symbol(ex[[1L]]) && as.character(ex[[1L]]) %in% names(known)) {
      name <- as.character(ex[[1L]])
      if (!trusted(ex[[1L]], name, known[[name]]))
        stop("Unqualified custom learner '", name,
          "' cannot be frozen as a standard learner; namespace-qualify it or precompute its predictor.", call. = FALSE)
    }
    for (name in names(known)) if (trusted(ex[[1L]], name, known[[name]])) {
      if (name == "scale") {
        ex <- normalize_scale(ex)
        parts <- as.list(ex)
      }
      fixed <- switch(name, scale = (is.numeric(ex$center) || identical(ex$center, FALSE)) &&
        (is.numeric(ex$scale) || identical(ex$scale, FALSE)),
        poly = isTRUE(ex$raw) || is.list(ex$coefs),
        ns = is.numeric(ex$knots) && is.numeric(ex$Boundary.knots),
        bs = is.numeric(ex$knots) && is.numeric(ex$Boundary.knots))
      if (fixed) {
        if (qualify) ex[[1L]] <- call("::", as.name(known[[name]]), as.name(name))
        return(ex)
      }
      if (!allow_learning)
        stop("Original training data are required to recover learned prediction metadata.", call. = FALSE)
      if (!all(vapply(seq_along(parts)[-1L], function(j)
          missing_arg(parts, j) || safe(parts[[j]]), logical(1))))
        stop(if (qualify)
          "Prediction metadata requires numeric inputs and supported pointwise arguments to scale/poly/ns/bs; precompute other learned transformations from original training data." else
          "Training-only prediction metadata requires numeric inputs and supported pointwise arguments to scale/poly/ns/bs; precompute other learned transformations within each analysis fold.", call. = FALSE)
      value <- eval(ex, envir = data, enclos = env)
      bare <- ex; bare[[1L]] <- as.name(name)
      # scale() can retain the class of its input (for example poly). Choose
      # the known learner's method rather than dispatching on that inner class.
      method <- utils::getS3method("makepredictcall", if (name == "scale") "default" else name,
        envir = asNamespace("stats"))
      learned <- method(value, bare)
      learned[[1L]] <- if (qualify) call("::", as.name(known[[name]]), as.name(name)) else ex[[1L]]
      return(if (name == "scale") normalize_scale(learned) else learned)
    }
    ex
  }
  prediction <- attr(terms, "predvars")
  if (is.null(prediction)) prediction <- attr(terms, "variables")
  attr(terms, "predvars") <- walk(prediction)
  terms
}

# Recover final training predictors from the fitted evaluator, never the initial
# scale or assessment rows. The analytic route retains its exact engine/row map.
.msr_cv_prediction_metadata <- function(fit, training, inputs = NULL) {
  contains <- function(ex) {
    if (!is.call(ex)) return(FALSE)
    h <- ex[[1L]]
    name <- if (is.symbol(h)) as.character(h) else
      if (is.call(h) && length(h) == 3L && identical(h[[1L]], as.name("::"))) as.character(h[[3L]]) else ""
    if (name %in% c("scale", "poly", "ns", "bs")) return(TRUE)
    parts <- as.list(ex)
    any(vapply(seq_along(parts)[-1L], function(j)
      if (identical(unname(parts[j]), unname(alist(x = )))) FALSE else contains(parts[[j]]), logical(1)))
  }
  model <- if (inherits(fit, "multiScaleR")) fit$opt_mod else fit
  tt <- stats::terms(model)
  if (!contains(attr(tt, "variables"))) return(fit)
  if (inherits(fit, "multiScaleR")) {
    contract <- fit$profile_fit$contract
    if (!is.null(contract)) {
      e <- contract$objective$engine
      if (inherits(e, "msr_surface_engine") && isTRUE(e$rebuild_from_inputs))
        e <- .msr_prepare_surface_engine(fit$kernel_inputs, e$covariates, e$rows,
          e$standardize, e$implementation)
      rows <- match(rownames(contract$frame), rownames(training))
      if (anyNA(rows) || anyDuplicated(rows) || length(rows) != nrow(training))
        stop("Cannot align final training predictors with CV analysis rows.", call. = FALSE)
      training <- training[rows, , drop = FALSE]
      values <- .msr_eval_kernel_engine(e, fit$profile_fit$par[seq_along(e$parameter_names)])$value
      if (nrow(values) != nrow(training))
        stop("Final training predictor rows do not match CV analysis data.", call. = FALSE)
      training[e$covariates] <- as.data.frame(values)
      fit$profile_fit$contract$terms <- .msr_cv_prediction_terms(contract$terms, training)
    } else training <- .msr_refit_prediction_data(fit, training, inputs, include_fixed = TRUE)
    fit$opt_mod <- .msr_cv_prediction_metadata(fit$opt_mod, training)
    return(fit)
  }
  fit$terms <- .msr_cv_prediction_terms(stats::terms(fit), training)
  if (is.data.frame(fit$model)) attr(fit$model, "terms") <- fit$terms
  fit
}

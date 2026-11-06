# One grammar for supported pointwise learner arguments. This preserves the
# existing CV whitelist and lets native frame recovery reject cohort functions.
.msr_prediction_pointwise_functions <- function() c(
  "+", "-", "*", "/", "^", "%%", "%/%", "(", "I", "c", "list", "log", "log1p",
  "log2", "log10", "exp", "expm1", "sqrt", "abs", "sin", "cos", "tan",
  "asin", "acos", "atan", "sinh", "cosh", "tanh", "round", "signif", "pmin", "pmax")

# terra can expose scale as an S4 generic whose numeric methods delegate to
# base's unchanged S3 generic. Verify that delegation instead of trusting its
# name or replacing an arbitrary custom learner with a base implementation.
.msr_standard_prediction_function <- function(fun, name, package) {
  expected <- get(name, asNamespace(package), inherits = FALSE)
  if (identical(fun, expected)) return(TRUE)
  if (!identical(name, "scale") || !identical(package, "base") ||
      !methods::is(fun, "standardGeneric") ||
      !identical(as.character(methods::slot(fun, "generic")), "scale") ||
      !identical(methods::slot(fun, "package"), "base") ||
      !identical(formals(fun), formals(expected)) ||
      !identical(body(fun), quote(standardGeneric("scale"))) ||
      !identical(get0("standardGeneric", environment(fun), mode = "function",
        inherits = TRUE), base::standardGeneric)) return(FALSE)
  same <- function(method) is.function(method) &&
    identical(formals(method), formals(expected)) &&
    identical(body(method), body(expected)) &&
    identical(environment(method), environment(expected))
  if (!same(methods::slot(fun, "default"))) return(FALSE)
  signatures <- expand.grid(x = c("numeric", "integer", "matrix", "array"),
    center = c("logical", "numeric"), scale = c("logical", "numeric"),
    stringsAsFactors = FALSE)
  all(vapply(seq_len(nrow(signatures)), function(i) {
    method <- tryCatch(methods::selectMethod("scale", as.character(signatures[i, ]),
      optional = TRUE, fdef = fun), error = function(e) NULL)
    same(method)
  }, logical(1)))
}

.msr_fixed_function_name <- function(head, env, learners_only = FALSE) {
  known <- list(scale = "base", poly = "stats", ns = "splines", bs = "splines")
  packages <- if (learners_only) known else c(known,
    stats::setNames(rep(list("base"), length(.msr_prediction_pointwise_functions())),
      .msr_prediction_pointwise_functions()), list(offset = "stats"))
  if (is.symbol(head)) {
    name <- as.character(head)
    if (!name %in% names(packages)) return(NULL)
    if (!.msr_standard_prediction_function(
        get0(name, env, mode = "function", inherits = TRUE), name, packages[[name]])) return(NULL)
    return(name)
  }
  if (!is.call(head) || length(head) != 3L || !identical(head[[1L]], as.name("::"))) return(NULL)
  name <- as.character(head[[3L]])
  if (length(name) != 1L || !name %in% names(packages) ||
      !identical(head[[2L]], as.name(packages[[name]]))) return(NULL)
  name
}

.msr_fixed_normalize_call <- function(ex, env) {
  tt <- stats::terms(~1)
  attr(tt, ".Environment") <- env
  attr(tt, "predvars") <- as.call(list(as.name("list"), ex))
  .msr_cv_prediction_terms(tt, list(), allow_learning = FALSE)
  # Validate completeness without renaming positional x arguments. glmmTMB's
  # native term matcher uses the original call prefix to find evaluated columns.
  # Keep the original tag/position when makepredictcall appends an exact name.
  # Native provider matching can depend on that original call prefix.
  walk <- function(x) {
    if (!is.call(x)) return(x)
    parts <- as.list(x)
    if (identical(.msr_fixed_function_name(x[[1L]], env, TRUE), "scale")) {
      for (parameter in c("center", "scale")) {
        nms <- names(parts)
        if (parameter %in% nms) {
          aliases <- which(nzchar(nms) & startsWith(parameter, nms) & nms != parameter)
          if (length(aliases)) {
            if (length(aliases) != 1L || sum(nms == parameter) != 1L)
              stop("Cannot align ambiguous duplicated scale arguments.", call. = FALSE)
            parts[[aliases]] <- parts[[which(nms == parameter)]]
            parts <- parts[-which(nms == parameter)]
          }
        }
      }
    }
    for (j in seq_along(parts)[-1L])
      if (!identical(unname(parts[j]), unname(alist(x = )))) parts[[j]] <- walk(parts[[j]])
    as.call(parts)
  }
  walk(ex)
}

# Recover a single learner from its own surviving original attributes. A sum
# of two learned bases can retain just one set, which cannot identify both.
.msr_fixed_nested_call <- function(ex, value, env) {
  complete <- tryCatch(list(call = .msr_fixed_normalize_call(ex, env)), error = function(e) NULL)
  if (!is.null(complete)) return(complete$call)
  missing_arg <- function(parts, j) identical(unname(parts[j]), unname(alist(x = )))
  collect <- function(x) {
    if (!is.call(x)) return(list())
    own <- if (!is.null(.msr_fixed_function_name(x[[1L]], env, TRUE))) list(x) else list()
    parts <- as.list(x)
    for (j in seq_along(parts)[-1L]) if (!missing_arg(parts, j)) own <- c(own, collect(parts[[j]]))
    own
  }
  calls <- collect(ex)
  if (length(calls) != 1L)
    stop("Nested learned prediction calls need unambiguous original training metadata; precompute multiple learners from training data.", call. = FALSE)
  safe <- function(x) {
    if (is.null(x) || is.atomic(x) || is.symbol(x)) return(TRUE)
    if (is.list(x)) return(all(vapply(x, safe, logical(1))))
    if (!is.call(x) || is.null(.msr_fixed_function_name(x[[1L]], env))) return(FALSE)
    parts <- as.list(x)
    all(vapply(seq_along(parts)[-1L], function(j)
      missing_arg(parts, j) || safe(parts[[j]]), logical(1)))
  }
  if (!safe(ex))
    stop("Nested learned prediction calls require supported pointwise expressions and original training metadata; precompute other transformations.", call. = FALSE)
  original <- calls[[1L]]; name <- .msr_fixed_function_name(original[[1L]], env, TRUE)
  bare <- original; bare[[1L]] <- as.name(name)
  method <- utils::getS3method("makepredictcall", if (name == "scale") "default" else name,
    envir = asNamespace("stats"))
  learned <- method(value, bare); learned[[1L]] <- original[[1L]]
  replace_call <- function(x) {
    if (identical(x, original)) return(learned)
    if (!is.call(x)) return(x)
    parts <- as.list(x)
    for (j in seq_along(parts)[-1L]) if (!missing_arg(parts, j)) parts[[j]] <- replace_call(parts[[j]])
    as.call(parts)
  }
  tryCatch(.msr_fixed_normalize_call(replace_call(ex), env), error = function(e)
    stop("Retained attributes cannot recover the original nested training learner; precompute its predictor or restore its original metadata.", call. = FALSE))
}

# Unlike the historical predicate, this considers predictor calls inside I()
# and offsets. Complete positional calls do not need an evaluated frame.
.msr_needs_nested_metadata <- function(terms) {
  if (!inherits(terms, "terms")) return(FALSE)
  tt <- stats::delete.response(terms); env <- environment(tt)
  prediction <- attr(tt, "predvars")
  if (is.null(prediction)) prediction <- attr(tt, "variables")
  failed <- tryCatch({.msr_cv_prediction_terms(tt, list(), allow_learning = FALSE); FALSE},
    error = function(e) TRUE)
  if (failed) return(TRUE)
  aliases <- function(x) {
    if (!is.call(x)) return(FALSE)
    if (identical(.msr_fixed_function_name(x[[1L]], env, TRUE), "scale")) {
      nms <- names(as.list(x))
      for (parameter in c("center", "scale")) if (parameter %in% nms &&
          any(nzchar(nms) & startsWith(parameter, nms) & nms != parameter)) return(TRUE)
    }
    parts <- as.list(x)
    any(vapply(seq_along(parts)[-1L], function(j)
      if (identical(unname(parts[j]), unname(alist(x = )))) FALSE else aliases(parts[[j]]), logical(1)))
  }
  aliases(prediction)
}

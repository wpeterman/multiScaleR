# Native GAM formulas can use a factor only inside an interaction. Some
# provider versions omit that factor when selecting saved contrasts. Scope
# the correction to one local prediction call, preserving the native method.
.msr_predict_gam <- function(object, ...) {
  if (!requireNamespace("mgcv", quietly = TRUE))
    stop("Native GAM prediction requires the mgcv package.", call. = FALSE)
  native <- utils::getS3method("predict", "gam", envir = asNamespace("mgcv"))
  method <- NULL
  for (kind in class(object)) {
    method <- utils::getS3method("predict", kind, optional = TRUE)
    if (!is.null(method)) break
  }
  # Retain native BAM branching, block controls and output assembly. Other
  # providers keep their own dispatch, without changing a fitting method.
  if (!identical(method, native)) {
    bam <- get("predict.bam", envir = asNamespace("mgcv"))
    if (identical(method, bam))
      return(.msr_predict_bam_saved(object, list(...), bam))
    return(stats::predict(object, ...))
  }
  args <- list(...)
  object <- .msr_gam_ziplss_response_metadata(object, args)
  if (!length(object$contrasts)) return(do.call(native, c(list(object = object), args)))
  .msr_predict_gam_saved(object, args, native)
}

# Keep the native BAM body and algorithms. Correct only its local delegation
# and ordinary worker-result assembly, without changing the provider namespace.
.msr_predict_bam_saved <- function(object, args, native) {
  gam <- get("predict.gam", envir = environment(native))
  scope <- new.env(parent = environment(native))
  scope$calls <- scope$parallel_calls <- 0L
  scope$warnings <- scope$warning_origins <- list()
  scope$constant_seen <- FALSE
  scope$constant <- NULL
  collect <- function(values) {
    lapply(values, function(value) {
      marker <- attr(value, ".msr_bam_delegate")
      if (!is.null(marker)) {
        attr(value, ".msr_bam_delegate") <- NULL
        seen <- any(vapply(scope$warning_origins, function(origin)
          identical(origin, marker$origin), logical(1)))
        if (length(marker$warnings) && !seen) {
          scope$warnings <- c(scope$warnings, marker$warnings)
          scope$warning_origins[[length(scope$warning_origins) + 1L]] <- marker$origin
        }
        if (isTRUE(marker$iterms)) {
          if (scope$constant_seen && !identical(scope$constant, marker$constant))
            stop("Cannot align the native BAM term intercept.", call. = FALSE)
          scope$constant <- marker$constant
          scope$constant_seen <- TRUE
        }
      }
      value
    })
  }
  scope$rbind <- function(...) do.call(base::rbind, collect(list(...)))
  scope$c <- function(...) {
    values <- collect(list(...))
    if (identical(args$type, "iterms") && scope$parallel_calls > 0L) {
      if (!all(vapply(values, is.matrix, logical(1))))
        stop("Cannot align native BAM term prediction matrices.", call. = FALSE)
      return(do.call(base::rbind, values))
    }
    do.call(base::c, values)
  }
  scope$predict.gam <- function(object, ...) {
    scope$calls <- scope$calls + 1L
    controls <- list(...)
    if (!length(object$contrasts))
      return(do.call(gam, c(list(object = object), controls)))
    .msr_predict_gam_saved(object, controls, gam)
  }
  # Native BAM drops general dots from its worker argument list. Bind those
  # controls to the portable callback while retaining explicit native args.
  controls <- args[setdiff(names(args), names(formals(native)))]
  worker <- NULL
  makeActiveBinding("pabapr", function(value) {
    if (!missing(value)) stop("Cannot replace the native BAM coding delegate.", call. = FALSE)
    scope$calls <- scope$calls + 1L
    scope$parallel_calls <- scope$parallel_calls + 1L
    if (is.null(worker)) worker <<- .msr_bam_saved_worker(gam, controls,
      coding = length(object$contrasts) > 0L)
    worker
  }, scope)
  # Workers keep original native warnings; also surface compact diagnostic copies
  # to the parent. Deduplicate only the same worker's fit/SE metadata bundle.
  on.exit(for (condition in scope$warnings) warning(condition), add = TRUE)
  adapted <- native
  environment(adapted) <- scope
  result <- do.call(adapted, c(list(object = object), args))
  discrete <- !is.null(object$dinfo) && if ("discrete" %in% names(args))
    isTRUE(as.logical(args$discrete)) else isTRUE(eval(formals(native)$discrete, environment(native)))
  rows <- if (!is.null(args$newdata)) nrow(args$newdata) else nrow(object$model)
  if (!discrete && (is.null(rows) || rows > 0L) && !scope$calls)
    stop("Cannot verify this native BAM prediction implementation.", call. = FALSE)
  if (identical(args$type, "iterms") && scope$constant_seen && !is.null(scope$constant))
    attr(result, "constant") <- scope$constant
  result
}

# All workers compact native call context; coding workers also need three helpers. The native method,
# controls and provider identity remain portable without a multiScaleR install.
.msr_bam_saved_worker <- function(native, controls = list(), coding = TRUE) {
  scope <- new.env(parent = baseenv())
  helpers <- c(".msr_bam_warning_context", if (coding)
    c(".msr_predict_gam_saved", ".msr_gam_factor_context", ".msr_response_frame_subset"))
  for (name in helpers) {
    fun <- get(name, envir = environment(.msr_predict_gam_saved))
    environment(fun) <- scope
    assign(name, fun, envir = scope)
  }
  .expected_version <- as.character(utils::packageVersion("mgcv"))
  .expected_body <- body(native)
  .expected_formals <- formals(native)
  .native_controls <- controls
  .needs_coding <- coding
  scope$.expected_version <- .expected_version
  scope$.expected_body <- .expected_body
  scope$.expected_formals <- .expected_formals
  scope$.native_controls <- .native_controls
  scope$.needs_coding <- .needs_coding
  worker <- function(arg) {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop("Cannot verify the native GAM provider on this BAM worker.", call. = FALSE)
    native <- get("predict.gam", envir = asNamespace("mgcv"))
    if (!identical(as.character(utils::packageVersion("mgcv")), .expected_version) ||
        !identical(body(native), .expected_body) || !identical(formals(native), .expected_formals))
      stop("Cannot verify the native GAM provider on this BAM worker.", call. = FALSE)
    object <- arg$object
    arg$object <- NULL
    shared <- intersect(names(arg), names(.native_controls))
    if (any(!vapply(shared, function(name) identical(arg[[name]], .native_controls[[name]]), logical(1))))
      stop("Cannot align the native BAM worker controls.", call. = FALSE)
    arg[names(.native_controls)] <- .native_controls
    warnings <- list()
    value <- withCallingHandlers(
      if (.needs_coding) .msr_predict_gam_saved(object, arg, native) else
        do.call(native, c(list(object = object), arg)),
      warning = function(w) warnings[[length(warnings) + 1L]] <<-
        .msr_bam_warning_context(w, native, .expected_version))
    iterms <- identical(arg$type, "iterms")
    if (length(warnings) || iterms) {
      marker <- list(warnings = warnings, iterms = iterms,
        origin = list(PID = Sys.getpid(), first_row = row.names(arg$newdata)[1L], rows = nrow(arg$newdata)),
        constant = if (iterms) attr(value, "constant") else NULL)
      mark <- function(x) { attr(x, ".msr_bam_delegate") <- marker; x }
      if (is.list(value) && all(c("fit", "se.fit") %in% names(value))) {
        value$fit <- mark(value$fit)
        value$se.fit <- mark(value$se.fit)
      } else value <- mark(value)
    }
    value
  }
  environment(worker) <- scope
  worker
}

# do.call can put an evaluated fit, rows and local function environment in a
# warning call. Compact only a verified native GAM call in the transported
# copy. The worker's original condition, message, class and controls stay intact.
.msr_bam_warning_context <- function(condition, native, version = NULL) {
  call <- conditionCall(condition)
  if (!is.list(condition) || !is.call(call)) return(condition)
  head <- call[[1L]]
  if (!is.function(head) || !identical(body(head), body(native)) ||
      !identical(formals(head), formals(native))) return(condition)
  call[[1L]] <- quote(mgcv::predict.gam)
  omitted <- character()
  if ("object" %in% names(call) && inherits(call[["object"]], "gam")) {
    call[["object"]] <- quote(.native_BAM_model)
    omitted <- c(omitted, "object")
  }
  if ("newdata" %in% names(call) &&
      (is.data.frame(call[["newdata"]]) || is.matrix(call[["newdata"]]))) {
    call[["newdata"]] <- quote(.native_prediction_rows)
    omitted <- c(omitted, "newdata")
  }
  condition$call <- call
  attr(condition, "multiScaleR_native_warning") <- list(provider = "mgcv",
    method = "predict.gam", version = version, symbolic_data_arguments = omitted)
  condition
}

.msr_predict_gam_saved <- function(object, args, native) {
  context <- .msr_gam_factor_context(object)
  saved <- context$saved; keys <- names(saved)
  retained <- context$levels; ordered <- context$ordered
  scope <- new.env(parent = environment(native))
  scope$calls <- 0L
  scope$model.matrix <- function(object, data, contrasts.arg = NULL, ...) {
    scope$calls <- scope$calls + 1L
    # An evaluated combined frame can contain factors from other predictors.
    # Select this component's actual variables before passing its coding.
    source_terms <- attr(data, "terms")
    if (!inherits(source_terms, "terms"))
      stop("Cannot align the native GAM coefficient frame.", call. = FALSE)
    data <- .msr_response_frame_subset(source_terms, data, object)
    active <- intersect(unique(keys), names(data))
    for (key in active) {
      value <- data[[key]]
      # Logical terms also carry saved contrasts. Keep their native logical
      # type; model.matrix supplies R's two levels independently of this block.
      if (key %in% context$logical) {
        if (!is.logical(value))
          stop("GAM logical input does not match its fitted type: ", key, call. = FALSE)
        next
      }
      target <- retained[[key]]
      if (is.null(target))
        stop("Saved GAM factor levels are unavailable: ", key, call. = FALSE)
      if (!is.factor(value) && !is.character(value))
        stop("GAM factor input does not match its fitted type: ", key, call. = FALSE)
      observed <- as.character(value[!is.na(value)])
      if (length(setdiff(observed, target)))
        stop("GAM factor input contains new levels: ", key, call. = FALSE)
      if (!is.factor(value) || !identical(levels(value), target) ||
          !identical(is.ordered(value), ordered[[key]])) {
        dims <- dim(value); dim_names <- dimnames(value)
        value <- factor(as.character(value), levels = target, ordered = ordered[[key]])
        if (!is.null(dims)) { dim(value) <- dims; dimnames(value) <- dim_names }
        data[[key]] <- value
      }
    }
    # Evaluate the unchanged native coefficient design with its saved coding.
    stats::model.matrix(object, data,
      contrasts.arg = saved[keys %in% names(data)], ...)
  }
  adapted <- native
  environment(adapted) <- scope
  result <- do.call(adapted, c(list(object = object), args))
  rows <- if (!is.null(args$newdata)) nrow(args$newdata) else nrow(object$model)
  if (!identical(args$type, "newdata") && (is.null(rows) || rows > 0L) && !scope$calls)
    stop("Cannot verify this native GAM prediction implementation.", call. = FALSE)
  result
}

# Saved levels and coding are training metadata, not prediction-row estimates.
.msr_gam_factor_context <- function(object) {
  saved <- object$contrasts
  keys <- names(saved)
  if (length(saved) && (is.null(keys) || anyNA(keys) || any(!nzchar(keys))))
    stop("Cannot align saved GAM factor contrasts.", call. = FALSE)
  if (is.null(keys)) keys <- character()
  for (key in unique(keys)) {
    entries <- saved[which(keys == key)]
    if (!all(vapply(entries, function(x) identical(x, entries[[1L]]), logical(1))))
      stop("Cannot align conflicting saved GAM factor contrasts: ", key, call. = FALSE)
  }
  # Native GAM families may store one named level list per linear predictor.
  # Collect those retained levels, rather than infer them from a new block.
  retained <- list()
  collect <- function(x) {
    if (!is.list(x)) return(invisible(NULL))
    for (i in seq_along(x)) {
      key <- names(x)[i]
      if (is.list(x[[i]])) collect(x[[i]]) else if (length(key) &&
          !is.na(key) && nzchar(key) && is.character(x[[i]])) {
        if (!is.null(retained[[key]]) && !identical(retained[[key]], x[[i]]))
          stop("Cannot align retained GAM factor levels: ", key, call. = FALSE)
        retained[[key]] <<- x[[i]]
      }
    }
    invisible(NULL)
  }
  collect(object$xlevels)
  frame_factors <- names(Filter(is.factor, as.list(object$model)))
  factor_names <- union(unique(keys), frame_factors)
  ordered <- stats::setNames(vector("list", length(factor_names)), factor_names)
  classes <- attr(object$terms, "dataClasses")
  for (key in factor_names) {
    value <- object$model[[key]]
    if (is.factor(value)) {
      if (!is.null(retained[[key]]) && !identical(retained[[key]], levels(value)))
        stop("Cannot align fitted GAM factor levels: ", key, call. = FALSE)
      retained[[key]] <- levels(value)
      ordered[[key]] <- is.ordered(value)
    } else {
      known <- if (key %in% names(classes)) classes[[key]] else NULL
      ordered[[key]] <- identical(known, "ordered")
    }
  }
  logical_keys <- unique(keys[vapply(keys, function(key) {
    is.logical(object$model[[key]]) ||
      (key %in% names(classes) && identical(classes[[key]], "logical"))
  }, logical(1))])
  # Logical columns are not factor inputs to raw-reference validation.
  retained[intersect(names(retained), logical_keys)] <- NULL
  list(saved = saved, levels = retained, ordered = ordered, logical = logical_keys)
}

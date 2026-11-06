# Keep the native component predictor and its bounded covariance chunks. Only
# construction of fixed designs changes, after checking the original fitted X.
.msr_unmarked_fixed_design <- function(basis, data) {
  absent <- setdiff(all.vars(basis$terms), names(data))
  if (length(absent)) stop("Prediction data must explicitly contain formula variables: ",
    paste(absent, collapse = ", "), call. = FALSE)
  # xlev reconstructs factors and discards their temporary contrast attributes.
  # The ORIGINAL fitted matrix's contrasts are applied explicitly below.
  for (name in intersect(names(basis$xlevels), names(data)))
    if (is.factor(data[[name]])) attr(data[[name]], "contrasts") <- NULL
  frame <- stats::model.frame(basis$terms, data, na.action = stats::na.pass,
    xlev = basis$xlevels, drop.unused.levels = FALSE)
  X <- stats::model.matrix(basis$terms, frame, contrasts.arg = basis$contrasts)
  if (!identical(colnames(X), basis$columns))
    stop("The prediction design does not match the fitted unmarked columns.", call. = FALSE)
  list(X = X, offset = stats::model.offset(frame))
}

# Subsetting a native fitted matrix can discard its contrasts attribute. For
# a bare main factor, its constant fitted columns identify the ORIGINAL code
# for each observed level. Recover those values, never guess from session options.
.msr_unmarked_training_contrasts <- function(terms, frame, fitted, keep) {
  out <- attr(fitted, "contrasts")
  factors <- names(frame)[vapply(frame, is.factor, logical(1))]
  missing <- setdiff(factors, names(out))
  if (!length(missing)) return(out)
  for (name in missing) if (!is.null(attr(frame[[name]], "contrasts")))
    out[[name]] <- attr(frame[[name]], "contrasts")
  missing <- setdiff(factors, names(out))
  if (!length(missing)) return(out)
  trial <- stats::model.matrix(terms, frame, contrasts.arg = out)
  labels <- attr(terms, "term.labels")
  syntax <- lapply(labels, str2lang)
  fail <- function() stop("Original factor contrasts cannot be recovered unambiguously from the fitted unmarked design; retain explicit training contrasts.", call. = FALSE)
  if (ncol(trial) != ncol(fitted)) fail()
  for (name in missing) {
    term <- which(vapply(syntax, identical, logical(1), as.name(name)))
    if (length(term) != 1L) fail()
    columns <- which(attr(trial, "assign") == term)
    factor_levels <- levels(frame[[name]])
    if (!length(columns) || !length(factor_levels)) fail()
    rows_by_level <- split(seq_len(sum(keep)), frame[[name]][keep], drop = FALSE)
    values <- lapply(rows_by_level, function(rows) {
      if (!length(rows)) fail()
      block <- fitted[rows, columns, drop = FALSE]
      first <- as.numeric(block[1L, ])
      if (any(!is.finite(block)) || !isTRUE(max(abs(sweep(block, 2L, first))) <= 1e-12)) fail()
      first
    })
    prefix <- labels[term]
    column_names <- colnames(fitted)[columns]
    if (any(!startsWith(column_names, prefix))) fail()
    suffix <- substring(column_names, nchar(prefix) + 1L)
    if (any(!nzchar(suffix)) || anyDuplicated(suffix)) fail()
    out[[name]] <- matrix(unlist(values, use.names = FALSE), nrow = length(factor_levels),
      byrow = TRUE, dimnames = list(factor_levels, suffix))
  }
  out
}

.msr_unmarked_prediction_basis <- function(object, type, verify_random = TRUE) {
  ns <- asNamespace("unmarked")
  raw <- get("get_covariates", ns)(object, type)
  formula <- get("get_formula", ns)(object, type)
  required <- if (verify_random) all.vars(formula) else
    .msr_marginal_variables(formula, offsets = TRUE)
  absent <- setdiff(required, c(names(raw), "."))
  if (length(absent)) stop("Original unmarked covariates must explicitly contain formula variables: ",
    paste(absent, collapse = ", "), call. = FALSE)
  # C/R fits retain the actual post-removal design in the likelihood closure.
  # TMB stores the same fixed matrices and offsets in its original data.
  dm <- if (is.list(object@TMB) && is.environment(object@TMB$env))
    object@TMB$env$data else get0("dm", environment(object@nllFun), inherits = FALSE)
  fitted <- dm[[paste0("X_", type)]]
  fitted_offset <- dm[[paste0("offset_", type)]]
  if (!is.matrix(fitted) || !is.numeric(fitted_offset))
    stop("Original fitted unmarked designs and offsets are required to freeze training transformations.", call. = FALSE)
  tt <- .msr_unmarked_curve_basis(formula, raw)$terms
  frame <- stats::model.frame(tt, raw, na.action = stats::na.pass,
    drop.unused.levels = FALSE)
  keep <- rep(TRUE, nrow(raw))
  removed <- object@sitesRemoved
  if (length(removed)) {
    removed <- if (type == "det") unlist(lapply(removed, function(i)
      (i - 1L) * ncol(object@data@y) + seq_len(ncol(object@data@y))), use.names = FALSE) else removed
    keep[removed] <- FALSE
  }
  basis <- list(terms = tt,
    xlevels = lapply(frame[vapply(frame, is.factor, logical(1))], levels),
    contrasts = .msr_unmarked_training_contrasts(tt, frame, fitted, keep),
    columns = colnames(fitted))
  original <- .msr_unmarked_fixed_design(basis, raw)
  offset <- original$offset
  if (is.null(offset)) offset <- rep(0, nrow(original$X))
  offset[is.na(offset)] <- 0
  same <- identical(dim(original$X[keep, , drop = FALSE]), dim(fitted)) &&
    isTRUE(all.equal(unname(original$X[keep, , drop = FALSE]), unname(fitted),
      tolerance = 1e-12, check.attributes = FALSE)) &&
    isTRUE(all.equal(as.numeric(offset[keep]), as.numeric(fitted_offset), tolerance = 1e-12))
  if (!same)
    stop("Original training transformations do not match the retained fitted unmarked design and offsets.", call. = FALSE)
  if (verify_random && get("has_random", ns)(formula)) {
    fitted_Z <- dm[[paste0("Z_", type)]]
    Z <- get("get_Z", ns)(formula, raw, raw)[keep, , drop = FALSE]
    if (is.null(dim(fitted_Z)) || !identical(dim(Z), dim(fitted_Z)) ||
        !isTRUE(max(abs(Z - fitted_Z)) <= 1e-12))
      stop("Original unmarked random designs do not match the retained fitted design.", call. = FALSE)
  }
  basis
}

# Native frame preparation would drop grouping columns for a subset of sites.
# Retain the ORIGINAL group's coordinate order, as native data-frame prediction
# does. Non-NULL controls request zero effects and do not need random-only inputs.
.msr_unmarked_random_design <- function(formula, raw, data, re.form = NULL) {
  ns <- asNamespace("unmarked")
  if (!is.null(re.form) || !get("has_random", ns)(formula)) return(get("get_Z", ns)(~1, data))
  absent <- setdiff(all.vars(formula), c(names(data), "."))
  if (length(absent)) stop("New unmarked random-effect rows must explicitly contain formula variables: ",
    paste(absent, collapse = ", "), call. = FALSE)
  get("get_Z", ns)(formula, raw, data)
}

# A fresh closure changes helper lookup locally. No provider namespace, method
# registration, original formula, fitted S4 slot or likelihood closure is edited.
.msr_unmarked_local_method <- function(method, bindings) {
  # S4's standard wrapper can embed an actual .local FUNCTION constant.
  # Rebinding the outer wrapper alone leaves that function in the namespace.
  expr <- body(method)
  if (is.call(expr) && identical(expr[[1L]], as.name("{")) && length(expr) == 3L &&
      is.call(expr[[2L]]) && identical(expr[[2L]][[1L]], as.name("<-")) &&
      identical(expr[[2L]][[2L]], as.name(".local")) &&
      is.function(expr[[2L]][[3L]]) && is.call(expr[[3L]]) &&
      identical(expr[[3L]][[1L]], as.name(".local"))) method <- expr[[2L]][[3L]]
  env <- list2env(bindings, parent = environment(method))
  eval(call("function", formals(method), body(method)), env)
}

.msr_unmarked_prediction_function <- function(object, re.form = NULL) {
  native <- unmarked::predict
  if (!class(object)[1L] %in% c("unmarkedFitOccu", "unmarkedFitPCount")) return(native)
  ns <- asNamespace("unmarked")
  formulas <- object@formlist
  zero_context <- !is.null(re.form)
  default_re_form <- re.form
  if (!all(vapply(formulas[c("state", "det")], inherits, logical(1), "formula"))) return(native)
  for (formula in formulas[c("state", "det")]) {
    labels <- attr(stats::terms(formula), "term.labels")
    random <- vapply(lapply(labels, str2lang), .msr_marginal_random_call, logical(1))
    if (!zero_context && any(random) && .msr_needs_nested_metadata(stats::terms(stats::reformulate(
        labels[random], env = environment(formula)))))
      stop("Learned random-effect transformations need precomputed original training predictors; trained random-design metadata are unavailable.", call. = FALSE)
  }
  raw <- lapply(c("state", "det"), function(type) get("get_covariates", ns)(object, type))
  needed <- any(vapply(formulas[c("state", "det")], function(f)
    .msr_needs_nested_metadata(stats::terms(f)) || get("has_random", ns)(f), logical(1))) ||
    any(vapply(seq_along(raw), function(i) {
      variables <- .msr_marginal_variables(formulas[[c("state", "det")[i]]], offsets = TRUE)
      if ("." %in% variables) variables <- names(raw[[i]])
      columns <- intersect(names(raw[[i]]), variables)
      any(vapply(raw[[i]][columns], is.factor, logical(1)))
    }, logical(1)))
  if (!needed) return(native)
  bases <- stats::setNames(lapply(c("state", "det"), function(type)
    .msr_unmarked_prediction_basis(object, type, verify_random = !zero_context)), c("state", "det"))
  generic <- methods::getGeneric("predict_internal", where = ns)
  selected <- methods::selectMethod("predict_internal", class(object)[1L], fdef = generic)
  common <- methods::selectMethod("predict_internal", "unmarkedFit", fdef = generic)
  if (!identical(body(selected), body(common)) || !identical(formals(selected), formals(common)))
    stop("This unmarked prediction method needs its own trained-design adapter.", call. = FALSE)
  function(object, type, newdata = NULL, backTransform = TRUE, na.rm = TRUE,
      appendData = FALSE, level = 0.95, re.form = default_re_form, ...) {
    if (zero_context && is.null(re.form))
      stop("This zero-effect prediction context cannot supply conditional random effects; prepare a conditional context with its original random-design metadata.", call. = FALSE)
    # Nonregression nuisance components keep their provider's original route.
    if (!type %in% names(bases)) return(native(object, type, newdata = newdata,
      backTransform = backTransform, na.rm = na.rm, appendData = appendData,
      level = level, re.form = re.form, ...))
    basis <- bases[[type]]
    make_matrix <- function(formula, data, newdata, re.form = NULL) {
      out <- .msr_unmarked_fixed_design(basis, newdata)
      if (is.null(re.form) && get("has_random", ns)(formula))
        out$X <- cbind(out$X, .msr_unmarked_random_design(formula, data, newdata, re.form))
      out
    }
    frame_inputs <- function(object, type, newdata, na.rm, re.form) {
      # Run the native site/observation removal algorithm with tagged LOCAL
      # formulas, so identical state and detection syntax can use distinct bases.
      tagged <- formulas
      for (name in names(bases)) attr(tagged[[name]], ".msr_prediction_component") <- name
      fixed <- function(formula, covs) {
        name <- attr(formula, ".msr_prediction_component")
        .msr_unmarked_fixed_design(bases[[name]], covs)
      }
      design <- methods::selectMethod("getDesign", class(newdata)[1L],
        fdef = methods::getGeneric("getDesign", where = ns))
      ordinary <- methods::selectMethod("getDesign", "unmarkedFrame",
        fdef = methods::getGeneric("getDesign", where = ns))
      if (!identical(body(design), body(ordinary)) || !identical(formals(design), formals(ordinary)))
        stop("This unmarked frame layout needs its own trained-design adapter.", call. = FALSE)
      get_design <- .msr_unmarked_local_method(design, list(
        get_Z = function(formula, covs) {
          name <- attr(formula, ".msr_prediction_component")
          original <- raw[[match(name, c("state", "det"))]]
          .msr_unmarked_random_design(formulas[[name]], original, covs, re.form)
        },
        get_model_matrix = function(formula, covs) fixed(formula, covs)$X,
        get_offset = function(formula, covs) {
          out <- fixed(formula, covs)$offset
          if (is.null(out)) return(rep(0, nrow(covs)))
          out[is.na(out)] <- 0
          as.vector(out)
        }))
      inputs <- .msr_unmarked_local_method(get("predict_inputs_from_umf", ns),
        list(getDesign = function(umf, formulas, na.rm) get_design(umf, tagged, na.rm = na.rm)))
      inputs(object, type, newdata, na.rm, re.form)
    }
    predict <- .msr_unmarked_local_method(selected,
      list(make_mod_matrix = make_matrix, predict_inputs_from_umf = frame_inputs))
    predict(object, type, newdata = newdata, backTransform = backTransform,
      na.rm = na.rm, appendData = appendData, level = level, re.form = re.form, ...)
  }
}

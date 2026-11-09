# Compact matrix diagnostics describe the stored design, never fit admission.
.msr_gdr_design_summary <- function(x) {
  if (ncol(x) > 256L || length(x) > 2e6)
    return(list(rows = nrow(x), columns = ncol(x), status = "work_limit_exceeded",
      maximum_columns = 256L, maximum_positions = 2e6))
  rank <- qr(x, tol = 1e-7)$rank
  list(rows = nrow(x), columns = ncol(x), status = "evaluated", column_names = colnames(x), rank = rank,
    rank_tolerance = 1e-7, column_rank_deficient = rank < ncol(x),
    condition_number = if (rank < ncol(x)) Inf else kappa(x, exact = TRUE),
    conditioning_scale = "native_columns_without_additional_rescaling",
    matrix_sha256 = .msr_gdr_hash(list(dimensions = dim(x), columns = colnames(x), values = as.numeric(x))))
}

.msr_gdr_design_audit <- function(gd, components, formula, site_data, projection_data = NULL) {
  if (is.null(site_data)) site_data <- data.frame(row.names = seq_len(nrow(gd$X_lambda)))
  capture <- function(expr) {
    warnings <- character(); error <- NULL
    value <- tryCatch(withCallingHandlers(force(expr),
      warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
      error = function(e) { error <<- conditionMessage(e); NULL })
    list(value = value, warnings = warnings, error = error)
  }
  summary <- stats::setNames(lapply(components, function(name)
    .msr_gdr_design_summary(gd[[paste0("X_", name)]])), components)
  result <- list(components = summary,
    state = list(status = "not_reconstructed", native_matrix_verified = FALSE),
    projection = list(status = if (is.null(projection_data)) "not_requested" else "not_evaluated"),
    uncertainty = "none_read_only_design_audit", fitting_admission_established = FALSE,
    prediction_context_established = FALSE, numerical_gates_changed = FALSE)
  if (identical(summary$lambda$status, "work_limit_exceeded")) {
    result$state$status <- "work_limit_exceeded"
    return(result)
  }
  reconstruction <- capture({
    frame <- stats::model.frame(formula, data = site_data, na.action = stats::na.fail)
    # Reuse the same training-only learner machinery as CV. Native predvars
    # alone can leave a nested scale() or poly() call incompletely frozen.
    terms <- .msr_cv_prediction_terms(attr(frame, "terms"), site_data, qualify = TRUE)
    frame <- stats::model.frame(terms, data = site_data, na.action = stats::na.fail)
    native <- gd$X_lambda
    matrix <- stats::model.matrix(terms, data = frame, contrasts.arg = attr(native, "contrasts"))
    if (!identical(dim(matrix), dim(native)) ||
        !identical(colnames(matrix), colnames(native)) ||
        !identical(as.numeric(matrix), as.numeric(native)))
      stop("Reconstructed state columns do not exactly reproduce the stored native design.")
    list(frame = frame, terms = terms, matrix = matrix)
  })
  result$state$warnings <- reconstruction$warnings
  result$state$error <- reconstruction$error
  if (is.null(reconstruction$value)) return(result)
  retained <- reconstruction$value
  frame <- retained$frame; terms <- retained$terms
  variables <- as.list(attr(terms, "variables"))[-1L]
  prediction <- as.list(attr(terms, "predvars"))[-1L]
  splines <- which(vapply(frame, function(x) inherits(x, "ns") || inherits(x, "bs"), logical(1)))
  support <- function(call, data, boundary, training_range = NULL) {
    # Do not reevaluate arbitrary transformed expressions or learn preprocessing.
    # Native matrix/recipe verification still applies to those expressions.
    if (!is.call(call) || length(call) < 2L || !is.symbol(call[[2L]]) ||
        !as.character(call[[2L]]) %in% names(data))
      return(list(status = "nonliteral_input_not_recomputed"))
    x <- data[[as.character(call[[2L]])]]
    if (!is.numeric(x) || any(!is.finite(x))) return(list(status = "nonfinite_or_nonnumeric_input"))
    list(status = "bare_numeric_input", input = as.character(call[[2L]]), range = range(x),
      zeros = sum(x == 0), below_boundary = sum(x < boundary[1L]),
      above_boundary = sum(x > boundary[2L]),
      outside_training_range = if (is.null(training_range)) NULL else
        sum(x < training_range[1L] | x > training_range[2L]))
  }
  recipes <- lapply(splines, function(i) {
    x <- frame[[i]]
    boundary <- attr(x, "Boundary.knots")
    list(term = names(frame)[i], basis = if (inherits(x, "ns")) "ns" else "bs",
      requested_call = paste(deparse(variables[[i]]), collapse = " "),
      retained_prediction_call = paste(deparse(prediction[[i]]), collapse = " "),
      knots = attr(x, "knots"), Boundary.knots = boundary,
      intercept = attr(x, "intercept"), degree = attr(x, "degree"),
      basis_columns = ncol(x), training = support(variables[[i]], site_data, boundary))
  })
  names(recipes) <- names(frame)[splines]
  result$state$status <- "exact_native_matrix_verified"
  result$state$native_matrix_verified <- TRUE
  result$state$formula <- paste(deparse(formula), collapse = " ")
  result$state$splines <- recipes
  if (!is.null(projection_data)) {
    projected <- capture({
      if (as.double(nrow(projection_data)) * ncol(retained$matrix) > 2e6)
        stop("Scenario state-design audit exceeds two million matrix positions.")
      xlevels <- lapply(frame[vapply(frame, function(x) is.factor(x) || is.character(x), logical(1))],
        function(x) if (is.factor(x)) levels(x) else levels(factor(x)))
      safe <- function(x) {
        if (!is.call(x)) return(TRUE)
        if (is.null(.msr_fixed_function_name(x[[1L]], environment(terms)))) {
          head <- x[[1L]]
          name <- if (is.symbol(head)) as.character(head) else if (is.call(head) &&
            length(head) == 3L && identical(head[[1L]], as.name("::")) &&
            identical(head[[2L]], as.name("base"))) as.character(head[[3L]]) else ""
          # A modeled factor call has retained native xlevels and contrasts.
          # A nested factor with no surviving levels is not interchangeable.
          if (!name %in% c("factor", "ordered", "as.factor") ||
              !paste(deparse(x), collapse = " ") %in% names(xlevels) ||
              !.msr_standard_prediction_function(get0(name, environment(terms),
                mode = "function", inherits = TRUE), name, "base")) return(FALSE)
        }
        parts <- as.list(x)
        all(vapply(seq_along(parts)[-1L], function(j)
          identical(unname(parts[j]), unname(alist(x = ))) || safe(parts[[j]]), logical(1)))
      }
      if (!safe(attr(terms, "predvars")))
        stop("Scenario audits require qualified fixed training transformations; precompute other transformations from training data.")
      # Training predvars retain actual native knots, including adjusted knots.
      mf <- stats::model.frame(terms, data = projection_data, na.action = stats::na.fail, xlev = xlevels)
      if (nrow(mf) != nrow(projection_data) || !identical(rownames(mf), rownames(projection_data)))
        stop("Scenario state rows differ from the supplied scenario identities.")
      matrix <- stats::model.matrix(terms, data = mf, contrasts.arg = attr(retained$matrix, "contrasts"))
      if (!identical(colnames(matrix), colnames(retained$matrix)))
        stop("Projected state columns differ from the verified native training design.")
      list(matrix = .msr_gdr_design_summary(matrix), splines = stats::setNames(lapply(seq_along(splines), function(j)
        support(variables[[splines[j]]], projection_data, recipes[[j]]$Boundary.knots,
          recipes[[j]]$training$range)), names(recipes)))
    })
    result$projection <- list(status = if (is.null(projected$value)) "failed" else "retained_state_basis_evaluated",
      warnings = projected$warnings, error = projected$error, summary = projected$value)
  }
  result
}

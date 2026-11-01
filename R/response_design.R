# Native multilinear response-design contract.
# Scalar kernel variables may appear in any native-expanded crossing/nesting
# term. Fixed numeric/matrix/factor blocks and contrast coding stay native.
.msr_compile_response_design <- function(terms, frame, covariates,
                                         columns = NULL, contrasts = NULL) {
  if (!inherits(terms, "terms") || !is.data.frame(frame) || !nrow(frame) ||
      !is.character(covariates) || !length(covariates) || anyNA(covariates) ||
      anyDuplicated(covariates)) stop("Invalid response-design inputs.")
  variables <- as.list(attr(terms, "variables"))[-1L]
  dependent <- lapply(variables, function(ex) intersect(all.vars(ex), covariates))
  factors <- attr(terms, "factors")
  active <- rowSums(factors != 0L) > 0L
  active[attr(terms, "offset")] <- TRUE
  if (attr(terms, "response") > 0L &&
      length(dependent[[attr(terms, "response")]]))
    stop("Optimized predictors cannot be the response variable.")
  dependent[!active] <- rep(list(character()), sum(!active))
  transformed <- vapply(seq_along(dependent), function(i) length(dependent[[i]]) > 0L &&
    !(length(dependent[[i]]) == 1L && identical(variables[[i]], as.name(dependent[[i]]))), logical(1))
  if (any(transformed))
    return(.msr_compile_expression_design(terms, frame, covariates, columns, contrasts))
  for (i in which(lengths(dependent) > 0L)) {
    if (length(dependent[[i]]) != 1L ||
        !identical(variables[[i]], as.name(dependent[[i]])))
      stop("Dynamic transformed expressions require their own derivative contract.")
  }
  rows <- vapply(covariates, function(nm) {
    i <- which(vapply(variables, identical, logical(1), y = as.name(nm)))
    if (length(i) != 1L) stop("Cannot identify each numeric kernel variable.")
    i
  }, integer(1))
  if (any(rows == attr(terms, "response")))
    stop("Optimized predictors cannot be the response variable.")
  factors <- attr(terms, "factors")
  if (!is.matrix(factors) || nrow(factors) != length(variables) ||
      ncol(frame) < length(variables)) stop("Inconsistent native formula metadata.")
  if (any(!vapply(frame[rows], function(x) is.numeric(x) && is.null(dim(x)), logical(1))))
    stop("Optimized kernel variables must be scalar numeric columns.")
  X <- stats::model.matrix(terms, frame, contrasts.arg = contrasts)
  if (is.null(columns)) columns <- colnames(X)
  if (!is.character(columns) || anyNA(columns) || anyDuplicated(columns) ||
      !all(columns %in% colnames(X)) || anyDuplicated(colnames(X)) ||
      any(!is.finite(X))) stop("Invalid native response-design columns.")
  # Replacing numeric kernel leaves by one extracts native fixed multipliers.
  # model.matrix retains full dummy coding for interaction-only factor terms.
  unit_frame <- frame
  for (i in rows) unit_frame[[i]] <- rep(1, nrow(frame))
  unit <- stats::model.matrix(terms, unit_frame, contrasts.arg = contrasts)
  if (!identical(colnames(unit), colnames(X))) stop("Native column layout changed.")
  position <- match(columns, colnames(unit))
  assignment <- attr(unit, "assign")[position]
  dependence <- lapply(assignment, function(j) {
    if (j == 0L) integer() else which(factors[rows, j] != 0L)
  })
  keys <- vapply(dependence, paste, character(1), collapse = ",")
  groups <- lapply(unique(keys), function(key) {
    selected <- which(keys == key)
    list(columns = selected, kernels = dependence[[selected[1L]]])
  })
  structure(list(multipliers = unit[, position, drop = FALSE],
    dependence = dependence, groups = groups, covariates = covariates, columns = columns,
    kernel_frame_columns = stats::setNames(names(frame)[rows], covariates),
    n = nrow(frame), p = length(columns), q = length(covariates),
    assignment = assignment, native_terms = terms, contrasts = contrasts),
    class = "msr_response_design")
}

.msr_response_design_matrix <- function(design, kernel) {
  if (!inherits(design, "msr_response_design") || !is.matrix(kernel) ||
      !is.numeric(kernel) || !identical(dim(kernel), c(design$n, design$q)) ||
      any(!is.finite(kernel)) || (!is.null(colnames(kernel)) &&
        !identical(colnames(kernel), design$covariates)))
    stop("Invalid ordered numeric kernel matrix for response design.")
  values <- if (is.null(design$expression_programs)) kernel else
    .msr_response_expression_values(design, kernel)
  out <- design$multipliers
  for (group in design$groups) {
    if (!length(group$kernels)) next
    value <- rep(1, design$n)
    for (k in group$kernels) value <- value *
      if (is.list(values)) values[[k]] else values[, k]
    out[, group$columns] <- out[, group$columns, drop = FALSE] * value
  }
  if (any(!is.finite(out))) stop("Dynamic response design is nonfinite.")
  out
}

# d eta / d kernel value, not d eta / d scale. Compose with the existing
# kernel scale/shape derivatives. Multiplication never divides by a predictor,
# so zero-valued kernels retain correct cross-term derivatives.
.msr_response_kernel_sensitivity <- function(design, kernel, beta) {
  if (length(beta) != design$p || any(!is.finite(beta))) stop("Invalid coefficients.")
  if (!is.null(design$expression_programs))
    return(.msr_response_expression_sensitivity(design, kernel, beta))
  out <- matrix(0, design$n, design$q,
    dimnames = list(rownames(design$multipliers), design$covariates))
  for (group in design$groups) {
    involved <- group$kernels
    if (!length(involved)) next
    coefficient <- drop(design$multipliers[, group$columns, drop = FALSE] %*% beta[group$columns])
    for (k in involved) {
      value <- coefficient
      for (other in setdiff(involved, k)) value <- value * kernel[, other]
      out[, k] <- out[, k] + value
    }
  }
  out
}

# AD-safe counterpart: no numeric matrix receives an AD value, and no native
# model.matrix call evaluates differentiated variables inside a tape.
.msr_response_design_eta <- function(design, kernel, beta) {
  if (!is.null(design$expression_programs)) {
    .msr_response_expression_ram(design)
    kernel <- .msr_response_expression_values(design, kernel)
  }
  eta <- rep(0, design$n)
  for (group in design$groups) {
    value <- drop(design$multipliers[, group$columns, drop = FALSE] %*% beta[group$columns])
    for (k in group$kernels) value <- value *
      if (is.list(kernel)) kernel[[k]] else kernel[, k]
    eta <- eta + value
  }
  eta
}

# Use native active terms, not removed symbols left in the formula expression.
.msr_active_formula_variables <- function(terms) {
  variables <- as.list(attr(terms, "variables"))[-1L]
  factors <- attr(terms, "factors")
  used <- rep(FALSE, length(variables))
  if (length(factors)) {
    if (!is.matrix(factors) || nrow(factors) != length(variables))
      stop("Cannot recover active native formula variables.")
    used <- rowSums(factors != 0) > 0L
  }
  used[attr(terms, "offset")] <- TRUE
  if (attr(terms, "response") > 0L) used[attr(terms, "response")] <- FALSE
  unique(unlist(lapply(variables[used], all.vars), use.names = FALSE))
}

.msr_validate_kernel_expressions <- function(terms, covariates) {
  variables <- as.list(attr(terms, "variables"))[-1L]
  response <- attr(terms, "response")
  if (response > 0L && length(intersect(covariates, all.vars(variables[[response]]))))
    stop("Optimized predictors cannot occur in the response expression.", call. = FALSE)
  indices <- setdiff(seq_along(variables), c(response, attr(terms, "offset")))
  factors <- attr(terms, "factors")
  indices <- intersect(indices, which(rowSums(factors != 0L) > 0L))
  for (i in indices) {
    involved <- intersect(all.vars(variables[[i]]), covariates)
    if (length(involved) && !(length(involved) == 1L &&
        identical(variables[[i]], as.name(involved))))
      .msr_kernel_expression(variables[[i]], covariates, attr(terms, ".Environment"))
  }
}

.msr_response_design_required <- function(terms, covariates, columns) {
  if (!all(covariates %in% columns)) return(TRUE)
  variables <- as.list(attr(terms, "variables"))[-1L]
  active <- rowSums(attr(terms, "factors") != 0L) > 0L
  for (i in which(active)) {
    involved <- intersect(all.vars(variables[[i]]), covariates)
    if (length(involved) && !(length(involved) == 1L &&
        identical(variables[[i]], as.name(involved)))) return(TRUE)
  }
  variables <- as.list(attr(terms, "variables"))[-1L]
  factors <- attr(terms, "factors")
  for (nm in covariates) {
    i <- which(vapply(variables, identical, logical(1), y = as.name(nm)))
    if (length(i) != 1L) stop("Cannot identify a modeled kernel variable.")
    active <- which(factors[i, ] != 0L)
    if (length(active) != 1L || any(colSums(factors[, active, drop = FALSE] != 0L) != 1L))
      return(TRUE)
  }
  FALSE
}

# Subset evaluated variables by expression identity, retaining learned metadata.
# In particular, never expand a conditional stratum into thousands of dummies.
.msr_response_frame_subset <- function(native_terms, frame, target_terms) {
  original <- as.list(attr(native_terms, "variables"))[-1L]
  selected <- as.list(attr(target_terms, "variables"))[-1L]
  indices <- vapply(selected, function(ex) {
    i <- which(vapply(original, identical, logical(1), y = ex))
    if (length(i) != 1L) stop("Cannot match native evaluated formula variables.")
    i
  }, integer(1))
  out <- frame[, indices, drop = FALSE]
  attr(out, "terms") <- target_terms
  out
}

# Prepared targets retain the legacy fixed-then-kernel coefficient order unless
# a native expanded response design is explicitly present.
.msr_scale_response_matrix <- function(spec, kernel) {
  if (is.null(spec$response_design))
    return(cbind(spec$X_fixed, do.call(cbind, kernel)))
  .msr_response_design_matrix(spec$response_design, do.call(cbind, kernel))
}

# Subset only observation-dependent multipliers. Term/contrast/column meanings
# remain frozen when independent groups are integrated separately.
.msr_response_design_rows <- function(design, rows) {
  design$multipliers <- design$multipliers[rows, , drop = FALSE]
  design$n <- length(rows)
  if (!is.null(design$reference_kernel))
    design$reference_kernel <- design$reference_kernel[rows, , drop = FALSE]
  design
}

# Remove only zero-factor variable rows from native terms. Rebuilding from
# labels could change factor coding; retaining the row order preserves it.
.msr_response_active_terms <- function(terms) {
  factors <- attr(terms, "factors")
  variables <- as.list(attr(terms, "variables"))
  if (!is.matrix(factors) || nrow(factors) != length(variables) - 1L)
    stop("Invalid native active-term metadata.")
  keep <- rowSums(factors != 0L) > 0L
  keep[c(attr(terms, "response"), attr(terms, "offset"))] <- TRUE
  if (all(keep)) return(terms)
  positions <- which(keep)
  for (field in c("variables", "predvars")) {
    value <- attr(terms, field)
    if (!is.null(value)) {
      value <- as.list(value)
      if (length(value) != length(variables)) stop("Native prediction metadata changed dimensions.")
      attr(terms, field) <- as.call(c(value[1L], value[-1L][keep]))
    }
  }
  attr(terms, "factors") <- factors[keep, , drop = FALSE]
  if (length(attr(terms, "offset")))
    attr(terms, "offset") <- match(attr(terms, "offset"), positions)
  if (attr(terms, "response") > 0L)
    attr(terms, "response") <- match(attr(terms, "response"), positions)
  classes <- attr(terms, "dataClasses")
  if (!is.null(classes)) {
    if (length(classes) != length(keep)) stop("Native data classes changed dimensions.")
    attr(terms, "dataClasses") <- classes[keep]
  }
  terms
}

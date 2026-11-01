# Restricted pointwise kernel expressions. Native formula expansion and learned
# fixed bases stay in model.matrix; no caller data or function is captured here.
.msr_kernel_expression <- function(expression, covariates, context) {
  unary <- c("I", "(", "exp", "log", "sqrt", "sin", "cos", "tan")
  binary <- c("+", "-", "*", "/", "^")
  nodes <- 0L
  normalize <- function(ex) {
    nodes <<- nodes + 1L
    if (is.numeric(ex) && !is.complex(ex) && length(ex) == 1L && is.finite(ex)) return(ex)
    if (is.symbol(ex)) {
      if (!as.character(ex) %in% covariates)
        stop("Dynamic expressions currently require kernel inputs and literal constants only.", call. = FALSE)
      return(ex)
    }
    if (!is.call(ex))
      stop("Unsupported dynamic expression call.", call. = FALSE)
    head <- ex[[1L]]
    qualified <- is.call(head) && identical(head[[1L]], as.name("::")) &&
      length(head) == 3L && identical(head[[2L]], as.name("base")) && is.symbol(head[[3L]])
    if (qualified) head <- head[[3L]]
    if (!is.symbol(head)) stop("Unsupported dynamic expression call.", call. = FALSE)
    name <- as.character(head)
    if (!name %in% c(unary, binary))
      stop("Unsupported dynamic function: ", name, ".", call. = FALSE)
    arguments <- as.list(ex)[-1L]
    if ((!is.null(names(arguments)) && any(nzchar(names(arguments)))) ||
        (name %in% unary && length(arguments) != 1L) ||
        (name %in% c("*", "/", "^") && length(arguments) != 2L) ||
        (name %in% c("+", "-") && !length(arguments) %in% 1:2))
      stop("Unsupported dynamic expression argument layout.", call. = FALSE)
    native <- if (qualified) get(name, envir = baseenv(), inherits = FALSE) else
      tryCatch(get(name, envir = context, mode = "function", inherits = TRUE),
        error = function(e) NULL)
    if (!identical(native, get(name, envir = baseenv(), inherits = FALSE)))
      stop("Dynamic expressions require the original base function binding: ", name, ".", call. = FALSE)
    arguments <- lapply(arguments, normalize)
    if (name %in% c("I", "(")) return(arguments[[1L]])
    as.call(c(list(as.name(name)), arguments))
  }
  normalized <- normalize(expression)
  program <- function(ex) {
    rewrite <- function(value) {
      if (is.symbol(value)) return(as.call(list(as.name("[["), as.name("kernel"),
        match(as.character(value), covariates))))
      if (is.call(value)) return(as.call(c(list(value[[1L]]), lapply(as.list(value)[-1L], rewrite))))
      value
    }
    # The closure has baseenv() only. No training frame, caller env, or tape.
    eval(as.call(list(as.name("function"), as.pairlist(alist(kernel = )), rewrite(ex))), baseenv())
  }
  derivatives <- lapply(covariates, function(name) stats::D(normalized, name))
  zero <- vapply(derivatives, function(ex) is.numeric(ex) && length(ex) == 1L && ex == 0, logical(1))
  list(value = program(normalized), derivative = lapply(derivatives, program),
       zero = zero, nodes = nodes, expression = normalized)
}

.msr_response_expression_values <- function(design, kernel) {
  if ((is.list(kernel) && (length(kernel) != design$q ||
       any(vapply(kernel, length, numeric(1)) != design$n))) ||
      (!is.list(kernel) && !identical(dim(kernel), c(design$n, design$q))))
    stop("Dynamic expressions require the ordered kernel rows and columns.")
  input <- if (is.list(kernel)) kernel else
    lapply(seq_len(design$q), function(j) kernel[, j])
  lapply(design$expression_programs, function(program) {
    value <- program$value(input)
    if (length(value) == 1L && design$n != 1L) value <- rep(value, design$n)
    if (length(value) != design$n) stop("Dynamic expression changed the response row count.")
    value
  })
}

.msr_response_expression_sensitivity <- function(design, kernel, beta) {
  values <- .msr_response_expression_values(design, kernel)
  input <- lapply(seq_len(design$q), function(j) kernel[, j])
  out <- matrix(0, design$n, design$q,
    dimnames = list(rownames(design$multipliers), design$covariates))
  for (group in design$groups) {
    involved <- group$kernels
    if (!length(involved)) next
    coefficient <- drop(design$multipliers[, group$columns, drop = FALSE] %*% beta[group$columns])
    for (leaf in involved) {
      value <- coefficient
      for (other in setdiff(involved, leaf)) value <- value * values[[other]]
      for (j in which(!design$expression_programs[[leaf]]$zero)) {
        derivative <- design$expression_programs[[leaf]]$derivative[[j]](input)
        if (any(!is.finite(derivative))) stop("Dynamic expression derivative is nonfinite.")
        out[, j] <- out[, j] + value * derivative
      }
    }
  }
  if (any(!is.finite(out))) stop("Dynamic response sensitivity is nonfinite.")
  out
}

# Planning allowance for the new expression tape, not a measured peak bound.
.msr_response_expression_ram <- function(design) {
  total <- .msr_total_ram_bytes()
  if (!is.finite(total) || total <= 0)
    stop("Cannot detect RAM for dynamic response differentiation.", call. = FALSE)
  estimate <- 1024^3 + 4 * as.numeric(utils::object.size(design)) +
    512 * as.double(design$n) * design$expression_nodes
  if (estimate > .75 * total)
    stop("Dynamic response exceeds the conservative differentiation RAM budget.", call. = FALSE)
  invisible(estimate)
}

# Refresh every evaluated expression without guessing raw kernels from powers.
.msr_response_update_frame <- function(design, frame, kernel) {
  values <- .msr_response_expression_values(design, kernel)
  for (j in seq_along(values)) {
    name <- design$expression_frame_columns[j]
    if (!name %in% names(frame)) stop("Missing evaluated dynamic expression column.")
    value <- values[[j]]
    if (any(!is.finite(value))) stop("Dynamic expression has nonfinite values.")
    if (inherits(frame[[name]], "AsIs")) value <- I(value)
    frame[[name]] <- value
  }
  for (j in seq_along(design$covariates)) {
    name <- design$kernel_frame_columns[j]
    if (!is.na(name) && name %in% names(frame)) frame[[name]] <- kernel[, j]
  }
  frame
}

.msr_compile_expression_design <- function(terms, frame, covariates, columns, contrasts) {
  variables <- as.list(attr(terms, "variables"))[-1L]
  dependent <- lapply(variables, function(ex) intersect(all.vars(ex), covariates))
  factors <- attr(terms, "factors")
  rows <- which(lengths(dependent) > 0L & rowSums(factors != 0L) > 0L)
  if (any(which(lengths(dependent) > 0L) %in%
          c(attr(terms, "response"), attr(terms, "offset"))))
    stop("Optimized predictors cannot occur in response expressions or offsets.")
  if (any(!vapply(frame[rows], function(x) is.numeric(x) && is.null(dim(x)), logical(1))))
    stop("Dynamic expressions must be scalar numeric columns.")
  programs <- lapply(variables[rows], .msr_kernel_expression,
    covariates = covariates, context = attr(terms, ".Environment"))
  X <- stats::model.matrix(terms, frame, contrasts.arg = contrasts)
  if (is.null(columns)) columns <- colnames(X)
  if (any(!is.finite(X)) || anyDuplicated(colnames(X)) ||
      anyNA(columns) || anyDuplicated(columns) || !all(columns %in% colnames(X)))
    stop("Invalid native response-design columns.")
  unit_frame <- frame
  for (j in rows) unit_frame[[j]] <- rep(1, nrow(frame))
  unit <- stats::model.matrix(terms, unit_frame, contrasts.arg = contrasts)
  if (!identical(colnames(unit), colnames(X))) stop("Native column layout changed.")
  position <- match(columns, colnames(unit))
  assignment <- attr(unit, "assign")[position]
  dependence <- lapply(assignment, function(j) if (j == 0L) integer() else which(factors[rows, j] != 0L))
  keys <- vapply(dependence, paste, character(1), collapse = ",")
  groups <- lapply(unique(keys), function(key) {
    selected <- which(keys == key)
    list(columns = selected, kernels = dependence[[selected[1L]]])
  })
  bare <- vapply(covariates, function(name) {
    j <- which(vapply(variables, identical, logical(1), y = as.name(name)))
    if (length(j) == 1L) names(frame)[j] else NA_character_
  }, character(1))
  structure(list(multipliers = unit[, position, drop = FALSE], dependence = dependence,
    groups = groups, covariates = covariates, columns = columns,
    kernel_frame_columns = bare, n = nrow(frame), p = length(columns), q = length(covariates),
    assignment = assignment, native_terms = terms, contrasts = contrasts,
    expression_programs = programs, expression_frame_columns = names(frame)[rows],
    expression_nodes = sum(vapply(programs, function(x) x$nodes, integer(1)))),
    class = "msr_response_design")
}

# Fixed-input cubic bases. Native splineDesign/QR run once, before AD;
# each candidate evaluates exact local polynomials, not fitted interpolation.
.msr_spline_call <- function(ex) {
  if (!is.call(ex)) return(NULL)
  head <- ex[[1L]]
  if (is.symbol(head) && as.character(head) %in% c("ns", "bs"))
    return(as.character(head))
  if (is.call(head) && length(head) == 3L &&
      identical(head[[1L]], as.name("::")) &&
      identical(head[[2L]], as.name("splines")) &&
      is.symbol(head[[3L]]) && as.character(head[[3L]]) %in% c("ns", "bs"))
    return(as.character(head[[3L]]))
  NULL
}

.msr_fixed_spline_program <- function(expression, covariates, context) {
  type <- .msr_spline_call(expression)
  native <- getExportedValue("splines", type)
  if (is.symbol(expression[[1L]]) && !identical(native,
      get0(type, envir = context, mode = "function", inherits = TRUE)))
    stop("Dynamic splines require the original splines function binding.", call. = FALSE)
  matched <- match.call(native, expression)
  # Evaluate only literal vectors. Caller variables/functions must not choose a
  # candidate-specific recipe or capture a frame in a serialized design.
  literal <- function(ex) {
    if (is.null(ex) || (is.atomic(ex) && !is.object(ex) &&
        (is.numeric(ex) || is.logical(ex)))) return(ex)
    if (is.call(ex) && is.symbol(ex[[1L]])) {
      name <- as.character(ex[[1L]])
      args <- as.list(ex)[-1L]
      if (name == "c" || (name %in% c("+", "-") && length(args) == 1L)) {
        if (!identical(get0(name, context, mode = "function", inherits = TRUE),
                       get(name, baseenv(), inherits = FALSE)))
          stop("Dynamic spline constants require original base functions.", call. = FALSE)
        return(do.call(get(name, baseenv()), lapply(args, literal)))
      }
    }
    stop("Dynamic spline recipes require literal fixed knots and boundaries; use the refit backend for adaptive recipes.", call. = FALSE)
  }
  if (!all(c("knots", "Boundary.knots") %in% names(matched)) ||
      (!is.null(matched$df) && !is.null(literal(matched$df))))
    stop("Dynamic splines require explicit fixed knots and Boundary.knots, without adaptive df; use the refit backend otherwise.", call. = FALSE)
  knots <- literal(matched$knots); boundary <- literal(matched$Boundary.knots)
  intercept <- if (is.null(matched$intercept)) FALSE else literal(matched$intercept)
  degree <- if (is.null(matched$degree)) 3L else literal(matched$degree)
  if (!is.null(matched$warn.outside)) validate_scalar_logical(literal(matched$warn.outside), "warn.outside")
  validate_scalar_logical(intercept, "intercept")
  if (!is.numeric(boundary) || length(boundary) != 2L || any(!is.finite(boundary)) ||
      boundary[1L] >= boundary[2L] ||
      (!is.null(knots) && (!is.numeric(knots) || any(!is.finite(knots)))) ||
      anyDuplicated(knots) || any(knots <= boundary[1L] | knots >= boundary[2L]))
    stop("Dynamic spline knots must be finite, distinct and strictly inside two increasing fixed boundary knots.", call. = FALSE)
  if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) || degree != 3)
    stop("Dynamic spline derivatives currently require cubic degree 3; use the refit backend otherwise.", call. = FALSE)
  knots <- sort(as.numeric(knots)); boundary <- as.numeric(boundary)
  input <- .msr_kernel_expression(matched$x, covariates, context)
  breaks <- c(boundary[1L], knots, boundary[2L])
  full_knots <- sort(c(rep(boundary, 4L), knots))
  transform <- function(value) if (!intercept) value[, -1L, drop = FALSE] else value
  if (type == "ns") {
    constraint <- transform(splines::splineDesign(full_knots, boundary,
      ord = 4L, derivs = c(2L, 2L)))
    constraint_qr <- qr(t(constraint))
    transform <- local({
      with_intercept <- intercept; Q <- constraint_qr
      function(value) {
        if (!with_intercept) value <- value[, -1L, drop = FALSE]
        t(qr.qty(Q, t(value)))[, -c(1L, 2L), drop = FALSE]
      }
    })
  }
  # Tail pivots match native bs extrapolation; ns tails use boundary value and
  # slope. Normalized local coordinates keep small knot intervals well scaled.
  centers <- c(if (type == "ns") boundary[1L] else boundary[1L] + diff(breaks)[1L]/4,
    utils::head(breaks, -1L) + diff(breaks)/2,
    if (type == "ns") boundary[2L] else boundary[2L] - utils::tail(diff(breaks), 1L)/4)
  scales <- c(diff(breaks)[1L], diff(breaks), utils::tail(diff(breaks), 1L))
  pieces <- lapply(seq_along(centers), function(j) {
    value <- transform(splines::splineDesign(full_knots, rep(centers[j], 4L),
      ord = 4L, derivs = 0:3))
    if (type == "ns" && j %in% c(1L, length(centers))) value[3:4, ] <- 0
    value * (scales[j]^(0:3)/factorial(0:3))
  })
  width <- ncol(pieces[[1L]])
  coefficients <- lapply(1:4, function(d) do.call(rbind, lapply(pieces, function(p) p[d, ])))
  recipe <- list(type = type, degree = 3L, knots = knots, Boundary.knots = boundary,
    intercept = intercept, width = width, breaks = breaks, centers = centers,
    scales = scales, coefficients = coefficients, knot_policy = "fixed_modeled_input",
    input_expression = input$expression)
  # A namespace-parent environment captures only numeric recipe and pure
  # scalar programs, never the caller frame, training data or an AD objective.
  value_env <- list2env(list(r = recipe, f = input$value),
    parent = environment(.msr_spline_columns))
  value <- eval(quote(function(kernel) .msr_spline_columns(r, f(kernel), 0L)), value_env)
  derivative <- lapply(input$derivative, function(df) {
    env <- list2env(list(r = recipe, f = input$value, d = df),
      parent = environment(.msr_spline_columns))
    eval(quote(function(kernel) {
      slope <- .msr_spline_columns(r, f(kernel), 1L); chain <- d(kernel)
      lapply(slope, function(x) x * chain)
    }), env)
  })
  list(value = value, derivative = derivative, zero = input$zero, width = width,
    nodes = as.integer(input$nodes + 16L * width * length(centers)),
    expression = expression, basis_recipe = recipe)
}

.msr_spline_columns <- function(recipe, x, derivative = 0L) {
  ad <- inherits(x, "advector")
  index <- if (ad) RTMB::findInterval(x, RTMB::AD(recipe$breaks)) + 1L else
    findInterval(x, recipe$breaks) + 1L
  at <- function(value) if (ad) RTMB::AD(value)[index] else value[index]
  u <- (x - at(recipe$centers))/at(recipe$scales)
  lapply(seq_len(recipe$width), function(j) {
    coefficient <- lapply(recipe$coefficients, function(value) at(value[, j]))
    if (derivative == 0L)
      ((coefficient[[4L]] * u + coefficient[[3L]]) * u + coefficient[[2L]]) * u + coefficient[[1L]]
    else ((3 * coefficient[[4L]] * u + 2 * coefficient[[3L]]) * u + coefficient[[2L]])/at(recipe$scales)
  })
}

# Retain matrix programs once; leaves identify native basis columns rather than
# inventing new kernel sources or independently optimized spatial scales.
.msr_compile_spline_design <- function(terms, frame, covariates, columns, contrasts) {
  variables <- as.list(attr(terms, "variables"))[-1L]
  factors <- attr(terms, "factors")
  dependent <- lengths(lapply(variables, function(ex) intersect(all.vars(ex), covariates))) > 0L
  if (any(which(dependent) %in% c(attr(terms, "response"), attr(terms, "offset"))))
    stop("Optimized predictors cannot occur in response expressions or offsets.")
  rows <- which(dependent & rowSums(factors != 0L) > 0L)
  programs <- lapply(variables[rows], .msr_kernel_expression,
    covariates = covariates, context = attr(terms, ".Environment"))
  widths <- vapply(programs, function(p) if (is.null(p$width)) 1L else p$width, integer(1))
  for (j in seq_along(rows)) {
    value <- frame[[rows[j]]]; recipe <- programs[[j]]$basis_recipe
    if (is.null(recipe)) {
      if (!is.numeric(value) || !is.null(dim(value))) stop("Dynamic expressions must be scalar numeric columns.")
    } else if (!inherits(value, recipe$type) || !is.matrix(value) ||
        ncol(value) != widths[j] ||
        !identical(as.numeric(attr(value, "knots")), recipe$knots) ||
        !identical(as.numeric(attr(value, "Boundary.knots")), recipe$Boundary.knots) ||
        !identical(attr(value, "intercept"), recipe$intercept))
      stop("Native dynamic spline frame does not match its fixed basis recipe.")
  }
  X <- stats::model.matrix(terms, frame, contrasts.arg = contrasts)
  if (is.null(columns)) columns <- colnames(X)
  if (!is.character(columns) || any(!is.finite(X)) || anyDuplicated(colnames(X)) ||
      anyNA(columns) || anyDuplicated(columns) || !all(columns %in% colnames(X)))
    stop("Invalid native response-design columns.")
  unit_frame <- frame
  for (j in rows) unit_frame[[j]][] <- 1
  unit <- stats::model.matrix(terms, unit_frame, contrasts.arg = contrasts)
  if (!identical(colnames(unit), colnames(X))) stop("Native column layout changed.")
  position <- match(columns, colnames(X)); assignment <- attr(unit, "assign")[position]
  dependence <- lapply(assignment, function(j) if (j == 0L) integer() else which(factors[rows, j] != 0L))
  offsets <- c(0L, cumsum(widths))
  spline_maps <- vector("list", length(rows))
  for (j in seq_along(rows)) {
    selected <- which(vapply(dependence, function(d) j %in% d, logical(1)))
    if (!length(selected)) next
    # Native one-hot probes recover matrix-column multiplication even with
    # arbitrary contrast coding/interactions. Names are never parsed heuristically.
    map <- if (is.null(programs[[j]]$basis_recipe)) rep(1L, length(selected)) else {
      incidence <- vapply(seq_len(widths[j]), function(k) {
        probe <- unit_frame; probe[[rows[j]]][] <- 0; probe[[rows[j]]][, k] <- 1
        M <- stats::model.matrix(terms, probe, contrasts.arg = contrasts)
        if (!identical(colnames(M), colnames(X))) stop("Native spline probe changed column layout.")
        colSums(abs(M[, position[selected], drop = FALSE])) > 0
      }, logical(length(selected)))
      incidence <- matrix(incidence, nrow = length(selected))
      if (any(rowSums(incidence) == 0L)) {
        # Assessment batches can legitimately omit a fitted factor level. Probe
        # native coding on synthetic level combinations, not assessment habitat
        # values or responses; this does not refit or change any coefficients.
        missing <- which(rowSums(incidence) == 0L)
        for (term in unique(assignment[selected[missing]])) {
          affected <- missing[assignment[selected[missing]] == term]
          probe <- .msr_spline_mapping_frame(unit_frame, factors, term, ncol(X))
          extra <- vapply(seq_len(widths[j]), function(k) {
            probe[[rows[j]]][] <- 0; probe[[rows[j]]][, k] <- 1
            M <- stats::model.matrix(terms, probe, contrasts.arg = contrasts)
            if (!identical(colnames(M), colnames(X))) stop("Native spline probe changed column layout.")
            colSums(abs(M[, position[selected[affected]], drop = FALSE])) > 0
          }, logical(length(affected)))
          incidence[affected, ] <- matrix(extra, nrow = length(affected))
        }
      }
      if (any(rowSums(incidence) != 1L))
        stop("Cannot identify a native spline column under the retained contrast coding.")
      max.col(incidence, ties.method = "first")
    }
    spline_maps[[j]] <- stats::setNames(as.integer(offsets[j] + map), as.character(selected))
  }
  # Replace program IDs in a separate pass: offsets may otherwise collide with
  # yet-unmapped program IDs when two bases occur in the same native term.
  dependence <- lapply(seq_along(assignment), function(k) {
    term <- assignment[k]
    if (term == 0L) return(integer())
    active <- which(factors[rows, term] != 0L)
    vapply(active, function(j) {
      unname(spline_maps[[j]][as.character(k)])
    }, integer(1), USE.NAMES = FALSE)
  })
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
    groups = groups, covariates = covariates, columns = columns, kernel_frame_columns = bare,
    n = nrow(frame), p = length(columns), q = length(covariates), assignment = assignment,
    native_terms = terms, contrasts = contrasts, expression_programs = programs,
    expression_frame_columns = names(frame)[rows], expression_widths = widths,
    expression_leaf_program = rep(seq_along(widths), widths),
    expression_leaf_column = sequence(widths),
    expression_nodes = sum(vapply(programs, function(p) p$nodes, integer(1)))),
    class = "msr_response_design")
}

# Bounded compile-time probes for empty assessment factor levels. Native
# model.matrix remains the authority on all column/contrast combinations.
.msr_spline_mapping_frame <- function(frame, factors, term, ncolumns) {
  active <- which(factors[, term] != 0L)
  categorical <- active[vapply(frame[active], function(x)
    is.factor(x) || is.character(x) || is.logical(x), logical(1))]
  levels <- lapply(frame[categorical], function(x) if (is.logical(x)) c(FALSE, TRUE) else
    if (is.factor(x)) levels(x) else levels(factor(x)))
  sizes <- lengths(levels); n <- prod(sizes)
  total <- .msr_total_ram_bytes()
  estimate <- 1024^3 + 8 * n * (6 * ncolumns + 4 * ncol(frame))
  if (!is.finite(total) || total <= 0 || !is.finite(estimate) || estimate > .75 * total)
    stop("Cannot safely allocate native spline contrast probes within the conservative RAM budget.")
  out <- frame[rep(1L, n), , drop = FALSE]
  for (i in active) if (!i %in% categorical) out[[i]][] <- 1
  for (j in seq_along(categorical)) {
    index <- rep(rep(seq_len(sizes[j]), each = prod(utils::head(sizes, j-1L))), length.out = n)
    i <- categorical[j]
    if (is.factor(out[[i]])) out[[i]][] <- levels[[j]][index] else
      out[[i]] <- if (is.logical(out[[i]])) levels[[j]][index] else
        factor(levels[[j]][index], levels = levels[[j]])
  }
  out
}

.msr_response_spline_sensitivity <- function(design, kernel, beta, values) {
  input <- lapply(seq_len(design$q), function(j) kernel[, j])
  # Evaluate each matrix derivative only once per original source/program.
  derivatives <- lapply(design$expression_programs, function(p)
    lapply(seq_len(design$q), function(j) if (p$zero[j]) NULL else p$derivative[[j]](input)))
  out <- matrix(0, design$n, design$q,
    dimnames = list(rownames(design$multipliers), design$covariates))
  for (group in design$groups) {
    involved <- group$kernels
    if (!length(involved)) next
    coefficient <- drop(design$multipliers[, group$columns, drop = FALSE] %*% beta[group$columns])
    for (leaf in involved) {
      value <- coefficient
      for (other in setdiff(involved, leaf)) value <- value * values[[other]]
      program <- design$expression_leaf_program[leaf]; column <- design$expression_leaf_column[leaf]
      for (j in which(!design$expression_programs[[program]]$zero)) {
        derivative <- derivatives[[program]][[j]]
        if (is.list(derivative)) derivative <- derivative[[column]]
        out[, j] <- out[, j] + value * derivative
      }
    }
  }
  if (any(!is.finite(out))) stop("Dynamic response sensitivity is nonfinite.")
  out
}

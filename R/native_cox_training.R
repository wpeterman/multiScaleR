# Recent survival versions can strip the attributes of matrix columns while
# applying na.omit, even when every row is retained. Compare known numeric
# learner output while allowing only this specific change of representation.
.msr_cox_plain_basis <- function(value, original) {
  learner <- c("poly", "ns", "bs")
  if (!is.matrix(value) || !is.numeric(value) || !is.matrix(original) ||
      !is.numeric(original) ||
      !all(class(value) %in% c("AsIs", "matrix", "array")) ||
      !all(class(original) %in% c(learner, "AsIs", "basis", "matrix", "array")) ||
      !any(class(original) %in% learner) ||
      !identical(inherits(value, "AsIs"), inherits(original, "AsIs"))) return(FALSE)
  dropped <- setdiff(names(attributes(original)), names(attributes(value)))
  all(dropped %in% c("class", "coefs", "degree", "knots", "Boundary.knots", "intercept"))
}

# Restore lost attributes from the raw data supplied to the package's final
# refit after reproducing the complete observed frame. Scaling loses the
# original units, so the fitted normalized values alone cannot verify
# externally supplied values.
.msr_cox_training_frame <- function(model, frame, terms, data) {
  call <- model$call
  required <- all.vars(terms)
  if (!is.call(call) || !is.data.frame(data) || !all(required %in% names(data)) ||
      length(model$na.action) ||
      !is.null(call$subset) || !is.null(call$weights) ||
      !is.null(call$offset) || !is.null(call$cluster)) return(frame)
  replay <- tryCatch(stats::model.frame(terms, data = data, xlev = model$xlevels,
    na.action = stats::na.pass), error = function(e) NULL, warning = function(w) NULL)
  if (!is.data.frame(replay) || !identical(names(frame), names(replay)) ||
      !identical(row.names(frame), row.names(replay))) return(frame)
  # Compare all columns, including the response, strata and formula offsets.
  # Only learners that can be identified by the existing parser are repaired.
  equivalent <- function(a, b) {
    if (!is.numeric(a) || !is.numeric(b) || inherits(a, "Surv") || inherits(b, "Surv"))
      return(identical(a, b))
    if (!identical(dimnames(a), dimnames(b)) || !identical(dim(a), dim(b))) return(FALSE)
    if (.msr_cox_plain_basis(a, b)) {
      a <- unclass(a)
      b <- unclass(b)
    }
    .msr_curve_numeric_equivalent(a, b)
  }
  if (!all(vapply(seq_along(frame), function(i) equivalent(frame[[i]], replay[[i]]),
      logical(1)))) return(frame)
  # Require the evaluated training data to preserve native factor coding,
  # coefficient order and every fitted design value before copying metadata.
  verified <- tryCatch({
    retained_design <- stats::model.matrix(model)
    for (name in intersect(names(model$contrasts), names(replay))) {
      if (is.factor(replay[[name]]) && is.null(attr(replay[[name]], "contrasts")))
        stats::contrasts(replay[[name]]) <- model$contrasts[[name]]
    }
    replay_design <- stats::model.matrix(model, data = replay)
    identical(dimnames(retained_design), dimnames(replay_design)) &&
      .msr_curve_numeric_equivalent(retained_design, replay_design)
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!verified) return(frame)
  # Restore only attributes after verification. Preserve the stored numbers,
  # contrasts and all other columns of the native fitted frame.
  for (i in seq_along(frame)) {
    original <- replay[[i]]
    restore <- intersect(setdiff(names(attributes(original)), names(attributes(frame[[i]]))),
      c("scaled:center", "scaled:scale", "coefs", "degree", "knots", "Boundary.knots", "intercept"))
    if (!length(restore)) next
    value <- frame[[i]]
    if (!is.numeric(value)) next
    for (name in restore) attr(value, name) <- attr(original, name)
    if (.msr_cox_plain_basis(value, original)) class(value) <- class(original)
    frame[[i]] <- value
  }
  frame
}

# Dynamic Cox splines already retain their original basis parameters in
# predvars. Supply just those attributes to the response-design compiler.
.msr_cox_spline_frame <- function(terms, frame, covariates) {
  variables <- as.list(attr(terms, "variables"))[-1L]
  prediction <- as.list(attr(terms, "predvars"))[-1L]
  for (j in seq_along(variables)) {
    ex <- variables[[j]]
    type <- .msr_spline_call(ex)
    value <- frame[[j]]
    if (is.null(type) || inherits(value, type) || !is.matrix(value) ||
        !is.numeric(value) || is.object(value) ||
        !length(intersect(all.vars(ex), covariates))) next
    saved <- prediction[[j]]
    if (!identical(.msr_spline_call(saved), type)) next
    # The compiler independently checks that these metadata match the formula
    # and rejects adaptive knots, arbitrary functions and unsupported degree.
    attr(value, "knots") <- saved$knots
    attr(value, "Boundary.knots") <- saved$Boundary.knots
    attr(value, "intercept") <- saved$intercept
    class(value) <- c(type, "basis", "matrix")
    frame[[j]] <- value
  }
  frame
}

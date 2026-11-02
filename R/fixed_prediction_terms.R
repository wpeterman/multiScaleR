# Recover known fixed prediction metadata from evaluated training columns.
# Keep native expression names and values; only prediction calls gain metadata.
.msr_fixed_prediction_terms <- function(terms, frame) {
  variables <- as.list(attr(terms, "variables"))
  prediction <- attr(terms, "predvars")
  if (is.null(prediction)) prediction <- attr(terms, "variables")
  prediction <- as.list(prediction)
  if (length(variables) != length(prediction) || length(variables) - 1L > ncol(frame))
    stop("Cannot align fixed prediction terms with the retained training frame.", call. = FALSE)
  for (j in seq_len(length(variables) - 1L)) {
    ex <- prediction[[j + 1L]]
    if (!is.call(ex) || !is.call(head <- ex[[1L]]) ||
        !identical(head[[1L]], as.name("::")) || length(head) != 3L) next
    known <- (identical(head[[2L]], as.name("stats")) && identical(head[[3L]], as.name("poly"))) ||
      (identical(head[[2L]], as.name("base")) && identical(head[[3L]], as.name("scale")))
    if (!known) next
    bare <- ex; bare[[1L]] <- head[[3L]]
    learned <- stats::makepredictcall(frame[[j]], bare)
    learned[[1L]] <- head
    prediction[[j + 1L]] <- learned
  }
  attr(terms, "predvars") <- as.call(prediction)
  terms
}

# Native prediction providers keep fixed training transformations in different
# terms fields. Restore only those fields, never a fitted numerical quantity.
.msr_needs_fixed_metadata <- function(terms) {
  prediction <- attr(terms, "predvars")
  if (is.null(prediction)) prediction <- attr(terms, "variables")
  variables <- as.list(prediction)[-1L]
  fixed <- function(value) is.numeric(value) || identical(value, FALSE)
  any(vapply(variables, function(ex) {
    if (!is.call(ex) || !is.call(head <- ex[[1L]]) || length(head) != 3L ||
        !identical(head[[1L]], as.name("::"))) return(FALSE)
    if (identical(head[[2L]], as.name("base")) && identical(head[[3L]], as.name("scale")))
      return(!fixed(ex$center) || !fixed(ex$scale))
    if (identical(head[[2L]], as.name("stats")) && identical(head[[3L]], as.name("poly")))
      return(!isTRUE(ex$raw) && !is.list(ex$coefs))
    FALSE
  }, logical(1)))
}

.msr_native_fixed_metadata <- function(model) {
  if (inherits(model, "merMod")) {
    tt <- stats::terms(model, fixed.only = TRUE)
    if (!.msr_needs_fixed_metadata(tt)) return(model)
    frame <- model@frame
  } else if (inherits(model, "glmmTMB")) {
    frame <- model$frame
    tt <- attr(frame, "terms")
    if (!.msr_needs_fixed_metadata(tt)) return(model)
  } else if (inherits(model, c("lm", "glm"))) {
    frame <- model$model
    tt <- model$terms
    parametric <- if (inherits(model, "gam")) model$pterms else NULL
    candidates <- if (is.list(parametric) && !inherits(parametric, "terms"))
      parametric else list(tt, parametric)
    if (!any(vapply(candidates, .msr_needs_fixed_metadata, logical(1)))) return(model)
  } else return(model)
  original <- attr(frame, "terms")
  if (!is.data.frame(frame) || !inherits(original, "terms"))
    stop("Fixed prediction transformations require the retained evaluated training frame.", call. = FALSE)
  repair <- function(terms) .msr_fixed_prediction_terms(terms,
    .msr_response_frame_subset(original, frame, terms))
  if (inherits(model, "gam") && !inherits(parametric, "terms")) {
    updated <- lapply(parametric, repair)
    if (!identical(updated, parametric))
      stop("Qualified fixed transformations in multi-component GAMs require separately checked component metadata.", call. = FALSE)
    return(model)
  }
  fixed <- repair(tt)
  if (inherits(model, "merMod")) {
    if (identical(fixed, tt)) return(model)
    attr(original, "predvars.fixed") <- attr(fixed, "predvars")
    attr(frame, "terms") <- original
    model@frame <- frame
  } else if (inherits(model, "glmmTMB")) {
    if (identical(fixed, tt)) return(model)
    attr(frame, "terms") <- fixed
    model$frame <- frame
  } else {
    if (inherits(model, "gam")) model$pterms <- repair(parametric)
    model$terms <- fixed
    attr(model$model, "terms") <- repair(original)
  }
  model
}

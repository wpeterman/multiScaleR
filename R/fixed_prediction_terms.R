# Recover known fixed prediction metadata from evaluated training columns.
# Keep native expression names and values; only prediction calls gain metadata.
.msr_fixed_prediction_terms <- function(terms, frame, nested = FALSE) {
  validate_scalar_logical(nested, "nested")
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
  if (nested) for (j in seq_len(length(variables)-1L)) {
    if (j == attr(terms, "response")) next
    prediction[[j+1L]] <- .msr_fixed_nested_call(prediction[[j+1L]], frame[[j]], environment(terms))
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

.msr_native_fixed_metadata <- function(model, nested = FALSE) {
  validate_scalar_logical(nested, "nested")
  needs <- function(tt) if (nested) .msr_needs_nested_metadata(tt) else
    .msr_needs_fixed_metadata(tt)
  if (inherits(model, "merMod")) {
    tt <- stats::terms(model, fixed.only = TRUE)
    if (!needs(tt)) return(model)
    frame <- model@frame
  } else if (inherits(model, "glmmTMB")) {
    frame <- model$frame
    tt <- attr(frame, "terms")
    if (!needs(tt)) return(model)
  } else if (inherits(model, "coxph")) {
    tt <- stats::terms(model)
    if (!needs(tt)) return(model)
    # Only retained time-fixed, unpenalized Cox prediction terms are qualified.
    # This repairs metadata without changing a native fit or its uncertainty.
    if (inherits(model, c("coxphms", "coxph.penal")) ||
        length(attr(tt, "specials")$tt))
      stop("Fixed prediction recovery requires a checked time-fixed, unpenalized native Cox structure.", call. = FALSE)
    frame <- model$model
  } else if (inherits(model, c("lm", "glm"))) {
    frame <- model$model
    tt <- model$terms
    parametric <- if (inherits(model, "gam")) model$pterms else NULL
    candidates <- if (is.list(parametric) && !inherits(parametric, "terms"))
      parametric else list(tt, parametric)
    if (!any(vapply(candidates, needs, logical(1)))) return(model)
  } else return(model)
  original <- attr(frame, "terms")
  if (!is.data.frame(frame) || !inherits(original, "terms"))
    stop("Fixed prediction transformations require the retained evaluated training frame.", call. = FALSE)
  repair <- function(terms) .msr_fixed_prediction_terms(terms,
    .msr_response_frame_subset(original, frame, terms), nested = nested)
  if (inherits(model, "gam") && !inherits(parametric, "terms")) {
    # mgcv stores combined term labels on the list itself. Preserve them,
    # component order and coefficient spans while repairing each terms object.
    labels <- attr(parametric, "term.labels")
    if (length(parametric) != length(model$nsdf) ||
        length(labels) != sum(vapply(parametric, function(t)
          length(attr(t, "term.labels")), integer(1))))
      stop("Cannot align native GAM component terms and coefficient spans.", call. = FALSE)
    updated <- parametric
    updated[] <- lapply(parametric, repair)
    model$pterms <- updated
    model$terms <- repair(tt)
    attr(model$model, "terms") <- repair(original)
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


# A native terms object can retain evaluated variables removed with '-'.
# Prune only those inactive frame rows in a LOCAL prediction copy. Preserve the
# original learned calls, term/column order, response, offsets and native specials.
.msr_native_active_metadata <- function(model) {
  if (!class(model)[1L] %in% c("lm", "glm", "negbin", "coxph", "clogit")) return(model)
  tt <- stats::terms(model)
  specials <- attr(tt, "specials")
  if (inherits(model, "coxph") && (inherits(model, c("coxphms", "coxph.penal")) ||
      length(specials$tt))) return(model)
  fail <- function() stop("Cannot align original native prediction-term metadata; retain the original fitted terms and frame.", call. = FALSE)
  variables <- attr(tt, "variables")
  if (!is.call(variables) || !identical(variables[[1L]], as.name("list"))) fail()
  n <- length(variables) - 1L
  factors <- attr(tt, "factors")
  used <- rep(FALSE, n)
  if (is.matrix(factors)) {
    if (nrow(factors) != n) fail()
    used <- rowSums(factors != 0L) > 0L
  } else if (length(factors)) fail()
  response <- attr(tt, "response")
  required <- c(if (response > 0L) response, attr(tt, "offset"),
    unlist(specials, use.names = FALSE))
  if (length(required) && (anyNA(required) || any(required < 1L | required > n) ||
      any(required != as.integer(required)))) fail()
  used[required] <- TRUE
  if (all(used)) return(model)
  keep <- which(used)
  subset_call <- function(x) {
    if (!is.call(x) || length(x) != n + 1L ||
        !identical(x[[1L]], as.name("list"))) fail()
    as.call(as.list(x)[c(1L, keep + 1L)])
  }
  reindex <- function(index) {
    if (is.null(index)) return(NULL)
    out <- match(index, keep)
    if (anyNA(out)) fail()
    out
  }
  attr(tt, "variables") <- subset_call(variables)
  predvars <- attr(tt, "predvars")
  if (!is.null(predvars)) attr(tt, "predvars") <- subset_call(predvars)
  if (is.matrix(factors)) attr(tt, "factors") <- factors[keep, , drop = FALSE]
  attr(tt, "offset") <- reindex(attr(tt, "offset"))
  if (!is.null(specials)) attr(tt, "specials") <- lapply(specials, reindex)
  if (response > 0L) attr(tt, "response") <- reindex(response)
  classes <- attr(tt, "dataClasses")
  if (!is.null(classes)) {
    if (length(classes) < n) fail()
    # model.frame may append native weights/offset classes after the variables.
    extra <- if (length(classes) > n) seq.int(n + 1L, length(classes)) else integer()
    attr(tt, "dataClasses") <- classes[c(keep, extra)]
  }
  model$terms <- tt
  model
}

# Compare retained likelihood outcomes, rather than a mutable formula call.
.msr_selection_frame <- function(model) {
  if (inherits(model, "merMod")) return(model@frame)
  if (inherits(model, "HLfit")) {
    if (!requireNamespace("spaMM", quietly = TRUE)) return(NULL)
    # spaMM builds this from its retained data, not a caller's mutable data name.
    frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    if (is.data.frame(frame)) return(frame)
    return(NULL)
  }
  if (!is.list(model)) return(NULL)
  if (is.data.frame(model$model)) return(model$model)
  if (is.data.frame(model$frame)) return(model$frame)
  NULL
}

.msr_selection_signature <- function(model) {
  model <- .analysis_model(model)
  frame <- .msr_selection_frame(model)
  unmarked <- methods::is(model, "unmarkedFit")
  if (unmarked) {
    response <- model@data@y
    ids <- rownames(model@data@siteCovs)
    rows <- setdiff(seq_len(nrow(response)), model@sitesRemoved)
    response <- response[rows, , drop = FALSE]
    ids <- ids[rows]
    weights <- rep(1, length(rows))
    family <- "unmarked"
  } else {
    response <- if (!is.null(frame)) tryCatch(stats::model.response(frame),
      error = function(e) NULL) else NULL
    if (is.null(response) && is.list(model)) response <- model$y
    if (is.null(response) && inherits(model, c("gls", "lme")))
      response <- stats::fitted(model) + stats::residuals(model, type = "response")
    if (is.null(response)) return(NULL)
    n <- NROW(response)
    ids <- if (!is.null(frame)) rownames(frame) else .msr_observation_ids(model)
    weights <- if (!is.null(frame)) stats::model.weights(frame) else NULL
    if (is.null(weights) && is.list(model) && !inherits(model, "glm"))
      weights <- model$weights
    if (is.null(weights)) weights <- rep(1, n)
    family <- tryCatch(stats::family(model)$family, error = function(e) "unknown")
    if (!is.character(family) || length(family) != 1L) family <- "unknown"
    if (inherits(model, c("lm", "gls", "lme")) && !inherits(model, "glm"))
      family <- "gaussian"
    if (inherits(model, c("zeroinfl", "hurdle"))) family <- "native_count"
    effective <- if (inherits(model, "glm")) model$prior.weights else weights
    keep <- if (length(effective) == n) effective != 0 else rep(TRUE, n)
    response <- if (is.matrix(response)) response[keep, , drop = FALSE] else response[keep]
    weights <- weights[keep]
    ids <- ids[keep]
  }
  if (is.null(ids) || length(ids) != NROW(response) || anyDuplicated(ids)) return(NULL)
  # Binomial counts and proportion/trial input are equivalent only with the
  # same native likelihood constants. Preserve case replication separately.
  trials <- NULL
  if (identical(family, "binomial") && !is.null(frame)) {
    if (is.factor(response)) response <- as.integer(response != levels(response)[1L])
    if (is.matrix(response) && ncol(response) == 2L) {
      trials <- rowSums(response)
      response <- response[, 1L]
    } else if (is.numeric(response) &&
        all(weights == round(weights)) &&
        all(response * weights == round(response * weights))) {
      trials <- weights
      response <- response * weights
      weights <- rep(1, length(response))
    }
  }
  reml <- if (inherits(model, "merMod")) lme4::isREML(model) else
    if (inherits(model, "glmmTMB")) isTRUE(model$modelInfo$REML) else
    if (inherits(model, c("gls", "lme"))) identical(model$method, "REML") else
      if (inherits(model, "HLfit")) {
        label <- tryCatch(names(stats::logLik(model)), error = function(e) NULL)
        isTRUE(tolower(label) %in% "p_bv")
      } else FALSE
  design <- if (reml) tryCatch(stats::model.matrix(model), error = function(e) NULL) else NULL
  if (reml && is.null(design)) {
    dat <- .model_data(model)
    design <- tryCatch(stats::model.matrix(stats::delete.response(stats::terms(model)),
      data = dat, contrasts.arg = model$contrasts), error = function(e) NULL)
  }
  if (reml && inherits(model, "HLfit") && !is.null(design) &&
      is.null(rownames(design)) && !is.null(frame) && nrow(design) == nrow(frame)) {
    # spaMM's native X omits row names; its rows follow the retained frame.
    rownames(design) <- rownames(frame)
  }
  partial <- inherits(model, "coxph")
  strata <- NULL
  if (partial && !is.null(frame)) {
    special <- attr(stats::terms(model), "specials")$strata
    if (length(special)) strata <- frame[special]
  }
  order <- order(ids)
  response <- if (is.matrix(response)) unname(response[order, , drop = FALSE]) else
    unname(response[order])
  measure <- if (partial) "partial_likelihood" else
    if (unmarked || family %in% c("poisson", "binomial", "native_count") ||
        grepl("Negative Binomial|nbinom", family)) "discrete" else
    if (tolower(family) %in% c("gaussian", "gamma", "inverse.gaussian", "beta", "betar"))
      "continuous" else "unknown"
  list(response = response, weights = unname(weights[order]), measure = measure,
    trials = if (is.null(trials)) NULL else unname(trials[order]), reml = reml,
    design = if (is.null(design)) NULL else design[match(sort(ids), rownames(design)), , drop = FALSE],
    partial_method = if (partial) model$method else NULL,
    strata = if (!partial) NULL else if (is.null(strata)) rep(1L, length(ids)) else {
      groups <- as.character(interaction(strata[order, , drop = FALSE], drop = TRUE))
      match(groups, unique(groups))
    })
}

.msr_assert_selection_targets <- function(models) {
  if (length(models) < 2L) return(invisible(NULL))
  signatures <- lapply(models, .msr_selection_signature)
  available <- Filter(Negate(is.null), signatures)
  if (length(available) < 2L) return(invisible(NULL))
  ref <- available[[1L]]
  equal <- function(x, y) isTRUE(all.equal(x, y, tolerance = 0, check.attributes = FALSE))
  for (other in available[-1L]) {
    if (!equal(ref$response, other$response))
      stop("Models have different retained response values or response layouts; their likelihoods are not comparable.", call. = FALSE)
    if (!equal(ref$weights, other$weights))
      stop("Models have different likelihood weights; compare the same weighted observations.", call. = FALSE)
    if (ref$measure != "unknown" && other$measure != "unknown" &&
        ref$measure != other$measure)
      stop("Models use incompatible likelihood measures (continuous, discrete, or partial).", call. = FALSE)
    if (!is.null(ref$trials) && !is.null(other$trials) && !equal(ref$trials, other$trials))
      stop("Binomial trial counts differ between models.", call. = FALSE)
    if (!identical(ref$reml, other$reml))
      stop("ML and REML likelihoods are not comparable.", call. = FALSE)
    if (ref$reml && (is.null(ref$design) || is.null(other$design) ||
        !equal(ref$design, other$design)))
      stop("REML comparisons require the same recoverable fixed-effect design; refit with ML to compare different fixed effects.", call. = FALSE)
    if (!identical(ref$partial_method, other$partial_method))
      stop("Partial-likelihood methods differ between models.", call. = FALSE)
    if (!is.null(ref$partial_method) &&
        !identical(as.integer(ref$strata), as.integer(other$strata)))
      stop("Partial-likelihood models condition on different strata.", call. = FALSE)
  }
  invisible(NULL)
}

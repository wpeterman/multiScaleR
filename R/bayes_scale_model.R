# Formula preparation is deliberately independent of RTMB. Freeze observation
# identity, contrasts and bin summaries before recording the differentiation tape.
.msr_bayes_formula <- function(formula, data, response = FALSE, active_frame = FALSE) {
  if (!inherits(formula, "formula") || length(formula) != if (response) 3L else 2L)
    stop("Use a two-sided response formula and one-sided scale formulas.", call. = FALSE)
  groups <- character()
  walk <- function(x, additive = TRUE) {
    if (!is.call(x)) return(x)
    op <- as.character(x[[1L]])[1L]
    if (op == "(") return(walk(x[[2L]], additive))
    if (op == "|") {
      if (!additive || !identical(x[[2L]], 1) || !is.symbol(x[[3L]]))
        stop("Random terms must be additive random intercepts, (1 | group).", call. = FALSE)
      groups <<- c(groups, as.character(x[[3L]]))
      return(NULL)
    }
    if (op == "+" && length(x) == 3L) {
      left <- walk(x[[2L]], additive)
      right <- walk(x[[3L]], additive)
      if (is.null(left)) return(right)
      if (is.null(right)) return(left)
      return(call("+", left, right))
    }
    for (i in seq.int(2L, length(x)))
      x[[i]] <- walk(x[[i]], additive && op == "+")
    x
  }
  fixed <- formula
  rhs <- walk(formula[[length(formula)]])
  fixed[[length(fixed)]] <- if (is.null(rhs)) 1 else rhs
  if (anyDuplicated(groups))
    stop("A grouping variable can occur only once in each formula.", call. = FALSE)
  if (length(setdiff(groups, names(data))))
    stop("Random-effect grouping variables must be columns in data.", call. = FALSE)
  # Newly expanded/removed prepared terms must not demand dead predictors at
  # prediction time. Keep native column order, offsets and learned metadata.
  evaluated <- if (active_frame)
    .msr_response_active_terms(stats::terms(fixed, data = data)) else fixed
  frame <- stats::model.frame(evaluated, data, na.action = stats::na.fail,
                              drop.unused.levels = TRUE)
  random <- lapply(groups, function(g) {
    v <- data[[g]]
    if (!is.atomic(v) || !is.null(dim(v)) || anyNA(v) || any(!nzchar(as.character(v))))
      stop("Grouping variables must have nonmissing, nonempty values.", call. = FALSE)
    f <- droplevels(as.factor(v))
    if (nlevels(f) < 2L)
      stop("Each random-effect grouping variable needs at least two levels.", call. = FALSE)
    list(group = g, levels = levels(f), index = as.integer(f))
  })
  names(random) <- groups
  out <- list(formula = formula, fixed = fixed, frame = frame,
       terms = stats::terms(frame), random = random)
  if (active_frame) out$active_frame <- TRUE
  out
}

.msr_bayes_priors <- function(priors, covariates, inputs) {
  defaults <- list(scale_median = sqrt(inputs$min_D * inputs$max_D),
    scale_intercept_sd = 1, scale_coef_sd = 0.5, scale_group_sd = 0.5,
    beta_sd = 2.5, response_group_sd = 1,
    dispersion_log_mean = 0, dispersion_log_sd = 1)
  if (!is.list(priors) || (length(priors) &&
      (is.null(names(priors)) || anyDuplicated(names(priors)) ||
       any(!names(priors) %in% names(defaults)))))
    stop("priors must be a named list of documented prior settings.", call. = FALSE)
  defaults[names(priors)] <- priors
  for (nm in setdiff(names(defaults), "scale_median"))
    validate_scalar_numeric(defaults[[nm]], nm,
      positive = nm != "dispersion_log_mean")
  med <- defaults$scale_median
  if (!is.numeric(med) || any(!is.finite(med) | med <= 0))
    stop("scale_median must contain positive distances in map units.", call. = FALSE)
  if (length(med) == 1L && is.null(names(med))) med <- rep(med, length(covariates))
  else if (!is.null(names(med)) && !anyDuplicated(names(med)) &&
           setequal(names(med), covariates)) med <- med[covariates]
  else stop("scale_median must be a scalar or name every kernel predictor.", call. = FALSE)
  defaults$scale_median <- stats::setNames(med, covariates)
  defaults
}

# Attach the existing proper priors after neutral design preparation. Keep the
# legacy field order and starts so saved objects and objective code retain their
# original contract. The likelihood preparer never calls this function.
.msr_bayes_attach_priors <- function(spec, priors) {
  for (nm in spec$covariates) {
    sm <- spec$scales[[nm]]
    mean <- ifelse(sm$distance_location,
      log(priors$scale_median[[nm]] / spec$unit_conv), 0)
    sd <- ifelse(sm$distance_location, priors$scale_intercept_sd, priors$scale_coef_sd)
    sm <- c(sm[setdiff(names(sm), c("prior_mean", "prior_sd", "representatives", "labels"))],
      list(prior_mean = mean, prior_sd = sd), sm[c("representatives", "labels")])
    spec$parameters$scale_beta[sm$beta_index] <- mean
    spec$scales[[nm]] <- sm
  }
  # Concatenation also preserves names on accepted scalar prior settings.
  # Subassignment into the neutral vector would silently discard those names.
  scale_sd <- response_sd <- numeric()
  for (i in seq_along(spec$parameters$log_scale_sd))
    scale_sd <- c(scale_sd, log(priors$scale_group_sd / 2))
  for (i in seq_along(spec$parameters$log_response_sd))
    response_sd <- c(response_sd, log(priors$response_group_sd / 2))
  spec$parameters$log_scale_sd <- scale_sd
  spec$parameters$log_response_sd <- response_sd
  spec$parameters$log_dispersion <- if (spec$family %in% c("gaussian", "nbinom2"))
    priors$dispersion_log_mean else numeric()
  c(spec[setdiff(names(spec), c("priors", "parameters"))],
    list(priors = priors, parameters = spec$parameters))
}

.msr_bayes_prepare <- function(formula, inputs, data, scales, family, row_id, priors) {
  prepared <- .msr_scale_prepare_inputs(formula, inputs, data, row_id)
  pr <- .msr_bayes_priors(priors, prepared$covariates, inputs)
  spec <- .msr_scale_prepare(formula, inputs, data, scales, family, row_id,
    reference_scale = pr$scale_median, prepared = prepared)
  .msr_bayes_attach_priors(spec, pr)
}

# Private value-only prerequisite. This consumes actual ordered raw feature
# matrices, not fabricated kernel_inputs. An external provider's geometry is
# still caller-declared, so public optimizer/AIC/CV admission remains closed.
#' Freeze a private raw-design binomial context
#'
#' @param model An ordinary retained binomial-logit GLM with fixed formula offsets.
#' @param covariates Ordered names of raw numeric prepared feature columns.
#' @param reference Numeric reference matrix with exact model row IDs and feature names.
#' @param source_id A declared source SHA256. Geometry authentication remains pending.
#' @param scale_parameters Names of future positive scale coordinates. Counts
#'   remain nominal until a separately qualified outer search is completed.
#' @return A private context with frozen design recipes, target and reference
#'   identities. No original model, mutable fit cache or inference eligibility is retained.
#' @noRd
.msr_prepare_raw_binomial <- function(model, covariates, reference, source_id,
    scale_parameters = character()) {
  if (!inherits(model, "glm") || inherits(model, "gam") ||
      !identical(model$family$family, "binomial") || !identical(model$family$link, "logit") ||
      !is.data.frame(model$model) || !is.null(model$call$offset))
    stop("Prepared raw binomial fitting requires an ordinary retained binomial-logit GLM with formula offsets.")
  if (!is.character(covariates) || !length(covariates) || anyNA(covariates) ||
      anyDuplicated(covariates) || any(!nzchar(covariates)) ||
      !is.character(source_id) || length(source_id) != 1L || is.na(source_id) ||
      !grepl("^[[:xdigit:]]{64}$", source_id) || !is.character(scale_parameters) ||
      anyNA(scale_parameters) || anyDuplicated(scale_parameters) || any(!nzchar(scale_parameters)))
    stop("Declare unique feature/scale names and one SHA256 source identifier.")
  frame <- model$model; ids <- rownames(frame)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) stop("Training point IDs must be explicit and unique.")
  valid_features <- function(x) is.matrix(x) && is.numeric(x) &&
    identical(dim(x), c(nrow(frame), length(covariates))) && all(is.finite(x)) &&
    identical(rownames(x), ids) && identical(colnames(x), covariates)
  if (!valid_features(reference)) stop("Raw reference rows, columns and point order must match the model.")
  terms <- .msr_fixed_prediction_terms(stats::terms(model), frame)
  for (index in attr(terms, "offset"))
    if (length(intersect(all.vars(attr(terms, "variables")[[index + 1L]]), covariates)))
      stop("Prepared raw features cannot occur in offsets.")
  X <- stats::model.matrix(terms, frame, contrasts.arg = model$contrasts)
  design <- .msr_compile_response_design(terms, frame, covariates, colnames(X), model$contrasts)
  checked <- .msr_response_design_matrix(design, reference)
  if (!isTRUE(all.equal(unname(checked), unname(X), tolerance = 1e-10, check.attributes = FALSE)))
    stop("The raw reference does not reproduce the original native design.")
  # The compiled programs retain numeric recipes. Do not retain the caller's
  # formula environment, model or mutable coefficient fitting state.
  attr(design$native_terms, ".Environment") <- baseenv()
  y <- stats::model.response(frame); trials <- rep(1, nrow(frame))
  if (is.matrix(y) && ncol(y) == 2L) {
    if (any(!is.finite(y) | y < 0 | y != floor(y))) stop("Use integer success/failure counts.")
    trials <- rowSums(y); y <- y[, 1L]
  } else {
    if (is.logical(y)) y <- as.integer(y)
    if (is.factor(y) && nlevels(y) == 2L) y <- as.integer(y) - 1L
  }
  if (!is.numeric(y) || length(y) != nrow(frame) || any(!is.finite(y) | y < 0 |
      y != floor(y) | y > trials | trials < 1)) stop("Use binary outcomes or positive integer binomial trials.")
  frequency <- stats::model.weights(frame)
  if (is.null(frequency)) frequency <- rep(1, nrow(frame))
  if (any(!is.finite(frequency) | frequency <= 0 | frequency != floor(frequency)))
    stop("Prepared binomial weights must be positive integer frequencies.")
  offset <- stats::model.offset(frame)
  if (is.null(offset)) offset <- rep(0, nrow(frame))
  if (any(!is.finite(offset))) stop("Offsets must be finite and fixed.")
  target <- list(ids = ids, successes = as.numeric(y), trials = as.numeric(trials),
    frequency = as.numeric(frequency), offset = as.numeric(offset))
  structure(list(design = design, target = target, source_id = tolower(source_id),
    design_id = digest::digest(design, algo = "sha256"),
    target_id = digest::digest(target, algo = "sha256"),
    reference_id = digest::digest(list(reference, colnames(X), model$contrasts), algo = "sha256"),
    covariates = covariates, columns = colnames(X), n = nrow(frame), p = ncol(X),
    scale_parameters = scale_parameters, transformation = "raw_identity",
    source_authentication = "caller_declared_geometry_unverified",
    selection_eligible = FALSE, uncertainty_available = FALSE), class = "msr_prepared_raw_binomial")
}

#' Evaluate private raw-design binomial coefficients
#'
#' @param context A context from .msr_prepare_raw_binomial().
#' @param features Finite raw matrix with the original ordered row and column IDs.
#' @param scale_values Named positive coordinates in the prepared scale order.
#' @return Value-only normalized NLL, coefficient estimates and diagnostics,
#'   finite-MLE status and nominal counts. No SEs, AICc or scale intervals are supplied.
#' @noRd
.msr_evaluate_raw_binomial <- function(context, features, scale_values = numeric()) {
  if (!inherits(context, "msr_prepared_raw_binomial") || !is.matrix(features) ||
      !is.numeric(features) || !identical(dim(features), c(context$n, length(context$covariates))) ||
      any(!is.finite(features)) || !identical(rownames(features), context$target$ids) ||
      !identical(colnames(features), context$covariates))
    stop("Prepared features must retain exact ordered point and variable identities.")
  if (!is.numeric(scale_values) || length(scale_values) != length(context$scale_parameters) ||
      any(!is.finite(scale_values) | scale_values <= 0) ||
      (length(scale_values) && !identical(names(scale_values), context$scale_parameters)))
    stop("Declare finite positive scale coordinates in the prepared order.")
  if (!identical(digest::digest(context$target, algo = "sha256"), context$target_id))
    stop("The prepared response/row/weight/offset target was modified.")
  if (!identical(digest::digest(context$design, algo = "sha256"), context$design_id))
    stop("The prepared raw design or recipe was modified.")
  X <- .msr_response_design_matrix(context$design, features)
  target <- context$target
  mle <- .msr_binomial_mle_status(X, target$successes, target$trials)
  result <- list(value = NA_real_, coefficients = NULL, coefficient_fit_accepted = FALSE,
    finite_mle = mle, coefficient_count = context$p,
    nominal_scale_count = length(context$scale_parameters),
    nominal_parameter_count = context$p + length(context$scale_parameters),
    nobs_rows = context$n, frequency_total = sum(target$frequency),
    target_id = context$target_id, source_id = context$source_id,
    point_id = digest::digest(list(features, scale_values), algo = "sha256"),
    selection_eligible = FALSE, uncertainty_available = FALSE,
    diagnostics = list(warnings = character(), error = NULL,
      qualification = "private_value_only_not_outer_optimization"))
  if (mle$status %in% c("nonexistent", "non_identifiable")) return(result)
  warnings <- character(); error <- NULL
  fit <- tryCatch(withCallingHandlers(stats::glm.fit(X, target$successes / target$trials,
    weights = target$frequency * target$trials, offset = target$offset,
    family = stats::binomial(), control = stats::glm.control(epsilon = 1e-11, maxit = 100L)),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
    error = function(e) { error <<- conditionMessage(e); NULL })
  result$diagnostics$warnings <- warnings; result$diagnostics$error <- error
  if (is.null(fit)) return(result)
  result$coefficients <- fit$coefficients
  result$diagnostics$converged <- isTRUE(fit$converged)
  result$diagnostics$rank <- fit$rank
  if (!isTRUE(fit$converged) || fit$rank != context$p || any(!is.finite(fit$coefficients))) return(result)
  result$value <- -sum(target$frequency * stats::dbinom(target$successes,
    target$trials, fit$fitted.values, log = TRUE))
  score <- drop(crossprod(X, target$frequency * (target$trials * fit$fitted.values - target$successes)))
  scaled_score <- max(abs(score) / (1 + colSums(abs(X) * target$frequency * target$trials)))
  result$diagnostics$maximum_coefficient_score <- max(abs(score))
  result$diagnostics$scaled_coefficient_score <- scaled_score
  result$coefficient_fit_accepted <- identical(mle$status, "finite") &&
    !length(warnings) && is.finite(result$value) && is.finite(scaled_score) && scaled_score <= 1e-8
  result
}

# Bind new locations to the actual retained numeric training representation.
# Source declarations are checked separately; matching numbers cannot certify
# the scientific provenance of an arbitrary user-supplied raster.
.msr_likelihood_prediction_inputs <- function(spec, training_inputs, inputs) {
  if (is.null(training_inputs) || is.null(training_inputs$binned))
    stop("New locations require the original prepared training_inputs.", call. = FALSE)
  for (x in list(training_inputs, inputs)) {
    if (!identical(x$kernel, spec$kernel) || is.null(x$binned))
      stop("Supply matching prepared kernel means for training and prediction.", call. = FALSE)
    validate_scalar_numeric(x$unit_conv, "unit_conv", positive = TRUE)
    if (!isTRUE(all.equal(c(min_D = x$min_D, max_D = x$max_D), spec$support, tolerance = 0)))
      stop("Prediction and training buffer support must match the retained fit.", call. = FALSE)
  }
  rows <- match(spec$observation_ids, .msr_scale_point_ids(training_inputs))
  if (anyNA(rows)) stop("Training point IDs do not match the retained specification.", call. = FALSE)
  engine <- .msr_prepare_kernel_engine(training_inputs$binned, spec$covariates, rows,
    spec$kernel, standardize = FALSE)
  distance <- engine$distance * training_inputs$unit_conv / spec$unit_conv
  for (nm in spec$covariates) {
    d <- matrix(rep(distance, each = nrow(engine$C[[nm]])), nrow(engine$C[[nm]]))
    nearest <- apply(replace(d, engine$C[[nm]] <= 0, Inf), 1L, min)
    delta <- pmax(if (spec$kernel == "gaussian") d^2 - nearest^2 else d - nearest, 0)
    same <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE))
    if (!same(engine$V[[nm]], spec$V[[nm]]) || !same(engine$C[[nm]], spec$C[[nm]]) ||
        !same(delta, spec$delta[[nm]]))
      stop("training_inputs do not reproduce the retained kernel representation.", call. = FALSE)
  }
  select <- function(x) {
    if (is.null(x$scale_vars)) return(NULL)
    z <- x$scale_vars[match(spec$covariates, x$scale_vars$covariate), , drop = FALSE]
    if (anyNA(z$covariate)) stop("Missing fitted source specifications.", call. = FALSE)
    as.data.frame(z)
  }
  reference <- select(training_inputs); actual <- select(inputs)
  if (!isTRUE(all.equal(reference, actual, tolerance = 0, check.attributes = FALSE)))
    stop("Prediction source and metric declarations must match training_inputs.", call. = FALSE)
  list(mean_specifications = reference,
    source_identity = if (is.null(reference)) "legacy_logical_names_and_verified_training_values" else
      "declared_sources_and_verified_training_values")
}

.msr_likelihood_parameters <- function(spec, parameters) {
  template <- spec$parameters
  if (!is.list(parameters) || is.null(names(parameters)) || anyDuplicated(names(parameters)) ||
      !setequal(names(parameters), names(template)))
    stop("Supply every retained parameter block with its original name.", call. = FALSE)
  parameters <- parameters[names(template)]
  for (nm in names(template)) if (!is.numeric(parameters[[nm]]) || is.object(parameters[[nm]]) || !is.null(dim(parameters[[nm]])) ||
      length(parameters[[nm]]) != length(template[[nm]]) || any(!is.finite(parameters[[nm]])))
    stop("Numeric finite parameter lengths must match the retained specification: ", nm, call. = FALSE)
  parameters
}

#' Prepare Private Conditional Likelihood Predictions
#'
#' Cache expected responses at one explicit population and realized-effect
#' vector. This private component does not fit a model or propagate uncertainty.
#'
#' @param spec A retained private marginal-likelihood specification, with no
#'   population priors. Its existing method/structure guards still apply.
#' @param parameters Complete finite numeric parameter blocks, including latent
#'   coordinates in the retained specification's order. Integrated coordinates
#'   must be reconstructed explicitly before calling this component.
#' @param newdata,kernel_inputs New observation data and matching prepared kernel
#'   means, or both NULL for retained training observations.
#' @param training_inputs Original prepared training inputs, required for new
#'   locations. Their values, counts, shifted distances and point IDs must
#'   reproduce the retained fit. Source declarations must match the new inputs.
#' @param row_id Optional point-ID column, with the existing shared prediction
#'   preparation semantics. Repeated observations at one point are allowed.
#' @param method Private likelihood method, "laplace" or "factorized".
#' @param scale_re.form,re.form NULL includes supplied scale/response effects;
#'   NA sets the corresponding effects to zero. Zero effects define a reference
#'   prediction, not integration over an effect distribution.
#' @param new_scale_levels,new_response_levels "error" rejects unseen groups;
#'   "zero" explicitly sets their effects to zero. No effects are sampled.
#' @param include_fixed Recompute verified fixed-radius predictors; default FALSE
#'   requires their explicit prepared values, as in the shared prediction API.
#' @return A private numeric snapshot with link/response vectors, kernel-scale
#'   and standardized-metric matrices, row/point IDs, compact provenance and
#'   separate prediction contracts. It holds no fit, raw buffers, response data,
#'   closures or AD tapes. Reuse it with .msr_likelihood_predict(); changed data,
#'   locations or parameter states require new preparation.
#' @details All population parameters and realized effects are fixed. Population
#'   estimation uncertainty, conditional latent variation, new-group variation
#'   and observation noise are omitted. No SE, interval, population-marginal mean
#'   or calibration is supplied. Supplied-state arithmetic grants no fit
#'   eligibility. Private spatial-field and weighted-SQ fitting guards remain.
#'   Binomial response output is success probability per trial; Poisson and NB2
#'   output is expected count. Gaussian output is the conditional response mean.
#'   Input checks verify numeric/source declarations, not scientific raster
#'   provenance. A legacy input without source definitions is labeled explicitly.
#' @rdname dot-msr_likelihood_point_snapshot
#' @keywords internal
.msr_likelihood_point_snapshot <- function(spec, parameters, newdata = NULL,
    kernel_inputs = NULL, training_inputs = NULL, row_id = NULL,
    method = c("laplace", "factorized"), scale_re.form = NULL, re.form = NULL,
    new_scale_levels = c("error", "zero"), new_response_levels = c("error", "zero"),
    include_fixed = FALSE) {
  method <- match.arg(method)
  .msr_scale_check_capability(spec, "marginal_likelihood", method)
  if (!is.null(spec$priors)) stop("A likelihood specification must not carry population priors.", call. = FALSE)
  parameters <- .msr_likelihood_parameters(spec, parameters)
  if (is.null(newdata) != is.null(kernel_inputs))
    stop("Supply newdata and kernel_inputs together.", call. = FALSE)
  if (is.null(newdata) && !is.null(training_inputs))
    stop("training_inputs applies only to new locations.", call. = FALSE)
  source <- if (is.null(newdata)) list(source_identity = "retained_training_specification") else
    .msr_likelihood_prediction_inputs(spec, training_inputs, kernel_inputs)
  scale_include <- .msr_bayes_prediction_re_form(scale_re.form, "scale_re.form")
  response_include <- .msr_bayes_prediction_re_form(re.form, "re.form")
  context <- .msr_bayes_prediction_prepare(list(specification = spec), newdata, kernel_inputs,
    row_id, match.arg(new_scale_levels), match.arg(new_response_levels), scale_include,
    response_include, include_fixed)
  p <- .msr_bayes_prediction_parameters(parameters, context)
  values <- .msr_bayes_components(context$spec, p, include_response = response_include)
  # A native one-column offset can carry matrix attributes into eta. Cached
  # response/link outputs retain the documented one-vector-per-row contract.
  values$eta <- as.numeric(values$eta)
  response <- switch(spec$family, gaussian = values$eta, binomial = stats::plogis(values$eta), exp(values$eta))
  scales <- do.call(cbind, values$scales); metrics <- do.call(cbind, values$kernel)
  if (any(!is.finite(c(values$eta, response, scales, metrics))) || any(scales <= 0))
    stop("Nonfinite conditional prediction or nonpositive scale; inspect the supplied state and inputs.", call. = FALSE)
  ids <- rownames(context$spec$data)
  if (is.null(ids)) ids <- as.character(seq_along(response))
  names(response) <- names(values$eta) <- ids
  rownames(scales) <- rownames(metrics) <- ids
  out <- list(linear_predictor = values$eta, expected_response = response, kernel_scale = scales,
    standardized_metric = metrics, row_ids = ids, point_ids = context$spec$observation_ids,
    new_groups = context$new_groups,
    provenance = c(list(method = method, family = spec$family, covariates = spec$covariates,
      unit_conv = spec$unit_conv, support = spec$support, center = spec$center, spread = spec$spread,
      fixed_specifications = spec$fixed_reference$specifications,
      parameter_state = parameters, origin = "supplied_parameter_vector"), source))
  out$prediction_contract <- lapply(c("expected_response", "linear_predictor", "kernel_scale"), function(quantity)
    .msr_prediction_contract("conditional_likelihood", method, quantity,
      if (quantity == "linear_predictor") "link" else if (quantity == "kernel_scale") "distance" else "response",
      population = "fixed_at_supplied_vector", kernel = "fixed_at_supplied_realized_effects",
      effects = list(scale = if (!length(spec$scales[[1L]]$random)) "no_active_scale_effects" else
        if (scale_include) "supplied_realized_effects" else "zero_effect_reference",
        response = if (!length(spec$response$random)) "no_active_response_effects" else
          if (response_include) "supplied_realized_effects" else "zero_effect_reference",
        new_groups = context$new_groups), conditioning = c("retained_training_transform", "supplied_population_and_latent_coordinates"),
      omitted = c("population_parameter_uncertainty", "conditional_latent_variation", "new_group_variation", "observation_noise"),
      coordinates = "original_retained_parameterization", fixed_inputs = TRUE))
  names(out$prediction_contract) <- c("expected_response", "linear_predictor", "kernel_scale")
  class(out) <- "msr_likelihood_prediction"
  out
}

#' @rdname dot-msr_likelihood_point_snapshot
#' @param object For preparation, a numerically admitted private marginal or
#'   factorized fit with its original selected start and zero-SD reduction.
#'   For prediction, a snapshot returned by either preparer.
#' @param ... New-location controls forwarded to the component preparer by the
#'   fitted adapter. Snapshot prediction rejects additional controls.
#' @param scale_effects,response_effects Optional named standardized latent
#'   coordinates z, keyed by exactly the retained group levels. The model's
#'   population SD multiplies these coordinates once. They are not already
#'   SD-scaled deviations. Marginal fits default to their retained modes;
#'   factorized fits with omitted integrated coordinates require explicit values,
#'   for example from their existing conditional reconstruction helper. This
#'   does not turn those values into a population-marginal mean or interval.
.msr_likelihood_prepare_prediction <- function(object, ..., scale_effects = NULL, response_effects = NULL) {
  method <- if (inherits(object, "msr_marginal_prototype")) "laplace" else
    if (inherits(object, "msr_factorized_prototype")) "factorized" else
      stop("Supply a private marginal or factorized fitted object.", call. = FALSE)
  .msr_scale_check_capability(object$specification, "marginal_likelihood", method)
  if (!isTRUE(object$diagnostics$passed) || is.null(object$parameters))
    stop("The original fit must pass all of its existing numerical and limiting gates.", call. = FALSE)
  selected <- object$selected
  if (!is.numeric(selected) || length(selected) != 1L || !is.finite(selected) ||
      selected < 1L || selected > length(object$attempts) || selected != as.integer(selected) ||
      !isTRUE(object$attempts[[selected]]$eligible) ||
      !identical(object$attempts[[selected]]$parameters, object$parameters) ||
      !identical(object$attempts[[selected]]$zero_sd, object$zero_sd))
    stop("The original selected eligible start and parameter state must be retained.", call. = FALSE)
  zero <- object$zero_sd
  if (!is.logical(zero) || !identical(sort(names(zero)), c("response", "scale")) || anyNA(zero))
    stop("The fitted zero-SD reduction must be retained explicitly.", call. = FALSE)
  spec <- .msr_marginal_reduce(object$specification, names(zero)[zero])
  p <- object$parameters
  for (component in c("scale", "response")) {
    name <- paste0(component, "_z"); effect <- if (component == "scale") scale_effects else response_effects
    terms <- if (component == "scale") spec$scales[[1L]]$random else spec$response$random
    if (!is.null(effect)) {
      if (length(terms) != 1L || !is.numeric(effect) || is.object(effect) || any(!is.finite(effect)) ||
          is.null(names(effect)) || anyDuplicated(names(effect)) || !setequal(names(effect), terms[[1L]]$levels))
        stop("Supply named realized ", component, " effects for exactly the retained group levels.", call. = FALSE)
      p[[name]] <- spec$parameters[[name]]
      p[[name]][terms[[1L]]$z_index] <- effect[terms[[1L]]$levels]
    } else if (length(p[[name]]) != length(spec$parameters[[name]]))
      stop("Integrated ", component, " effects require explicit reconstructed coordinates.", call. = FALSE)
  }
  out <- .msr_likelihood_point_snapshot(spec, p, ..., method = method)
  out$provenance$origin <- "numerically_admitted_private_fit_conditional_state"
  for (nm in names(out$prediction_contract)) {
    out$prediction_contract[[nm]]$diagnostics_status <- "passed"
    out$prediction_contract[[nm]]$uncertainty$population_parameters <- "fixed_at_fitted_estimate"
    out$prediction_contract[[nm]]$conditioning <- c(out$prediction_contract[[nm]]$conditioning,
      "observed_training_responses_through_fitted_state")
    out$prediction_contract[[nm]]$effects$coordinate_origin <- list(
      scale = if (!is.null(scale_effects)) "explicit_named_coordinates" else "retained_fitted_state",
      response = if (!is.null(response_effects)) "explicit_named_coordinates" else "retained_fitted_state")
  }
  out
}

#' @rdname dot-msr_likelihood_point_snapshot
#' @param quantity Output: "expected_response", "linear_predictor" or
#'   "kernel_scale". Scales are in the retained model's map units.
#' @param se.fit Must be FALSE. This snapshot has no qualified SE estimator.
#' @param interval Must be "none". Fixed-state outputs are not intervals.
.msr_likelihood_predict <- function(object, quantity = c("expected_response", "linear_predictor", "kernel_scale"),
    se.fit = FALSE, interval = "none", ...) {
  if (!inherits(object, "msr_likelihood_prediction")) stop("Supply a prepared private likelihood snapshot.", call. = FALSE)
  validate_scalar_logical(se.fit, "se.fit")
  if (isTRUE(se.fit) || !identical(interval, "none"))
    stop("This conditional snapshot supplies no parameter/latent SEs or intervals.", call. = FALSE)
  if (length(list(...))) stop("Snapshots do not accept new data or changed controls; prepare again.", call. = FALSE)
  quantity <- match.arg(quantity)
  out <- object[[quantity]]
  attr(out, "prediction_contract") <- object$prediction_contract[[quantity]]
  out
}

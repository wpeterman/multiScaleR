#' Read Prediction Contracts
#'
#' Shared metadata describes which quantity a prediction represents, which
#' effects it includes, and which uncertainty it propagates.
#'
#' @name prediction_contract
#' @details
#' Bayesian prediction summaries and posterior predictive simulations store this
#' list in \code{attr(result, "prediction_contract")}. Raw posterior draw matrices
#' and lists retain their existing attributes. Prepared prediction snapshots
#' instead store a named \code{prediction_contract} list, one entry per prepared
#' output type. Their native numerical outputs retain native attributes.
#'
#' The fields describe the calculation:
#' \describe{
#'   \item{schema_version}{Version of the metadata structure.}
#'   \item{inference_target, method}{The posterior, likelihood, conditional
#'     likelihood, or native target, and its retained computation label.}
#'   \item{quantity, output_scale}{Expected response, linear predictor, kernel
#'     scale, replicate response, choice probability, or named native component,
#'     on the stated output scale. Native component meanings remain model specific.}
#'   \item{target_population, effects}{The supplied prediction rows and policies
#'     for known and unseen groups. Sampled new effects are shared by each group's
#'     rows within a complete joint vector. Zero deviations are reference values,
#'     not integration over the group's distribution.}
#'   \item{conditioning, integration}{The quantities held fixed and the recorded
#'     treatment of uncertainty. Missing historical solver metadata is marked
#'     unrecorded rather than inferred.}
#'   \item{uncertainty}{Population and kernel parameter treatment, inclusion of
#'     observation noise, interval meaning, and omitted uncertainty. Analytic
#'     snapshots also record whether joint covariance was available at preparation.
#'     This does not bypass the prediction method's covariance checks.}
#'   \item{original_coordinates, joint_parameter_map_complete}{The parameter
#'     representation used by the adapter. For Bayesian predictions, completeness
#'     checks stored block dimensions against the original prepared parameter map.
#'     This structural check does not validate a coordinate transformation,
#'     reconstruct integrated coefficients, or establish posterior accuracy.}
#'   \item{diagnostics_status}{The existing fit status when available. It can be
#'     \code{"passed"}, \code{"failed"}, or \code{"unavailable"}. Joint covariance
#'     availability is separate from this status.}
#'   \item{metadata_establishes_eligibility}{Always \code{FALSE}. Read the original
#'     fit diagnostics before interpreting predictions.}
#'   \item{fixed_inputs, component}{Whether the rows and controls are bound to a
#'     snapshot, and its native component when applicable.}
#' }
#' Contracts describe the original prediction call. Ordinary data-frame
#' subsetting can retain its attributes; it does not recompute the contract for
#' a smaller batch. Use the existing row-level new-group flags for individual rows.
#'
#' Posterior intervals for expected responses exclude observation noise.
#' Posterior predictive replicates include it. New-group predictions also carry
#' the explicitly requested group variation. Native confidence intervals and
#' standard errors condition on selected kernel parameters. Likelihood summaries
#' conditional on one population parameter vector omit its estimation uncertainty.
#' These uncertainty targets are distinct. No contract establishes interval
#' calibration or upgrades a rejected fit.
#' @seealso \code{\link{predict.multiScaleR_bayes}},
#'   \code{\link{predict_scales}}, \code{\link{simulate.multiScaleR_bayes}},
#'   \code{\link{prepare_prediction}}
NULL

# Prediction metadata describes the calculation and its conditioning. It does
# not grant fit eligibility, add uncertainty or change the original diagnostics.
.msr_prediction_contract <- function(target, method, quantity, output_scale,
    population, kernel, effects = list(), conditioning = character(),
    integration = list(), intervals = "not_provided", omitted = character(),
    coordinates = "not_applicable", joint_map = NA, diagnostics = NULL,
    fixed_inputs = FALSE, component = NULL) {
  target <- match.arg(target, c("posterior", "likelihood", "conditional_likelihood", "native"))
  quantity <- match.arg(quantity, c("expected_response", "linear_predictor",
    "replicate_response", "kernel_scale", "choice_probability", "native_component"))
  status <- if (isTRUE(diagnostics$passed)) "passed" else
    if (isFALSE(diagnostics$passed)) "failed" else "unavailable"
  list(schema_version = 1L, inference_target = target, method = method,
    quantity = quantity, output_scale = output_scale,
    target_population = "specified_prediction_rows", effects = effects,
    conditioning = conditioning, integration = integration,
    uncertainty = list(population_parameters = population,
      kernel_parameters = kernel, observation_noise = quantity == "replicate_response",
      intervals = intervals, omitted = omitted, calibration = "not_established"),
    original_coordinates = coordinates, joint_parameter_map_complete = joint_map,
    diagnostics_status = status, metadata_establishes_eligibility = FALSE,
    fixed_inputs = fixed_inputs, component = component)
}

# Keep new-group sampling distinct from setting effects to zero. The existing
# new_groups attribute supplies row/group identities without duplicating them.
.msr_bayes_prediction_contract <- function(object, context, quantity,
    output_scale, include_scale = TRUE, include_response = TRUE) {
  spec <- object$specification
  group_contract <- function(present, include, policy, count) {
    if (!present) return(list(known = "not_applicable", unseen = "not_applicable",
                             unseen_present = FALSE))
    list(known = if (include) "joint_posterior_effects" else "zero_deviation_reference",
      unseen = if (!include) "zero_deviation_reference" else switch(policy,
        error = "reject", sample = "sample_once_per_group_and_joint_vector",
        zero = "zero_deviation_reference"), unseen_present = count > 0L)
  }
  scale_groups <- any(vapply(spec$scales, function(sm) length(sm$random) > 0L, TRUE))
  response_groups <- quantity != "kernel_scale" && length(spec$response$random) > 0L
  effects <- list(
    scale = group_contract(scale_groups, include_scale, context$scale_policy, context$n_scale),
    response = group_contract(response_groups, include_response,
      context$response_policy, context$n_response),
    dependence = "complete_joint_vectors_and_shared_effects_within_each_call")
  n <- nrow(object$draws$beta)
  complete <- !is.null(n) && all(vapply(names(spec$parameters), function(nm) {
    x <- object$draws[[nm]]
    is.matrix(x) && nrow(x) == n && ncol(x) == length(spec$parameters[[nm]])
  }, TRUE))
  info <- object$inference_metadata
  sampled_new <- (include_scale && context$n_scale > 0L && context$scale_policy == "sample") ||
    (response_groups && include_response && context$n_response > 0L && context$response_policy == "sample")
  reference <- (scale_groups && !include_scale) || (response_groups && !include_response) ||
    (context$n_scale > 0L && context$scale_policy == "zero") ||
    (response_groups && context$n_response > 0L && context$response_policy == "zero")
  out <- .msr_prediction_contract("posterior", object$inference, quantity, output_scale,
    population = "joint_posterior_vectors", kernel = "joint_posterior_vectors",
    effects = effects, conditioning = c("observed_training_responses",
      "frozen_training_design_and_normalization", "declared_prediction_covariates"),
    integration = if (is.null(info)) list(recorded = FALSE) else info$integration,
    intervals = if (quantity == "replicate_response") "posterior_predictive_replicates" else
      if (sampled_new) "new_group_predictive_quantiles" else
        if (reference) "posterior_reference_quantiles" else "conditional_posterior_quantiles",
    omitted = if (quantity == "replicate_response") character() else "observation_noise",
    coordinates = "original_stored_parameter_blocks", joint_map = complete,
    diagnostics = object$diagnostics)
  out$integration$prediction_vectors <- if (!is.null(object$importance$resample_ids))
    "stored_joint_resamples" else "stored_joint_draws"
  out
}

# Analytic snapshots bind a design, while uncertainty is selected later by
# predict(). Describe available controls without treating availability as a gate.
.msr_analytic_prediction_contract <- function(object, type) {
  conditional <- identical(object$backend, "conditional_logit")
  out <- .msr_prediction_contract(if (conditional) "conditional_likelihood" else "likelihood",
    object$backend, if (type == "link") "linear_predictor" else
      if (conditional) "choice_probability" else "expected_response", type,
    population = "estimated_coefficients_with_requested_covariance",
    kernel = if (conditional) "fixed_at_estimate" else "joint_or_conditional_by_request",
    effects = list(choice_sets = if (conditional) "supplied_complete_choice_sets" else "not_applicable"),
    conditioning = c("frozen_training_design_and_normalization", "fixed_prediction_rows",
      if (conditional) "supplied_choice_sets" else character()),
    integration = list(uncertainty = "delta_method_when_requested"),
    intervals = if (conditional) "not_provided" else "mean_confidence_when_eligible",
    omitted = c("observation_noise", if (conditional) "kernel_parameter_uncertainty" else character()),
    diagnostics = object$diagnostics,
    fixed_inputs = TRUE)
  out$uncertainty$joint_available <- !conditional && isTRUE(object$uncertainty$available)
  out
}

# Native tables/vectors retain their exact attributes and component meaning.
# Store this description on the snapshot, never on the cached native values.
.msr_native_prediction_contract <- function(object, kind, type, native_args) {
  model <- .analysis_model(object$opt_mod)
  zero <- .msr_native_zero_random_effects(model, native_args, kind)
  selected <- .msr_native_selected_random_effects(model, native_args, kind)
  re_policy <- if (zero) "zero_deviation_reference" else
    if (selected) "selected_conditional_effects" else
    if ("re.form" %in% names(native_args) && is.null(native_args$re.form)) "native_conditional_effects" else
      if (kind == "glmmTMB" && !"re.form" %in% names(native_args)) "native_conditional_effects" else
        "native_re_form_or_default"
  out <- .msr_prediction_contract("native", kind, if (kind == "unmarked" ||
      !type %in% c("response", "link")) "native_component" else
      if (type == "link") "linear_predictor" else "expected_response",
    output_scale = if (kind == "unmarked")
      if (isFALSE(native_args$backTransform)) "native_link" else "native_component_scale" else type,
    population = "native_parameter_uncertainty_conditional_on_kernel_estimate",
    kernel = "fixed_at_estimate",
    effects = list(response = re_policy, new_groups = if (kind == "glmmTMB")
      if (isTRUE(native_args$allow.new.levels)) "native_new_level_handling" else "reject" else
        "native_component_handling"),
    conditioning = c("fixed_prediction_rows_and_native_controls",
      "selected_kernel_parameters", "native_component_and_observation_year_layout"),
    integration = list(random_effects = if (zero) "set_to_zero_not_integrated" else "native_method"),
    intervals = if (kind == "unmarked") "native_conditional_confidence" else "conditional_standard_errors_only",
    omitted = c("kernel_parameter_uncertainty", "observation_noise"),
    fixed_inputs = TRUE, component = type)
  if (selected) out$effects$selection <- paste(deparse(native_args$re.form), collapse = " ")
  if (inherits(model, "gam") || identical(kind, "gam")) {
    out$uncertainty$smoothing_parameters <- if (isTRUE(native_args$unconditional))
      "native_unconditional_requested" else "native_default_conditional"
    out$term_selection <- list(terms = native_args$terms, exclude = native_args$exclude,
      iterms.type = native_args$iterms.type)
  }
  out
}

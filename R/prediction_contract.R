#' Read prediction contracts
#'
#' @description
#' A prediction contract is a list of metadata that travels with a prediction.
#' It names the quantity that was predicted, its scale, how random effects
#' were treated, and which uncertainty the standard errors or intervals
#' include and leave out. Read it when you are unsure what a result
#' describes. A contract describes a calculation. It does not certify that a
#' fit is eligible for prediction or that an interval is calibrated.
#' [multiScaleR-prediction] explains the quantities and uncertainties that the
#' fields name.
#'
#' @name prediction_contract
#' @details
#' # Find the contract
#'
#' | Result | Where the contract is |
#' |:-------|:----------------------|
#' | Snapshot from [prepare_prediction()] | `snapshot$prediction_contract`, a named list with one contract per prepared output type, such as `$response` and `$link` |
#' | Native refit prediction from [predict.multiScaleR()] | `attr(result, "prediction_contract")` |
#' | Raster prediction from [predict.multiScaleR()] | `attr(map, "prediction_contract")`, beside the `projection_contract` attribute |
#' | Effect curves from [plot_marginal_effects()] for `unmarked`, glmmTMB, lme4, multi-component GAM, and Bayesian fits | `attr(curves, "prediction_contract")` |
#' | Bayesian summaries and simulations from [predict.multiScaleR_bayes()], [simulate.multiScaleR_bayes()], and [predict_scales()] | `attr(result, "prediction_contract")` |
#'
#' Point predictions from the analytic backends return plain vectors or
#' lists without a contract. Prepare a snapshot to get one. Raw posterior
#' draw matrices and lists keep their existing attributes, and native
#' numerical outputs inside a snapshot keep their native attributes.
#'
#' # Read the fields
#'
#' * `schema_version`: version of the metadata structure.
#' * `inference_target`, `method`: the posterior, likelihood, conditional
#'   likelihood, or native target, and the label of the computation that
#'   produced it.
#' * `quantity`, `output_scale`: what the numbers are, on the stated scale.
#'   The quantity is an expected response, a linear predictor, a kernel
#'   scale, a contrast between kernel scales, a replicate response, a choice
#'   probability, or a named native component. Native component meanings
#'   remain model specific.
#' * `target_population`: the supplied prediction rows, or the specified
#'   location pairs for private spatial scale contrasts.
#' * `effects`: the policy for known and unseen groups. Zero deviations
#'   (`"zero_deviation_reference"`) are reference values, not integration
#'   over the group distribution. Sampled new effects are shared by each
#'   group's rows within a complete joint vector.
#' * `conditioning`, `integration`: the quantities held fixed and the
#'   recorded treatment of uncertainty. Missing historical solver metadata is
#'   marked as unrecorded rather than inferred.
#' * `uncertainty`: the treatment of the population and kernel parameters,
#'   whether observation noise is included (`observation_noise`), what the
#'   intervals mean, what is omitted (`omitted`), and `calibration`, which is
#'   always `"not_established"`. Analytic snapshots also record whether joint
#'   covariance was available at preparation (`joint_available`), and
#'   `joint_conditional_on` names any kernel parameters at a search bound that
#'   it holds fixed. This does not bypass the covariance checks in
#'   [predict.multiScaleR()].
#' * `original_coordinates`, `joint_parameter_map_complete`: the parameter
#'   representation used by the adapter. For Bayesian predictions,
#'   completeness checks the stored block dimensions against the original
#'   prepared parameter map. This structural check does not validate a
#'   coordinate transformation, reconstruct integrated coefficients, or
#'   establish posterior accuracy.
#' * `diagnostics_status`: the existing fit status when available,
#'   `"passed"`, `"failed"`, or `"unavailable"`. Joint covariance
#'   availability is separate from this status.
#' * `kernel_support`: for native predictions, any recoverable table of
#'   optimization bounds and the existing extraction-window diagnostic. These
#'   flags qualify the selected kernel. They do not add kernel uncertainty,
#'   alter conditional coefficient standard errors, or establish prediction
#'   eligibility.
#' * `metadata_establishes_eligibility`: always `FALSE`. Read the original
#'   fit diagnostics before you interpret predictions.
#' * `fixed_inputs`, `component`: whether the rows and controls are bound to
#'   a snapshot, and the named component when one applies.
#'
#' # Read the uncertainty record
#'
#' Expected-response intervals exclude observation noise, and replicates
#' include it. New-group predictions also carry the group variation you
#' requested. Native confidence intervals and standard errors condition on
#' the selected kernel parameters. Likelihood summaries conditional on one
#' population parameter vector omit its estimation uncertainty.
#'
#' * Native point and raster contracts record `uncertainty$standard_errors`
#'   as `"provided"` or `"not_requested"`, and separately keep
#'   `se.fit_requested`. Provider controls such as glmmTMB's `cov.fit` can
#'   return standard errors even when the wrapper's `se.fit` is `FALSE`.
#'   Returned native covariance is conditional on the selected kernel
#'   parameters.
#' * spaMM variance attributes keep their native meaning and are listed in
#'   `native_variances`, with their numerical attributes unchanged.
#'   `fixefVar` and `predVar` are on the linear-predictor scale. `residVar`
#'   and `respVar` are on the response scale, and the response variance
#'   includes observation noise, which is distinct from uncertainty in the
#'   mean.
#' * Native interval attributes are listed in `native_intervals` and keep the
#'   provider's target. Native `lm` interval columns record
#'   `"native_mean_confidence"` or `"native_prediction_interval"` and the
#'   bound `interval_level`. Only the prediction interval includes
#'   future-observation noise. Both include native parameter uncertainty and
#'   hold the selected kernel fixed.
#' * Native `mgcv::ziplss()` response means combine both linear predictors.
#'   Their standard errors propagate the full native coefficient covariance
#'   through that mean, keeping signed cross-covariance. The contract records
#'   this as `response_standard_error_method`, distinct from the provider's
#'   own response standard-error formula. Kernel uncertainty and count noise
#'   remain omitted, and native smoothing-covariance requests keep their
#'   meaning.
#' * Refit snapshot contracts describe their cached point output, and
#'   `conditional_se_available` records whether a cached standard-error
#'   request succeeded. `unmarked` snapshots describe the retained native
#'   confidence tables. Mixed native snapshots keep both point and standard
#'   error choices, and `"selection_by_request"` identifies that capability.
#'   None of these fields certifies finite standard errors, covariance
#'   eligibility, or calibration.
#'
#' # Random effects in unmarked contracts
#'
#' `unmarked` treats any non-`NULL` `re.form` as zero random effects. The
#' contract records zero deviations and no integration for that request,
#' including formula controls that other providers use to select some
#' effects. Native argument abbreviations are resolved before this policy is
#' recorded. Zero-effect contexts do not evaluate unused random
#' transformations or add latent variation. Conditional requests still need
#' their original random-design metadata.
#'
#' # Choice sets and selection scores
#'
#' Conditional logistic linear predictors and exponentiated relative
#' selection scores need no prediction choice sets. Choice probabilities need
#' complete sets. A snapshot can keep the reason probabilities are
#' unavailable while its link and risk outputs remain available. Relative
#' scores are not occurrence probabilities or utilization distributions.
#'
#' # Spatial response fields
#'
#' Native glmmTMB response and link predictions with recognized exponential,
#' Gaussian, or Matern response covariance add `response_spatial_field`. It
#' records the provider and layout status, the conditional-response
#' component, and each field's native term, covariance, coordinate and
#' grouping variable names, number of locations, independent fields, and
#' covariance-parameter count. Conditioning records either the native
#' spatial prediction or all response random effects set to zero.
#' New-location controls keep their native meaning. Coordinate units are not
#' inferred from factor labels. The field range is distinct from the
#' ecological kernel scale, and kernel uncertainty is not included.
#' Identification checks the provider's registered covariance codes and the
#' fitted block layout. Coordinate-aware field maps also record `geometry`:
#' the projected coordinate reference system and units, native column names
#' and shared group, reference and context hashes, serial batch size, and
#' memory budget. Cell centers supply the coordinates, and the working-memory
#' estimate is not an operating-system memory cap.
#'
#' Identified native `gdistremoval()` response fields record a fixed
#' correlation range, conditional physical abundance effects, and a distinct
#' residual variance at new locations. Their standard errors keep the full
#' fitted fixed and latent covariance. Density point predictions are
#' inverse-link conditional predictors, not lognormal field averages. These
#' native field calculations do not enable the private RTMB fitting route,
#' field simulation, spatial cross-validation, or ecological-scale
#' uncertainty.
#'
#' Spatial-field cross-validation records a joint assessment-response
#' component. It uses exact Gaussian conditioning for Gaussian responses and
#' the Laplace marginal likelihood ratio for other responses. Both include
#' observation noise and conditional field variation at the fitted training
#' parameters. Log-link response moments use exponentially tilted Laplace
#' ratios, and Bernoulli moments use a Gaussian approximation. The predictive
#' standard deviation describes observation variation. For non-Gaussian
#' responses, the contract identifies field variance on the
#' linear-predictor scale and observation variance on the response scale.
#'
#' # Kernel-scale contrasts
#'
#' Private spatial contrasts use `quantity = "kernel_scale_contrast"` and
#' `output_scale = "log_ratio"` or `"ratio"`. Both locations use the same
#' complete joint vector. Their intervals summarize the contrasts rather
#' than differences between marginal interval endpoints. Observation noise
#' is omitted, and calibration remains unestablished.
#'
#' # Subsets keep the original contract
#'
#' A contract describes the original prediction call. Ordinary data-frame
#' subsetting can keep the attribute, but it does not recompute the contract
#' for a smaller batch. Use the existing row-level new-group flags for
#' individual rows. No contract establishes interval calibration or upgrades
#' a rejected fit.
#'
#' @seealso [multiScaleR-prediction], [prepare_prediction()],
#'   [predict.multiScaleR()], [predict.multiScaleR_bayes()],
#'   [predict_scales()], [simulate.multiScaleR_bayes()], and
#'   `vignette("v09-estimates-and-intervals", package = "multiScaleR")`.
#'
#' @examples
#' ## A stored N-mixture fit used in vignette("v02-user-guide")
#' nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
#'                             package = "multiScaleR"))
#'
#' ## A snapshot stores one contract per prepared component
#' snap <- prepare_prediction(nmix, component = "state")
#' contract <- snap$prediction_contract$state
#' contract$quantity                       # a named native component
#' contract$uncertainty$kernel_parameters  # "fixed_at_estimate": sigma held fixed
#' contract$uncertainty$omitted            # what the limits leave out
#'
#' ## A native refit prediction carries its contract as an attribute
#' abundance <- predict(nmix, type = "state", se.fit = TRUE,
#'                      uncertainty = "conditional")
#' head(abundance)                         # expected abundance with 95% limits
#' attr(abundance, "prediction_contract")$uncertainty$intervals
NULL

# Prediction metadata describes the calculation and its conditioning. It does
# not grant fit eligibility, add uncertainty or change the original diagnostics.
.msr_prediction_contract <- function(target, method, quantity, output_scale,
    population, kernel, effects = list(), conditioning = character(),
    integration = list(), intervals = "not_provided", omitted = character(),
    coordinates = "not_applicable", joint_map = NA, diagnostics = NULL,
    fixed_inputs = FALSE, component = NULL) {
  target <- match.arg(target, c("posterior", "likelihood", "conditional_likelihood", "native"))
  # Admit the new exact label without changing older abbreviation matching.
  if (!identical(quantity, "kernel_scale_contrast"))
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
  response_groups <- !(quantity %in% c("kernel_scale", "kernel_scale_contrast")) &&
    length(spec$response$random) > 0L
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
  choice_probability <- conditional && identical(type, "response")
  relative_score <- conditional && identical(type, "risk")
  out <- .msr_prediction_contract(if (conditional) "conditional_likelihood" else "likelihood",
    object$backend, if (type == "link") "linear_predictor" else
      if (relative_score) "native_component" else
      if (conditional) "choice_probability" else "expected_response", type,
    population = "estimated_coefficients_with_requested_covariance",
    kernel = if (conditional) "fixed_at_estimate" else "joint_or_conditional_by_request",
    effects = list(choice_sets = if (choice_probability) "supplied_complete_choice_sets" else
      if (conditional) "not_required_for_selection_score" else "not_applicable"),
    conditioning = c("frozen_training_design_and_normalization", "fixed_prediction_rows",
      if (choice_probability) "supplied_choice_sets" else character(),
      if (conditional && !choice_probability) "unshifted_formula_offsets" else character()),
    integration = list(uncertainty = "delta_method_when_requested"),
    intervals = if (conditional) "not_provided" else "mean_confidence_when_eligible",
    omitted = c("observation_noise", if (conditional) "kernel_parameter_uncertainty" else character()),
    diagnostics = object$diagnostics,
    fixed_inputs = TRUE)
  if (relative_score) {
    out$component <- "relative_selection_score"
    out$uncertainty$risk_standard_errors <- "delta_method_risk_times_link_SE"
  }
  u <- object$uncertainty
  held <- !isTRUE(u$available) && isTRUE(u$joint_available) && length(u$joint_conditional_on) > 0L
  out$uncertainty$joint_available <- !conditional && (isTRUE(u$available) || held)
  # Record bounded kernel parameters that joint uncertainty holds fixed.
  if (!conditional && held) out$uncertainty$joint_conditional_on <- u$joint_conditional_on
  out
}

# A snapshot can store several uncertainty choices; returned values describe
# only the choice actually requested. This metadata never supplies a missing SE.
.msr_returned_prediction_uncertainty <- function(contract, se.fit = NULL,
    se_returned = se.fit, native_variances = contract$uncertainty$native_variances$components,
    native_intervals = contract$uncertainty$native_intervals,
    covariance = FALSE) {
  contract$uncertainty$se.fit_requested <- se.fit
  contract$uncertainty$standard_errors <- if (is.null(se_returned)) "selection_by_request" else
    if (se_returned) "provided" else if (isTRUE(se.fit)) "not_provided" else "not_requested"
  if (isFALSE(se_returned)) {
    contract$uncertainty$intervals <- "not_provided"
    if (!any(c("fixefVar", "predVar", "respVar") %in%
        c(native_variances, native_intervals)) && !length(setdiff(native_intervals, "residVar"))) {
      contract$uncertainty$population_parameters <- "fixed_at_estimate"
      contract$uncertainty$omitted <- unique(c(contract$uncertainty$omitted,
        "population_parameter_uncertainty"))
    } else contract$uncertainty$population_parameters <-
      if (any(c("confidence", "prediction") %in% native_intervals))
        "native_parameter_uncertainty_conditional_on_kernel_estimate" else
          "native_uncertainty_attributes_conditional_on_kernel_estimate"
  }
  if (length(native_variances)) contract$uncertainty$native_variances <- list(
    components = native_variances,
    linear_predictor_scale = intersect(native_variances, c("fixefVar", "predVar")),
    response_scale = intersect(native_variances, c("residVar", "respVar")))
  if (length(native_intervals)) {
    contract$uncertainty$intervals <- if ("prediction" %in% native_intervals)
      "native_prediction_interval" else if ("confidence" %in% native_intervals)
        "native_mean_confidence" else "native_attribute_intervals"
    contract$uncertainty$native_intervals <- native_intervals
  }
  if (any(c("residVar", "respVar", "prediction") %in% c(native_variances, native_intervals))) {
    contract$uncertainty$observation_noise <- TRUE
    contract$uncertainty$omitted <- setdiff(contract$uncertainty$omitted, "observation_noise")
  }
  if (covariance) contract$uncertainty$covariance <- "native_conditional_covariance"
  contract
}

# Describe native components without changing their numerical values. Raw
# native snapshots keep their original attributes and store contracts separately.
.msr_native_prediction_contract <- function(object, kind, type, native_args, se.fit = NULL, output = NULL) {
  object <- .msr_refit_bound_view(object)
  model <- .analysis_model(object$opt_mod)
  zero <- .msr_native_zero_random_effects(model, native_args, kind)
  selected <- .msr_native_selected_random_effects(model, native_args, kind)
  re_policy <- if (zero) "zero_deviation_reference" else
    if (selected) "selected_conditional_effects" else
    if ("re.form" %in% names(native_args) && is.null(native_args[["re.form"]])) "native_conditional_effects" else
      if (kind == "glmmTMB" && !"re.form" %in% names(native_args)) "native_conditional_effects" else
        "native_re_form_or_default"
  out <- .msr_prediction_contract("native", kind, if (kind == "unmarked" ||
      !type %in% c("response", "link")) "native_component" else
      if (type == "link") "linear_predictor" else "expected_response",
    output_scale = if (kind == "unmarked")
      if (isFALSE(native_args[["backTransform"]])) "native_link" else "native_component_scale" else type,
    population = "native_parameter_uncertainty_conditional_on_kernel_estimate",
    kernel = "fixed_at_estimate",
    effects = list(response = re_policy, new_groups = if (kind == "glmmTMB")
      if (isTRUE(native_args[["allow.new.levels"]])) "native_new_level_handling" else "reject" else
        "native_component_handling"),
    conditioning = c("fixed_prediction_rows_and_native_controls",
      "selected_kernel_parameters", "native_component_and_observation_year_layout"),
    integration = list(random_effects = if (zero) "set_to_zero_not_integrated" else "native_method"),
    intervals = if (kind == "unmarked") "native_conditional_confidence" else "conditional_standard_errors_only",
    omitted = c("kernel_parameter_uncertainty", "observation_noise"),
    fixed_inputs = TRUE, component = type)
  if (methods::is(model, "unmarkedFitGDRField")) {
    out$response_spatial_field <- .msr_gdr_field_contract(model, zero)
    out$uncertainty$population_parameters <- "joint_native_coefficients_and_response_field_conditional_on_kernel_and_fixed_range"
  }
  if (kind == "glmmTMB" && type %in% c("response", "link", "conditional")) {
    fields <- .msr_glmmTMB_response_fields(model, zero,
      isTRUE(native_args[["allow.new.levels"]]))
    if (!is.null(fields)) out$response_spatial_field <- fields
  }
  out$kernel_support <- list(optimization_bounds = object$diagnostics$refit_hessian$bounds,
    extraction_window = object$diagnostics$max_distance,
    establishes_eligibility = FALSE)
  if (selected) out$effects$selection <- paste(deparse(native_args[["re.form"]]), collapse = " ")
  if (inherits(model, "gam") || identical(kind, "gam")) {
    if (is.list(model$pterms) && !inherits(model$pterms, "terms")) {
      # Multiple response columns need their family's native interpretation.
      if (type == "response") out$quantity <- "native_component"
      out$native_components <- list(provider = "mgcv", family = model$family$family,
        linear_predictors = length(model$pterms), interpretation = "native_family_specific")
    }
    out$uncertainty$smoothing_parameters <- if (isTRUE(native_args[["unconditional"]]))
      "native_unconditional_requested" else "native_default_conditional"
    if (type == "response" && .msr_gam_ziplss_model(model)) {
      out$native_components$response_outputs <- 1L
      out$native_components$response_quantities <- "overall_expected_count"
      out$uncertainty$response_standard_error_method <- "delta_full_native_coefficient_covariance"
      out$uncertainty$native_provider_response_SE_used <- FALSE
    }
    out$term_selection <- list(terms = native_args[["terms"]], exclude = native_args[["exclude"]],
      iterms.type = native_args[["iterms.type"]])
    if (inherits(model, "bam") && "discrete" %in% names(native_args))
      out$term_selection$discrete_requested <- native_args[["discrete"]]
  }
  provided <- if (is.null(output)) se.fit else if (kind == "unmarked")
    is.data.frame(output) && "SE" %in% names(output) else
      is.list(output) && !is.null(output$se.fit)
  native_variances <- native_intervals <- character()
  if (inherits(model, "HLfit") && !is.null(output)) {
    native_variances <- Filter(function(nm) length(attr(output, nm)) > 0L,
      c("fixefVar", "predVar", "residVar", "respVar"))
    limits <- attr(output, "intervals")
    if (!is.null(limits)) native_intervals <- if (length(names(limits))) names(limits) else
      if (length(native_args[["intervals"]])) native_args[["intervals"]] else "native"
  }
  interval <- .msr_native_lm_interval(model, native_args, type)
  bounds <- if (is.list(output) && !is.null(output$fit)) output$fit else output
  if (interval != "none" && (is.null(output) ||
      (is.matrix(bounds) && all(c("fit", "lwr", "upr") %in% colnames(bounds))))) {
    native_intervals <- interval
    out$uncertainty$interval_level <- if (is.null(native_args[["level"]])) .95 else native_args[["level"]]
  }
  .msr_returned_prediction_uncertainty(out, se.fit, se_returned = provided,
    native_variances = native_variances, native_intervals = native_intervals,
    covariance = is.list(output) && !is.null(output$cov.fit))
}

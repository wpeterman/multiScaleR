#' Prepare a fixed set of multiscale predictions for repeated use
#'
#' @description Prepare a spatial prediction design or native prediction output
#'   once and retain a compact snapshot for repeated summaries of fixed rows.
#' @param object A fitted \code{multiScaleR} object using \code{profile_glm},
#'   \code{joint_nb2}, \code{conditional_logit}, or \code{joint_mixed}, or
#'   containing a supported native-refit model in \code{opt_mod}.
#' @param newdata,kernel_inputs Prediction data and matching prepared spatial
#'   inputs, with the same requirements as \code{\link{predict.multiScaleR}}.
#'   The default \code{NULL} for both uses the retained fitting rows.
#' @param max_bytes A positive finite byte limit for the returned snapshot.
#'   The default is 64 MiB. Exceeding the limit raises an error. This limit
#'   applies to retained size, not temporary memory used during preparation.
#' @param component For unmarked models, the required native prediction
#'   \code{type}, such as \code{"state"}, \code{"psi"}, \code{"det"}, \code{"col"}, or
#'   \code{"ext"}, as supported by the fitted model. For \code{joint_mixed},
#'   optionally select one native glmmTMB prediction type, such as
#'   \code{"zprob"}; \code{NULL} prepares both response and link predictions.
#'   Leave \code{NULL} for the analytic GLM and conditional-logistic backends.
#'   General refit snapshots default to response and link for GLM-like models,
#'   or linear predictor and risk for Cox models. Other native providers
#'   require an explicit component.
#' @param native_args A named list of native prediction controls to bind during
#'   preparation. For \code{joint_mixed}, supports \code{re.form},
#'   \code{allow.new.levels}, and \code{na.action}. For unmarked, supports
#'   \code{backTransform}, \code{na.rm}, \code{appendData}, \code{level}, and
#'   \code{re.form}. Defaults to an empty list, preserving native defaults.
#'   General refits additionally support \code{reference} and nlme
#'   \code{level}, when accepted by the native method. For mgcv GAMs,
#'   also supports \code{exclude}, \code{terms}, \code{unconditional},
#'   \code{iterms.type}, and \code{block.size}. The controls and native
#'   warnings are bound to the snapshot. Unconditional GAM SEs request native
#'   smoothing-parameter covariance; they do not add kernel-scale uncertainty.
#'   Ordinary lm refits also accept native \code{interval} and \code{level}
#'   controls for mean confidence or future-observation prediction intervals.
#'   spaMM refits accept \code{variances} and \code{intervals}; cross-row
#'   covariance attributes require point predictions instead of a snapshot.
#'   GAM rows retain all formula-input requirements, including excluded-term
#'   inputs. The unchecked \code{newdata.guaranteed} shortcut is unsupported.
#'   Selected lme4/spaMM effects and GAM terms are recorded in prediction contracts.
#'   Supported native S3 control abbreviations are resolved to their full names
#'   before checking the snapshot allowlist and are recorded by those names.
#'   spaMM refits use retained native frames to verify fixed features; zero or
#'   selected effects follow the provider's row requirements. Unsupported native
#'   SE requests stay explicit errors and are recorded separately from point output.
#' @param join_by For native-refit spatial predictions, explicit site keys with
#'   one row per prepared prediction site, in spatial input order. Requires
#'   both \code{newdata} and \code{kernel_inputs}. See
#'   \code{\link{prepare_unmarked_data}}. The default \code{NULL} uses strict
#'   row-name alignment for general refits. Unmarked spatial snapshots require
#'   explicit keys; without spatial inputs their native data pass unchanged.
#'   Unsupported for analytic and
#'   joint mixed-model snapshots.
#' @param include_fixed Logical. The default \code{FALSE} keeps fixed-radius
#'   predictors explicit. For supported spatial snapshots, set \code{TRUE} to
#'   use the verified fixed-radius preparation in
#'   \code{\link{prepare_unmarked_data}}. Requires spatial inputs and verified
#'   row alignment, using explicit site keys for repeated sites and unmarked.
#'   Analytic GLM/NB2, conditional-logistic, and joint mixed snapshots also
#'   support this verification, with their existing strict row-ID rules.
#' @details
#' The snapshot binds coefficients, covariance, spatial scales, training
#' standardization, row order, and any choice sets to the supplied model and
#' inputs. It does not retain raw spatial buffers or the fitted native model.
#' Later changes to the original model or data do not update the snapshot.
#' Call \code{prepare_prediction()} again after changing the model, locations,
#' choice sets, metric definitions, or distance units. New data cannot be
#' supplied to \code{predict()} on a snapshot.
#'
#' Reuse avoids rebuilding the spatial engine when requesting several output
#' types for the same rows. Uncertainty restrictions remain those of the fitted
#' model. Conditional-logistic snapshots without valid choice-set IDs support
#' linear predictors and relative scores, but reject response probabilities and
#' retain the reason in their response contract. In particular, their SEs require
#' \code{uncertainty = "conditional"}; scale uncertainty is not added.
#' Snapshots can be saved with \code{saveRDS()} and restored with \code{readRDS()}.
#'
#' Mixed and unmarked snapshots cache the native numerical prediction outputs
#' for the fixed controls, rather than retaining a native model or rebuilding
#' its prediction machinery. Native warnings are recorded and repeated when
#' the corresponding cached output is requested. Native errors abort preparation.
#' Mixed standard errors remain conditional on selected scales. Changing
#' random-effect settings requires a new snapshot. Zero random-effect controls
#' use formula syntax, independently of its environment; grouping columns may
#' be omitted when the native provider permits zero-effect predictions. Native glmmTMB \code{re.form = NA}
#' sets random effects to zero and does not integrate over their distribution.
#'
#' For unmarked without \code{kernel_inputs}, \code{newdata} is passed unchanged
#' to native prediction. Supply spatial covariates on the fitted standardized
#' scale. With \code{kernel_inputs}, supply explicit \code{join_by} keys and
#' native prediction rows. \code{\link{prepare_unmarked_data}} computes spatial
#' columns and applies training standardization without changing native
#' observation/year ordering or nonspatial factor/offset inputs. Spatial
#' snapshots retain prediction row names and site keys in their provenance.
#'
#' General refit snapshots also cache numerical outputs without retaining the
#' fitted native model, formula environments, or raw cell buffers. If a native
#' provider lacks standard errors, point predictions remain available and a
#' later SE request replays the recorded error. Errors in point predictions
#' still abort preparation. These snapshots retain native uncertainty
#' distinctions and do not add kernel-scale uncertainty. Native prediction
#' intervals and response variances retain the provider's observation noise.
#' Cached native intervals or variance attributes require
#' \code{uncertainty = "conditional"}, including when \code{se.fit = FALSE}.
#' Native controls remain bound to the snapshot; leave the wrapper's
#' \code{interval = "none"} when retrieving those native outputs.
#' The returned
#' table preserves native columns, intervals, missing-value handling, and
#' component interpretation. Native standard errors and intervals condition on
#' selected spatial scales. Their native confidence level and transformation
#' are bound by \code{native_args}, not changed at prediction time.
#' For dynamic occupancy, \code{"psi"} describes initial occupancy and
#' \code{"col"}/\code{"ext"} describe transitions. These native component
#' predictions are distinct from smoothed or projected occupancy trajectories.
#' @return An object of class \code{multiScaleR_prediction}. Its \code{predict()}
#'   method returns the same output as \code{predict()} on the original model
#'   and fixed prediction inputs. \code{max_bytes} counts its R object size.
#'   Unmarked snapshots return the native \code{predict(object$opt_mod, ...)}
#'   table, including its \code{Predicted}, \code{SE}, \code{lower}, and
#'   \code{upper} columns. Mixed snapshots return native glmmTMB vectors or
#'   standard-error lists for the bound prediction type and controls.
#'   The named \code{prediction_contract} list describes each prepared output
#'   type: its target, conditioning, random-effect policy, and uncertainty
#'   restrictions. This descriptive metadata does not establish fit eligibility
#'   or interval calibration. Numerical prediction outputs retain their usual
#'   structure, including native attributes.
#' @examples
#' \dontrun{
#' # Prepare once for several summaries at the same locations.
#' prepared <- prepare_prediction(fit, newdata = sites, kernel_inputs = inputs)
#' predict(prepared)
#' predict(prepared, type = "link", se.fit = TRUE,
#'         uncertainty = "conditional")
#'
#' # Bind a native zero-inflation probability with random effects set to zero.
#' zero_prob <- prepare_prediction(mixed_fit, component = "zprob",
#'   native_args = list(re.form = NA))
#' predict(zero_prob)
#'
#' # Preserve native unmarked state predictions and their conditional intervals.
#' state <- prepare_prediction(unmarked_fit, component = "state")
#' predict(state)
#' }
#' @seealso \code{\link{predict.multiScaleR}}, \code{\link{prediction_contract}}
#' @export
prepare_prediction <- function(object, newdata = NULL, kernel_inputs = NULL,
                               max_bytes = 64 * 1024^2,
                               component = NULL, native_args = list(), join_by = NULL,
                               include_fixed = FALSE) {
  validate_scalar_numeric(max_bytes, "max_bytes", positive = TRUE)
  validate_scalar_logical(include_fixed, "include_fixed")
  refit <- inherits(object, "multiScaleR") &&
    !.msr_analytic_backend(object) &&
    !isTRUE(object$backend %in% c("joint_mixed", "conditional_logit"))
  if (include_fixed && !refit && !methods::is(object$opt_mod, "unmarkedFit") &&
      !isTRUE(object$backend %in% c("profile_glm", "joint_nb2", "conditional_logit", "joint_mixed")))
    stop("include_fixed requires a supported spatial prediction model.", call. = FALSE)
  if (!is.null(join_by) && !refit && !methods::is(object$opt_mod, "unmarkedFit"))
    stop("join_by applies only to unmarked spatial predictions.", call. = FALSE)
  if (isTRUE(object$backend %in% c("profile_glm", "joint_nb2", "conditional_logit")))
    newdata <- .msr_fixed_point_data(object, newdata, kernel_inputs, include_fixed)
  if (inherits(object, "multiScaleR") &&
      (identical(object$backend, "joint_mixed") || methods::is(object$opt_mod, "unmarkedFit")))
    return(.msr_prepare_native_prediction(object, newdata, kernel_inputs,
                                         max_bytes, component, native_args, join_by,
                                         include_fixed))
  if (refit)
    return(.msr_prepare_refit_prediction(object, newdata, kernel_inputs,
      max_bytes, component, native_args, join_by, include_fixed))
  if (!is.null(component) || length(native_args))
    stop("component and native_args apply only to native mixed or unmarked snapshots.", call. = FALSE)
  if (!inherits(object, "multiScaleR") ||
      !isTRUE(object$backend %in% c("profile_glm", "joint_nb2", "conditional_logit")))
    stop("Preparation requires a profile_glm, joint_nb2, or conditional_logit fit.", call. = FALSE)
  contract <- object$profile_fit$contract
  groups <- contract$groups
  choice_groups_error <- NULL
  if (identical(object$backend, "conditional_logit")) {
    groups <- tryCatch(.msr_clogit_prediction_groups(object, newdata, kernel_inputs),
      error = function(e) { choice_groups_error <<- conditionMessage(e); NULL })
  }
  design <- .msr_profile_prediction_design(object, newdata, kernel_inputs)
  out <- list(backend = object$backend, design = design,
    profile_fit = list(contract = list(covariates = contract$covariates,
      groups = groups, strata_name = contract$strata_name)),
    kernel_inputs = list(kernel = object$kernel_inputs$kernel),
    provenance = list(scale_parameters = object$profile_fit$par,
      include_fixed = include_fixed,
      distance_unit_conversion = object$kernel_inputs$unit_conv,
      surface_specs = contract$objective$engine$surface_specs))
  if (identical(object$backend, "conditional_logit")) {
    model <- .analysis_model(object$opt_mod)
    out$coefficients <- stats::coef(model)[colnames(design$X)]
    out$conditional_vcov <- stats::vcov(model)[colnames(design$X), colnames(design$X), drop = FALSE]
  } else {
    out$uncertainty <- object$uncertainty
    out$opt_mod <- list(family = object$opt_mod$family)
  }
  types <- c("response", "link", if (identical(object$backend, "conditional_logit")) "risk")
  out$prediction_contract <- stats::setNames(lapply(types,
    function(type) .msr_analytic_prediction_contract(object, type)), types)
  if (!is.null(choice_groups_error)) {
    out$choice_groups_error <- choice_groups_error
    out$prediction_contract$response$available <- FALSE
    out$prediction_contract$response$unavailable_reason <- choice_groups_error
  }
  class(out) <- "multiScaleR_prediction"
  if (as.numeric(utils::object.size(out)) > max_bytes)
    stop("Prepared predictions exceed max_bytes; use fewer rows or an explicit larger budget.", call. = FALSE)
  out
}

#' @rdname prepare_prediction
#' @param type,se.fit,uncertainty,interval,level Output and uncertainty controls
#'   passed to \code{\link{predict.multiScaleR}}, with the same defaults and
#'   safeguards. They change the summary of the fixed prediction design.
#' @param ... Reserved. Additional arguments, including new prediction data,
#'   raise an error; create a new snapshot for changed data.
#' @exportS3Method stats::predict
predict.multiScaleR_prediction <- function(object, type = c("response", "link"),
    se.fit = FALSE, uncertainty = c("joint", "conditional"),
    interval = c("none", "confidence"), level = .95, ...) {
  if (length(list(...)))
    stop("Prepared predictions have fixed inputs; call prepare_prediction() again for new data.", call. = FALSE)
  if (!is.null(object$native_kind)) {
    if (identical(object$native_kind, "unmarked") &&
        (!missing(se.fit) || !missing(uncertainty) || !missing(interval) || !missing(level)))
      stop("Unmarked snapshots return native prediction tables; bind native controls in prepare_prediction().", call. = FALSE)
    return(.msr_predict_native_snapshot(object,
      if (missing(type)) object$default_type else type, se.fit, uncertainty, interval))
  }
  predict.multiScaleR(object, type = type, se.fit = se.fit,
    uncertainty = uncertainty, interval = interval, level = level)
}

#' @rdname prepare_prediction
#' @param x A prepared prediction snapshot.
#' @export
print.multiScaleR_prediction <- function(x, ...) {
  cat("Prepared multiscale predictions (", x$backend, "): ",
      if (is.null(x$native_kind)) nrow(x$design$X) else x$n_rows, " fixed predictions, ",
      format(utils::object.size(x), units = "auto"), ".\n", sep = "")
  invisible(x)
}

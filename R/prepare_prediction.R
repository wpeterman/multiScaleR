#' Prepare a fixed set of predictions for repeated use
#'
#' @description
#' Compute the predictions for one fixed set of rows once and keep them in a
#' compact snapshot. For the analytic and conditional logistic backends, the
#' snapshot holds the prepared spatial design, so `predict()` on it can switch
#' between scales and uncertainty choices without rebuilding the kernels. For
#' joint mixed, `unmarked`, and other native refit models, it holds the native
#' numerical outputs for the controls you bind. A snapshot does not add
#' uncertainty, update the fit, or predict new locations.
#'
#' @param object A fitted `multiScaleR` object from [multiScale_optim()], using
#'   `backend = "profile_glm"`, `"joint_nb2"`, `"conditional_logit"`, or
#'   `"joint_mixed"`, or the default refit backend with a supported native
#'   model in `object$opt_mod`. For `predict()`, a snapshot returned by
#'   `prepare_prediction()`. A Bayesian fit from [multiScale_bayes()] is
#'   refused; use [predict.multiScaleR_bayes()] with new kernel inputs.
#' @param newdata,kernel_inputs Prediction rows and matching spatial inputs
#'   from [kernel_prep()], with the same requirements as in
#'   [predict.multiScaleR()]. The default `NULL` for both uses the rows
#'   retained from fitting.
#' @param max_bytes Positive number. The largest allowed size of the returned
#'   snapshot, in bytes. The default, `64 * 1024^2`, is 64 MiB. A larger
#'   snapshot raises an error, so use fewer rows or raise the limit. It limits
#'   the stored size, not the temporary memory used while preparing.
#' @param component Character scalar naming the native prediction type to
#'   keep.
#'   * `unmarked` models: required. Name one component, such as `"state"`,
#'     `"det"`, `"psi"`, `"col"`, or `"ext"`, as the fitted model supports.
#'   * `"joint_mixed"`: optional. Name one native glmmTMB type, such as
#'     `"zprob"`. The default `NULL` prepares both `"response"` and `"link"`.
#'   * Other refits: the default `NULL` prepares `"response"` and `"link"` for
#'     `lm`, `glm`, glmmTMB, lme4, `gls`, and `lme` models, and `"lp"` and
#'     `"risk"` for Cox models. Other native providers need an explicit
#'     component.
#'   * Analytic and conditional logistic backends: leave `NULL`. Their
#'     snapshots keep every type (`"response"`, `"link"`, and, for
#'     conditional logistic fits, `"risk"`).
#' @param native_args A named list of native prediction controls to bind to
#'   the snapshot. Default `list()`, which keeps the native defaults. Changing
#'   any control later needs a new snapshot. Supported controls:
#'   * `"joint_mixed"`: `re.form`, `allow.new.levels`, and `na.action`.
#'   * `unmarked`: `backTransform`, `na.rm`, `appendData`, `level`, and
#'     `re.form`.
#'   * Other refits: `re.form`, `allow.new.levels`, `na.action`, `reference`,
#'     and `level` (such as nlme `level = 0`), when the native method accepts
#'     them. Ordinary `lm` refits also accept `interval` (`"confidence"` for
#'     mean limits, `"prediction"` for future-observation limits). spaMM
#'     refits accept `variances` and `intervals`, but cross-row covariance
#'     needs point predictions instead of a snapshot. mgcv GAMs accept
#'     `exclude`, `terms`, `unconditional`, `iterms.type`, and `block.size`,
#'     and BAMs also accept `discrete`. `unconditional = TRUE` requests the
#'     native smoothing-parameter covariance, which adds no kernel-scale
#'     uncertainty.
#'
#'   Native argument abbreviations are resolved to their full names before
#'   the list is checked and are recorded by those names. Selected lme4 and
#'   spaMM effects and GAM terms are recorded in the prediction contracts.
#' @param join_by For native refit and `unmarked` spatial predictions, a data
#'   frame of site keys with one row per prepared site, in `kernel_inputs`
#'   order. It needs both `newdata` and `kernel_inputs`. `unmarked` spatial
#'   snapshots require it (see [prepare_unmarked_data()]). The default `NULL`
#'   aligns general refit rows by row name, and `unmarked` data without
#'   spatial inputs pass unchanged. Not used by the analytic, conditional
#'   logistic, and joint mixed backends.
#' @param include_fixed Logical. The default `FALSE` keeps fixed-radius
#'   predictors as explicit inputs in `newdata`. Set `TRUE` to recompute them
#'   from the spatial inputs after verifying their training preparation, as
#'   described in [prepare_unmarked_data()] and [predict.multiScaleR()]. It
#'   needs spatial inputs and verified row alignment, with explicit site keys
#'   for repeated sites and for `unmarked` models.
#'
#' @details
#' # What a snapshot keeps
#'
#' A snapshot binds the coefficients, covariance, selected scales, training
#' standardization, row order, choice sets, and native controls of one call.
#' It keeps neither raw spatial buffers nor the fitted native model, so it
#' can be saved with [saveRDS()] and restored with [readRDS()]. Later changes
#' to the original model or data do not update it. Call
#' `prepare_prediction()` again after changing the model, locations, choice
#' sets, metric definitions, distance units, or native controls. New data
#' cannot be passed to `predict()` on a snapshot.
#'
#' Supported native and analytic models keep their trained standard scaling,
#' polynomial, and spline calls, including those inside nested expressions
#' and formula offsets. Recovery needs one unambiguous learner's attributes
#' per evaluated predictor. Missing or ambiguous metadata raises an error
#' rather than learning a transformation from the prediction rows. Ordinary
#' native `lm`, `glm`, and time-fixed unpenalized Cox snapshots use the
#' original fitted terms for dot shortcuts and active variables, so removed
#' variables add no prediction inputs and new data columns never enter the
#' design.
#'
#' # Analytic and conditional logistic snapshots
#'
#' These snapshots keep the spatial design, so `predict()` on them accepts
#' the same `type`, `se.fit`, `uncertainty`, `interval`, and `level` choices
#' as [predict.multiScaleR()], with the same restrictions. Reuse avoids
#' rebuilding the spatial engine when you need several summaries of the same
#' rows. Conditional logistic snapshots without valid choice-set IDs still
#' give linear predictors and relative scores, but they reject response
#' probabilities and keep the reason in their response contract. Their
#' standard errors need `uncertainty = "conditional"`, and no scale
#' uncertainty is added.
#'
#' # Native snapshots
#'
#' Joint mixed, `unmarked`, and other refit snapshots cache the native
#' numerical outputs for the bound controls, rather than keeping a native
#' model or rebuilding its prediction machinery. Native warnings are recorded
#' and repeated each time the cached output is requested. A native error in
#' the point predictions stops preparation.
#'
#' * **Joint mixed models.** Point and standard-error outputs are both cached.
#'   Standard errors stay conditional on the selected scales. Zero
#'   random-effect controls use formula syntax independently of its
#'   environment, and grouping columns may be omitted when the provider
#'   permits zero-effect predictions. glmmTMB `re.form = NA` sets random
#'   effects to zero and does not integrate over their distribution.
#' * **General refits.** If the provider cannot supply standard errors, the
#'   point predictions remain available and a later standard-error request
#'   repeats the recorded error. Native prediction intervals and response
#'   variances keep the provider's observation noise. Cached native intervals
#'   or variance attributes need `uncertainty = "conditional"`, even with
#'   `se.fit = FALSE`. Leave the wrapper's `interval = "none"` when you
#'   retrieve those native outputs. GAM rows keep every formula-input
#'   requirement, including inputs to excluded terms, and the unchecked
#'   `newdata.guaranteed` shortcut is not supported. spaMM refits use their
#'   retained native frames to verify fixed features. An unsupported native
#'   standard-error request stays an explicit error, recorded separately from
#'   the point output.
#' * **`unmarked` models.** `predict()` returns the cached native table. Its
#'   confidence level and back-transformation are fixed by `native_args` and
#'   cannot be changed at prediction time.
#'
#' # Unmarked snapshots
#'
#' Without `kernel_inputs`, supply native prediction rows with spatial
#' covariates already on their fitted scale, including raw units when the
#' native refit used them. With `kernel_inputs`, supply explicit `join_by`
#' keys and native prediction rows. [prepare_unmarked_data()] then computes
#' the spatial columns and applies the training standardization without
#' changing native observation or year ordering or non-spatial factor and
#' offset inputs. Spatial snapshots keep the prediction row names and site
#' keys in their provenance.
#'
#' Ordinary occupancy and `pcount` state and detection terms use the
#' original raw training transformations, fitted factor contrasts, and native
#' component covariance. Missing contrasts can be recovered from unambiguous
#' original main-factor columns. Other missing coding needs explicit training
#' contrasts. Preparation verifies the fixed designs and offsets against the
#' retained native fit. Explicit data-frame rows prepare only the requested
#' component, while native frame inputs or an omitted `newdata` keep both
#' components for native missing-row removal. ZIP outputs keep their native
#' zero-inflation estimate and uncertainty, and every output keeps the full
#' model's parameter count and diagnostics. Keep the original covariates and
#' fitted designs. Frame subsets keep the original known-group columns.
#'
#' `gdistremoval()` snapshots likewise freeze the active `"lambda"`,
#' `"phi"`, `"dist"`, or `"rem"` component's original training terms and
#' native covariance, with distinct site, primary-period, and
#' removal-interval recipes. Data-frame rows limit preparation to the
#' requested component, and native GDR frames keep every active design and
#' their original missing-input checks. No model or executable TMB tape is
#' saved in the snapshot. GDR models with active formula offsets or TMB zero
#' inflation are rejected before fitting, model selection, or prediction, as
#' [aic_tab()] explains.
#'
#' Conditional learned random-effect terms need precomputed training
#' predictors. Zero-effect contexts need only their fixed design and
#' offsets, and they reject a later conditional-effect request. Any
#' non-`NULL` native `re.form` sets the effects to zero. It neither selects
#' some effects nor integrates over them. For dynamic occupancy, `"psi"`
#' describes initial occupancy, and `"col"` and `"ext"` describe transitions.
#' These component predictions are distinct from smoothed or projected
#' occupancy trajectories.
#'
#' @return `prepare_prediction()` returns an object of class
#'   `multiScaleR_prediction`, printed as a one-line summary of the backend,
#'   the number of fixed predictions, and the stored size. Its elements are
#'   internal, except `prediction_contract`, a named list with one
#'   [prediction_contract] per prepared output type. Each records the
#'   quantity, its scale and conditioning, the random-effect policy, and the
#'   uncertainty restrictions. This metadata does not establish fit
#'   eligibility or interval calibration.
#'
#'   `predict()` on a snapshot returns the same output as
#'   [predict.multiScaleR()] would return for the original model and fixed
#'   rows. `unmarked` snapshots return the native table with `Predicted`,
#'   `SE`, `lower`, and `upper` columns. Joint mixed snapshots return native
#'   glmmTMB vectors or standard-error lists for the bound type and controls.
#'   Numerical outputs keep their usual structure, including native
#'   attributes.
#'
#'   `print()` returns the snapshot invisibly.
#'
#' @seealso [predict.multiScaleR()], [prediction_contract],
#'   [multiScaleR-prediction], [prepare_unmarked_data()]
#'
#' @examples
#' ## Analytic fit from the running example: counts and two landscape layers
#' data("landscape_counts")
#' data("surv_pts")
#' pts <- terra::vect(surv_pts)
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))[[c("land1", "land2")]]
#' inputs <- kernel_prep(pts, land, max_D = 1700, kernel = "gaussian",
#'                       bin = TRUE, verbose = FALSE)
#' dat <- data.frame(landscape_counts, inputs$kernel_dat)
#' start <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
#' fit <- multiScale_optim(start, inputs, backend = "profile_glm",
#'                         n_cores = 1, verbose = FALSE)
#'
#' ## Prepare three new survey points once. Row names of new_sites match the
#' ## prepared point IDs.
#' new_inputs <- kernel_prep(pts[1:3], land, max_D = 1700, kernel = "gaussian",
#'                           bin = TRUE, verbose = FALSE)
#' new_sites <- data.frame(site = c(-1, 0, 1))
#' prepared <- prepare_prediction(fit, new_sites, new_inputs)
#' prepared
#'
#' ## Several summaries of the same rows, without rebuilding the kernels
#' predict(prepared)                                  # expected counts
#' predict(prepared, type = "link", se.fit = TRUE)    # log scale, joint SEs
#' predict(prepared, interval = "confidence",
#'         uncertainty = "conditional")$upper         # upper 95% limits
#'
#' ## A native unmarked snapshot binds one component. This stored N-mixture
#' ## fit is used in vignette("v02-user-guide").
#' nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
#'                             package = "multiScaleR"))
#' detection <- prepare_prediction(nmix, component = "det")
#' head(predict(detection))   # detection probability per survey, 95% limits
#' @export
prepare_prediction <- function(object, newdata = NULL, kernel_inputs = NULL,
                               max_bytes = 64 * 1024^2,
                               component = NULL, native_args = list(), join_by = NULL,
                               include_fixed = FALSE) {
  # A Bayesian fit gets a pointer to the Bayesian tools, not a class error.
  .msr_stop_bayes_fit(object, "prepare_prediction")
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
#' @param type,se.fit,uncertainty,interval,level Output and uncertainty
#'   controls for `predict()` on a snapshot, with the defaults and checks of
#'   [predict.multiScaleR()]. For analytic and conditional logistic snapshots
#'   they choose the scale and the uncertainty summary of the fixed design.
#'   For joint mixed and refit snapshots, `type` selects one prepared
#'   component (the first one prepared by default), `se.fit` selects the
#'   cached point or standard-error output, and `uncertainty` must be
#'   `"conditional"` when standard errors or native interval attributes are
#'   returned. `unmarked` snapshots accept only `type`.
#' @param ... For `predict()`, any additional argument, including new
#'   prediction data, raises an error. Prepare a new snapshot for new rows.
#'   Ignored by `print()`.
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
#' @param x A `multiScaleR_prediction` snapshot from `prepare_prediction()`.
#' @export
print.multiScaleR_prediction <- function(x, ...) {
  cat("Prepared multiscale predictions (", x$backend, "): ",
      if (is.null(x$native_kind)) nrow(x$design$X) else x$n_rows, " fixed predictions, ",
      format(utils::object.size(x), units = "auto"), ".\n", sep = "")
  invisible(x)
}

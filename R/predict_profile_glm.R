# New locations use the fitting sample's center and spread. Their derivatives
# must also include how those training quantities change with kernel scale.
.msr_profile_prediction_design <- function(object, newdata, kernel_inputs) {
  # Prepared predictions bind the numerical design to one model and row set.
  if (inherits(object, "multiScaleR_prediction")) return(object$design)
  contract <- object$profile_fit$contract
  obj <- contract$objective
  e <- obj$engine
  if (inherits(e, "msr_surface_engine") && isTRUE(e$rebuild_from_inputs))
    e <- .msr_prepare_surface_engine(object$kernel_inputs, e$covariates, e$rows,
                                     e$standardize, e$implementation)
  nk <- length(e$parameter_names)
  kp <- object$profile_fit$par[seq_len(nk)]
  if (is.null(newdata) && is.null(kernel_inputs)) {
    kd <- .msr_eval_kernel_engine(e, kp)
    out <- list(X = if (is.null(obj$response_design)) cbind(obj$X_fixed, kd$value) else
      .msr_response_design_matrix(obj$response_design, kd$value), offset = obj$offset,
      ds = kd$d_scale, dh = kd$d_shape, rows = rownames(contract$frame))
    if (!is.null(obj$response_design) && !identical(object$backend, "conditional_logit"))
      out$kernel_sensitivity <- .msr_response_kernel_sensitivity(obj$response_design, kd$value,
        object$uncertainty$par[object$uncertainty$beta_index])
    return(out)
  }
  if (is.null(kernel_inputs) || (is.null(kernel_inputs$binned) && !inherits(e, "msr_surface_engine")) ||
      !identical(kernel_inputs$kernel, e$kernel))
    stop("New predictions require matching smooth kernel_inputs.", call. = FALSE)
  validate_scalar_numeric(kernel_inputs$unit_conv, "kernel_inputs$unit_conv", positive = TRUE)
  b <- kernel_inputs$binned
  if (is.null(newdata)) {
    ids <- if (is.null(b)) rownames(kernel_inputs$kernel_dat) else b$point_ids
    if (is.null(ids) && is.null(b)) stop("Surface-only predictions require preparation row IDs.", call. = FALSE)
    if (is.null(ids)) ids <- seq_len(nrow(b$vsum[[1L]]))
    newdata <- data.frame(row.names = as.character(ids))
  }
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.", call. = FALSE)
  required <- setdiff(all.vars(stats::delete.response(contract$terms)), e$covariates)
  missing <- setdiff(required, names(newdata))
  if (length(missing))
    stop("newdata is missing formula variables: ", paste(missing, collapse = ", "),
         ". Supply them explicitly; formula-environment fallback is not used.", call. = FALSE)
  rows <- .msr_profile_input_rows(newdata, kernel_inputs)
  # Convert prepared distance units, not fitted scales or sample normalization.
  if (!is.null(b)) b$bin_rep <- b$bin_rep * kernel_inputs$unit_conv / object$kernel_inputs$unit_conv
  if (inherits(e, "msr_surface_engine")) {
    specs <- kernel_inputs$scale_vars
    expected <- e$surface_specs
    if (is.null(specs)) stop("New surface specifications are required.", call. = FALSE)
    actual <- specs[match(expected$covariate, specs$covariate), , drop = FALSE]
    fields <- c("covariate", "source", "type", "metric", "weighted", "optimize")
    if (is.null(specs) || !isTRUE(all.equal(actual[, fields], expected[, fields],
                                           check.attributes = FALSE)))
      stop("New surface source and metric specifications must match the fitted model.", call. = FALSE)
    surface_inputs <- kernel_inputs
    surface_inputs$binned <- b
    surface_inputs$d_list <- lapply(kernel_inputs$d_list, function(d)
      d * kernel_inputs$unit_conv / object$kernel_inputs$unit_conv)
    en <- .msr_prepare_surface_engine(surface_inputs, e$covariates, rows, FALSE,
                                      e$implementation)
  } else en <- .msr_prepare_kernel_engine(b, e$covariates, rows, e$kernel, FALSE,
    implementation = if (is.null(e$implementation)) "reference" else e$implementation)
  er <- e
  er$standardize <- FALSE
  train <- .msr_eval_kernel_engine(er, kp)
  kd <- .msr_eval_kernel_engine(en, kp)
  center <- colMeans(train$raw)
  spread <- apply(train$raw, 2L, stats::sd)
  value <- sweep(sweep(kd$raw, 2L, center), 2L, spread, "/")
  centered <- sweep(train$raw, 2L, center)
  derivative <- function(dnew, dtrain) {
    dc <- colMeans(dtrain)
    dsp <- colSums(centered * sweep(dtrain, 2L, dc)) / ((e$n - 1) * spread)
    sweep(sweep(dnew, 2L, dc), 2L, spread, "/") -
      sweep(value, 2L, dsp / spread, "*")
  }
  ds <- derivative(kd$d_scale, train$d_scale)
  dh <- if (e$kernel == "expow") derivative(kd$d_shape, train$d_shape) else kd$d_shape
  for (nm in e$covariates) newdata[[nm]] <- value[, nm]
  tt <- contract$terms
  # Older full fits can retain incomplete prediction terms while their native
  # evaluated frame still contains the original fixed training transformation.
  if (.msr_needs_fixed_metadata(tt)) {
    native <- .analysis_model(object$opt_mod)
    tt <- .msr_fixed_prediction_terms(tt, .msr_response_frame_subset(
      stats::terms(native), native$model, tt))
  }
  tt <- stats::delete.response(tt)
  mf <- stats::model.frame(tt, data = newdata, xlev = contract$xlevels,
                          na.action = stats::na.fail)
  X <- stats::model.matrix(tt, mf, contrasts.arg = contract$contrasts)
  columns <- if (is.null(obj$response_design)) c(contract$fixed_columns, e$covariates) else
    obj$response_design$columns
  X <- X[, columns, drop = FALSE]
  if (any(!is.finite(X))) stop("New design values must be finite.", call. = FALSE)
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  if (any(!is.finite(offset))) stop("New offsets must be finite.", call. = FALSE)
  out <- list(X = X, offset = offset, ds = ds, dh = dh, rows = rownames(mf))
  if (!is.null(obj$response_design) && !identical(object$backend, "conditional_logit")) {
    response_design <- .msr_compile_response_design(tt, mf, e$covariates, columns, contract$contrasts)
    out$kernel_sensitivity <- .msr_response_kernel_sensitivity(response_design, value,
      object$uncertainty$par[object$uncertainty$beta_index])
  }
  out
}

#' Predict from an analytic multiscale model
#'
#' @description Predict the expected response, with optional uncertainty from
#'   joint estimation of regression coefficients and spatial kernel parameters.
#' @param object A \code{multiScaleR} fit using \code{backend = "profile_glm"}
#'   or \code{backend = "joint_nb2"}, or an experimental
#'   \code{backend = "joint_mixed"} fit.
#'   Also accepts \code{backend = "conditional_logit"} and native-refit fits.
#' @param newdata A data frame containing fixed covariates and formula offset
#'   variables. Row names must match point IDs in \code{kernel_inputs}. Kernel
#'   columns are computed from the prepared data, replacing any supplied values.
#'   Alternatively, supply checked \code{SpatRaster} output from
#'   \code{\link{kernel_scale.raster}}; see Raster prediction below.
#' @param kernel_inputs Prepared inputs from \code{\link{kernel_prep}} for new
#'   locations, with the same kernel and covariate names as the fit. If both
#'   \code{newdata} and \code{kernel_inputs} are omitted, predict the retained
#'   fitting rows. For analytic GLM fits, if only \code{newdata} is omitted,
#'   an empty frame is constructed
#'   from prepared point IDs; this works only without additional covariates.
#'   Conditional-logistic fits with optimized weighted surface SQ also require
#'   matching raw buffers and source/metric specifications in these inputs.
#' @param type Predict on the \code{"response"} or \code{"link"} scale for
#'   analytic fits. Native-refit fits also accept their provider's prediction
#'   components. Their defaults are response, unmarked state, or Cox linear
#'   predictor, as appropriate. Unmarked raster prediction requires an explicit
#'   component. Cox response probabilities are not supplied by this wrapper.
#' @param se.fit Return standard errors of the estimated mean.
#' @param uncertainty \code{"joint"} includes estimated kernel scales, shapes,
#'   and their covariance with coefficients. \code{"conditional"} holds kernel
#'   parameters fixed. Both account for estimated dispersion in the coefficient
#'   information matrix, when applicable.
#' @param interval \code{"none"} or \code{"confidence"}. Confidence limits
#'   describe the expected response, not future observations.
#' @param level Confidence level strictly between zero and one.
#' @param ... Native prediction controls for refit and \code{joint_mixed}
#'   models. For refit points, also accepts \code{join_by} and
#'   \code{include_fixed} for verified spatial preparation. Use a named
#'   \code{native_args} list for controls that share a wrapper argument name,
#'   such as \code{native_args = list(level = 0)} for nlme population predictions.
#'   Native S3 argument abbreviations are resolved before row and uncertainty
#'   checks using the provider's exact-first argument matching. They cannot
#'   replace the wrapper's verified prediction inputs or bypass its guards.
#'   Zero-effect predictions may omit random-only grouping or slope inputs,
#'   while all fixed predictors and offsets remain required. Setting effects
#'   to zero does not integrate over their distribution. For lme4 partial
#'   \code{re.form} formulas, supply the fixed predictors and the grouping/slope
#'   variables in the selected effects. Native validation still checks that the
#'   selected effects belong to the fitted model.
#'   Raster controls go in a named \code{raster_args} list; see Details.
#' @details
#' For \code{conditional_logit}, \code{type = "response"} returns probabilities
#' normalized within each supplied choice set, not occurrence probabilities or
#' a utilization distribution. \code{type = "link"} returns the uncentered
#' linear predictor; differences within a set describe relative selection.
#' The exact type \code{"risk"} returns \code{exp(eta)}, a relative selection
#' score, not an occurrence probability or utilization distribution. New link
#' and risk predictions require matching spatial inputs and all fixed formula
#' inputs; they do not require choice-set IDs. Response predictions additionally
#' require explicit choice-set IDs, at least two alternatives per set, and
#' complete prediction choice sets defined by the caller. \code{se.fit = TRUE}
#' requires \code{uncertainty = "conditional"}, holding scales fixed; intervals
#' and joint uncertainty are unavailable. Formula offsets and training
#' standardization are preserved. Link and risk outputs leave formula offsets
#' unshifted; risk SEs use \code{risk * SE(link)}. Native
#' \code{opt_mod} prediction methods use
#' their own survival conventions and require already standardized covariates.
#'
#' New predictors use the training sample's center and sample standard deviation;
#' they are never standardized within the prediction sample. Fixed transforms,
#' factor contrasts, and offsets follow the fitted terms object. New factor
#' levels, missing values, and unmatched point IDs are errors for the GLM backends.
#' Every non-kernel variable named in the formula must be supplied in
#' \code{newdata}; values are not borrowed from the formula environment.
#' Use literal constants in transform arguments such as \code{poly(z, 2)}.
#'
#' Standard errors use the delta method. Confidence limits are computed on the
#' link scale and then transformed. Joint uncertainty is refused for boundary
#' or weakly identified fits; point predictions remain available. The
#' \code{joint_nb2} exact-Poisson boundary supports point predictions and
#' explicitly conditional intervals, holding both scales and the Poisson
#' dispersion boundary fixed. Joint intervals remain unavailable there.
#' The ordinary
#' \code{object$opt_mod} is a conditional GLM whose own prediction methods expect
#' already standardized kernel predictors and do not include scale uncertainty.
#' @return Without standard errors or intervals, a named numeric vector.
#'   Otherwise a list containing \code{fit}, \code{se.fit}, \code{type}, and
#'   \code{uncertainty}; confidence intervals add \code{lower}, \code{upper},
#'   and \code{level}. Binomial response predictions are probabilities.
#'   For \code{joint_mixed}, returns the native glmmTMB prediction vector or,
#'   with \code{se.fit = TRUE}, its list containing \code{fit} and \code{se.fit}.
#'   Native-refit results follow the provider's output. Computed-raster inputs
#'   return a \code{SpatRaster} of predictions, with SE layers when requested
#'   and available. Unmarked SE output retains the native \code{Predicted},
#'   \code{SE}, \code{lower}, and \code{upper} columns as layers, with native
#'   confidence-level semantics. No observation noise is added.
#' @seealso \code{\link{multiScale_optim}}, \code{\link{vcov.multiScaleR}},
#'   \code{\link{prepare_prediction}} for repeated summaries of fixed rows.
#' @details For \code{joint_mixed}, prediction delegates to the final native
#'   glmmTMB model. Standard errors require \code{uncertainty = "conditional"}
#'   and exclude scale estimation uncertainty; confidence intervals are not
#'   provided by this wrapper. Supply both \code{newdata} and matching
#'   \code{kernel_inputs} for new rows. Include all fixed-effect, offset, and
#'   grouping variables. Kernel means need binned summaries; weighted surface
#'   SQ needs matching source specifications, raw cell buffers, and row IDs.
#'   Spatial columns are recomputed using the
#'   training means and standard deviations. Additional arguments, including
#'   \code{re.form} and \code{allow.new.levels}, follow native glmmTMB behavior.
#'   The default includes fitted random effects; \code{re.form = NA} sets them
#'   to zero and does not integrate over their distribution. Use native
#'   \code{predict(object$opt_mod, ...)} for other glmmTMB prediction types.
#'
#' \strong{Fixed-radius point predictors.} For analytic GLM/NB2, conditional
#' logistic, joint mixed, and native-refit models, \code{include_fixed = TRUE}
#' recomputes declared fixed-radius predictors from matching spatial inputs.
#' Their original prepared numeric columns and means/sample SDs must be verified
#' against the retained native model frame. Source, metric, radius, kernel,
#' resolution, buffer extent (\code{max_D}), and observation-to-site alignment
#' must match. The default keeps
#' these columns explicit. Fixed radii add no scale parameters or uncertainty.
#'
#' \strong{Native-refit points.} Spatial means and metrics are recomputed from
#' matching source definitions, prepared rows, physical distance units, and
#' selected scales. Unused prepared variables are ignored. Training means and
#' sample standard deviations are retained. Supply both \code{newdata} and
#' \code{kernel_inputs}, or omit both for fitting rows. Explicit \code{join_by}
#' keys support repeated site scenarios. Set \code{include_fixed = TRUE} to
#' verify and reproduce modeled fixed-radius predictors from the original
#' preparation. See \code{\link{prepare_unmarked_data}} for the key and fixed
#' normalization requirements. Standard errors require
#' \code{uncertainty = "conditional"} and native support. Unsupported native
#' SE requests raise an error. Wrapper confidence intervals are unavailable.
#' Fixed qualified \code{base::scale()} and \code{stats::poly()} calls use
#' their retained training transformation for supported LM/GLM, single-component
#' GAM, lme4, and glmmTMB predictions. New rows do not relearn the center or
#' spread. Recovery requires an aligned evaluated training frame. Legacy full
#' analytic fits can recover this metadata; numerical-only snapshots must be
#' prepared again. Nested learned calls and additional multi-component layouts
#' require separate preparation.
#' Random-effect controls retain native meaning; setting effects to zero is
#' not integration over the random-effect distribution.
#' For ordinary lm response predictions, \code{native_args = list(interval =
#' "confidence", level = 0.9)} requests native mean confidence limits;
#' \code{interval = "prediction"} within that list requests limits for a future
#' observation, including residual noise. Both require
#' \code{uncertainty = "conditional"} and retain the native fit/lower/upper
#' columns for points and maps. They hold kernel parameters fixed. The wrapper's
#' own \code{interval} remains \code{"none"}.
#'
#' \strong{Raster prediction.} First call \code{kernel_scale.raster()} with
#' the fitted object, \code{scale_center = TRUE}, and
#' \code{match_buffer = TRUE}. Modeled fixed metrics additionally require
#' \code{include_fixed = TRUE}. Supply categorical scenarios through
#' \code{site_values}; logical predictors are restored to their fitted type
#' during prediction. The retained descriptor checks scales, source and metric
#' definitions, normalization, and support against this fit. It does not
#' reconstruct a lost descriptor or establish that supplied maps have the
#' correct scientific provenance. Recompute after saving and reloading cells
#' through file formats that discard R attributes.
#' spaMM refits use retained native model frames for fixed-feature verification.
#' Zero and selected \code{re.form} predictions require only fixed inputs and
#' selected random-effect inputs. Native variance attributes retain their own
#' scales and targets; the generic \code{se.fit} interface remains unavailable
#' if the native method does not supply a standard-error list.
#' Cross-row \code{cov.fit = TRUE} and spaMM variance/interval attributes require
#' point predictions. Raster chunks cannot retain those objects and reject
#' these controls explicitly. Point predictions preserve native covariance
#' outputs, including SEs triggered by native controls.
#' Native covariance, variance and interval controls also require
#' \code{uncertainty = "conditional"}, even when \code{se.fit = FALSE};
#' they cannot add uncertainty in the selected kernel parameters.
#'
#' Raster outputs retain \code{prediction_contract} and
#' \code{projection_contract} attributes, including requested clamping,
#' smoothing, and offset controls. Recompute the preparation after editing cells.
#'
#' Raster SEs are conditional on the selected kernel parameters, with
#' \code{interval = "none"}. Analytic GLM/NB2 maps use the same retained
#' coefficient information as point predictions. Other maps use native
#' prediction. Native matrix components return one band per column; when native
#' SEs are supplied, fit and SE bands retain their corresponding column order.
#' GAM term and basis-matrix outputs keep their native component meaning.
#' Native smoothing-covariance requests do not add kernel-scale uncertainty. Cox maps default to linear predictors; conditional-logistic
#' maps require \code{type = "link"} or \code{"risk"}, since choice
#' probabilities need complete point choice sets. Supply native controls in
#' \code{native_args}. The \code{raster_args} list accepts \code{filename},
#' \code{overwrite}, \code{wopt}, \code{const}, \code{cores}, and \code{na.rm}.
#' Internally created worker clusters are stopped on exit. Rasterization,
#' distance binning, and partial source coverage remain separate approximation
#' choices. A single selected-scale raster is not a Bayesian posterior map.
#' @exportS3Method stats::predict
predict.multiScaleR <- function(object, newdata = NULL, kernel_inputs = NULL,
                                type = c("response", "link"), se.fit = FALSE,
                                uncertainty = c("joint", "conditional"),
                                interval = c("none", "confidence"), level = .95, ...) {
  if (inherits(newdata, "SpatRaster")) {
    if (!is.null(kernel_inputs)) stop("Computed raster predictions do not take kernel_inputs.", call. = FALSE)
    return(.msr_predict_raster(object, newdata, if (missing(type)) NULL else type,
      se.fit, match.arg(uncertainty), match.arg(interval), ...))
  }
  if (identical(object$backend, "conditional_logit")) {
    dots <- list(...)
    include_fixed <- if (is.null(dots$include_fixed)) FALSE else dots$include_fixed
    dots$include_fixed <- NULL
    if (length(dots)) stop("Additional prediction arguments are unsupported for conditional_logit.", call. = FALSE)
    newdata <- .msr_fixed_point_data(object, newdata, kernel_inputs, include_fixed)
    # Admit the exact new name without changing old type abbreviations.
    conditional_type <- if (identical(type, "risk")) "risk" else match.arg(type)
    return(.msr_predict_clogit(object, newdata, kernel_inputs, conditional_type,
      se.fit, match.arg(uncertainty), match.arg(interval)))
  }
  if (identical(object$backend, "joint_mixed"))
    return(.msr_predict_mixed(object, newdata, kernel_inputs, match.arg(type),
      se.fit, match.arg(uncertainty), match.arg(interval), ...))
  if (!.msr_analytic_backend(object))
    return(.msr_predict_refit(object, newdata, kernel_inputs,
      if (missing(type)) NULL else type, se.fit, match.arg(uncertainty),
      match.arg(interval), ...))
  type <- match.arg(type)
  dots <- list(...)
  include_fixed <- if (is.null(dots$include_fixed)) FALSE else dots$include_fixed
  newdata <- .msr_fixed_point_data(object, newdata, kernel_inputs, include_fixed)
  uncertainty <- match.arg(uncertainty)
  interval <- match.arg(interval)
  if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit) ||
      !is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1)
    stop("se.fit must be logical and level must be strictly between zero and one.", call. = FALSE)
  design <- .msr_profile_prediction_design(object, newdata, kernel_inputs)
  info <- object$uncertainty
  beta <- info$par[info$beta_index]
  eta <- drop(design$X %*% beta) + design$offset
  fam <- object$opt_mod$family
  value <- if (type == "link") eta else fam$linkinv(eta)
  names(value) <- design$rows
  if (!se.fit && interval == "none") return(value)
  V <- vcov.multiScaleR(object, type = uncertainty)
  if (uncertainty == "joint") {
    J <- matrix(0, nrow(design$X), nrow(V))
    p <- length(object$profile_fit$contract$covariates)
    if (is.null(design$kernel_sensitivity)) {
      J[, seq_len(p)] <- sweep(design$ds, 2L, utils::tail(beta, p), "*")
      if (object$kernel_inputs$kernel == "expow")
        J[, p + seq_len(p)] <- sweep(design$dh, 2L, utils::tail(beta, p), "*")
    } else {
      J[, seq_len(p)] <- design$ds * design$kernel_sensitivity
      if (object$kernel_inputs$kernel == "expow")
        J[, p + seq_len(p)] <- design$dh * design$kernel_sensitivity
    }
    J[, info$beta_index] <- design$X
  } else J <- design$X
  se_link <- sqrt(pmax(0, rowSums((J %*% V) * J)))
  out <- list(fit = value, se.fit = stats::setNames(if (type == "link") se_link else
    abs(fam$mu.eta(eta)) * se_link, design$rows), type = type, uncertainty = uncertainty)
  if (interval == "confidence") {
    z <- stats::qnorm((1 + level) / 2)
    out$lower <- stats::setNames(if (type == "link") eta - z * se_link else
      fam$linkinv(eta - z * se_link), design$rows)
    out$upper <- stats::setNames(if (type == "link") eta + z * se_link else
      fam$linkinv(eta + z * se_link), design$rows)
    out$level <- level
  }
  out
}

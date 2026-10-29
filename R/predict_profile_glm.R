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
    return(list(X = cbind(obj$X_fixed, kd$value), offset = obj$offset,
                ds = kd$d_scale, dh = kd$d_shape, rows = rownames(contract$frame)))
  }
  if (is.null(kernel_inputs) || is.null(kernel_inputs$binned) ||
      !identical(kernel_inputs$kernel, e$kernel))
    stop("New predictions require matching smooth, binned kernel_inputs.", call. = FALSE)
  validate_scalar_numeric(kernel_inputs$unit_conv, "kernel_inputs$unit_conv", positive = TRUE)
  b <- kernel_inputs$binned
  if (is.null(newdata)) {
    ids <- b$point_ids
    if (is.null(ids)) ids <- seq_len(nrow(b$vsum[[1L]]))
    newdata <- data.frame(row.names = as.character(ids))
  }
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.", call. = FALSE)
  required <- setdiff(all.vars(stats::delete.response(contract$terms)), e$covariates)
  missing <- setdiff(required, names(newdata))
  if (length(missing))
    stop("newdata is missing formula variables: ", paste(missing, collapse = ", "),
         ". Supply them explicitly; formula-environment fallback is not used.", call. = FALSE)
  rows <- .msr_profile_rows(newdata, b)
  # Convert prepared distance units, not fitted scales or sample normalization.
  b$bin_rep <- b$bin_rep * kernel_inputs$unit_conv / object$kernel_inputs$unit_conv
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
  tt <- stats::delete.response(contract$terms)
  mf <- stats::model.frame(tt, data = newdata, xlev = contract$xlevels,
                          na.action = stats::na.fail)
  X <- stats::model.matrix(tt, mf, contrasts.arg = contract$contrasts)
  X <- X[, c(contract$fixed_columns, e$covariates), drop = FALSE]
  if (any(!is.finite(X))) stop("New design values must be finite.", call. = FALSE)
  offset <- stats::model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nrow(X))
  if (any(!is.finite(offset))) stop("New offsets must be finite.", call. = FALSE)
  list(X = X, offset = offset, ds = ds, dh = dh, rows = rownames(mf))
}

#' Predict from an analytic multiscale model
#'
#' @description Predict the expected response, with optional uncertainty from
#'   joint estimation of regression coefficients and spatial kernel parameters.
#' @param object A \code{multiScaleR} fit using \code{backend = "profile_glm"}
#'   or \code{backend = "joint_nb2"}, or an experimental
#'   \code{backend = "joint_mixed"} fit.
#'   Also accepts \code{backend = "conditional_logit"} fits.
#' @param newdata A data frame containing fixed covariates and formula offset
#'   variables. Row names must match point IDs in \code{kernel_inputs}. Kernel
#'   columns are computed from the prepared data, replacing any supplied values.
#' @param kernel_inputs Binned inputs from \code{\link{kernel_prep}} for new
#'   locations, with the same kernel and covariate names as the fit. If both
#'   \code{newdata} and \code{kernel_inputs} are omitted, predict the retained
#'   fitting rows. If only \code{newdata} is omitted, an empty frame is constructed
#'   from prepared point IDs; this works only without additional covariates.
#'   Conditional-logistic fits with optimized weighted surface SQ also require
#'   matching raw buffers and source/metric specifications in these inputs.
#' @param type Predict on the \code{"response"} or \code{"link"} scale.
#' @param se.fit Return standard errors of the estimated mean.
#' @param uncertainty \code{"joint"} includes estimated kernel scales, shapes,
#'   and their covariance with coefficients. \code{"conditional"} holds kernel
#'   parameters fixed. Both account for estimated dispersion in the coefficient
#'   information matrix, when applicable.
#' @param interval \code{"none"} or \code{"confidence"}. Confidence limits
#'   describe the expected response, not future observations.
#' @param level Confidence level strictly between zero and one.
#' @param ... Passed to native glmmTMB prediction for \code{joint_mixed};
#'   otherwise reserved for future use.
#' @details
#' For \code{conditional_logit}, \code{type = "response"} returns probabilities
#' normalized within each supplied choice set, not occurrence probabilities or
#' a utilization distribution. \code{type = "link"} returns the uncentered
#' linear predictor; differences within a set describe relative selection.
#' New predictions require both \code{newdata} with explicit choice-set IDs and
#' matching \code{kernel_inputs}; each set must have at least two alternatives.
#' The caller defines complete prediction choice sets. \code{se.fit = TRUE}
#' requires \code{uncertainty = "conditional"}, holding scales fixed; intervals
#' and joint uncertainty are unavailable. Formula offsets and training
#' standardization are preserved. Native \code{opt_mod} prediction methods use
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
#' @seealso \code{\link{multiScale_optim}}, \code{\link{vcov.multiScaleR}},
#'   \code{\link{prepare_prediction}} for repeated summaries of fixed rows.
#' @details For \code{joint_mixed}, prediction delegates to the final native
#'   glmmTMB model. Standard errors require \code{uncertainty = "conditional"}
#'   and exclude scale estimation uncertainty; confidence intervals are not
#'   provided by this wrapper. Supply both \code{newdata} and matching binned
#'   \code{kernel_inputs} for new rows. Include all fixed-effect, offset, and
#'   grouping variables explicitly. Spatial columns are recomputed using the
#'   training means and standard deviations. Additional arguments, including
#'   \code{re.form} and \code{allow.new.levels}, follow native glmmTMB behavior.
#'   The default includes fitted random effects; \code{re.form = NA} sets them
#'   to zero and does not integrate over their distribution. Use native
#'   \code{predict(object$opt_mod, ...)} for other glmmTMB prediction types.
#' @exportS3Method stats::predict
predict.multiScaleR <- function(object, newdata = NULL, kernel_inputs = NULL,
                                type = c("response", "link"), se.fit = FALSE,
                                uncertainty = c("joint", "conditional"),
                                interval = c("none", "confidence"), level = .95, ...) {
  if (identical(object$backend, "conditional_logit")) {
    if (length(list(...))) stop("Additional prediction arguments are unsupported for conditional_logit.", call. = FALSE)
    return(.msr_predict_clogit(object, newdata, kernel_inputs, match.arg(type),
      se.fit, match.arg(uncertainty), match.arg(interval)))
  }
  if (identical(object$backend, "joint_mixed"))
    return(.msr_predict_mixed(object, newdata, kernel_inputs, match.arg(type),
      se.fit, match.arg(uncertainty), match.arg(interval), ...))
  if (!.msr_analytic_backend(object))
    stop("This method requires backend = 'profile_glm' or 'joint_nb2'; use predict(object$opt_mod) for other backends.",
         call. = FALSE)
  type <- match.arg(type)
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
    J[, seq_len(p)] <- sweep(design$ds, 2L, utils::tail(beta, p), "*")
    if (object$kernel_inputs$kernel == "expow")
      J[, p + seq_len(p)] <- sweep(design$dh, 2L, utils::tail(beta, p), "*")
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

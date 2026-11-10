# Numeric prediction never executes or restores the native TMB tape.
.msr_gdr_field_prediction_function <- function(model) {
  field <- .msr_gdr_field_guard(model)
  forms <- field$forms
  force(field); force(forms)
  predictor <- function(object, type = "lambda", newdata = NULL, backTransform = TRUE,
      na.rm = TRUE, appendData = FALSE, level = .95, re.form = NULL,
      coords = NULL, crs = NULL, ...) {
    if (length(list(...))) stop("Unrecognized GDR field prediction controls.", call. = FALSE)
    if (!methods::is(object, "unmarkedFitGDRField") ||
        !identical(object@response_field$prediction_sha256, field$prediction_sha256) ||
        !identical(.msr_gdr_hash(object@data), field$response_sha256))
      stop("This field prediction context belongs to a different fitted target.", call. = FALSE)
    if (identical(type, "state")) type <- "lambda"
    if (!is.character(type) || length(type) != 1L || !type %in% names(field$bases))
      stop("Predict an active GDR field component: lambda, dist or rem.", call. = FALSE)
    if (methods::is(newdata, "unmarkedFrameGDR")) {
      if (newdata@numPrimary != 1L ||
          !identical(newdata@dist.breaks, object@data@dist.breaks) ||
          !identical(newdata@period.lengths, object@data@period.lengths))
        stop("Prediction frames must retain the declared GDR survey protocol.", call. = FALSE)
      temporary <- object; temporary@data <- newdata
      newdata <- get("get_covariates", asNamespace("unmarked"))(temporary, type)
    }
    if (is.null(newdata)) newdata <- field$raw[[type]]
    .msr_gdr_field_predict(field, forms, newdata, type, backTransform, na.rm,
      appendData, level, re.form, coords, crs)
  }
  # Retain only the numeric packet and fixed formulas, never the input fit
  # or its native TMB constructor/tape through an unused caller binding.
  scope <- list2env(list(field = field, forms = forms),
    parent = environment(.msr_gdr_field_prediction_function))
  environment(predictor) <- scope
  predictor
}

.msr_gdr_field_predict <- function(field, forms, data, type, backTransform,
    na.rm, appendData, level, re.form, coords, crs) {
  valid_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!all(vapply(list(backTransform, na.rm, appendData), valid_flag, logical(1))) ||
      (!is.null(re.form) && !(is.logical(re.form) && length(re.form) == 1L && is.na(re.form))))
    stop("Use scalar logical prediction flags and re.form = NULL or NA.", call. = FALSE)
  if (!is.null(level) && (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1))
    stop("level must be NULL or a confidence level strictly between zero and one.", call. = FALSE)
  if (!is.data.frame(data)) stop("Supply explicit GDR field prediction rows as a data frame or native GDR frame.", call. = FALSE)
  design <- .msr_unmarked_fixed_design(field$bases[[type]], data)
  X <- design$X; rows <- nrow(X); point <- field$point
  bytes <- 8 * (4 * as.double(rows) * length(point) + 4 * as.double(length(point))^2)
  if (!is.finite(bytes) || bytes > field$max_memory_mb * 2^20)
    stop("GDR field prediction work exceeds max_memory_mb; use smaller batches.", call. = FALSE)
  if (any(!is.finite(X)))
    stop("GDR field prediction rows must have complete finite fixed predictors.", call. = FALSE)
  beta <- which(names(point) == paste0("beta_", type))
  latent <- which(names(point) == paste0("b_", type))
  if (length(beta) != ncol(X)) stop("Prediction component coordinates changed.", call. = FALSE)
  gradient <- matrix(0, rows, length(point)); gradient[, beta] <- X
  new_variance <- numeric(rows)
  if (is.null(re.form) && length(latent)) {
    if (type == "lambda") {
      internal <- c(".msr_field_x", ".msr_field_y")
      if (is.null(coords) && all(internal %in% names(data))) coords <- as.matrix(data[internal])
      if (is.null(coords)) {
        Z <- .msr_unmarked_random_design(forms$lambda, field$raw$lambda, data)
        if (ncol(Z) != nrow(field$lower))
          stop("New abundance states require explicit projected coords and crs.", call. = FALSE)
        index <- .msr_gdr_field_map(Z, rows)
        locations <- field$coords[match(seq_len(nrow(field$lower)), field$block), , drop = FALSE]
        coords <- locations[index, , drop = FALSE]; crs <- field$crs
      } else if (is.null(crs))
        stop("Explicit new field coordinates require their declared crs.", call. = FALSE)
      projection <- field$project(coords, as.character(seq_len(rows)), crs)
      gradient[, latent] <- projection$loading
      sd <- exp(point[names(point) == "lsigma_lambda"])
      new_variance <- sd^2 * projection$conditional_unit_variance
    } else {
      Z <- .msr_unmarked_random_design(forms[[type]], field$raw[[type]], data)
      if (ncol(Z) != length(latent) || any(!is.finite(Z)))
        stop("Retain original known observer states or request re.form = NA.", call. = FALSE)
      gradient[, latent] <- as.matrix(Z)
    }
  }
  fitted_variance <- rowSums((gradient %*% field$joint_covariance) * gradient)
  allowance <- 64 * .Machine$double.eps * length(point) * pmax(1, abs(fitted_variance))
  if (any(fitted_variance < -allowance) || any(!is.finite(fitted_variance)))
    stop("Field prediction covariance is not numerically valid.", call. = FALSE)
  fitted_variance <- pmax(0, fitted_variance)
  eta <- drop(gradient %*% point); se <- sqrt(fitted_variance + new_variance)
  inverse <- if (type == "rem") stats::plogis else exp
  if (backTransform) {
    predicted <- inverse(eta)
    derivative <- if (type == "rem") predicted * (1 - predicted) else predicted
    reported_se <- derivative * se
  } else { predicted <- eta; reported_se <- se }
  quantile <- if (is.null(level)) NA_real_ else stats::qnorm((1 + level)/2)
  lower <- eta - quantile * se; upper <- eta + quantile * se
  if (backTransform) { lower <- inverse(lower); upper <- inverse(upper) }
  if (any(!is.finite(predicted)) || any(!is.finite(reported_se)))
    stop("GDR field predictions are not representable on the requested scale.", call. = FALSE)
  out <- data.frame(Predicted = predicted, SE = reported_se, lower = lower, upper = upper)
  if (appendData) out <- cbind(data, out)
  attr(out, "response_field_uncertainty") <- list(field_id = field$field_id,
    component = type, output_scale = if (backTransform) "native_component" else "link",
    point = if (is.null(re.form)) "inverse_link_conditional_predictor_not_field_average" else
      "zero_effect_reference_not_integrated",
    fitted_joint_predictor_variance = fitted_variance,
    new_field_residual_predictor_variance = new_variance,
    SE_includes_new_field_residual = any(new_variance > 0),
    omitted = c("observation_noise", "ecological_scale_uncertainty", "range_uncertainty"),
    range_fixed = TRUE, kernel_fixed = TRUE)
  out
}

#' Predict Components of a Fixed-Range GDR Response Field
#'
#' @param object An [unmarkedFitGDRField-class] returned by
#'   [fit_gdr_response_field()].
#' @param type Active native component, `"lambda"`, `"dist"` or `"rem"`.
#'   Default `"lambda"`; `"state"` is an alias for abundance/density.
#' @param newdata Explicit prediction data frame or native GDR frame. Default
#'   `NULL` uses retained training rows for the selected component. Fixed
#'   transformations and factors retain their training definitions.
#' @param backTransform Logical; default `TRUE` returns the native density,
#'   distance-key or removal-probability scale. `FALSE` returns the link scale.
#' @param na.rm Logical native interface flag. Default `TRUE`. This field route
#'   requires complete prediction rows and rejects missing fixed predictors.
#' @param appendData Logical; default `FALSE`. If `TRUE`, append prediction
#'   columns to the supplied component rows.
#' @param level Confidence level between zero and one, default 0.95. `NULL`
#'   omits confidence limits while retaining point predictions and SEs.
#' @param re.form `NULL` (default) retains conditional response effects. `NA`
#'   sets every response effect to zero; this does not integrate the field.
#' @param coords Optional two-column projected coordinates for abundance
#'   predictions at new locations. Default `NULL` uses original known abundance
#'   states named in the new rows. New states require explicit coordinates.
#' @param crs CRS of explicit new coordinates. Required with `coords`; it must
#'   equal the fitted projected CRS. Coordinates are never transformed silently.
#' @param ... Reserved. Unrecognized controls stop.
#' @return A data frame with `Predicted`, `SE`, `lower` and `upper` on the
#'   requested component scale. Confidence limits are constructed on the link
#'   scale before back-transformation. The `response_field_uncertainty`
#'   attribute identifies fitted full joint predictor variance separately from
#'   new-location residual field variance. Both enter the conditional predictor
#'   SE; neither is observation noise or ecological-scale uncertainty. Density
#'   point predictions are inverse-link conditional predictors, not lognormal
#'   field averages. Fixed range has no uncertainty in this model.
#' @seealso [fit_gdr_response_field()], [make_gdr_refitter()]
#' @importMethodsFrom unmarked predict
#' @md
#' @export
methods::setMethod("predict", "unmarkedFitGDRField", function(object, type = "lambda",
    newdata = NULL, backTransform = TRUE, na.rm = TRUE, appendData = FALSE,
    level = .95, re.form = NULL, coords = NULL, crs = NULL, ...) {
  .msr_gdr_field_prediction_function(object)(object, type, newdata, backTransform,
    na.rm, appendData, level, re.form, coords, crs, ...)
})

#' @rdname unmarkedFitGDRField-class
#' @importMethodsFrom unmarked update simulate profile ranef
#' @export
methods::setMethod("update", "unmarkedFitGDRField", function(object, ...) {
  stop("Use make_gdr_refitter() to preserve the audited GDR response field.", call. = FALSE)
})

#' @rdname unmarkedFitGDRField-class
#' @export
methods::setMethod("simulate", "unmarkedFitGDRField", function(object, nsim = 1, seed = NULL, ...) {
  stop("Spatial GDR simulation requires separate field qualification; iid simulation is not a substitute.", call. = FALSE)
})

#' @rdname unmarkedFitGDRField-class
#' @export
methods::setMethod("profile", "unmarkedFitGDRField", function(fitted, ...) {
  stop("Spatial GDR profile intervals require separate field qualification; use the retained joint diagnostics.", call. = FALSE)
})

#' @rdname unmarkedFitGDRField-class
#' @export
methods::setMethod("ranef", "unmarkedFitGDRField", function(object, ...) {
  stop("Joint latent-abundance posterior summaries need separate spatial GDR qualification.", call. = FALSE)
})

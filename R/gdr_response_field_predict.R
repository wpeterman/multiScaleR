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

#' Predict from a Distance-Removal Model with a Fixed-Range Abundance Field
#'
#' Predict density, the distance scale or the removal probability from a
#' model fitted by [fit_gdr_response_field()]. Abundance predictions include
#' the fitted field, at fitted blocks or at new projected locations, and their
#' standard errors include the field's uncertainty.
#'
#' @param object An [unmarkedFitGDRField][unmarkedFitGDRField-class] returned
#'   by [fit_gdr_response_field()].
#' @param type Character. The model component to predict: `"lambda"`
#'   (density, the default), `"dist"` (distance scale) or `"rem"` (removal
#'   probability). `"state"` is accepted as another name for `"lambda"`.
#' @param newdata A data frame of covariate values, or an `unmarkedFrameGDR`
#'   with the fitted survey protocol (one visit, the same distance bands and
#'   the same interval lengths). The default `NULL` predicts for the rows used
#'   in fitting. Transformations, factor levels and spline bases keep their
#'   definitions from fitting. With the default, removal predictions have one
#'   row per survey and removal interval.
#' @param backTransform Logical, default `TRUE`. `TRUE` returns density per
#'   square kilometer, the distance scale in meters or the removal probability
#'   per unit of interval time. `FALSE` returns the link scale (log, log or
#'   logit).
#' @param na.rm Logical, default `TRUE`. Kept for compatibility with
#'   `unmarked`. Rows with missing predictor values are not dropped. They
#'   stop the prediction, so remove them first.
#' @param appendData Logical, default `FALSE`. If `TRUE`, the prediction
#'   columns are added to the right of the prediction rows.
#' @param level Confidence level between 0 and 1, default 0.95. `NULL` skips
#'   the limits, which are then returned as `NA`.
#' @param re.form `NULL` (default) includes the fitted field, and any observer
#'   intercepts in the distance component. `NA` sets all of them to zero. A
#'   zero-field prediction describes a block with no deviation. It is not an
#'   average over blocks.
#' @param coords Optional two-column matrix or data frame of projected
#'   coordinates for
#'   abundance predictions at new locations, one row per row of `newdata`.
#'   The default `NULL` locates each row by its level of the grouping factor
#'   (for example `block = "b3"`), which must be a fitted level.
#' @param crs The CRS of `coords`, as an EPSG code or other input accepted by
#'   `sf::st_crs()`. Required with `coords`. It must equal the CRS used in
#'   fitting, because coordinates are never transformed.
#' @param ... Not used. Unrecognized arguments stop the prediction.
#' @return A data frame with one row per prediction row:
#'   * `Predicted`: the prediction on the scale set by `backTransform`.
#'   * `SE`: its standard error on the same scale. Back-transformed SEs use the
#'     delta method.
#'   * `lower`, `upper`: confidence limits at `level`, computed on the link
#'     scale and then back-transformed, so they are not symmetric around
#'     `Predicted`.
#'
#'   The attribute `response_field_uncertainty` explains the SE. Its element
#'   `fitted_joint_predictor_variance` is the link-scale variance from the
#'   fitted coefficients and field. `new_field_residual_predictor_variance` is
#'   the extra link-scale variance at new locations (zero at fitted blocks).
#'   `SE_includes_new_field_residual` is `TRUE` when any row has that extra
#'   variance. `omitted` lists what the SE leaves out: observation noise and
#'   uncertainty in the ecological scales and the range.
#' @details
#' # Choose where to predict
#'
#' For abundance, each prediction needs a field value. A row whose grouping
#' factor names a fitted block uses that block's fitted deviation. A row with
#' a location in `coords` uses the field value predicted at that location
#' from the fitted block deviations. At a fitted block location this equals
#' the block's own deviation. Far from all blocks, compared with the range, it
#' shrinks toward zero. For raster maps, use [response_field_geometry()] with
#' [predict.multiScaleR()] instead of calling this method cell by cell.
#'
#' For `type = "dist"` with an observer random intercept, rows must name
#' fitted observers unless `re.form = NA`.
#'
#' # What the standard errors include
#'
#' The SE combines two link-scale variances. The first comes from the joint
#' covariance of the fitted coefficients and field values, including their
#' cross-covariances. The second applies only at new locations. It is the
#' field variance left after conditioning on the fitted blocks, evaluated at
#' the estimated field SD. It is zero at a fitted block and grows toward the
#' full field variance far from all blocks. Neither
#' variance is observation noise, so these are confidence intervals for the
#' expected value, not prediction intervals for a new count.
#'
#' A density prediction is `exp()` of the linear predictor at the fitted field
#' values. It is not the average density over possible field values. The
#' intervals treat the ecological scales of effect and the correlation range
#' as known.
#' @seealso [fit_gdr_response_field()], [make_gdr_refitter()]
#' @examples
#' \donttest{
#' # The field needs unmarked 1.5.2 and TMB 1.9.25.
#' if (requireNamespace("TMB", quietly = TRUE) &&
#'     utils::packageVersion("unmarked") == "1.5.2" &&
#'     utils::packageVersion("TMB") == "1.9.25") {
#'   # Eight survey blocks 1 km apart, six point counts in each
#'   block <- factor(rep(paste0("b", 1:8), each = 6))
#'   hab <- rep(seq(-1, 1, length.out = 6), 8)
#'   shift <- c(-0.4, -0.3, -0.1, 0, 0.1, 0.2, 0.35, 0.4)[block]
#'   sim <- sim_dat_gdr(density = exp(log(400) + 0.5 * hab + shift),
#'     distance_scale = 40, removal = 0.2, dist_breaks = c(0, 25, 50, 75),
#'     interval_lengths = c(2, 3, 5), interval_units = "minutes",
#'     site_covs = data.frame(hab = hab, block = block), seed = 1)
#'   template <- unmarked::gdistremoval(lambdaformula = ~ hab + (1 | block),
#'     removalformula = ~ 1, distanceformula = ~ 1, data = sim$observed,
#'     keyfun = "halfnorm", mixture = "P", output = "density",
#'     unitsOut = "kmsq", K = 60, engine = "TMB",
#'     starts = c(log(400), 0, log(40), qlogis(0.2), log(0.3)))
#'   block_xy <- cbind(x = 500000 + 1000 * (as.integer(block) - 1),
#'                     y = 4400000)
#'   rownames(block_xy) <- rownames(unmarked::siteCovs(sim$observed))
#'   cold <- c("beta_lambda[1]" = log(300), "beta_lambda[2]" = 0,
#'             "beta_dist[1]" = log(30), "beta_rem[1]" = qlogis(0.4),
#'             "lsigma_lambda[1]" = log(0.5))
#'   settings <- list(method = "BFGS", se = TRUE, threads = 1L,
#'                    control = list(maxit = 2500L, reltol = 1e-10))
#'   field_fit <- fit_gdr_response_field(template, block_xy, crs = 32617,
#'     range = 2000, cold_starts = list(cold = cold),
#'     audit_dir = tempfile("gdr-field-"), tmb_settings = settings)
#'
#'   # Density (birds per square km) at average habitat: in block b3, at a
#'   # point halfway between blocks b3 and b4, and with the field set to zero
#'   in_block <- predict(field_fit, type = "lambda",
#'     newdata = data.frame(hab = 0, block = factor("b3", levels(block))))
#'   between <- predict(field_fit, type = "lambda",
#'     newdata = data.frame(hab = 0), coords = cbind(502500, 4400000),
#'     crs = 32617)
#'   no_field <- predict(field_fit, type = "lambda",
#'     newdata = data.frame(hab = 0), re.form = NA)
#'   dens <- rbind(in_block, between, no_field)
#'   rownames(dens) <- c("block b3", "between b3 and b4", "no field")
#'   print(dens)
#'   # The SE between blocks adds field variance left after conditioning on
#'   # the fitted blocks (both variances on the log scale)
#'   u <- attr(between, "response_field_uncertainty")
#'   print(c(fitted = unname(u$fitted_joint_predictor_variance),
#'           between_blocks = unname(u$new_field_residual_predictor_variance)))
#'   # Removal probability per minute, one row per survey and interval
#'   print(head(predict(field_fit, type = "rem"), 3))
#' }
#' }
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

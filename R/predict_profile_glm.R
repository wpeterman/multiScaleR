# New locations use the fitting sample's center and spread. Their derivatives
# must also include how those training quantities change with kernel scale.
.msr_profile_prediction_terms <- function(object) {
  tt <- object$profile_fit$contract$terms
  # Legacy full fits retain the evaluated training columns needed for recovery.
  # This same check runs before both point designs and analytic raster chunks.
  if (.msr_needs_nested_metadata(tt)) {
    native <- .analysis_model(object$opt_mod)
    tt <- .msr_fixed_prediction_terms(tt, .msr_response_frame_subset(
      stats::terms(native), native$model, tt), nested = TRUE)
  }
  tt
}

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
  tt <- stats::delete.response(.msr_profile_prediction_terms(object))
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

#' Predict from a fitted multiscale model
#'
#' @description
#' Predict the expected response, a linear predictor, or a native model
#' component from a [multiScale_optim()] fit, at the fitted rows, at new
#' points, or across a raster. Kernel-weighted covariates are recomputed at the
#' fitted scales of effect and standardized with the mean and standard
#' deviation of the fitted data. Fits from the analytic backends
#' (`"profile_glm"` and `"joint_nb2"`) can also carry joint uncertainty in the
#' kernel scales and the regression coefficients.
#'
#' [multiScaleR-prediction] explains which quantity each route returns and
#' what its uncertainty includes. This page gives the arguments and the rules
#' for each backend. `vignette("v09-estimates-and-intervals", package =
#' "multiScaleR")` compares conditional and joint intervals, and
#' `vignette("v03-map-predictions", package = "multiScaleR")` maps predictions.
#'
#' @param object A fitted `multiScaleR` object from [multiScale_optim()], from
#'   any backend: `"profile_glm"`, `"joint_nb2"`, `"conditional_logit"`, the
#'   experimental `"joint_mixed"`, or the default `"refit"`.
#' @param newdata Where to predict. One of:
#'   * `NULL` (the default), with `kernel_inputs = NULL`, to predict the rows
#'     used to fit the model.
#'   * A data frame of new rows, supplied together with `kernel_inputs`. It
#'     must hold every non-spatial variable in the model formula: fixed
#'     covariates, offset inputs, and grouping variables. Row names must match
#'     the point IDs in `kernel_inputs`, unless you supply `join_by` keys (refit
#'     backend). Kernel-weighted columns are computed for you and replace any
#'     columns of the same name.
#'   * A `SpatRaster` from [kernel_scale.raster()], to predict every cell.
#'     Leave `kernel_inputs` as `NULL`. See "Predict across a raster".
#' @param kernel_inputs A `multiScaleR_data` object from [kernel_prep()] for
#'   the new points. Prepare it with the fitted kernel, the same covariates and
#'   source or metric definitions, and the same binned or exact evaluation as
#'   the fit. Kernel means in the analytic, conditional logistic, and joint
#'   mixed backends need binned inputs (`bin = TRUE`). Its own distance unit
#'   conversion is applied, and the
#'   reported scales stay in map units. Default `NULL`. For the analytic
#'   backends, a model with no non-spatial covariates can omit `newdata`, and a
#'   frame is built from the prepared point IDs. Conditional logistic fits with
#'   an optimized weighted surface SQ metric also need matching raw cell
#'   buffers and source or metric specifications here.
#' @param type Scale or component to predict. For the analytic and joint mixed
#'   backends, `"response"` (the default, the model's expected response) or
#'   `"link"` (the linear predictor). For `"conditional_logit"`, `"response"`
#'   (the default, a probability within each choice set), `"link"`, or
#'   `"risk"`. For native refits, any type the provider accepts. When `type`
#'   is not supplied, a refit predicts the response for most models, `"state"`
#'   for `unmarked` occupancy and count models, and the linear predictor
#'   (`"lp"`) for Cox models, which have no response probability here. Raster
#'   prediction from an `unmarked` model needs an explicit `type`.
#' @param se.fit Logical. If `TRUE`, also return standard errors of the
#'   predicted mean, on the scale of `type`. Default `FALSE`.
#' @param uncertainty Which parameters the standard errors and intervals let
#'   vary. `"joint"` (the default) includes the estimated kernel scales and
#'   shapes and their covariance with the coefficients. It is available only
#'   for the analytic backends at points. `"conditional"` holds the kernel
#'   parameters at their estimates. Both include an estimated dispersion in
#'   the coefficient information when the family has one. Every other route
#'   needs `uncertainty = "conditional"` whenever it returns a standard error
#'   and stops otherwise.
#' @param interval `"none"` (the default) or `"confidence"`. Confidence limits
#'   describe the expected response, not a future observation. They are
#'   available only for the analytic backends at points.
#' @param level Confidence level for `interval = "confidence"`, strictly
#'   between zero and one. Default `0.95`.
#' @param ... Further controls. Which ones apply depends on the backend:
#'   * `include_fixed`: logical, default `FALSE`. Set `TRUE` to recompute
#'     declared fixed-radius predictors at new points. See "Recompute
#'     fixed-radius predictors". This is the only control that the analytic
#'     and conditional logistic backends use.
#'   * `join_by`: refit backend only. A data frame of site keys, one row per
#'     prepared point in `kernel_inputs` order, that lets several `newdata`
#'     rows share one site, for example one row per site and year. See
#'     [prepare_unmarked_data()] for the key rules.
#'   * Native prediction controls for `"joint_mixed"` and refit models, such
#'     as `re.form` or `allow.new.levels`. Put a control whose name matches a
#'     wrapper argument in a named `native_args` list, for example
#'     `native_args = list(level = 0)` for nlme population predictions.
#'     Native argument abbreviations are resolved with the provider's
#'     exact-first matching before any other check. They cannot replace the
#'     wrapper's verified inputs or bypass its checks.
#'   * `raster_args`: a named list of raster controls. See "Predict across a
#'     raster".
#'
#' @details
#' # Choose the rows to predict
#'
#' * **Fitted rows.** Omit `newdata` and `kernel_inputs`. The analytic
#'   backends rebuild the fitted design. Native refits return the provider's
#'   prediction for its fitted rows.
#' * **New points.** Run [kernel_prep()] at the new points with the fitted
#'   kernel, `max_D`, and covariate definitions, then supply the matching
#'   `newdata`. For the refit backend, supply both or neither.
#' * **A map.** Smooth and standardize the raster with
#'   [kernel_scale.raster()], then pass the result as `newdata`.
#'
#' New kernel-weighted values are standardized with the fitted data's mean and
#' sample standard deviation. They are never standardized within the
#' prediction sample. Every non-spatial variable named in the formula must be
#' supplied in `newdata`. Values are not borrowed from the formula
#' environment. For the analytic backends, new factor levels, missing values,
#' and unmatched point IDs are errors. Use literal constants in transform
#' arguments, such as `poly(z, 2)`. Preparing distance bins again at new
#' points can shift their representative distances slightly, so a prediction
#' at a fitted location can differ slightly from the fitted value.
#'
#' # Match the fitted formula
#'
#' Fixed transformations, factor contrasts, and offsets follow the fitted
#' terms. For ordinary `lm`, `glm`, `MASS::glm.nb`, and time-fixed,
#' unpenalized `coxph` and `clogit` refits, prediction uses the original
#' fitted terms. A dot shortcut such as `~ .` keeps its fitted meaning, and
#' new data columns never enter the design. Variables removed with a minus
#' term need no prediction column, while active interactions and offsets keep
#' their inputs. The original coefficients, covariance, factor coding, and
#' learned transformations are retained. Other native subclasses use their
#' provider's own prediction path. `unmarked::pcount()` (version 1.5.2)
#' rejects dot formulas when fitting, so write `unmarked` component formulas
#' explicitly.
#'
#' Standard `scale()`, `poly()`, `splines::ns()`, and `splines::bs()` calls,
#' including namespace-qualified calls, keep their training constants (center
#' and spread, polynomial coefficients, spline knots) for supported LM, GLM,
#' GAM, lme4, glmmTMB, time-fixed unpenalized Cox, and ordinary `unmarked`
#' occupancy and count predictions. New rows never relearn them, and a single
#' new value or a constant raster does not define a new standard deviation.
#' Nested expressions such as `I(scale(z)^2)` and transformed formula offsets
#' recover one standard learner per evaluated predictor from its retained
#' attributes, using an aligned evaluated training frame. Missing or ambiguous
#' attributes, several learners in one predictor, and unsupported
#' transformations need precomputed training predictors or complete saved
#' prediction calls. Scaling arguments keep their original order and standard
#' abbreviations. Custom functions keep their own semantics. Legacy analytic
#' fits can recover this metadata, but numerical snapshots from
#' [prepare_prediction()] must be prepared again.
#'
#' # Analytic GLM backends
#'
#' For `"profile_glm"` and `"joint_nb2"`, standard errors use the delta
#' method. With `uncertainty = "joint"`, they include the derivative of each
#' prediction with respect to the kernel scales and shapes, including how the
#' training mean and standard deviation change with scale. Confidence limits
#' are computed on the link scale and then transformed, and response-scale
#' standard errors are `abs(d mu / d eta)` times the link-scale standard
#' error. Binomial response predictions are probabilities.
#'
#' Joint uncertainty is refused for fits at a parameter bound or with weak
#' identifiability. `fit$uncertainty$reasons` says why, and point predictions
#' and conditional uncertainty remain available. Kernel parameters at a
#' search bound do not block joint uncertainty when they are the only failed
#' check. Joint standard errors and limits then hold those
#' parameters fixed at their bounds and let the other kernel parameters vary.
#' The result names them in `conditional_on`, and
#' `fit$uncertainty$joint_conditional_on` records them. These limits are
#' narrower than limits that also let the bounded scale vary, and the true
#' scale may lie beyond the bound. Widen the search, for example with a larger
#' `max_D` in [kernel_prep()], before you rely on them. At the exact Poisson
#' boundary of a `"joint_nb2"` fit, point predictions and explicitly
#' conditional intervals hold both the scales and the Poisson dispersion
#' boundary fixed. Use [vcov.multiScaleR()] for the covariance matrices.
#'
#' The native model in `fit$opt_mod` is a conditional GLM. Its own prediction
#' methods expect kernel predictors that are already standardized and do not
#' include scale uncertainty.
#'
#' # Conditional logistic backend
#'
#' For `"conditional_logit"`, `type = "response"` returns probabilities
#' normalized within each supplied choice set. They are not occurrence
#' probabilities or a utilization distribution. `type = "link"` returns the
#' linear predictor, including formula offsets, and differences within a set
#' describe relative selection. `type = "risk"` returns `exp(link)`, a
#' relative selection score. Link and risk predictions at new points need
#' matching spatial inputs and every fixed formula input, but no choice-set
#' IDs. Response predictions also need explicit choice-set IDs, at least two
#' alternatives per set, and complete choice sets that you define.
#'
#' `se.fit = TRUE` requires `uncertainty = "conditional"`. Risk standard
#' errors are `risk x SE(link)`, and response standard errors use the delta
#' method through the within-set normalization. Intervals and joint
#' uncertainty are unavailable. Formula offsets and training standardization
#' are preserved. The native `fit$opt_mod` prediction methods use their own
#' survival conventions and require already standardized covariates.
#'
#' # Joint mixed backend
#'
#' For `"joint_mixed"`, prediction delegates to the final native glmmTMB
#' model at the selected scales, with `type = "response"` or `"link"`.
#' Standard errors require `uncertainty = "conditional"` and exclude scale
#' uncertainty. This wrapper provides no confidence intervals. For new rows,
#' supply both `newdata` and matching `kernel_inputs`, including every
#' fixed-effect, offset, and grouping variable. Kernel means need binned
#' summaries, and weighted surface SQ needs matching source specifications,
#' raw cell buffers, and row IDs. Spatial columns are recomputed with the
#' training mean and standard deviation.
#'
#' Additional arguments, such as `re.form` and `allow.new.levels`, follow
#' native glmmTMB behavior. The default includes the fitted random effects,
#' and `re.form = NA` sets them to zero without integrating over their
#' distribution. Use native `predict(fit$opt_mod, ...)` for other glmmTMB
#' prediction types, or bind one with [prepare_prediction()].
#'
#' # Native refits
#'
#' The default `"refit"` backend returns the provider's own prediction for
#' the model in `fit$opt_mod`, with spatial means and metrics recomputed from
#' matching source definitions, prepared rows, physical distance units, and
#' the selected scales. Unused prepared variables are ignored. Supply both
#' `newdata` and `kernel_inputs`, or omit both for the fitted rows. Explicit
#' `join_by` keys let several prediction rows share one site.
#'
#' Standard errors require `uncertainty = "conditional"` and a provider that
#' supplies them. An unsupported native standard-error request raises an
#' error rather than returning zero. The wrapper's own `interval` must stay
#' `"none"`. Native covariance, variance, and interval controls also require
#' `uncertainty = "conditional"`, even with `se.fit = FALSE`, and cannot add
#' uncertainty in the kernel parameters. Cross-row native covariance, such as
#' glmmTMB's `cov.fit = TRUE`, needs point predictions, which keep the native
#' covariance output.
#' * **Linear models.** For an ordinary `lm` refit,
#'   `native_args = list(interval = "confidence", level = 0.9)` returns
#'   native confidence limits for the mean, and `interval = "prediction"`
#'   within that list returns limits for a future observation, including the
#'   residual variance. Both keep the native `fit`, `lwr`, and `upr` columns
#'   for points and maps.
#' * **Random effects.** Native random-effect controls keep their meaning.
#'   For glmmTMB and lme4, `re.form = NA` sets the random effects to zero.
#'   For nlme, use `native_args = list(level = 0)`. Zero-effect predictions can
#'   omit random-only grouping and slope inputs, but every fixed predictor and
#'   offset remains required. For partial lme4 `re.form` formulas, supply the
#'   fixed predictors and the grouping and slope variables of the selected
#'   effects. Native validation still checks that the selected effects belong
#'   to the fitted model.
#' * **spaMM.** spaMM refits use their retained native model frames to verify
#'   fixed features. Zero and selected `re.form` predictions need only the
#'   fixed inputs and the selected random-effect inputs. The generic `se.fit`
#'   interface is unavailable when the native method does not return a
#'   standard-error list. Native `variances` and `intervals` attributes keep
#'   their own scales and targets and require point predictions. `fixefVar`
#'   and `predVar` are on the linear-predictor scale. `residVar` and
#'   `respVar` are on the response scale and include observation noise.
#' * **Cox models.** Predictions are linear predictors (`"lp"`, also selected
#'   by `type = "link"`) or relative risks (`"risk"`). Neither is a
#'   probability. Use the `"conditional_logit"` backend for choice-set
#'   probabilities.
#'
#' # Predict from unmarked models
#'
#' Name the component with `type`, such as `"state"` or `"det"`.
#' [multiScaleR-prediction] lists what each component means. Supply
#' `newdata` in the row layout that native `unmarked` prediction expects for
#' that component, with row names that match the prepared point IDs or with
#' explicit `join_by` keys. See [prepare_unmarked_data()] for the key and
#' fixed-radius rules. Without
#' `se.fit`, the result is the vector of native `Predicted` values. With
#' `se.fit = TRUE` and `uncertainty = "conditional"`, it is the native table
#' with confidence limits at the native level.
#'
#' Ordinary occupancy and `pcount` state and detection predictions freeze
#' supported fixed scaling, polynomial, and spline terms from the original raw
#' covariates before site removal, including nested calls and formula
#' offsets. Keep those original covariates and the fitted designs. The
#' reconstructed designs and offsets must match the native fit before new rows
#' are evaluated. This check preserves the fitted target, but it cannot
#' detect edits to raw data that leave the evaluated design unchanged. Factor
#' contrasts use retained metadata, or original columns for factor main
#' effects when those recover the coding exactly. Other missing coding needs
#' explicit saved training contrasts.
#'
#' For explicit data-frame rows, only the requested component's training
#' metadata is prepared. Raster chunks and effect curves reuse that
#' component's checked basis. A native `unmarkedFrame` input, or an omitted
#' `newdata`, keeps both components, because native missing-row removal
#' depends on both designs. Derived ZIP means and limits still use the native
#' zero-inflation estimate and covariance, and the full model's parameter
#' count and diagnostics are unchanged. Native observation rows, missing-site
#' handling, and confidence levels stay native, and frame subsets keep their
#' original known-group columns.
#'
#' Conditional predictions with learned random-effect terms need precomputed
#' original training predictors. Zero-effect requests do not evaluate those
#' unused terms or need new random-only grouping or slope inputs. `unmarked`
#' treats every non-`NULL` `re.form` as zero random effects. It does not
#' select partial effects or integrate over them. Other `unmarked` classes
#' keep their provider's prediction requirements.
#'
#' Qualified `gdistremoval()` predictions keep the original `"lambda"`,
#' active `"phi"`, `"dist"`, and `"rem"` component recipes, with their native
#' site, primary-period, and removal-interval rows. Native C and TMB designs
#' are checked before frozen scaling, polynomial, spline, or factor
#' interaction terms are used. Points, prepared outputs, effect curves, and
#' raster chunks share these terms. Rows for one component need only that
#' component's inputs, and native GDR frames keep their native layout and
#' missing-covariate checks. Native component covariance and conditioning are
#' unchanged, and zero effects do not integrate over unknown groups.
#'
#' # Recompute fixed-radius predictors
#'
#' By default, declared fixed-radius predictors are inputs you supply in
#' `newdata`. For the analytic, conditional logistic, joint mixed, and refit
#' backends, `include_fixed = TRUE` recomputes them from matching spatial
#' inputs instead. Their original prepared values and training means and
#' sample standard deviations must first be verified against the retained
#' native model frame. Source, metric, radius, kernel, resolution, buffer
#' extent (`max_D`), and observation-to-site alignment must match. Fixed radii
#' add no scale parameters or uncertainty.
#'
#' # Predict across a raster
#'
#' First call [kernel_scale.raster()] with the fitted object,
#' `scale_center = TRUE`, and `match_buffer = TRUE`. `match_buffer = TRUE`
#' needs square cells in a projected or local planar coordinate reference
#' system, at the prepared resolution.
#' Modeled fixed-radius metrics also require `include_fixed = TRUE` there.
#' Set categorical and other non-spatial scenario values with `site_values`,
#' and offset inputs with `offset_values`. A simple `log()` offset defaults to
#' a raw value of 1 and a direct offset to 0, and a compound offset such as
#' `offset(log(effort * duration))` needs explicit values. Logical predictors
#' are coded 0 and 1 in the raster and regain their fitted type during
#' prediction. The descriptor attached to the raster checks scales, source and
#' metric definitions, normalization, and support against this fit. It does
#' not rebuild a lost descriptor or establish that the maps have the right
#' scientific provenance. Recompute the raster after editing cells or saving
#' it in a file format that drops R attributes.
#'
#' Raster standard errors are conditional on the selected kernel parameters,
#' with `interval = "none"`. Analytic GLM maps use the same coefficient
#' information as point predictions. Other maps use native prediction, and
#' native matrix components return one layer per column, with standard-error
#' layers in the same column order. Cox maps default to linear predictors.
#' Conditional logistic maps need `type = "link"` or `"risk"`, because choice
#' probabilities need complete point choice sets. `unmarked` maps need an
#' explicit component. Native cross-row covariance (`cov.fit = TRUE`) and
#' spaMM variance or interval attributes are rejected, because raster chunks
#' cannot keep them. Supply native controls in `native_args`. A raster
#' prediction from one set of selected scales is not a Bayesian posterior map.
#'
#' The `raster_args` list accepts `filename`, `overwrite`, `wopt`, `const`,
#' `cores`, `na.rm`, and `field_geometry`. `cores` defaults to one worker,
#' and clusters created internally are stopped on exit. `field_geometry`
#' takes a checked [response_field_geometry()] descriptor, as described in
#' "Spatial response fields". Rasterization, distance binning, and partial source
#' coverage remain separate approximation choices.
#'
#' # Spatial response fields
#'
#' Native glmmTMB exponential, Gaussian, and Matern response fields keep their
#' provider predictions and mean-prediction standard errors at the selected
#' kernel scales. Conformance checks cover Gaussian, log-Poisson,
#' logit-binomial, log-NB2, and log-Gamma responses with one shared projected
#' field without a year effect. Supply the native coordinate factor and
#' grouping columns explicitly, including new [glmmTMB::numFactor()]
#' positions and `allow.new.levels = TRUE` when appropriate. A field's
#' covariance range is distinct from the ecological kernel scale.
#' `re.form = NA` sets the response random effects to zero without
#' integrating over their distribution. These checks do not qualify scale
#' recovery, interval coverage, or a spatial field in the `"joint_mixed"`
#' likelihood. A full raster needs matching field coordinates in each cell,
#' which a constant site scenario does not provide. Use explicit native point
#' rows or `raster_args = list(field_geometry = ...)` for coordinate-aware
#' maps. That route requires a matching projected coordinate reference
#' system, uses cell-center coordinates and one worker, limits native batch
#' expansion to an explicit working-memory budget, and needs
#' `allow.new.levels = TRUE` in `native_args` for new positions.
#'
#' Identified native TMB `gdistremoval()` response-field models keep their
#' conditional component predictor and full joint covariance through the
#' `unmarked` adapter. New abundance coordinates require the fitted projected
#' coordinate reference system. Field maps need a checked geometry descriptor
#' and `type = "lambda"`, and they run on one worker, including zero-effect
#' references. The residual field variance at new locations is distinct from
#' observation noise and from ecological-scale uncertainty.
#'
#' # GAM and BAM predictions
#'
#' GAMs with several linear predictors keep each component's terms,
#' coefficient spans, and native output columns. Response columns keep the
#' fitted family's meaning. For example, the second `mgcv::gaulss()` response
#' column is a reciprocal standard deviation, and multinomial response columns
#' are category probabilities. For `mgcv::ziplss()` response predictions, the
#' native mean is unchanged, and its standard error propagates the full native
#' coefficient covariance through both predictors, keeping signed
#' cross-covariance. It omits count noise and ecological-scale uncertainty.
#' Custom family prediction hooks remain guarded.
#'
#' Native term selection, exclusions, and smoothing-parameter covariance
#' requests (`unconditional`) keep their native behavior, and none adds
#' uncertainty in the selected scales. Package predictions keep the saved
#' factor coding of interaction-only terms when session contrast options
#' change. BAM predictions keep the native discrete and block algorithms,
#' restore fitted coding in serial and native worker calls, and pass explicit
#' term exclusions and controls such as `iterms.type` and `unconditional` to
#' those workers. Worker mgcv versions must match the parent. Optimized
#' discrete predictions for binomial and Poisson BAMs keep their native
#' values, warnings, and output geometry across point, snapshot, and raster
#' routes, including the native warning that resets `iterms` to `terms`. Set
#' `discrete = FALSE` for ordinary `iterms`. Native warnings, such as one for
#' selecting terms from a model without parametric coefficients, reach the
#' calling session with their native meaning.
#'
#' @return The class depends on the route.
#' * **Analytic backends, no standard errors or interval:** a numeric vector
#'   of predictions on the scale of `type`, named by row ID.
#' * **Analytic backends with `se.fit = TRUE` or `interval = "confidence"`:**
#'   a list with `fit` (the predictions), `se.fit` (standard errors on the
#'   scale of `type`), `type`, and `uncertainty` (`"joint"` or
#'   `"conditional"`). Confidence intervals add `lower` and `upper`, on the
#'   scale of `type`, and `level`. Joint uncertainty from a fit with kernel
#'   parameters at a search bound adds `conditional_on`, which names the
#'   parameters held fixed.
#' * **`"conditional_logit"`:** a named numeric vector, or with
#'   `se.fit = TRUE` a list with `fit`, `se.fit`, `type`, and
#'   `uncertainty = "conditional"`.
#' * **`"joint_mixed"`:** the native glmmTMB prediction vector, or with
#'   `se.fit = TRUE` its list with `fit` and `se.fit`.
#' * **Native refits:** the provider's own output, with a
#'   [prediction_contract] attribute that records the quantity and its
#'   uncertainty. For an `unmarked` model this is the vector of native
#'   `Predicted` values, or with `se.fit = TRUE` the native data frame with
#'   `Predicted` (the component on its natural scale), `SE`, and the
#'   confidence limits `lower` and `upper`.
#' * **Raster input:** a `SpatRaster` of predictions, with a standard-error
#'   layer when requested and available, and `prediction_contract` and
#'   `projection_contract` attributes that record the controls used.
#'   `unmarked` output keeps the native `Predicted`, `SE`, `lower`, and
#'   `upper` columns as layers.
#'
#' No route adds observation noise, except the native `lm` prediction
#' interval and the spaMM response variances (see "Native refits").
#'
#' @seealso [multiScaleR-prediction] for the quantity and uncertainty of each
#'   route, [vcov.multiScaleR()] for the covariance matrices,
#'   [prepare_prediction()] for repeated summaries of fixed rows,
#'   [prepare_unmarked_data()] for keyed `unmarked` prediction rows,
#'   [kernel_scale.raster()] for prediction rasters, [plot_marginal_effects()]
#'   for effect curves, and [multiScale_optim()].
#'
#' @examples
#' ## Running example: counts at 100 survey points and two landscape layers
#' data("landscape_counts")
#' data("surv_pts")
#' pts <- terra::vect(surv_pts)
#' land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                 package = "multiScaleR"))[[c("land1", "land2")]]
#'
#' ## The simulated layers are in meters but carry no coordinate reference
#' ## system. Declare a local planar one so raster prediction can check cells.
#' terra::crs(land) <- "local"
#' terra::crs(pts) <- "local"
#'
#' ## bin = TRUE stores the distance bins that the analytic backend needs
#' inputs <- kernel_prep(pts, land, max_D = 1700, kernel = "gaussian",
#'                       bin = TRUE, verbose = FALSE)
#' dat <- data.frame(landscape_counts, inputs$kernel_dat)
#' start <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
#' fit <- multiScale_optim(start, inputs, backend = "profile_glm",
#'                         n_cores = 1, verbose = FALSE)
#'
#' ## Expected counts at the first three survey points
#' predict(fit)[1:3]
#'
#' ## 95% limits for the expected count. Joint limits let the scales vary,
#' ## conditional limits hold them at their estimates.
#' joint <- predict(fit, interval = "confidence")
#' cond <- predict(fit, interval = "confidence", uncertainty = "conditional")
#' round(cbind(expected = joint$fit,
#'             joint_lower = joint$lower, joint_upper = joint$upper,
#'             cond_lower = cond$lower, cond_upper = cond$upper)[1:3, ], 2)
#'
#' ## Back-transform a link-scale prediction yourself: exponentiate the
#' ## estimate and its limits, never the standard error. This reproduces the
#' ## joint limits for the first point.
#' eta <- predict(fit, type = "link", se.fit = TRUE)
#' exp(eta$fit[1] + c(-1, 1) * qnorm(0.975) * eta$se.fit[1])
#'
#' ## New points: prepare their landscape with the fitted settings, then give
#' ## site values whose row names ("1", "2", "3") match the prepared point IDs
#' new_inputs <- kernel_prep(pts[1:3], land, max_D = 1700, kernel = "gaussian",
#'                           bin = TRUE, verbose = FALSE)
#' new_sites <- data.frame(site = c(-1, 0, 1))
#' predict(fit, new_sites, new_inputs, se.fit = TRUE)
#'
#' ## Map the expected count and its conditional standard error. The site
#' ## layer holds the mean of site, so the map describes an average site.
#' layers <- kernel_scale.raster(land, multiScaleR = fit, scale_center = TRUE,
#'                               match_buffer = TRUE, verbose = FALSE)
#' count_map <- predict(fit, layers, se.fit = TRUE, uncertainty = "conditional")
#' terra::plot(count_map)
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
  # The analytic backends accept only include_fixed through `...`. Stop on
  # anything else (for example re.form or a misspelled name) rather than
  # silently ignoring it, as the other backends do.
  extra <- setdiff(names(dots), "include_fixed")
  if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) || length(extra)))
    stop("Unsupported prediction arguments for the analytic backends: ",
         if (length(extra)) paste(extra, collapse = ", ") else "unnamed arguments",
         ". Only include_fixed is accepted through `...`.", call. = FALSE)
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
  # Joint uncertainty from a profile_glm fit at a search bound holds the
  # bounded kernel parameters fixed; name them with the result.
  held <- attr(V, "conditional_on")
  if (length(held)) out$conditional_on <- held
  out
}

# Fixed radial bases use projected coordinate units. All centering/projection
# constants are learned on training rows once and reused at new locations.
.msr_spatial_coordinates <- function(data, coordinates) {
  if (!is.character(coordinates) || length(coordinates) != 2L || anyNA(coordinates) ||
      anyDuplicated(coordinates) || !all(coordinates %in% names(data)))
    stop("coordinates must name two projected numeric columns in data.", call. = FALSE)
  x <- as.matrix(data[coordinates])
  if (!is.numeric(x) || any(!is.finite(x)))
    stop("Spatial coordinates must be finite numeric projected coordinates.", call. = FALSE)
  x
}

.msr_spatial_radial <- function(coordinates, knots, range) {
  squared <- outer(coordinates[, 1L], knots[, 1L], "-")^2 / range^2 +
    outer(coordinates[, 2L], knots[, 2L], "-")^2 / range^2
  exp(-.5 * squared)
}

.msr_spatial_basis <- function(coordinates, knots, range, fixed) {
  validate_scalar_numeric(range, "range", positive = TRUE)
  if (!is.matrix(knots) || !is.numeric(knots) || ncol(knots) != 2L || !nrow(knots) ||
      any(!is.finite(knots)) || anyDuplicated(as.data.frame(knots)))
    stop("knots must be a finite matrix of distinct projected coordinate pairs.", call. = FALSE)
  if (nrow(fixed) != nrow(coordinates) || qr(fixed)$rank < ncol(fixed))
    stop("Spatial training rows need a full-rank fixed scale design.", call. = FALSE)
  radial <- .msr_spatial_radial(coordinates, knots, range)
  projection <- qr.solve(fixed, radial, tol = 1e-10)
  basis <- radial - fixed %*% projection
  if (qr(basis, tol = 1e-10)$rank != ncol(basis))
    stop("Spatial basis is rank deficient after projecting out the fixed scale design; use fewer or different knots.", call. = FALSE)
  normalizer <- sqrt(mean(rowSums(basis^2)))
  if (!is.finite(normalizer) || normalizer <= 1e-12)
    stop("Spatial basis has no remaining variation.", call. = FALSE)
  list(design = basis / normalizer, knots = knots, range = range,
    fixed_projection = projection, normalizer = normalizer)
}

.msr_spatial_basis_predict <- function(basis, coordinates, fixed) {
  (.msr_spatial_radial(coordinates, basis$knots, basis$range) -
    fixed %*% basis$fixed_projection) / basis$normalizer
}

#' Prepare an Experimental Spatially Varying Scale Model
#'
#' Add one smooth Gaussian basis field to a kernel's log scale. This private
#' prototype preserves the existing public formula interfaces and priors.
#'
#' @param formula,inputs,data,scales,family,row_id,priors Existing Bayesian or
#'   private likelihood preparation arguments. Existing scale formulas may
#'   include independent annual intercepts.
#' @param covariate One kernel predictor receiving the spatial scale field.
#' @param coordinates Names of two projected coordinate columns in data.
#' @param knots Numeric matrix with two coordinate columns and one row per basis.
#' @param range Positive fixed smoothing distance in the coordinate units.
#' @param field_sd_prior Positive half-normal SD for the field's RMS pointwise
#'   log-scale SD, required explicitly for the posterior target. Use NULL for ML.
#' @param target Posterior or unpenalized marginal-likelihood target.
#' @param reference_scale Positive distance used by private likelihood preparation.
#' @return A private \code{msr_spatial_scale_prototype} prepared object. Its
#'   specification contains one Gaussian basis term and appended whitened latent
#'   coefficients/field SD. Range and knots are fixed, not estimated. Bayesian
#'   field SD has the explicitly declared prior; ML adds no population prior.
#' @details Basis columns are projected off the fixed scale design and normalized
#'   so the mean pointwise prior variance is one before multiplying by field SD.
#'   Prediction retains these training constants. Spatial variation is variation
#'   in ecological kernel scale, not a spatially varying response coefficient.
#'   The scale-field covariance range is distinct from ecological kernel scale.
#'   This is a finite basis model with no additional new-location innovation.
#'   The first prototype requires an intercept-only fixed scale design, with
#'   optional annual random intercepts, so the field is a function of location
#'   alone. Other fixed scale formulas remain available in the existing models.
#'   It provides no claim of interval calibration or detection of nonstationarity.
#'   Use the dedicated private prediction adapter; public summaries/prediction
#'   and scalar/one-SD integrators do not support Gaussian-block specifications.
#' @keywords internal
.msr_spatial_scale_prepare <- function(formula, inputs, data, scales = NULL,
    family = c("nbinom2", "poisson", "binomial", "gaussian"), row_id = NULL,
    priors = list(), covariate, coordinates, knots, range, field_sd_prior = NULL,
    target = c("posterior", "marginal_likelihood"), reference_scale = NULL) {
  target <- match.arg(target); family <- match.arg(family)
  if (target == "posterior") {
    validate_scalar_numeric(field_sd_prior, "field_sd_prior", positive = TRUE)
    spec <- .msr_bayes_prepare(formula, inputs, data, scales, family, row_id, priors)
  } else {
    if (length(priors) || !is.null(field_sd_prior))
      stop("The likelihood target does not use population or field-SD priors.", call. = FALSE)
    spec <- .msr_marginal_prepare(formula, inputs, data, scales, family, row_id, reference_scale)
  }
  if (length(covariate) != 1L || !covariate %in% spec$covariates)
    stop("covariate must select one prepared kernel predictor.", call. = FALSE)
  fixed <- spec$scales[[covariate]]$X
  if (ncol(fixed) != 1L || any(fixed[, 1L] != 1))
    stop("The first spatial-field prototype requires an intercept-only fixed scale formula; annual random intercepts are supported.", call. = FALSE)
  xy <- .msr_spatial_coordinates(data, coordinates)
  basis <- .msr_spatial_basis(xy, knots, range, spec$scales[[covariate]]$X)
  block <- .msr_gaussian_block(basis$design, label = "spatial log-scale basis")
  block$basis <- basis[names(basis) != "design"]
  block$coordinates <- coordinates
  q <- block$dimension
  old_z <- length(spec$parameters$scale_z); old_sd <- length(spec$parameters$log_scale_sd)
  term <- list(group = ".spatial_field", levels = paste0("basis", seq_len(q)),
    z_index = old_z + seq_len(q), sd_index = old_sd + 1L, gaussian_block = block)
  if (".spatial_field" %in% names(spec$scales[[covariate]]$random))
    stop("The reserved spatial-field term already exists.", call. = FALSE)
  spec$scales[[covariate]]$random[[".spatial_field"]] <- term
  spec$parameters$scale_z <- c(spec$parameters$scale_z, rep(0, q))
  spec$parameters$log_scale_sd <- c(spec$parameters$log_scale_sd, log(.2))
  if (target == "posterior")
    spec$priors$scale_group_sd <- c(rep(spec$priors$scale_group_sd, length.out = old_sd), field_sd_prior)
  structure(list(specification = spec, parameters = spec$parameters,
    inference = "prepared spatial scale prototype", target = target,
    prototype = list(covariate = covariate, fixed_range = TRUE,
      field_sd_prior = field_sd_prior, calibration = "not_established",
      field_sd_meaning = "RMS pointwise log-scale SD on training rows")),
    class = "msr_spatial_scale_prototype")
}

.msr_spatial_scale_objective <- function(object, ad = FALSE) {
  if (!inherits(object, "msr_spatial_scale_prototype") || length(object$target) != 1L ||
      is.na(object$target) || !object$target %in% c("posterior", "marginal_likelihood"))
    stop("Unknown spatial prototype target; use its prepared posterior or likelihood specification.", call. = FALSE)
  if ((object$target == "posterior" && is.null(object$specification$priors)) ||
      (object$target == "marginal_likelihood" && length(object$specification$priors)))
    stop("Spatial prototype target and population-prior specification do not agree.", call. = FALSE)
  if (object$target == "posterior") .msr_bayes_objective(object$specification, ad) else
    .msr_marginal_objective(object$specification, ad)
}

# A normalized target adapter, not an ML fitting/interval/eligibility interface.
.msr_spatial_scale_adfun <- function(object, integrate_latent = FALSE) {
  if (!requireNamespace("RTMB", quietly = TRUE)) stop("Spatial prototype requires RTMB.", call. = FALSE)
  random <- intersect(c("scale_z", "response_z"), names(object$parameters)[lengths(object$parameters) > 0L])
  RTMB::MakeADFun(.msr_spatial_scale_objective(object, TRUE), object$parameters,
    random = if (integrate_latent && length(random)) random else NULL, silent = TRUE)
}

.msr_spatial_scale_laplace <- function(object, draws = 2000L, seed = NULL, control = list()) {
  if (object$target != "posterior")
    stop("This prototype fitter supports the posterior target; likelihood has an objective adapter only.", call. = FALSE)
  .msr_bayes_laplace(object, draws, seed, control)
}

# All observation rows are retained. New-group effects use the existing policy
# and layout, with new coefficients appended after the spatial latent block.
.msr_spatial_scale_predict <- function(object, newdata, parameters = NULL,
    scale = c("distance", "log"), new_levels = c("error", "sample", "zero"), draws = NULL) {
  scale <- match.arg(scale); new_levels <- match.arg(new_levels)
  if (is.null(parameters)) {
    ids <- .msr_bayes_draw_ids(object, draws)
    values <- lapply(ids, function(i) .msr_spatial_scale_predict(object, newdata,
      parameters = .msr_bayes_draw(object, i), scale = scale, new_levels = new_levels))
    return(do.call(rbind, values))
  }
  nm <- object$prototype$covariate; spec <- object$specification
  term <- spec$scales[[nm]]$random[[".spatial_field"]]
  spec$scales[[nm]]$random[[".spatial_field"]] <- NULL
  context <- .msr_bayes_prediction_scales(spec, newdata, new_levels, TRUE)
  context$n_response <- 0L
  p <- .msr_bayes_prediction_parameters(parameters, context)
  sm <- context$spec$scales[[nm]]
  xy <- .msr_spatial_coordinates(newdata, term$gaussian_block$coordinates)
  term$gaussian_block$design <- .msr_spatial_basis_predict(term$gaussian_block$basis, xy, sm$X)
  sm$random[[".spatial_field"]] <- term
  log_s <- drop(sm$X %*% p$scale_beta[sm$beta_index])
  for (effect in sm$random) {
    log_s <- log_s + if (is.null(effect$gaussian_block))
      exp(p$log_scale_sd[effect$sd_index]) * p$scale_z[effect$z_index][effect$index] else
      .msr_gaussian_effect(effect$gaussian_block, p$scale_z[effect$z_index], p$log_scale_sd[effect$sd_index])
  }
  log_s <- log_s + log(spec$unit_conv)
  if (scale == "log") log_s else exp(log_s)
}

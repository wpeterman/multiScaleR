#' Evaluate Private Spatial Conditional Response Log Densities
#'
#' Return one normalized response log density per original training observation,
#' conditional on each supplied complete parameter vector and its shared field.
#'
#' @param object A prepared private \code{msr_spatial_scale_prototype}. The
#'   existing posterior or unpenalized likelihood target must remain consistent
#'   with its prior specification. Ordinary scale-group effects are excluded.
#' @param parameters Optional named list containing one complete finite original
#'   parameter vector. NULL uses stored complete posterior vectors. Supplying a
#'   fixed vector does not propagate parameter uncertainty.
#' @param draws Optional positive number of evenly spaced stored joint vectors.
#'   NULL uses all. This controls the returned matrix size, not chain precision.
#'   It cannot be combined with \code{parameters}.
#' @return A numeric joint-vectors-by-training-observations matrix in natural
#'   log density units, including each response distribution's normalizing
#'   constants. Columns retain original training row order. Attributes
#'   \code{observation_ids} and \code{training_rows} retain the original buffer
#'   identities and buffer-row mapping, including repetitions; \code{draw_ids}
#'   identifies stored vectors, or is NULL for a supplied fixed vector.
#'   \code{likelihood_contract} records:
#'   \describe{
#'     \item{target, source_target}{Conditional response density and the object's
#'       unchanged posterior or likelihood preparation target, respectively.}
#'     \item{parameter_source, method, vector_integration}{Supplied fixed versus
#'       stored joint/resampled vectors, the original computation label and its
#'       original recorded integration metadata. These do not integrate this
#'       conditional density over effects.}
#'     \item{population_priors, latent_densities, marginalization}{All omitted
#'       from this response-only quantity.}
#'     \item{conditioning, spatial}{Original responses, frozen designs/reference
#'       transforms and shared realized effects; fixed field range/basis and no
#'       independent new-location innovation.}
#'     \item{parameter_blocks, population_dimension, latent_dimension}{Original
#'       parameter block sizes and their population/latent dimensions. The
#'       realized scale/response latent coordinates are excluded from the
#'       population dimension. Fixed range and knots add no fitted coordinates.}
#'     \item{uncertainty}{Fixed or joint-vector evaluation, with no intervals,
#'       MCSE or simulated observation noise. Observation dispersion still enters
#'       the response density.}
#'     \item{raw_diagnostics}{The object's original diagnostics, unchanged.}
#'     \item{zero_density_count}{Number of negative-infinite log densities in
#'       this matrix. These remain explicit zero-density evaluations.}
#'     \item{metadata_establishes_eligibility, IC_eligible}{Always FALSE. Neither
#'       the density nor the dimension report grants likelihood-criterion or
#'       model eligibility.}
#'   }
#' @details This private training-row component includes response offsets,
#'   binomial trial totals, dispersion and original mean/weighted-SQ reference
#'   transforms. It excludes all population priors, Jacobians and latent-normal
#'   density terms. It is neither the posterior objective nor the normalized
#'   full marginal likelihood. No observation or vector is dropped.
#'
#'   Negative infinity is retained as zero response density. Nonfinite original
#'   parameters, nonfinite predictors, NaN or positive-infinite log densities
#'   cause explicit errors. Stored vectors must preserve every original block,
#'   including unselected rows. The return size is the number of selected
#'   vectors times the number of training observations.
#'
#'   Rows share a realized field and any response-group coefficients. These
#'   conditional row densities are not held-out-location or new-group predictive
#'   densities. This function does not compute AIC, BIC, WAIC, LOO or intervals,
#'   qualify a failed fit, or introduce spatial fitting admission. Dimension
#'   accounting alone does not supply a maximized likelihood or an IC penalty.
#' @seealso \code{\link{.msr_spatial_response_predict}},
#'   \code{\link{.msr_spatial_scale_probabilities}}
#' @keywords internal
.msr_spatial_response_log_lik <- function(object, parameters = NULL, draws = NULL) {
  if (!inherits(object, "msr_spatial_scale_prototype"))
    stop("Use a prepared private spatial scale prototype.", call. = FALSE)
  .msr_spatial_scale_objective(object)
  spec <- object$specification; nm <- object$prototype$covariate
  field <- spec$scales[[nm]]$random[[".spatial_field"]]$gaussian_block
  allowed <- vapply(names(spec$scales), function(covariate) {
    effects <- spec$scales[[covariate]]$random
    if (identical(covariate, nm))
      identical(names(effects), ".spatial_field") && length(effects) == 1L else !length(effects)
  }, TRUE)
  if (is.null(field$basis) || !all(allowed))
    stop("Spatial-only response densities require the frozen field and no ordinary scale-group effects.", call. = FALSE)
  blocks <- lengths(spec$parameters)
  if (!is.null(parameters) && !is.null(draws))
    stop("Supply either parameters or draws, not both.", call. = FALSE)
  fixed <- !is.null(parameters)
  ids <- if (fixed) NULL else .msr_spatial_response_joint_ids(object, draws)
  n <- if (fixed) 1L else length(ids)
  values <- matrix(NA_real_, n, length(spec$y))
  for (i in seq_len(n)) {
    p <- if (fixed) parameters else .msr_bayes_draw(object, ids[i])
    if (!is.list(p) || anyDuplicated(names(p)) || !setequal(names(p), names(spec$parameters)) ||
        !all(vapply(names(spec$parameters), function(block) {
          x <- p[[block]]
          is.numeric(x) && !is.complex(x) && is.null(dim(x)) &&
            length(x) == blocks[[block]] && all(is.finite(x))
        }, TRUE)))
      stop("Retain every complete finite original parameter block.", call. = FALSE)
    eta <- .msr_bayes_components(spec, p)$eta
    if (length(eta) != length(spec$y) || any(!is.finite(eta)))
      stop("Nonfinite conditional response predictor; no training rows were dropped.", call. = FALSE)
    density <- .msr_scale_response_loglik(spec$family, spec$y, eta,
      p$log_dispersion, spec$trials)
    if (length(density) != length(spec$y) || any(is.na(density) | density == Inf))
      stop("Invalid conditional response log density; no training rows were dropped.", call. = FALSE)
    values[i, ] <- density
  }
  latent <- names(blocks) %in% c("scale_z", "response_z")
  contract <- list(target = "conditional_training_response_log_density", source_target = object$target,
    parameter_source = if (fixed) "supplied_fixed_vector" else
      if (!is.null(object$importance$resample_ids)) "stored_joint_resamples" else "stored_joint_draws",
    method = object$inference, vector_integration = object$inference_metadata$integration,
    population_priors = "omitted", latent_densities = "omitted", marginalization = "none",
    conditioning = c("observed_training_responses", "frozen_training_design_and_normalization",
      "shared_realized_field_and_response_groups"),
    spatial = list(covariate = nm, basis = "frozen_training_projection_and_normalization",
      range = "fixed", new_location_innovation = FALSE),
    parameter_blocks = blocks, population_dimension = sum(blocks[!latent]),
    latent_dimension = sum(blocks[latent]),
    uncertainty = list(parameters = if (fixed) "fixed_at_supplied_vector" else "joint_stored_vectors",
      intervals = "not_provided", mcse = "not_computed", observation_noise = "density_evaluated_not_simulated"),
    raw_diagnostics = object$diagnostics, zero_density_count = sum(values == -Inf),
    metadata_establishes_eligibility = FALSE, IC_eligible = FALSE)
  attr(values, "observation_ids") <- spec$observation_ids
  attr(values, "training_rows") <- spec$rows
  attr(values, "draw_ids") <- ids
  attr(values, "likelihood_contract") <- contract
  values
}

# Private prerequisites for a variable response correlation range. This does
# not alter the published native fixed-range field or admit a new fitted class.
#' Prepare a private GDR correlation-range likelihood
#'
#' @param data,parameters Complete original finite-K Poisson/half-normal numeric
#'   layout and native parameter blocks. Abundance states are physical field
#'   values here, rather than correlation-whitened coordinates.
#' @param coords,crs Original named survey coordinates and projected CRS.
#' @param range Positive initial exponential correlation distance in CRS units.
#'   The prepared program adds lrange_lambda = log(range) to its parameters.
#' @param max_memory_mb Finite conservative dense-work planning budget in MiB.
#' @return An internal context with a prepare function, initial parameter list,
#'   range coordinate metadata and explicit closed fit/prediction/CV gates.
#' @details The original normalized observation likelihood and observer SD are
#'   unchanged. The abundance density is replaced by its normalized correlated
#'   Gaussian density, including the correlation determinant. One log-range
#'   likelihood coordinate adds no prior or hyperparameter Jacobian. Physical
#'   and whitened joint densities differ by their state-coordinate Jacobian;
#'   their integrated Laplace likelihoods agree. No jitter or range penalty is
#'   added to conceal an unidentified or singular covariance.
#' @noRd
.msr_gdr_range_context <- function(data, parameters, coords, crs, range,
    max_memory_mb = 512) {
  geometry <- .msr_gdr_response_field_context(data, parameters, coords, crs,
    range, max_memory_mb, .geometry_only = TRUE)
  m <- length(parameters$b_lambda)
  if (m < 2L)
    stop("Correlation-range qualification requires at least two distinct abundance states.", call. = FALSE)
  n <- nrow(data$Xlambda); dimension <- sum(lengths(parameters)) + 1
  planned <- 8 * (64 * as.double(m)^2 + 32 * as.double(n) * m +
    32 * as.double(dimension)^2)
  if (!is.finite(planned) || planned > max_memory_mb * 2^20)
    stop("Private GDR range derivative planning exceeds max_memory_mb.", call. = FALSE)
  locations <- unname(as.matrix(coords))[match(seq_len(m), data$block), , drop = FALSE]
  distance <- as.matrix(stats::dist(locations))
  base <- .msr_gdr_gaussian_context(data, parameters)
  original_names <- names(parameters)
  parameters$lrange_lambda <- log(range)
  definition <- geometry$field
  definition$range <- NULL
  definition$loading_sha256 <- NULL
  definition$initial_range <- range
  definition$range_estimated <- TRUE
  definition$range_coordinate <- "lrange_lambda"
  definition$additional_estimated_range_parameters <- 1L
  definition$state_coordinates <- "physical_abundance_field"
  definition$prior <- "original_normalized_Gaussian_SD_with_variable_exponential_correlation"
  implementation <- list(deparse(formals(.msr_gdr_range_context)),
    deparse(body(.msr_gdr_range_context)))
  field_id <- .msr_gdr_hash(list(definition, distance, implementation,
    base$response_fingerprint))
  force(distance); force(original_names); force(field_id); force(definition)
  prepare <- function(Xlambda, row_ids) {
    program <- base$prepare(Xlambda, row_ids)
    evaluate_base <- program$evaluate
    evaluate <- function(p) {
      # The original observation law is reused exactly. Its iid abundance
      # normalizer must be removed before adding the correlated normalizer.
      iid <- sum(RTMB::dnorm(p$b_lambda, 0, exp(p$lsigma_lambda), log = TRUE))
      correlated <- RTMB::dmvnorm(p$b_lambda, mu = 0,
        Sigma = exp(-distance / exp(p$lrange_lambda)),
        scale = exp(p$lsigma_lambda), log = TRUE)
      evaluate_base(p[original_names]) + iid - correlated
    }
    layout <- program$parameter_map
    layout <- rbind(layout, data.frame(position = nrow(layout) + 1L,
      block = "lrange_lambda", within_block = 1L, kind = "log_range"))
    list(evaluate = evaluate, parameter_map = layout,
      fixed_coordinate_count = program$fixed_coordinate_count + 1L,
      latent_coordinate_count = program$latent_coordinate_count,
      range_coordinate_count = 1L, range_estimated = TRUE,
      state_coordinates = "physical_abundance_field", field = definition,
      field_id = field_id, response_fingerprint = base$response_fingerprint,
      fingerprint = .msr_gdr_hash(list(program$fingerprint, field_id)),
      target = "P_finite_K_native_margins_normalized_variable_range_spatial_Gaussian_joint",
      hyperparameter_jacobian_included = FALSE, added_range_prior = FALSE,
      planned_bytes = planned, coordinate_count_is_not_AIC_admission = TRUE,
      fit_admission_established = FALSE, prediction_recipe_qualified = FALSE,
      CV_admission_established = FALSE, requires_serial_R_callbacks = TRUE,
      portable_objective = FALSE, format_version = 1L,
      qualification = "private_variable_range_derivative_prerequisite_only")
  }
  list(prepare = prepare, parameters = parameters, field = definition,
    field_id = field_id, planned_bytes = planned,
    fixed_coordinate_count_change = 1L, annual_field_coordinates = 0L,
    fit_admission_established = FALSE, prediction_recipe_qualified = FALSE,
    CV_admission_established = FALSE, requires_serial_R_callbacks = TRUE,
    portable_objective = FALSE, fingerprint_is_not_full_provider_provenance = TRUE,
    qualification = "private_variable_range_context_not_fitted_model")
}

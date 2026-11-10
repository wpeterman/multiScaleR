# A separate private response-field target. It does not admit a native fitted
# model, replace the observation likelihood, or vary ecological kernel scales.
#' Prepare a private fixed-range GDR abundance response field
#'
#' @param data,parameters Qualified private finite-K GDR numeric inputs.
#' @param coords Two finite projected coordinate columns in the exact survey row
#'   order, with row names equal to the retained survey identities. Each abundance
#'   state must have one distinct location and at least one observed survey.
#' @param crs Known projected two-dimensional CRS describing coords without an
#'   implicit transformation or change of units.
#' @param range Positive fixed exponential correlation distance in those units.
#'   This supplied distance is not an estimated parameter or ecological scale.
#' @param max_memory_mb Positive finite planning budget in MiB. The default is
#'   512. Dense correlation/loading work must fit a conservative estimate. This
#'   estimate is not an operating-system or automatic-differentiation tape cap.
#' @return A private immutable observation context with a prepare function and
#'   bounded new-coordinate projection. Abundance states are correlation-whitened
#'   coordinates with the original normalized Gaussian SD. Conditional new-field
#'   residual variance is reported separately from parameter uncertainty. Public
#'   fitting, prediction, AIC and CV admission remain false.
#' @noRd
.msr_gdr_response_field_context <- function(data, parameters, coords, crs, range,
    max_memory_mb = 512) {
  .msr_gdr_gaussian_data(data, parameters)
  m <- length(parameters$b_lambda); n <- nrow(data$Xlambda)
  if (!m || !length(data$block))
    stop("The private response field needs an explicit abundance-state row map.", call. = FALSE)
  validate_scalar_numeric(range, "range", positive = TRUE)
  validate_scalar_numeric(max_memory_mb, "max_memory_mb", positive = TRUE)
  if (!is.finite(range) || !is.finite(max_memory_mb))
    stop("Use a finite fixed range and memory budget.", call. = FALSE)
  budget <- function(rows) {
    bytes <- 8 * (32 * as.double(m)^2 + 16 * as.double(n + rows) * m)
    if (!is.finite(bytes) || bytes > max_memory_mb * 2^20)
      stop("Private dense GDR response-field work exceeds max_memory_mb.", call. = FALSE)
    bytes
  }
  budget(0)
  declared <- .msr_response_field_crs(crs)
  if ((!is.matrix(coords) && !is.data.frame(coords)) || ncol(coords) != 2L ||
      !is.numeric(as.matrix(coords)) || nrow(coords) != n ||
      any(!is.finite(as.matrix(coords))) || !identical(rownames(coords), data$row_ids$survey))
    stop("Retain finite projected coordinate pairs in exact named GDR survey order.", call. = FALSE)
  coords <- unname(as.matrix(coords))
  first <- match(seq_len(m), data$block)
  if (anyNA(first))
    stop("Every spatial abundance state needs an observed location; unused states need separate qualification.", call. = FALSE)
  locations <- coords[first, , drop = FALSE]
  if (anyDuplicated(locations, MARGIN = 1L) ||
      !all(coords == locations[data$block, , drop = FALSE]))
    stop("Each abundance state must map to one distinct consistent physical location.", call. = FALSE)
  distance <- as.matrix(stats::dist(locations))
  if (any(!is.finite(distance)))
    stop("Response-field distances are not numerically representable.", call. = FALSE)
  correlation <- exp(-distance / range)
  lower <- tryCatch(t(chol(correlation)), error = function(e)
    stop("The fixed response correlation is not numerically positive definite; no jitter is added.", call. = FALSE))
  loading <- lower[data$block, , drop = FALSE]
  base <- .msr_gdr_gaussian_context(data, parameters, .lambda_loading = loading)
  definition <- list(covariance = "exponential", range = range, range_estimated = FALSE,
    crs = declared$wkt, units = declared$units_gdal, state_count = m, survey_count = n,
    coordinate_sha256 = .msr_gdr_hash(list(coords, data$block, data$row_ids$survey)),
    loading_sha256 = .msr_gdr_hash(loading),
    state_coordinates = "correlation_whitened_abundance",
    prior = "original_normalized_Gaussian_SD_in_whitened_coordinates",
    annual_field_coordinates = 0L, additional_estimated_range_parameters = 0L,
    geometry_authentication = "explicit_caller_coordinates_and_exact_survey_state_map")
  field_id <- .msr_gdr_hash(definition)
  response_id <- .msr_gdr_hash(list(base$response_fingerprint, field_id))
  prepare <- function(Xlambda, row_ids) {
    program <- base$prepare(Xlambda, row_ids)
    program$field <- definition; program$field_id <- field_id
    program$response_fingerprint <- response_id
    program$fingerprint <- .msr_gdr_hash(list(program$fingerprint, field_id))
    program$target <- "P_finite_K_native_margins_normalized_fixed_range_spatial_Gaussian_joint"
    program$format_version <- 3L
    program$qualification <- "private_fixed_range_spatial_derivative_prerequisite_only"
    program
  }
  project <- function(newcoords, row_ids, crs) {
    if ((!is.matrix(newcoords) && !is.data.frame(newcoords)) || ncol(newcoords) != 2L ||
        !is.numeric(as.matrix(newcoords)) || any(!is.finite(as.matrix(newcoords))) ||
        !is.character(row_ids) || length(row_ids) != nrow(newcoords) ||
        anyNA(row_ids) || any(!nzchar(row_ids)))
      stop("Supply finite projected new coordinates and explicit aligned row identities.", call. = FALSE)
    if (!isTRUE(declared == .msr_response_field_crs(crs)))
      stop("New field coordinates must retain the declared CRS and linear units.", call. = FALSE)
    budget(nrow(newcoords))
    xy <- unname(as.matrix(newcoords)); rows <- nrow(xy)
    if (!rows) return(list(row_ids = row_ids, loading = matrix(numeric(), 0L, m),
      conditional_unit_variance = numeric(), variance_roundoff_corrections = 0L,
      field_id = field_id, full_prediction_covariance_retained = FALSE,
      variance_scale = "multiply_by_exp(2*lsigma_lambda)",
      uncertainty = "new_field_residual_conditional_on_parameters_and_training_field",
      omitted = c("parameter_uncertainty", "training_field_uncertainty",
        "ecological_scale_uncertainty", "observation_noise", "cross_row_field_covariance"),
      fit_admission_established = FALSE, public_prediction_admission = FALSE))
    cross <- sqrt(outer(xy[, 1L], locations[, 1L], "-")^2 +
      outer(xy[, 2L], locations[, 2L], "-")^2)
    if (any(!is.finite(cross)))
      stop("New response-field distances are not numerically representable.", call. = FALSE)
    weights <- t(forwardsolve(lower, t(exp(-cross / range))))
    variance <- 1 - rowSums(weights^2)
    known <- vapply(seq_len(rows), function(i) {
      found <- which(locations[, 1L] == xy[i, 1L] & locations[, 2L] == xy[i, 2L])
      if (length(found)) found[[1L]] else NA_integer_
    }, 1L)
    if (any(!is.na(known))) {
      weights[!is.na(known), ] <- lower[known[!is.na(known)], , drop = FALSE]
      variance[!is.na(known)] <- 0
    }
    allowance <- 64 * .Machine$double.eps * m * pmax(1, rowSums(weights^2))
    if (any(variance < -allowance))
      stop("Conditional response-field variance is numerically inconsistent.", call. = FALSE)
    clipped <- sum(variance < 0)
    list(row_ids = row_ids, loading = weights,
      conditional_unit_variance = pmax(0, variance), variance_roundoff_corrections = clipped,
      field_id = field_id, full_prediction_covariance_retained = FALSE,
      variance_scale = "multiply_by_exp(2*lsigma_lambda)",
      uncertainty = "new_field_residual_conditional_on_parameters_and_training_field",
      omitted = c("parameter_uncertainty", "training_field_uncertainty",
        "ecological_scale_uncertainty", "observation_noise", "cross_row_field_covariance"),
      fit_admission_established = FALSE, public_prediction_admission = FALSE)
  }
  # A prediction worker needs numeric geometry only, not the observation
  # preparation, RTMB callbacks or a constructor's complete call frame.
  projection_env <- list2env(list(m = m, n = n, max_memory_mb = max_memory_mb,
    declared = declared, locations = locations, range = range, lower = lower,
    field_id = field_id), parent = environment(.msr_gdr_response_field_context))
  numeric_budget <- budget; environment(numeric_budget) <- projection_env
  projection_env$budget <- numeric_budget
  environment(project) <- projection_env
  list(prepare = prepare, project = project, field = definition, field_id = field_id,
    response_fingerprint = response_id,
    fixed_coordinate_count_change = 0L, annual_field_coordinates = 0L,
    shared_observation_preparation = TRUE, requires_serial_R_callbacks = TRUE,
    numeric_projection_portable = TRUE,
    portable_objective = FALSE, fingerprint_is_not_full_provider_provenance = TRUE,
    fit_admission_established = FALSE, prediction_recipe_qualified = FALSE,
    qualification = "private_fixed_range_response_field_context")
}

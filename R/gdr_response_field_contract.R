.msr_gdr_field_contract <- function(model, zero = FALSE) {
  # This describes an already-bound prediction; it is not fit admission and
  # must not rehash provider DLLs once per raster sub-batch.
  field <- model@response_field
  list(provider = "native_unmarked_TMB", provider_version = field$provider$TMB_version,
    fields = list(field$definition), field_id = field$field_id,
    conditioning = if (zero) "all_response_random_effects_set_to_zero" else
      "conditional_fitted_field_with_identified_new_location_residual",
    target = "native_finite_K_P_halfnormal_distance_removal_with_normalized_Gaussian_field",
    fixed_range = field$definition$range, range_estimated = FALSE,
    additional_estimated_range_parameters = 0L, annual_field_coordinates = 0L,
    physical_field_coefficients = TRUE, kernel_uncertainty_included = FALSE,
    new_location_residual_variance = if (zero) "excluded_in_zero_effect_reference" else
      "separate_component_in_conditional_predictor_SE",
    response_point = "inverse_link_predictor_not_lognormal_field_average",
    observation_noise_included = FALSE, full_grid_covariance_returned = FALSE,
    descriptor_is_not_new_fit_admission = TRUE)
}

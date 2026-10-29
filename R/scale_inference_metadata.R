# Describe the inferential target independently of the solver and its gates.
# This is descriptive metadata, never a replacement for stored diagnostics.
.msr_scale_inference_metadata <- function(target, method) {
  .msr_scale_capability(target, method)
  posterior <- target == "posterior"
  prepared <- method == "prepare"
  integration <- switch(method,
    prepare = c(scale_effects = "not_computed", response_effects = "not_computed",
      population_parameters = "not_computed"),
    nuts = c(scale_effects = "joint_sampling", response_effects = "joint_sampling",
      population_parameters = "joint_sampling"),
    joint_laplace = c(scale_effects = "joint_gaussian", response_effects = "joint_gaussian",
      population_parameters = "joint_gaussian"),
    marginal_laplace = c(scale_effects = "conditional_gaussian", response_effects = "conditional_gaussian",
      population_parameters = "scale_sd_quadrature_and_conditional_gaussian"),
    importance_corrected = c(scale_effects = "joint_importance_sampling",
      response_effects = "joint_importance_sampling", population_parameters = "joint_importance_sampling"),
    laplace = c(scale_effects = "laplace", response_effects = "laplace",
      population_parameters = "optimized"),
    factorized = c(scale_effects = "independent_annual_quadrature",
      response_effects = "exact_gaussian_when_present", population_parameters = "optimized"))
  list(schema_version = 1L, target = target, method = method,
    objective = if (posterior) "negative_log_posterior" else "negative_marginal_log_likelihood",
    population_priors = if (posterior) "proper" else "none",
    response_normalizers = "retained", integration = as.list(integration),
    solver = if (prepared) "none" else if (method == "nuts") "nuts" else
      if (method == "importance_corrected") "importance_sampling" else "optimization",
    uncertainty = list(
      intervals = if (posterior && !prepared) "marginal_credible" else "not_provided",
      population_parameters = if (prepared) "not_computed" else
        if (posterior) "propagated_in_joint_posterior" else "not_propagated",
      annual_scales = if (prepared) "not_computed" else if (posterior) "joint_posterior" else
        if (method == "factorized") "conditional_log_scale_modes_and_mean_logs_at_fitted_population" else
          "conditional_log_scale_modes_at_fitted_population",
      calibration = "not_established"),
    diagnostics_field = "diagnostics", metadata_establishes_eligibility = FALSE)
}

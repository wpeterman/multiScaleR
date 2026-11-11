# A normalized Gaussian field with SD zero has all its mass at zero. Build
# its exact parent by removing only the selected latent blocks from the
# conditional design and omitting their integrated Gaussian factors.
.msr_gdr_gaussian_zero_sd_parent <- function(data, parameters, fields) {
  .msr_gdr_gaussian_data(data, parameters)
  maps <- c(lambda = "block", dist = "observer")
  if (!is.character(fields) || !length(fields) || anyNA(fields) ||
      anyDuplicated(fields) || !all(fields %in% names(maps)) ||
      any(lengths(data[maps[fields]]) == 0L))
    stop("Choose existing abundance or distance Gaussian fields for the exact zero-SD parent.", call. = FALSE)
  for (field in fields) {
    data[[maps[[field]]]] <- integer()
    parameters[[paste0("b_", field)]] <- numeric()
    parameters[[paste0("lsigma_", field)]] <- numeric()
  }
  list(program = .msr_gdr_gaussian_program(data, parameters), parameters = parameters,
    zero_sd_fields = fields, boundary = "normalized_Gaussian_mass_at_zero",
    fit_admission_established = FALSE, uncertainty_evaluated = FALSE)
}

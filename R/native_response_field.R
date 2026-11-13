# Identify native response covariance without copying its dense distance matrix.
# This describes prediction conditioning; it does not admit joint scale fitting.
.msr_glmmTMB_response_fields <- function(model, zero = FALSE, allow_new = FALSE) {
  if (!inherits(model, "glmmTMB")) return(NULL)
  blocks <- model$modelInfo$reStruc$condReStruc
  if (!length(blocks)) return(NULL)
  code <- vapply(blocks, function(x) if (is.numeric(x$blockCode) &&
    length(x$blockCode) == 1L && is.finite(x$blockCode)) as.integer(x$blockCode) else NA_integer_, 1L)
  provider <- asNamespace("glmmTMB")
  version <- as.character(getNamespaceVersion("glmmTMB"))
  registered <- get0(".valid_covstruct", envir = provider, inherits = FALSE)
  covariance_names <- c(exp = "exponential", gau = "gaussian", mat = "matern")
  if (!is.numeric(registered) || is.null(names(registered)) ||
      anyDuplicated(names(registered)) || anyDuplicated(registered) ||
      !all(names(covariance_names) %in% names(registered)))
    return(list(provider = "glmmTMB", provider_version = version,
      status = "unqualified_provider_layout", fit_eligibility_established = FALSE))
  field_code <- registered[names(covariance_names)]
  if (any(!is.finite(field_code) | field_code != floor(field_code)))
    return(list(provider = "glmmTMB", provider_version = version,
      status = "unqualified_provider_layout", fit_eligibility_established = FALSE))
  structure_index <- match(code, unname(field_code))
  selected <- which(!is.na(structure_index))
  if (!length(selected)) return(NULL)
  fields <- lapply(selected, function(i) {
    block <- blocks[[i]]
    covariance <- unname(covariance_names[structure_index[i]])
    term <- names(blocks)[i]
    expression <- tryCatch(str2lang(term), error = function(e) NULL)
    scalar_count <- function(x) is.numeric(x) && length(x) == 1L &&
      is.finite(x) && x >= 1 && x == floor(x)
    if (!scalar_count(block$blockSize) || !scalar_count(block$blockReps) ||
        !scalar_count(block$blockNumTheta) ||
        block$blockNumTheta != (if (covariance == "matern") 3L else 2L) ||
        !identical(dim(block$dist), rep(as.integer(block$blockSize), 2L)) ||
        !is.call(expression) || !identical(expression[[1L]], as.name("|")) ||
        length(expression) != 3L)
      stop("Native response-field covariance metadata is incomplete; retain the original provider context.", call. = FALSE)
    list(term = term, covariance = covariance,
      coordinate_variables = all.vars(expression[[2L]]),
      grouping_variables = all.vars(expression[[3L]]),
      locations_per_field = as.integer(block$blockSize),
      independent_fields = as.integer(block$blockReps),
      covariance_parameter_count = as.integer(block$blockNumTheta))
  })
  list(provider = "glmmTMB", provider_version = version, status = "native_layout_identified",
    component = "conditional_response", fields = fields,
    conditioning = if (zero) "all_response_random_effects_set_to_zero" else
      "native_conditional_spatial_prediction",
    new_locations = if (allow_new) "native_new_location_handling_requested" else
      "native_default_new_location_warning_policy",
    coordinate_units = "native_coordinate_units_not_inferred",
    range_is_ecological_kernel_scale = FALSE, kernel_uncertainty_included = FALSE,
    fit_eligibility_established = FALSE)
}

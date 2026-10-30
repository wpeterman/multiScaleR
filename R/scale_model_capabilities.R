# Structural support is separate from package availability, numerical eligibility
# and interval calibration. These rows describe the existing implementations;
# they never select a fallback method or change a diagnostic decision.
.msr_scale_capabilities <- function() {
  data.frame(
    target = c(rep("posterior", 5L), rep("marginal_likelihood", 2L)),
    method = c("prepare", "nuts", "joint_laplace", "marginal_laplace",
      "importance_corrected", "laplace", "factorized"),
    max_kernel_predictors = c(rep(Inf, 5L), 1, 1),
    max_scale_groups_per_predictor = c(rep(Inf, 5L), 1, 1),
    max_response_groups = c(rep(Inf, 5L), 1, 1),
    required_scale_sd_count = c(NA_integer_, NA_integer_, NA_integer_, 1L, 1L,
      NA_integer_, NA_integer_),
    response_effects = c(rep("joint", 6L), "gaussian_within_scale_group"),
    public_interface = c(rep(TRUE, 5L), FALSE, FALSE),
    stringsAsFactors = FALSE)
}

.msr_scale_capability <- function(target, method) {
  if (length(target) != 1L || is.na(target) ||
      !target %in% c("posterior", "marginal_likelihood"))
    stop("Unknown scale inference target.", call. = FALSE)
  if (length(method) != 1L || is.na(method))
    stop("Integration method is incompatible with the scale inference target.", call. = FALSE)
  table <- .msr_scale_capabilities()
  row <- which(table$target == target & table$method == method)
  if (length(row) != 1L)
    stop("Integration method is incompatible with the scale inference target.", call. = FALSE)
  table[row, , drop = FALSE]
}

# Inspect only the frozen model. Optional packages and all numerical gates stay
# with their original solvers. Return the specification unchanged for wrappers.
.msr_scale_check_capability <- function(spec, target, method) {
  capability <- .msr_scale_capability(target, method)
  if (target == "marginal_likelihood" && length(spec$surface$buffers))
    stop("Weighted SQ is not yet supported by private marginal-likelihood fit and limiting adapters.", call. = FALSE)
  if (.msr_has_gaussian_blocks(spec) &&
      (target == "marginal_likelihood" || method %in% c("marginal_laplace", "importance_corrected")))
    stop("Gaussian projection blocks require the private joint spatial adapter; current scalar and likelihood fit interfaces support grouping intercepts only.", call. = FALSE)
  # Family, kernel, data and formula validation belongs to shared preparation.
  # Optimizer adapters also use this structural guard with reduced test inputs.
  if (!length(spec$covariates) || length(spec$covariates) > capability$max_kernel_predictors ||
      any(vapply(spec$scales, function(sm) length(sm$random), 0L) >
        capability$max_scale_groups_per_predictor) ||
      length(spec$response$random) > capability$max_response_groups) {
    message <- if (method == "factorized")
        "Factorized likelihood requires one kernel and at most one scale and one response grouping term." else
        "The private marginal prototype supports one kernel and at most one scale and one response grouping term."
    stop(message, call. = FALSE)
  }
  required <- capability$required_scale_sd_count
  if (!is.na(required) && length(spec$parameters$log_scale_sd) != required)
    stop("Marginal Laplace requires exactly one scale-group SD; use one scale random-intercept term.", call. = FALSE)
  if (capability$response_effects == "gaussian_within_scale_group") {
    response <- spec$response$random
    if (length(response)) {
      if (spec$family != "gaussian")
        stop("Factorized response random intercepts currently require a Gaussian response.", call. = FALSE)
      scale <- spec$scales[[1L]]$random
      if (length(scale) && any(vapply(split(scale[[1L]]$index, response[[1L]]$index),
        function(i) length(unique(i)) > 1L, TRUE)))
        stop("Response random effects coupling different scale groups invalidate factorization.", call. = FALSE)
    } else if (length(spec$parameters$response_z))
      stop("Response latent parameters have no corresponding grouping term.", call. = FALSE)
  }
  invisible(spec)
}

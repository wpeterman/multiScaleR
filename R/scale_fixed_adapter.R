# Fixed spatial features are ordinary response/scale-design columns. This proof
# binds their original prepared values and normalization without adding a free
# radius, derivative, prior, or native fitted model.
.msr_fixed_reference_contract <- function(inputs, data, rows, covariates) {
  specs <- inputs$scale_vars
  empty <- list(covariates = character(), center = numeric(), spread = numeric())
  # Legacy mean-only inputs predate scale_vars and have no declared fixed metrics.
  if (is.null(specs)) return(empty)
  fixed <- specs$covariate[!specs$optimize & specs$covariate %in% covariates]
  if (!length(fixed)) return(empty)
  validate_scalar_numeric(inputs$max_D, "training max_D", positive = TRUE)
  if (!is.data.frame(data) || !all(fixed %in% names(data)) ||
      !all(fixed %in% names(inputs$kernel_dat)))
    stop("Automatic fixed predictors require retained prepared and native numeric covariates.", call. = FALSE)
  if (length(rows) != nrow(data) || anyNA(rows) ||
      any(rows < 1 | rows > nrow(inputs$kernel_dat)))
    stop("Fixed predictors require verified observation-to-site rows.", call. = FALSE)
  center <- inputs$scl_params$mean[fixed]
  spread <- inputs$scl_params$sd[fixed]
  if (length(center) != length(fixed) || any(!is.finite(center)) ||
      length(spread) != length(fixed) || any(!is.finite(spread) | spread <= 0))
    stop("Fixed predictors require finite retained training means and positive standard deviations.", call. = FALSE)
  context <- list(covs = character(), n_covs = 0L, scale_vars = specs,
    unit_conv = inputs$unit_conv, resolution = inputs$resolution,
    n_cols = inputs$n_cols, binned = inputs$binned)
  raw <- .msr_kernel_cov_w(numeric(), inputs$d_list, inputs$raw_cov,
    inputs$kernel, context, covariates = fixed)
  scaled <- sweep(sweep(raw, 2L, center), 2L, spread, "/")
  for (nm in fixed) {
    if (!is.numeric(data[[nm]]) || !is.numeric(inputs$kernel_dat[[nm]]) ||
        !isTRUE(all.equal(as.numeric(data[[nm]]),
          as.numeric(inputs$kernel_dat[[nm]][rows]), tolerance = 0,
          check.attributes = FALSE)) ||
        !isTRUE(all.equal(as.numeric(inputs$kernel_dat[[nm]]),
          as.numeric(scaled[, nm]), tolerance = 0, check.attributes = FALSE)))
      stop("Fixed predictor '", nm,
        "' does not match its retained training preparation; supply it explicitly.", call. = FALSE)
  }
  list(covariates = fixed, center = center, spread = spread,
    specifications = specs[match(fixed, specs$covariate), , drop = FALSE],
    kernel = inputs$kernel, resolution = inputs$resolution, max_D = inputs$max_D)
}

# Evaluate each prepared site once, then map it to explicit prediction rows.
# A new distance normalizer changes numerical distances, not physical radii.
.msr_fixed_reference_values <- function(reference, inputs, rows) {
  fixed <- reference$covariates
  if (!length(fixed)) return(matrix(numeric(), length(rows), 0L))
  if (!identical(inputs$kernel, reference$kernel))
    stop("Fixed prediction inputs must use the fitted kernel.", call. = FALSE)
  if (is.null(reference$max_D))
    stop("Fixed prediction requires the original training buffer extent.", call. = FALSE)
  validate_scalar_numeric(inputs$max_D, "kernel_inputs$max_D", positive = TRUE)
  if (!isTRUE(all.equal(inputs$max_D, reference$max_D, tolerance = 0)))
    stop("Fixed prediction buffer extent must match training max_D.", call. = FALSE)
  validate_scalar_numeric(inputs$unit_conv, "kernel_inputs$unit_conv", positive = TRUE)
  if (!isTRUE(all.equal(inputs$resolution, reference$resolution, tolerance = 0)))
    stop("Fixed-radius prediction requires the prepared training resolution.", call. = FALSE)
  if (is.null(inputs$scale_vars))
    stop("Fixed-radius source and metric definitions must match the fit.", call. = FALSE)
  actual <- inputs$scale_vars[match(fixed, inputs$scale_vars$covariate), , drop = FALSE]
  if (anyNA(actual$covariate) || !isTRUE(all.equal(as.data.frame(actual),
      as.data.frame(reference$specifications), tolerance = 0, check.attributes = FALSE)))
    stop("Fixed-radius source and metric definitions must match the fit.", call. = FALSE)
  context <- list(covs = character(), n_covs = 0L,
    scale_vars = reference$specifications, unit_conv = inputs$unit_conv,
    resolution = inputs$resolution, n_cols = inputs$n_cols, binned = inputs$binned)
  raw <- .msr_kernel_cov_w(numeric(), inputs$d_list, inputs$raw_cov,
    inputs$kernel, context, covariates = fixed)
  if (any(is.infinite(raw))) stop("Spatial prediction values cannot be infinite.", call. = FALSE)
  scaled <- sweep(sweep(raw, 2L, reference$center), 2L, reference$spread, "/")
  scaled[rows, , drop = FALSE]
}

.msr_fixed_point_data <- function(object, newdata, kernel_inputs, include_fixed) {
  validate_scalar_logical(include_fixed, "include_fixed")
  if (!include_fixed || (is.null(newdata) && is.null(kernel_inputs))) return(newdata)
  if (is.null(kernel_inputs))
    stop("Automatic fixed prediction requires explicit kernel_inputs.", call. = FALSE)
  if (is.null(newdata)) newdata <- data.frame(row.names = rownames(kernel_inputs$kernel_dat))
  if (!is.data.frame(newdata)) stop("newdata must be a data frame.", call. = FALSE)
  rows <- .msr_profile_input_rows(newdata, kernel_inputs)
  reference <- .msr_fixed_training_contract(object)
  if (length(reference$covariates))
    newdata[reference$covariates] <- as.data.frame(.msr_fixed_reference_values(reference, kernel_inputs, rows))
  newdata
}

.msr_attach_fixed_reference <- function(object) {
  reference <- .msr_fixed_training_contract(object)
  if (length(reference$covariates)) {
    object$fixed_reference <- reference
    object$scl_params$mean[reference$covariates] <- reference$center
    object$scl_params$sd[reference$covariates] <- reference$spread
  }
  object
}

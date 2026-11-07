# Retain point estimates when outer curvature cannot support Wald uncertainty.
# Hard-radius positions and inactive shapes have no smooth Wald interpretation.
.msr_refit_hessian_se <- function(hessian, active, multiplier) {
  se <- rep(NA_real_, length(active))
  reason <- "no_active_smooth_parameters"
  if (any(active)) {
    reason <- "invalid_hessian"
    if (is.matrix(hessian) && is.numeric(hessian) &&
        identical(dim(hessian), rep(length(active), 2L))) {
      h <- hessian[active, active, drop = FALSE]
      if (all(is.finite(h)) && isTRUE(all.equal(h, t(h), tolerance = 1e-10))) {
        ch <- tryCatch(chol(h), error = function(e) NULL)
        if (!is.null(ch)) {
          covariance <- chol2inv(ch)
          candidate <- sqrt(diag(covariance)) * multiplier[active]
          if (all(is.finite(candidate)) && all(candidate > 0)) {
            se[active] <- candidate
            reason <- "available"
          }
        }
      }
    }
  }
  list(SE = se, diagnostic = list(code = "refit_hessian",
    triggered = any(active) && reason != "available", status = reason,
    active = active, uncertainty = "outer_refit_wald"))
}

# Numerical log-gap tolerance is independent of the distance unit.
.msr_refit_bounds <- function(parameter, lower, upper, active, multiplier,
    parameter_names, tolerance = 1e-5) {
  n <- length(parameter)
  numeric_vector <- function(x) is.numeric(x) && is.null(dim(x)) &&
    length(x) == n && all(is.finite(x))
  if (!n || !numeric_vector(parameter) || !numeric_vector(lower) ||
      !numeric_vector(upper) || !numeric_vector(multiplier) ||
      any(parameter < 0 | lower < 0 | upper <= lower | multiplier <= 0) ||
      !is.logical(active) || length(active) != n || anyNA(active) ||
      !is.character(parameter_names) || length(parameter_names) != n ||
      anyNA(parameter_names) || any(!nzchar(parameter_names)) || anyDuplicated(parameter_names))
    stop("Finite aligned refit parameters, bounds, activity and names are required.", call. = FALSE)
  validate_scalar_numeric(tolerance, "boundary log-gap tolerance", positive = TRUE)
  # A zero lower bound cannot be approached on the log scale. Flag it only
  # when reached; positive bounds use the same log-gap scale as analytic fits.
  at_lower <- parameter <= lower
  positive <- parameter > 0 & lower > 0
  at_lower[positive] <- log(parameter[positive]) - log(lower[positive]) <= tolerance
  at_upper <- parameter >= upper
  positive <- parameter > 0
  at_upper[positive] <- log(upper[positive]) - log(parameter[positive]) <= tolerance
  list(active_boundary = any(active & (at_lower | at_upper)),
    tolerance = tolerance, tolerance_scale = "absolute_log_parameter_gap",
    parameters = data.frame(parameter = parameter_names,
      estimate = parameter * multiplier, lower = lower * multiplier,
      upper = upper * multiplier, at_lower = at_lower, at_upper = at_upper,
      smooth = active, row.names = parameter_names))
}

# Preserve the original curvature result separately. A bound invalidates the
# unconstrained joint outer Wald block, including correlated interior scales.
.msr_refit_bound_uncertainty <- function(uncertainty, bounds) {
  if (is.null(uncertainty$diagnostic$curvature_status))
    uncertainty$diagnostic$curvature_status <- uncertainty$diagnostic$status
  uncertainty$diagnostic$bounds <- bounds
  if (isTRUE(bounds$active_boundary)) {
    uncertainty$SE[uncertainty$diagnostic$active] <- NA_real_
    uncertainty$diagnostic$status <- "parameter_at_optimization_bound"
    uncertainty$diagnostic$triggered <- TRUE
  }
  uncertainty
}

# Reporting copy for older saved refit fits with recoverable optimizer bounds.
# Other backends and incomplete legacy metadata retain their own contracts.
.msr_refit_bound_view <- function(object) {
  if (!is.null(object$backend) && !identical(object$backend, "refit")) return(object)
  kernel <- object$kernel_inputs$kernel
  if (is.null(kernel) || !kernel %in% c("gaussian", "exp", "expow", "fixed")) return(object)
  covariates <- rownames(object$scale_est)
  n <- length(covariates)
  p <- n * if (kernel == "expow") 2L else 1L
  parameter <- object$optim_results$par
  unit <- object$kernel_inputs$unit_conv
  valid_scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!n || !is.numeric(parameter) || !is.null(dim(parameter)) || length(parameter) != p ||
      any(!is.finite(parameter)) || any(parameter < 0) || !valid_scalar(unit) || unit <= 0 ||
      !valid_scalar(object$min_D) || !valid_scalar(object$max_D) ||
      object$min_D < 0 || object$max_D <= object$min_D ||
      (kernel == "expow" && (!is.data.frame(object$shape_est) ||
        nrow(object$shape_est) != n))) return(object)
  active <- rep(!covariates %in% .msr_hard_radius_covariates(
    object$kernel_inputs$scale_vars, kernel, covariates), length.out = p)
  multiplier <- c(rep(unit, n), if (kernel == "expow") rep(1, n))
  bounds <- .msr_refit_bounds(parameter,
    c(rep(object$min_D / unit, n), if (kernel == "expow") rep(.75, n)),
    c(rep(object$max_D / unit, n), if (kernel == "expow") rep(50, n)),
    active, multiplier, c(paste0("sigma:", covariates),
      if (kernel == "expow") paste0("shape:", covariates)))
  if (!isTRUE(bounds$active_boundary)) return(object)
  # Preserve point estimates, raw Hessians, native models and parameter counts.
  se <- c(object$scale_est$SE, if (kernel == "expow") object$shape_est$SE)
  diagnostic <- object$diagnostics$refit_hessian
  if (is.null(diagnostic)) diagnostic <- list(code = "refit_hessian", triggered = FALSE,
    status = "legacy_stored_uncertainty", active = active, uncertainty = "outer_refit_wald")
  result <- .msr_refit_bound_uncertainty(list(SE = se, diagnostic = diagnostic), bounds)
  object$scale_est$SE <- result$SE[seq_len(n)]
  if (kernel == "expow") object$shape_est$SE <- result$SE[n + seq_len(n)]
  object$diagnostics$refit_hessian <- result$diagnostic
  if (!is.null(object$profile_scale_est)) {
    tab <- object$profile_scale_est
    smooth <- intersect(covariates[active[seq_len(n)]], rownames(tab))
    methods <- attr(tab, "interval_method")
    keep <- if (is.null(methods)) character() else smooth[
      !is.na(methods[smooth]) & grepl("^profile", methods[smooth])]
    if ("SE" %in% names(tab)) tab[smooth, "SE"] <- NA_real_
    withheld <- setdiff(smooth, keep)
    if (length(withheld)) tab[withheld, c("2.5%", "97.5%")] <- NA_real_
    object$profile_scale_est <- tab
  }
  object
}

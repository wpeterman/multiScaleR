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

# Per-scale Wald SEs for the joint_mixed backend.
#
# The selected start's TMB sdreport holds `cov.fixed`, the inverse of the
# Hessian of the Laplace-approximated negative log marginal likelihood over
# every fixed parameter (kernel_par, beta, gamma, log_sd, log_theta). The
# kernel_par entries are log(sigma / unit_conv), so their SEs are log-scale
# SEs and the delta method gives SE(sigma) = sigma x SE(log sigma).
#
# The withholding rule follows .msr_profile_scale_se():
# * A missing sdreport, a Hessian that is not positive definite, or unusable
#   covariance withholds every scale SE.
# * A log scale within `tolerance` of a search bound loses its own SE. The
#   other SEs come from the information matrix with the bounded scales
#   dropped, which holds them fixed at their bounds. `conditional_on` names
#   those bounded scales when some other scale keeps an SE.
# * An interior scale also loses its SE when its log-scale SE exceeds one, or
#   when another eligible start within the 95% likelihood threshold sits more
#   than 0.5 log units away on that scale.
# Coefficient and prediction uncertainty are untouched; they still condition
# on the selected scales.
.msr_mixed_scale_se <- function(sdreport, log_scale, lower, upper, outcomes,
                                covariates, tolerance = 1e-5) {
  p <- length(log_scale)
  at_bound <- log_scale - lower < tolerance | upper - log_scale < tolerance
  out <- list(status = "available", pd_hessian = NA,
    log_se = stats::setNames(rep(NA_real_, p), covariates),
    reasons = stats::setNames(rep(NA_character_, p), covariates),
    at_bound = stats::setNames(at_bound, covariates),
    conditional_on = character())
  withhold_all <- function(status, reason) {
    out$status <- status
    out$reasons[] <- reason
    out
  }
  V <- if (is.list(sdreport)) sdreport$cov.fixed else NULL
  if (is.null(V))
    return(withhold_all("no_sdreport", "no joint sdreport for the selected start"))
  out$pd_hessian <- isTRUE(sdreport$pdHess)
  if (!out$pd_hessian)
    return(withhold_all("hessian_not_positive_definite",
                        "joint Hessian is not positive definite"))
  ki <- if (is.matrix(V)) which(rownames(V) == "kernel_par") else integer()
  if (!is.numeric(V) || nrow(V) != ncol(V) || length(ki) != p || any(!is.finite(V)))
    return(withhold_all("invalid_covariance", "joint covariance is unusable"))
  if (all(at_bound))
    return(withhold_all("all_scales_at_bound", "parameter at an optimization bound"))
  # Rebuild the observed information and drop the bounded scales. With no
  # bounded scale this returns cov.fixed; the shared condition check applies
  # either way.
  V <- (V + t(V)) / 2
  H <- tryCatch(solve(V), error = function(e) NULL)
  keep <- setdiff(seq_len(nrow(V)), ki[at_bound])
  vcov <- if (is.null(H)) NULL else
    .msr_information_inverse(H[keep, keep, drop = FALSE])$vcov
  if (is.null(vcov))
    return(withhold_all("ill_conditioned_information",
                        "nonpositive or ill-conditioned observed information"))
  free <- which(!at_bound)
  variance <- diag(vcov)[match(ki[free], keep)]
  log_se <- rep(NA_real_, p)
  log_se[free] <- ifelse(is.finite(variance) & variance > 0, sqrt(pmax(variance, 0)), NA_real_)
  # Competing starts: eligible starts within the 95% likelihood threshold of
  # the selected objective that moved a scale by more than 0.5 log units.
  eligible <- vapply(outcomes, function(x) isTRUE(x$eligible), logical(1))
  objective <- vapply(outcomes, function(x)
    if (isTRUE(x$eligible)) as.numeric(x$objective) else NA_real_, numeric(1))
  good <- if (any(eligible)) which(eligible &
    objective <= min(objective, na.rm = TRUE) + stats::qchisq(.95, 1) / 2) else integer()
  moved <- vapply(seq_len(p), function(k) any(vapply(good, function(i)
    abs(outcomes[[i]]$log_scale[k] - log_scale[k]) > 0.5, logical(1))), logical(1))
  reasons <- rep("available", p)
  reasons[!is.finite(log_se)] <- "nonpositive or ill-conditioned observed information"
  reasons[is.finite(log_se) & log_se > 1] <- "log-kernel standard error exceeds one"
  reasons[!at_bound & moved] <- "competing scale modes within the 95% likelihood threshold"
  reasons[at_bound] <- "parameter at an optimization bound"
  log_se[reasons != "available"] <- NA_real_
  out$log_se[] <- log_se
  out$reasons[] <- reasons
  reported <- is.finite(log_se)
  out$status <- if (all(reported)) "available" else if (any(reported)) "partial" else "withheld"
  if (any(reported) && any(at_bound))
    out$conditional_on <- paste0("sigma:", covariates[at_bound])
  out
}

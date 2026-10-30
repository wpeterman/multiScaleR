# Fresh physical-vector audits distinguish conditional nuisance stationarity
# from the unrestricted full-parameter stationarity requirement.
.msr_response_sd_statistic_audit <- function(spec, parameters, zero, objective,
    fixed = FALSE) {
  tryCatch({
    reduced <- if (zero) .msr_marginal_reduce(spec, "response") else spec
    parameters$scale_z <- parameters$response_z <- numeric()
    ordinary_object <- .msr_factorized_adfun(reduced, parameters, step = 1)
    par <- ordinary_object$par
    fixed_indices <- if (fixed) which(names(par) == "log_response_sd") else integer()
    free <- setdiff(seq_along(par), fixed_indices)
    value <- as.numeric(ordinary_object$fn(par))
    gradient <- as.numeric(ordinary_object$gr(par))
    hessian <- ordinary_object$he(par)
    cold <- .msr_factorized_adfun(reduced, parameters, rel.tol = 1e-10, step = .5)
    if (!identical(names(par), names(cold$par)) ||
        !identical(as.numeric(par), as.numeric(cold$par)))
      stop("Cold audit changed the physical parameter vector.", call. = FALSE)
    cold_value <- as.numeric(cold$fn(par))
    cold_gradient <- as.numeric(cold$gr(par))
    ordinary <- .msr_factorized_objective(reduced, FALSE, rel.tol = 1e-10, step = .5)(parameters)
    geometry <- .msr_marginal_curvature(hessian[free, free, drop = FALSE])
    value_error <- max(abs(c(value - cold_value, value - ordinary, objective - ordinary)))
    score_error <- max(abs(gradient[free] - cold_gradient[free]))
    score_max <- max(abs(c(gradient[free], cold_gradient[free])))
    passed <- all(is.finite(c(objective, value, cold_value, ordinary, gradient,
      cold_gradient, hessian))) && value_error <= 1e-6 && score_error <= 1e-4 &&
      score_max <= 1e-4 && geometry$passed
    list(passed = passed, stored_objective = objective, ordinary_objective = as.numeric(ordinary),
      AD_objective = value, cold_objective = cold_value, value_error = value_error,
      score_error = score_error, score_max = score_max, condition = geometry$condition,
      minimum_eigenvalue = geometry$minimum, full_gradient = gradient,
      full_hessian = hessian, parameter_blocks = names(par), physical_vector = par,
      free_indices = free, fixed_indices = fixed_indices, zero_sd_response = zero,
      stationarity = if (fixed) "nuisance only" else "all retained population coordinates")
  }, error = function(e) list(passed = FALSE, stored_objective = objective,
    ordinary_objective = NA_real_, error = conditionMessage(e)))
}

# Preserve every audited feasible objective, including an ineligible nuisance
# attempt or an exact limiting model that prevents finite identification.
.msr_response_sd_statistic_witnesses <- function(profile) {
  rows <- list()
  for (k in seq_along(profile$points)) {
    point <- profile$points[[k]]
    groups <- c(list(profile = point$attempts), lapply(point$limits$references, `[[`, "attempts"))
    for (kind in names(groups)) for (j in seq_along(groups[[kind]])) {
      a <- groups[[kind]][[j]]
      if (isTRUE(a$audited_value) && is.finite(a$objective))
        rows[[length(rows) + 1L]] <- data.frame(point = k, sd = point$sd,
          kind = kind, attempt = j, objective = a$objective)
    }
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(point = integer(), sd = numeric(), kind = character(),
      attempt = integer(), objective = numeric())
}

#' Compare Private Response-SD Profiles With an Unrestricted Likelihood Fit
#'
#' Fit the same Gaussian model at fixed response-group standard deviations and
#' without that constraint. Retain the likelihood statistic and its numerical
#' qualification separately from confidence calibration.
#'
#' @param spec A private prepared unpenalized Gaussian likelihood specification
#'   accepted by \code{.msr_factorized_response_sd_profile}. It has one response
#'   random intercept and a shared intercept-only kernel scale.
#' @param sd Distinct finite nonnegative response-group SDs, including exact zero,
#'   in response linear-predictor units. Their order is retained.
#' @param control List with \code{starts} (default three) and \code{iter.max}
#'   (default 250), used for both fits. All numerical gates remain fixed.
#' @param initial_parameters Additional original parameter-block lists, as in the
#'   response-SD profile. Selected profile solutions supplement these starts in
#'   the unrestricted fit. Its original face mapping and deduplication are kept.
#' @return A private \code{msr_response_sd_statistic} list. \code{profile} and
#'   \code{unrestricted} retain both complete original fit records, including
#'   diagnostics, attempts and limits. \code{audits} holds fresh physical-vector
#'   checks of the full fit and selected profile solutions. \code{warm_starts}
#'   retains supplied/profile starts and their mapping to original full-fit
#'   attempts. \code{reference} records the unrestricted ordinary negative log
#'   likelihood, qualification and all audited feasible \code{witnesses}.
#'   \code{curve} keeps SD, raw objectives, \code{raw_statistic}, \code{available},
#'   qualified \code{statistic} and a reason for every requested point. Negative
#'   raw discrepancies remain negative; unavailable statistics are \code{NA}.
#'   \code{target}, \code{parameter} and \code{units} identify the comparison.
#'   \code{global_fit_eligible} and \code{uncertainty$available} remain false.
#' @details The normalized full marginal likelihood contains no population priors
#'   or REML adjustment. Both fits use the same specification. If P is the
#'   fixed-SD negative log likelihood and M is the qualified unrestricted candidate
#'   value, the reported diagnostic is \code{2 * (P - M)}. A lower audited
#'   feasible profile/limiting witness makes the reference unresolved. Any negative
#'   discrepancy is retained without clamping, even within numerical tolerance.
#'
#'   A valid nuisance block does not require a zero fixed-SD score. Positive-SD
#'   scores retain log-SD coordinates, and exact zero removes that coordinate.
#'   Neither full Hessian entries nor these statistics supply profile curvature,
#'   a confidence cutoff, p-value, calibrated interval or global-optimum proof.
#'   Rejected fits and limiting witnesses cannot be promoted by this comparison.
#' @seealso \code{\link{.msr_factorized_response_sd_profile}}
#' @keywords internal
.msr_factorized_response_sd_statistic <- function(spec, sd, control = list(),
    initial_parameters = list()) {
  profile <- .msr_factorized_response_sd_profile(spec, sd, control, initial_parameters)
  warm <- lapply(profile$points, function(point) point$parameters)
  warm <- warm[!vapply(warm, is.null, TRUE)]
  completed <- list()
  fit <- tryCatch(.msr_factorized_fit(spec, control,
    initial_parameters = c(initial_parameters, warm),
    checkpoint = function(attempts) completed <<- attempts), error = function(e)
      list(attempts = completed, diagnostics = list(passed = FALSE,
        reason = paste("Unrestricted fit failed:", conditionMessage(e))), error = conditionMessage(e)))
  full <- if (is.null(fit$parameters)) list(passed = FALSE,
    ordinary_objective = NA_real_, reason = "No unrestricted candidate.") else
      .msr_response_sd_statistic_audit(spec, fit$parameters,
        isTRUE(fit$zero_sd[["response"]]), -fit$log_likelihood)
  audits <- lapply(profile$points, function(point) {
    if (is.null(point$parameters)) return(list(passed = FALSE,
      ordinary_objective = NA_real_, reason = "No selected profile candidate."))
    .msr_response_sd_statistic_audit(spec, point$parameters, point$sd == 0,
      point$objective, fixed = TRUE)
  })
  witnesses <- .msr_response_sd_statistic_witnesses(profile)
  for (k in seq_along(audits)) if (isTRUE(audits[[k]]$passed) &&
      is.finite(audits[[k]]$ordinary_objective)) witnesses <- rbind(witnesses,
        data.frame(point = k, sd = profile$points[[k]]$sd, kind = "fresh_profile",
          attempt = if (is.null(profile$points[[k]]$selected)) NA_integer_ else profile$points[[k]]$selected,
          objective = audits[[k]]$ordinary_objective))
  M <- full$ordinary_objective
  witnesses$lower_than_reference <- is.finite(M) & witnesses$objective < M
  lower <- any(witnesses$lower_than_reference)
  qualified <- isTRUE(fit$diagnostics$passed) && isTRUE(full$passed) && !lower
  reference_reason <- if (!isTRUE(fit$diagnostics$passed)) fit$diagnostics$reason else
    if (!isTRUE(full$passed)) "Fresh unrestricted physical-vector audit failed." else
      if (lower) "A lower audited feasible profile or limiting objective was retained." else
        "Original full-fit gates and fresh audits passed; global optimality is unproved."
  curve <- data.frame(sd = profile$curve$sd, profile_objective = NA_real_,
    unrestricted_objective = M, raw_statistic = NA_real_, available = FALSE,
    statistic = NA_real_, reason = "")
  for (k in seq_along(audits)) {
    P <- audits[[k]]$ordinary_objective
    curve$profile_objective[k] <- P
    raw <- 2 * (P - M)
    curve$raw_statistic[k] <- raw
    curve$available[k] <- qualified && isTRUE(profile$points[[k]]$available) &&
      isTRUE(audits[[k]]$passed) && is.finite(raw) && raw >= 0
    curve$reason[k] <- if (is.finite(raw) && raw < 0)
      "Negative profile-minus-unrestricted discrepancy retained; statistic unavailable." else
        if (!qualified) paste(reference_reason, collapse = " ") else
          if (!isTRUE(profile$points[[k]]$available)) paste(profile$points[[k]]$reason, collapse = " ") else
            if (!isTRUE(audits[[k]]$passed)) "Fresh profile physical-vector audit failed." else
              "Qualified numerical-candidate statistic; confidence calibration is unavailable."
    if (curve$available[k]) curve$statistic[k] <- raw
  }
  mapping <- lapply(warm, function(p) lapply(c(FALSE, TRUE), function(zero) {
    p$scale_z <- p$response_z <- numeric()
    if (zero) p$log_response_sd <- numeric()
    flat <- unlist(p, use.names = FALSE); names(flat) <- rep(names(p), lengths(p))
    hits <- which(vapply(fit$attempts, function(a)
      isTRUE(a$zero_sd[["response"]]) == zero &&
        identical(names(a$initial), names(flat)) &&
        identical(as.numeric(a$initial), as.numeric(flat)), TRUE))
    list(zero_sd_response = zero, present_after_original_deduplication = length(hits) > 0L,
      attempt_indices = hits)
  }))
  out <- list(target = profile$target, parameter = profile$parameter, units = profile$units,
    profile = profile, unrestricted = fit, audits = list(unrestricted = full, points = audits),
    warm_starts = list(supplied = initial_parameters, from_profile = warm, mapping = mapping),
    reference = list(qualified = qualified, ordinary_objective = M,
      lower_feasible_profile_witness = lower, witnesses = witnesses, reason = reference_reason),
    curve = curve, global_fit_eligible = FALSE, uncertainty = list(available = FALSE,
      reason = "No confidence cutoff, test calibration or set inversion is implemented."))
  class(out) <- "msr_response_sd_statistic"
  out
}

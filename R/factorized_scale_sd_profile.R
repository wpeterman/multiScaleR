# Exact broad/nearest kernels no longer depend on any annual scale. Their
# normalized latent prior integrates to one, so the SD coordinate is omitted.
.msr_scale_profile_limits <- function(spec, sd, parameters, candidate_nll) {
  references <- list(); witnesses <- character()
  for (kernel in c("broad", "nearest")) {
    limit <- .msr_factorized_limit_spec(spec, kernel)
    p <- parameters; p$scale_beta <- p$scale_z <- p$log_scale_sd <- numeric()
    attempts <- lapply(seq_len(3L), function(i) {
      initial <- p
      if (i > 1L) {
        sizes <- lengths(initial); flat <- unlist(initial, use.names = FALSE)
        flat <- flat + .05 * sin(seq_along(flat) + i)
        at <- 0L
        for (j in seq_along(initial)) {
          n <- sizes[j]
          initial[[j]] <- if (n) flat[at + seq_len(n)] else numeric()
          at <- at + n
        }
      }
      .msr_fixed_sd_profile_attempt(limit, sd, initial, 250L, "scale", allow_irrelevant = TRUE)
    })
    usable <- which(vapply(attempts, function(a) isTRUE(a$audited_value), TRUE))
    best <- if (length(usable)) min(vapply(attempts[usable], `[[`, 0, "objective")) else Inf
    references[[kernel]] <- list(completed = length(usable) > 0L, objective = best, attempts = attempts,
      sd_coordinate = "structurally_omitted_kernel_limit")
    if (is.finite(best) && best <= candidate_nll + 1e-5) witnesses <- c(witnesses, kernel)
  }
  list(completed = all(vapply(references, function(x) isTRUE(x$completed), TRUE)),
    witnesses = witnesses, references = references, tolerance = 1e-5, fixed_sd = sd,
    interpretation = "Kernel limits are independent of annual SD because the normalized latent prior integrates to one; feasible exclusions, never accepted finite fits.")
}

#' Trace a Private Gaussian Annual Log-Scale SD Profile
#'
#' Fix the SD of one independent kernel-scale group, including exact zero, and
#' audit the remaining likelihood parameters. This private trace supplies no CI.
#'
#' @param spec An unpenalized specification from \code{.msr_marginal_prepare}
#'   with Gaussian response, one kernel, an intercept-only fixed scale and one
#'   independent scale random intercept. Response random effects, population
#'   priors and Gaussian projection blocks are unsupported.
#' @param sd Distinct finite nonnegative SDs, including zero, in log kernel-scale
#'   deviations. Requested order is retained. These are not distance radii.
#' @param control List with \code{starts} (default three) and \code{iter.max}
#'   (default 250). Numerical audit and integration thresholds remain fixed.
#' @param initial_parameters Optional additional parameter lists retaining the
#'   original block names. All starts remain explicit, including duplicates and
#'   failures; integrated latent coordinates are omitted and the requested SD fixed.
#' @return A private \code{msr_scale_sd_profile} list. \code{curve} records each
#'   requested SD, exact-zero flag, availability, objective/log likelihood and
#'   objective relative only to the best requested valid point. \code{points}
#'   retain all attempts, full/nuisance gradients and curvature, selected parameters,
#'   failures and required limiting audits. \code{control}, target, parameter,
#'   units, group/covariate and observation/group counts describe the calculation.
#'   \code{diagnostics} reports trace completeness. \code{global_fit_eligible}
#'   and \code{uncertainty$available} remain false. No interval is returned.
#' @details Positive annual SD stays fixed in the normalized full marginal
#'   likelihood while independent latent scale effects are integrated. Exact zero
#'   removes their block, without a floor. Cold integration uses a fresh, stricter
#'   partition and ordinary-value audit. The existing alternate integration
#'   coordinate is also checked for SDs of eight or more. Population priors and
#'   a REML adjustment are absent. Nuisance stationarity/curvature does not imply
#'   a global fit; the full fixed-direction score is retained.
#'
#'   Broad/nearest kernel limits remove annual scale dependence exactly. Their
#'   normalized latent density integrates to one at any requested SD, so its
#'   coordinate is recorded as structurally omitted rather than reoptimized or
#'   given artificial curvature. Both limits retain three starts and cold/ordinary
#'   audits. An audited feasible limit can exclude a finite kernel even if its
#'   curvature is singular; it is never accepted as a finite fit.
#'
#'   Groups are independent and exchangeable, without temporal correlation.
#'   Annual log-scale SD is distinct from response-group SD and spatial-field
#'   amplitude. Every requested point, lower ineligible candidate and competing
#'   nuisance optimum remains explicit. A later boundary-aware bootstrap or
#'   test-inversion design is required before reporting confidence intervals.
#' @keywords internal
.msr_factorized_scale_sd_profile <- function(spec, sd, control = list(), initial_parameters = list()) {
  if (!is.list(spec) || !is.null(spec$priors) || is.null(spec$zero_sd))
    stop("A prepared unpenalized marginal-likelihood specification is required.", call. = FALSE)
  .msr_factorized_validate(spec)
  sm <- spec$scales[[1L]]
  if (spec$family != "gaussian" || length(spec$response$random) ||
      length(sm$random) != 1L || ncol(sm$X) != 1L || any(sm$X[, 1L] != 1) ||
      .msr_has_gaussian_blocks(spec))
    stop("Annual-SD profiles require Gaussian response, no response group and one intercept-only kernel with one independent scale group.", call. = FALSE)
  if (!is.numeric(sd) || is.complex(sd) || !is.null(dim(sd)) || !length(sd) ||
      any(!is.finite(sd) | sd < 0) || anyDuplicated(sd) || !any(sd == 0))
    stop("sd must contain distinct finite nonnegative values including exact zero.", call. = FALSE)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% c("starts", "iter.max")))))
    stop("Profile control supports starts and iter.max only.", call. = FALSE)
  ctl <- utils::modifyList(list(starts = 3L, iter.max = 250L), control)
  for (nm in names(ctl)) {
    if (!is.null(dim(ctl[[nm]]))) stop("Profile controls must be numeric scalars.", call. = FALSE)
    validate_scalar_numeric(ctl[[nm]], nm, integerish = TRUE, lower = 1)
  }
  if (!is.list(initial_parameters)) stop("initial_parameters must be a list.", call. = FALSE)
  for (p in initial_parameters) {
    if (!is.list(p) || !identical(names(p), names(spec$parameters)) ||
        any(!vapply(p, function(x) is.numeric(x) && !is.complex(x) && is.null(dim(x)) && all(is.finite(x)), TRUE)))
      stop("Each supplied start must retain finite numeric original parameter blocks.", call. = FALSE)
    population <- setdiff(names(p), c("scale_z", "response_z", "log_scale_sd"))
    if (any(lengths(p[population]) != lengths(spec$parameters[population])) ||
        !length(p$log_scale_sd) %in% 0:1 ||
        !length(p$scale_z) %in% c(0L, length(spec$parameters$scale_z)) || length(p$response_z))
      stop("Supplied starts must preserve original population and latent block dimensions.", call. = FALSE)
  }
  reason <- spec$identification
  index <- sm$random[[1L]]$index
  out <- list(target = "unpenalized full marginal likelihood", parameter = "annual_log_scale_sd",
    units = "log kernel scale deviations", scale_group = sm$random[[1L]]$group,
    covariate = spec$covariates, observations = length(spec$y),
    scale_groups = length(unique(index)), control = ctl, points = list(), global_fit_eligible = FALSE,
    uncertainty = list(available = FALSE, reason = "Profile trace only; no interval cutoff or boundary calibration."))
  class(out) <- "msr_scale_sd_profile"
  curve <- data.frame(sd = unname(sd), zero_sd_scale = sd == 0, available = FALSE,
    objective = NA_real_, log_likelihood = NA_real_, relative_to_best_requested = NA_real_, reason = "")
  if (length(reason)) {
    curve$reason <- paste(reason, collapse = " "); out$curve <- curve
    out$points <- lapply(sd, function(value) list(sd = value, available = FALSE, attempts = list(), reason = reason))
    out$diagnostics <- list(trace_complete = FALSE, reason = reason)
    return(out)
  }
  for (k in seq_along(sd)) {
    value <- sd[k]; p <- spec$parameters
    p$scale_z <- p$response_z <- numeric()
    p$log_scale_sd <- if (value == 0) numeric() else log(value)
    initials <- lapply(seq_len(ctl$starts), function(i) {
      sizes <- lengths(p); flat <- unlist(p, use.names = FALSE); names(flat) <- rep(names(p), sizes)
      if (i > 1L) {
        flat <- flat + .05 * sin(seq_along(flat) + i)
        flat[names(flat) == "scale_beta"] <- flat[names(flat) == "scale_beta"] + if (i %% 2L == 0L) -.6 else .6
      }
      answer <- p; at <- 0L
      for (j in seq_along(answer)) {
        n <- sizes[j]; answer[[j]] <- if (n) unname(flat[at + seq_len(n)]) else numeric(); at <- at + n
      }
      answer$log_scale_sd <- if (value == 0) numeric() else log(value)
      answer
    })
    initials <- c(initials, initial_parameters) # Keep duplicates and failed starts explicitly.
    attempts <- lapply(initials, function(initial) .msr_fixed_sd_profile_attempt(spec, value, initial, ctl$iter.max, "scale"))
    eligible <- which(vapply(attempts, function(a) isTRUE(a$eligible_nuisance), TRUE))
    point <- list(sd = value, zero_sd_scale = value == 0, available = FALSE,
      attempts = attempts, selected = NA_integer_, global_fit_eligible = FALSE)
    if (!length(eligible)) point$reason <- "No nuisance start passed the fixed numerical gates." else {
      best <- eligible[which.min(vapply(attempts[eligible], `[[`, 0, "objective"))]
      selected <- attempts[[best]]
      better <- any(vapply(attempts, function(a) is.finite(a$objective) && a$objective < selected$objective - 1e-5, TRUE))
      competing <- any(vapply(attempts[eligible], function(a) {
        d <- as.numeric(a$par - selected$par)
        a$objective - selected$objective <= 2 && sqrt(max(0, drop(crossprod(d, selected$nuisance_hessian %*% d)))) > 1
      }, TRUE))
      limits <- tryCatch(.msr_scale_profile_limits(spec, value, selected$parameters, selected$objective),
        error = function(e) list(completed = FALSE, witnesses = character(), error = conditionMessage(e)))
      point$selected <- best; point$parameters <- selected$parameters
      point$objective <- selected$objective; point$limits <- limits
      point$better_ineligible_candidate <- better; point$competing_nuisance_optima <- competing
      point$available <- !better && !competing && isTRUE(limits$completed) && !length(limits$witnesses)
      point$reason <- if (better) "A lower ineligible nuisance candidate was retained." else if (competing)
        "Competing nuisance optima remain unresolved." else if (!isTRUE(limits$completed))
          "Required SD-independent kernel-limit audit failed." else if (length(limits$witnesses))
            "SD-independent kernel-limit witnesses reject finite-kernel identification." else "Conditional nuisance and limiting audits passed; not a global fit or interval."
    }
    out$points[[k]] <- point; curve$available[k] <- point$available; curve$reason[k] <- point$reason
    if (isTRUE(point$available)) {
      curve$objective[k] <- point$objective; curve$log_likelihood[k] <- -point$objective
    }
    gc(FALSE)
  }
  if (any(curve$available)) curve$relative_to_best_requested[curve$available] <-
    curve$objective[curve$available] - min(curve$objective[curve$available])
  out$curve <- curve
  out$diagnostics <- list(trace_complete = all(curve$available), attempted_points = length(sd),
    available_points = sum(curve$available), interval_available = FALSE,
    reason = "Constrained profile trace; existing global-fit gates and interval-calibration requirements are unchanged.")
  out
}

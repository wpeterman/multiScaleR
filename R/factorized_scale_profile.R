# Preserve the original response-profile interfaces and output while sharing
# fixed-SD optimization/audits with the annual log-scale profile.
.msr_response_profile_problem <- function(spec, sd, parameters) {
  .msr_fixed_sd_profile_problem(spec, sd, parameters, "response")
}

.msr_response_profile_attempt <- function(spec, sd, initial, iter.max) {
  .msr_fixed_sd_profile_attempt(spec, sd, initial, iter.max, "response")
}

# Kernel-limit witnesses retain the same response SD. The unrestricted fitter's
# witness optimizer cannot be reused because it would reoptimize this coordinate.
.msr_response_profile_limits <- function(spec, sd, parameters, candidate_nll) {
  references <- list(); witnesses <- character()
  for (kernel in c("broad", "nearest")) {
    reduced <- if (sd == 0) .msr_marginal_reduce(spec, "response") else spec
    limit <- .msr_factorized_limit_spec(reduced, kernel)
    p <- parameters; p$scale_beta <- numeric()
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
      .msr_response_profile_attempt(limit, sd, initial, 250L)
    })
    usable <- which(vapply(attempts, function(a) isTRUE(a$audited_value), TRUE))
    best <- if (length(usable)) min(vapply(attempts[usable], `[[`, 0, "objective")) else Inf
    references[[kernel]] <- list(completed = length(usable) > 0L, objective = best, attempts = attempts)
    if (is.finite(best) && best <= candidate_nll + 1e-5) witnesses <- c(witnesses, kernel)
  }
  list(completed = all(vapply(references, function(x) isTRUE(x$completed), TRUE)),
    witnesses = witnesses, references = references, tolerance = 1e-5,
    fixed_sd = sd, interpretation = "Same-SD feasible limiting witnesses, not accepted finite fits.")
}

#' Trace a Private Gaussian Response-Group SD Profile
#'
#' Hold a response random-intercept SD fixed, including exact zero, and audit
#' nuisance likelihood optima. This diagnostic trace supplies no interval cutoff.
#'
#' @param spec A private unpenalized likelihood specification from
#'   \code{.msr_marginal_prepare}: Gaussian response, one shared intercept-only
#'   kernel scale and one indexed response random intercept. Population priors,
#'   scale random effects and Gaussian projection blocks are unsupported.
#' @param sd Distinct finite nonnegative response-group SDs, including zero.
#'   Values use response linear-predictor units and retain their supplied order.
#' @param control List with \code{starts} (default three) and \code{iter.max}
#'   (default 250). Numerical audit thresholds are fixed.
#' @param initial_parameters Optional list of additional population-parameter
#'   lists with the original block names. Every supplied start is retained;
#'   integrated latent coordinates are omitted and the requested SD is imposed.
#' @return A private \code{msr_response_sd_profile} list with a numeric
#'   \code{curve}, all \code{points}/attempts/limiting witnesses, controls and
#'   diagnostics. Rejected points remain explicit with unavailable curve values.
#'   Relative objectives use the best requested valid point, not a global MLE.
#'   \code{uncertainty$available} and \code{global_fit_eligible} remain false.
#' @details The target is normalized full marginal likelihood, without population
#'   priors or a REML adjustment. Positive SDs remain fixed during nuisance and
#'   limiting-model optimization. Zero removes the response latent/SD block
#'   exactly. No boundary floor or Hessian repair is applied.
#'
#'   Nuisance stationarity and positive nuisance curvature do not imply a global
#'   optimum: full gradients and Hessians remain in each attempt. Every original
#'   start, cold/ordinary audit, lower ineligible candidate, competing nuisance
#'   optimum and same-SD limiting witness is retained. Singular limiting models
#'   can be feasible exclusion witnesses but are never accepted finite fits.
#'
#'   This first trace concerns response-group SD, not annual kernel-scale SD or
#'   spatial-field amplitude. It provides no confidence level, Wilks cutoff,
#'   p-value, predictive interval or calibration claim. A later boundary-aware
#'   bootstrap or test-inversion design is required before reporting intervals.
#' @keywords internal
.msr_factorized_response_sd_profile <- function(spec, sd, control = list(), initial_parameters = list()) {
  if (!is.list(spec) || !is.null(spec$priors) || is.null(spec$zero_sd))
    stop("A prepared unpenalized marginal-likelihood specification is required.", call. = FALSE)
  .msr_factorized_validate(spec)
  sm <- spec$scales[[1L]]
  if (spec$family != "gaussian" || length(spec$response$random) != 1L ||
      length(sm$random) || ncol(sm$X) != 1L || any(sm$X[, 1L] != 1) ||
      .msr_has_gaussian_blocks(spec))
    stop("Response-SD profiles require Gaussian response, one indexed response intercept and one shared intercept-only kernel scale.", call. = FALSE)
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
    population <- setdiff(names(p), c("scale_z", "response_z", "log_response_sd"))
    if (any(lengths(p[population]) != lengths(spec$parameters[population])) ||
        !length(p$log_response_sd) %in% 0:1 ||
        !length(p$response_z) %in% c(0L, length(spec$parameters$response_z)) || length(p$scale_z))
      stop("Supplied starts must preserve original population and latent block dimensions.", call. = FALSE)
  }
  reason <- spec$identification
  index <- spec$response$random[[1L]]$index
  if (max(tabulate(index)) <= 1L) reason <- c(reason, "Singleton response groups cannot identify separate residual and group variances.")
  if (max(tabulate(index)) > 1L && all(vapply(split(spec$y - spec$offset, index), function(y) length(unique(y)) == 1L, TRUE)))
    reason <- c(reason, "Gaussian likelihood is unbounded for constant responses within every response group.")
  out <- list(target = "unpenalized full marginal likelihood", parameter = "response_group_sd",
    units = "response linear predictor", response_group = spec$response$random[[1L]]$group,
    covariate = spec$covariates, observations = length(spec$y),
    response_groups = length(unique(index)), control = ctl, points = list(), global_fit_eligible = FALSE,
    uncertainty = list(available = FALSE, reason = "Profile trace only; no interval cutoff or boundary calibration."))
  class(out) <- "msr_response_sd_profile"
  curve <- data.frame(sd = unname(sd), zero_sd_response = sd == 0, available = FALSE,
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
    p$log_response_sd <- if (value == 0) numeric() else log(value)
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
      answer$log_response_sd <- if (value == 0) numeric() else log(value)
      answer
    })
    initials <- c(initials, initial_parameters) # Keep duplicates and failed starts explicitly.
    attempts <- lapply(initials, function(initial) .msr_response_profile_attempt(spec, value, initial, ctl$iter.max))
    eligible <- which(vapply(attempts, function(a) isTRUE(a$eligible_nuisance), TRUE))
    point <- list(sd = value, zero_sd_response = value == 0, available = FALSE,
      attempts = attempts, selected = NA_integer_, global_fit_eligible = FALSE)
    if (!length(eligible)) point$reason <- "No nuisance start passed the fixed numerical gates." else {
      best <- eligible[which.min(vapply(attempts[eligible], `[[`, 0, "objective"))]
      selected <- attempts[[best]]
      better <- any(vapply(attempts, function(a) is.finite(a$objective) && a$objective < selected$objective - 1e-5, TRUE))
      competing <- any(vapply(attempts[eligible], function(a) {
        d <- as.numeric(a$par - selected$par)
        a$objective - selected$objective <= 2 && sqrt(max(0, drop(crossprod(d, selected$nuisance_hessian %*% d)))) > 1
      }, TRUE))
      limits <- tryCatch(.msr_response_profile_limits(spec, value, selected$parameters, selected$objective),
        error = function(e) list(completed = FALSE, witnesses = character(), error = conditionMessage(e)))
      point$selected <- best; point$parameters <- selected$parameters
      point$objective <- selected$objective; point$limits <- limits
      point$better_ineligible_candidate <- better; point$competing_nuisance_optima <- competing
      point$available <- !better && !competing && isTRUE(limits$completed) && !length(limits$witnesses)
      point$reason <- if (better) "A lower ineligible nuisance candidate was retained." else if (competing)
        "Competing nuisance optima remain unresolved." else if (!isTRUE(limits$completed))
          "Required same-SD limiting audit failed." else if (length(limits$witnesses))
            "Same-SD limiting witnesses reject finite-kernel identification." else "Conditional nuisance and limiting audits passed; not a global fit or interval."
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

# Fixed-kernel limits give feasible likelihood witnesses. They do not provide
# finite annual-scale estimates, confidence intervals, or a global certificate.
.msr_factorized_limit_spec <- function(spec, kernel, poisson = FALSE) {
  .msr_factorized_validate(spec)
  .msr_factorized_kernel_limit(spec, kernel, poisson)
}

# Pure fixed-kernel reduction. This does not admit a finite SQ fit: the wrapper
# above retains structural validation, and all numerical eligibility gates stay.
.msr_factorized_kernel_limit <- function(spec, kernel, poisson = FALSE) {
  if (length(spec$covariates) != 1L || length(spec$scales) != 1L)
    stop("Fixed-kernel reduction requires one predictor.", call. = FALSE)
  if (!kernel %in% c("broad", "nearest")) stop("Unknown kernel limit.", call. = FALSE)
  if (poisson && spec$family != "nbinom2") stop("Poisson limit requires NB2.", call. = FALSE)
  # Shared or independently fixed scales have no annual SD block to remove.
  # Their fixed-kernel witnesses still need exactly the same likelihood audit.
  out <- if (length(spec$scales[[1L]]$random))
    .msr_marginal_reduce(spec, "scale") else spec
  nm <- spec$covariates[[1L]]
  if (nm %in% names(spec$surface$buffers)) {
    raw <- .msr_scale_surface_limit(spec$surface$buffers[[nm]], kernel)
    if (length(raw) != length(spec$y))
      stop("Exact SQ limit rows do not match the response.", call. = FALSE)
    # SQ-only specifications have no binned V/C names. Name every new entry.
    out$V[[nm]] <- matrix(raw, ncol = 1L)
    out$C[[nm]] <- matrix(1, length(raw), 1L)
    out$delta[[nm]] <- matrix(0, length(raw), 1L)
    out$surface$buffers[[nm]] <- NULL
    out$surface$memory <- .msr_scale_surface_size(out$surface)
    out$fixed_kernel_limit <- list(kernel = kernel, covariate = nm, metric = "sq",
      source_specification = spec$surface$specs,
      normalization = "original_preparation_reference",
      finite_scale_estimate = FALSE)
  } else {
    keep <- if (kernel == "nearest") spec$delta[[1L]] == 0 else
      matrix(TRUE, nrow(spec$V[[1L]]), ncol(spec$V[[1L]]))
    out$V[[1L]] <- matrix(rowSums(spec$V[[1L]] * keep), ncol = 1L)
    out$C[[1L]] <- matrix(rowSums(spec$C[[1L]] * keep), ncol = 1L)
    out$delta[[1L]] <- matrix(0, length(spec$y), 1L)
  }
  out$scales[[1L]]$X <- matrix(numeric(), length(spec$y), 0L)
  out$scales[[1L]]$beta_index <- integer()
  out$parameters$scale_beta <- numeric()
  if (poisson) {
    out$family <- "poisson"
    out$parameters$log_dispersion <- numeric()
  }
  out
}

.msr_factorized_limit_audit <- function(spec, parameters, candidate_nll) {
  .msr_factorized_limit_audit_engine(spec, parameters, candidate_nll,
    .msr_factorized_limit_spec, .msr_factorized_path_limits,
    .msr_factorized_objective, .msr_factorized_scale_objective)
}

# Share the original audit algorithm without changing strict fitting admission.
# Factories supply exact reference specifications and normalized point values;
# all starts, failed attempts, thresholds and curvature diagnostics stay here.
.msr_factorized_limit_audit_engine <- function(spec, parameters, candidate_nll,
    kernel_factory, path_factory, objective_factory, scale_objective_factory) {
  if (length(spec$response$random) && !length(parameters$log_response_sd))
    spec <- .msr_marginal_reduce(spec, "response")
  references <- list(); witnesses <- character()
  for (kernel in c("broad", "nearest")) for (poisson in
    if (spec$family == "nbinom2") c(FALSE, TRUE) else FALSE) {
    label <- paste(kernel, if (poisson) "poisson" else spec$family, sep = "-")
    s <- kernel_factory(spec, kernel, poisson)
    p <- parameters
    p$scale_beta <- p$scale_z <- p$log_scale_sd <- numeric()
    if (poisson) p$log_dispersion <- numeric()
    attempts <- list()
    obj <- tryCatch({
      # Explicit component factories retain partial endpoint failures as
      # conditions. The strict default factory still raises before this point.
      if (inherits(s, "error")) stop(s)
      .msr_factorized_adfun(s, p)
    }, error = function(e)e)
    if (inherits(obj, "error")) {
      references[[label]] <- list(completed = FALSE, error = conditionMessage(obj))
      next
    }
    for (i in 1:3) {
      initial <- obj$par
      if (i > 1L) initial <- initial + .05 * sin(seq_along(initial) + i)
      a <- tryCatch({
        fit <- stats::nlminb(initial, obj$fn, obj$gr,
          control = list(iter.max = 250L, eval.max = 750L, rel.tol = 1e-12))
        refined <- .msr_marginal_polish(fit$par, obj$fn, obj$gr, obj$he, maxit = 6L)
        if (isTRUE(refined$passed) && refined$objective <= fit$objective + 1e-10)
          fit$par <- refined$par
        value <- as.numeric(obj$fn(fit$par)); gradient <- as.numeric(obj$gr(fit$par))
        q <- obj$env$parList(x = fit$par, par = obj$env$last.par)
        ordinary <- objective_factory(s, FALSE)(q)
        audited <- all(is.finite(c(value, ordinary, gradient))) && abs(value - ordinary) <= 1e-6
        list(audited = audited, objective = ordinary, ad_objective = value,
          audit_error = abs(value - ordinary), parameters = q, initial = initial,
          gradient_max = max(abs(gradient)), optimizer = fit, refinement = refined,
          curvature = .msr_marginal_curvature(obj$he(fit$par)))
      }, error = function(e) list(audited = FALSE, initial = initial, error = conditionMessage(e)))
      attempts[[i]] <- a
    }
    checked <- Filter(function(a)isTRUE(a$audited), attempts)
    best <- if (length(checked)) min(vapply(checked, `[[`, 0, "objective")) else Inf
    references[[label]] <- list(completed = length(checked) > 0L, objective = best, attempts = attempts)
    if (is.finite(best) && best <= candidate_nll + 1e-5) witnesses <- c(witnesses, label)
  }
  if (spec$family == "nbinom2") {
    s <- if (length(spec$scales[[1L]]$random) && !length(parameters$log_scale_sd))
      .msr_marginal_reduce(spec, "scale") else spec
    s$family <- "poisson"; q <- parameters; q$log_dispersion <- numeric()
    a <- tryCatch({
      value <- objective_factory(s, FALSE, 1e-10, .5)(q)
      coordinate <- if (length(q$log_scale_sd) && exp(q$log_scale_sd) >= 8)
        scale_objective_factory(s, FALSE)(q) else value
      list(completed = all(is.finite(c(value, coordinate))) && abs(value - coordinate) <= 1e-6,
        objective = value, coordinate_error = abs(value - coordinate))
    }, error = function(e) list(completed = FALSE, error = conditionMessage(e)))
    references[["fitted-scales-poisson"]] <- a
    if (isTRUE(a$completed) && a$objective <= candidate_nll + 1e-5)
      witnesses <- c(witnesses, "fitted-scales-poisson")
  }
  # The input model, rather than its selected variance face, defines whether
  # an unbounded annual-SD direction exists. These are fixed-nuisance feasible
  # references; all original optimized references and failure records remain.
  path <- tryCatch(path_factory(spec, parameters),
    error = function(e) list(applicable = TRUE, reason = conditionMessage(e),
      references = list(`annual-sd-infinite-error` = list(completed = FALSE,
        error = conditionMessage(e)))))
  references <- c(references, path$references)
  for (label in names(path$references)) {
    a <- path$references[[label]]
    if (isTRUE(a$completed) && is.finite(a$objective) &&
        a$objective <= candidate_nll + 1e-5) witnesses <- c(witnesses, label)
  }
  list(completed = all(vapply(references, function(a)isTRUE(a$completed), TRUE)),
    witnesses = witnesses, references = references, tolerance = 1e-5,
    infinite_sd_path = path[setdiff(names(path), "references")],
    interpretation = "Feasible limiting-model exclusion witnesses; absence does not prove global identification.")
}

# Negative log likelihood for independent groups, each sharing one latent choice
# of the exact broad or nearest kernel. Inputs are GROUP joint-density NLLs.
# This arithmetic does not validate a scale-intercept/iid-group limiting path,
# admit a finite fit, add a prior/Jacobian, or change the existing limiting audit.
.msr_factorized_mixture_nll <- function(broad, nearest, weight, ad = FALSE) {
  if (!is.numeric(weight) || length(weight) != 1L || !is.finite(weight) ||
      weight < 0 || weight > 1)
    stop("Mixture weight must be one finite number between zero and one.", call. = FALSE)
  if (!length(broad) || length(broad) != length(nearest) ||
      !is.null(dim(broad)) || !is.null(dim(nearest)))
    stop("Mixture inputs must be matching nonempty group NLL vectors.", call. = FALSE)
  if (!ad && (!is.numeric(broad) || !is.numeric(nearest) || anyNA(c(broad, nearest)) ||
      any(c(broad, nearest) == -Inf)))
    stop("Group NLLs must be numeric without missing or negative-infinite values.", call. = FALSE)
  if (weight == 0) return(sum(nearest))
  if (weight == 1) return(sum(broad))
  a <- log(weight) - broad
  b <- log1p(-weight) - nearest
  if (ad) return(-sum(RTMB::logspace_add(a, b)))
  # Both components can have zero density. Retain Inf NLL rather than NaN.
  zero <- is.infinite(broad) & is.infinite(nearest)
  value <- pmax(a, b) + log1p(exp(-abs(a - b)))
  value[zero] <- -Inf
  -sum(value)
}

# With one iid annual intercept, tau -> Inf and mu/tau -> ratio select the
# broad kernel when Z > -ratio and the nearest kernel otherwise. The choice
# belongs to the whole group. Fixed Gaussian response covariance is already
# integrated by the response leaf; the annual normal mass supplies Phi(ratio).
# This neutral likelihood component does not certify full-model factorization,
# posterior limits, finite identification, or fitting eligibility.
.msr_factorized_group_mixture_limit <- function(group, parameters, ratio) {
  .msr_factorized_group_check(group)
  if (!is.null(group$priors))
    stop("Group likelihood limits require a neutral specification, without population priors.", call. = FALSE)
  sm <- group$scales[[1L]]
  if (is.null(.msr_factorized_location_direction(sm, group$parameters)) ||
      nrow(sm$X) != length(group$y))
    stop("Group mixture paths require a canonical scale intercept or supported common location direction, and one annual SD.", call. = FALSE)
  sizes <- lengths(group$parameters)
  sizes[c("scale_z", "response_z")] <- 0L
  if (!is.list(parameters) || is.null(names(parameters)) || anyDuplicated(names(parameters)) ||
      !identical(lengths(parameters), sizes) ||
      !all(vapply(parameters, function(x) is.numeric(x) && !is.complex(x), TRUE)) ||
      any(!is.finite(unlist(parameters, use.names = FALSE))))
    stop("Supply a matching finite population layout without realized latent effects.", call. = FALSE)
  if (!is.numeric(ratio) || is.complex(ratio) || length(ratio) != 1L || !is.finite(ratio))
    stop("Supply one finite real intercept-to-SD limiting ratio.", call. = FALSE)
  p <- parameters
  p$scale_beta <- p$log_scale_sd <- p$scale_z <- numeric()
  component <- vapply(c("broad", "nearest"), function(kernel) {
    endpoint <- .msr_factorized_kernel_limit(group, kernel)
    value <- sum(.msr_factorized_response_nll(endpoint, p))
    if (!is.finite(value))
      stop("Finite normalized endpoint likelihoods are required for a group limiting witness.", call. = FALSE)
    value
  }, 0)
  weight <- stats::pnorm(ratio)
  list(nll = .msr_factorized_mixture_nll(component[[1L]], component[[2L]], weight),
    component_nll = component, broad_weight = weight, ratio = ratio,
    target = "marginal_likelihood_limit", finite_scale_estimate = FALSE,
    interpretation = "Fixed-ratio large-SD group likelihood limit; not a finite fit, posterior limit or identification certificate.")
}

# The unrestricted prepared model retains positive-SD directions even when its
# maximizing candidate is on the exact zero face. A structurally zero-fixed
# model has no such direction. The finite seed below is explicitly a path
# origin that the group limit removes, never a fitted/restored SD or epsilon.
.msr_factorized_path_limits <- function(spec, parameters) {
  .msr_factorized_validate(spec)
  .msr_factorized_path_component(spec, parameters)
}

# Component-only operation: exact normalized limiting densities on an intact
# neutral model. Existing fitting and constructor wrappers retain admission.
.msr_factorized_path_component <- function(spec, parameters) {
  .msr_factorized_limit_component_validate(spec)
  if (!is.null(spec$priors))
    stop("Annual-SD likelihood paths require a neutral specification.", call. = FALSE)
  sm <- spec$scales[[1L]]
  if (!length(sm$random))
    return(list(applicable = FALSE, reason = "Prepared model has no annual-SD direction.", references = list()))
  if (is.null(.msr_factorized_location_direction(sm, spec$parameters)) ||
      nrow(sm$X) != length(spec$y))
    return(list(applicable = FALSE, reason = "Fixed-ratio paths require a supported common scale-location direction and annual SD.",
      references = list()))
  if (!is.list(parameters) || !is.numeric(parameters$log_scale_sd) ||
      is.complex(parameters$log_scale_sd) || length(parameters$log_scale_sd) > 1L)
    stop("Supply the candidate population layout and an explicit annual-SD coordinate or zero face.", call. = FALSE)
  p <- parameters
  candidate_zero <- !length(p$log_scale_sd)
  if (candidate_zero) p$log_scale_sd <- spec$parameters$log_scale_sd
  references <- list()
  ratios <- c(-.75, 0, .75)
  for (poisson in if (spec$family == "nbinom2") c(FALSE, TRUE) else FALSE) {
    s <- spec; q <- p
    if (poisson) {
      s$family <- "poisson"
      s$parameters$log_dispersion <- q$log_dispersion <- numeric()
    }
    for (j in seq_along(ratios)) {
      label <- paste("annual-sd-infinite", s$family, c("negative", "zero", "positive")[j], sep = "-")
      references[[label]] <- tryCatch({
        values <- vapply(seq_along(sm$random[[1L]]$levels), function(i)
          .msr_factorized_group_mixture_limit(.msr_factorized_group(s, i), q, ratios[j])$nll, 0)
        value <- sum(values)
        if (!length(values) || !is.finite(value)) stop("Finite group path likelihoods are required.", call. = FALSE)
        list(completed = TRUE, objective = value, ratio = ratios[j], broad_weight = stats::pnorm(ratios[j]),
          groups = length(values), target = "marginal_likelihood_limit", finite_scale_estimate = FALSE)
      }, error = function(e) list(completed = FALSE, ratio = ratios[j], error = conditionMessage(e)))
    }
  }
  list(applicable = TRUE, reason = "Feasible original-model annual-SD directions at fixed nuisance parameters.",
    candidate_zero_sd = candidate_zero, direction = "original_unrestricted_model",
    path_origin = if (candidate_zero) "prepared_finite_SD_seed_removed_by_limit" else "candidate_finite_SD_removed_by_limit",
    references = references)
}

# A deterministic feasible direction X*a=1 is enough for a common-ratio path.
# Shared preparation permits an intercept or disjoint factor cell means. Fixed
# contrasts remain finite while the location moves with the annual SD. Group
# subsets keep all original columns even when some are locally unobserved.
.msr_factorized_location_direction <- function(sm, parameters) {
  X <- sm$X
  if (!is.matrix(X) || !is.numeric(X) || is.complex(X) || !nrow(X) || !ncol(X) ||
      any(!is.finite(X)) || length(sm$random) != 1L ||
      !identical(as.integer(sm$random[[1L]]$sd_index), 1L) ||
      !identical(as.integer(sm$beta_index), seq_len(ncol(X))) ||
      length(parameters$scale_beta) != ncol(X) || length(parameters$log_scale_sd) != 1L)
    return(NULL)
  intercept <- which(colSums(X == 1) == nrow(X))
  if (length(intercept)) {
    direction <- numeric(ncol(X)); direction[intercept[1L]] <- 1
  } else if (all(X %in% c(0, 1)) && all(rowSums(X) == 1)) {
    direction <- rep(1, ncol(X))
  } else return(NULL)
  direction
}

# Separate exact-limit structure from model admission. This never constructs a
# tape, optimizer, finite estimate, interval, or a surrogate mean specification.
.msr_factorized_limit_component_validate <- function(spec) {
  if (!is.null(spec$priors))
    stop("Annual-SD likelihood paths require a neutral specification.", call. = FALSE)
  n <- length(spec$y)
  if ((length(spec$rows) && (length(spec$rows) != n || anyNA(spec$rows))) ||
      (length(spec$observation_ids) && (length(spec$observation_ids) != n || anyNA(spec$observation_ids))) ||
      (!is.null(spec$data) && nrow(spec$data) != n))
    stop("Component observation metadata does not match its density payload.", call. = FALSE)
  .msr_scale_check_structure(spec, "marginal_likelihood", "factorized",
    .msr_scale_capability("marginal_likelihood", "factorized"))
}

# Exact endpoint operation, not finite fitting admission. The intact source is
# validated before its mathematical reduction; fixed references are retained.
.msr_factorized_kernel_component <- function(spec, kernel, poisson = FALSE) {
  .msr_factorized_limit_component_validate(spec)
  out <- .msr_factorized_kernel_limit(spec, kernel, poisson)
  out$exact_limit_component <- list(operation = "exact_kernel_limit",
    target = "marginal_likelihood_limit", source_family = spec$family,
    kernel = kernel, poisson = poisson, finite_scale_estimate = FALSE)
  out
}

# Numeric callbacks for the explicit limiting-component audit only. Exact
# Poisson reduction omits dispersion from the prepared template as well as the
# evaluated parameters; the general density component keeps strict layout.
.msr_factorized_audit_density <- function(spec, ad = FALSE, rel.tol = NULL,
    step = NULL, coordinate = c("standard", "physical")) {
  if (!identical(ad, FALSE))
    stop("Audit density components support ordinary numeric evaluation only.", call. = FALSE)
  if (spec$family == "poisson") spec$parameters$log_dispersion <- numeric()
  .msr_factorized_density_component(spec, coordinate, rel.tol, step)
}

#' Exact Limiting-Likelihood Component Audit
#'
#' Compare a supplied candidate negative log likelihood with feasible exact
#' kernel and annual-SD limits. This private operation reuses the original
#' audit algorithm without admitting a finite weighted-SQ fit.
#'
#' @param spec Intact neutral prepared specification with one supported scale
#'   predictor and independent annual groups. Original surface reference
#'   normalization, observation rows, and dependence restrictions are required.
#' @param parameters Named finite population parameter blocks for the candidate,
#'   without realized latent effects. An empty annual or response SD block
#'   denotes its explicitly selected zero face in the original prepared model.
#' @param candidate_nll One finite numeric negative log likelihood at the
#'   candidate. A reference is an exclusion witness when its value is no larger
#'   than this value plus the original audit tolerance, \code{1e-5}.
#' @return A list retaining the original audit fields: \code{completed} records
#'   whether every required reference completed; \code{witnesses} names feasible
#'   limits meeting the exclusion rule; \code{references} retains every endpoint
#'   attempt, initial values, objective/ordinary-value discrepancy, optimizer,
#'   refinement, curvature, and failure, together with the fitted-scale Poisson
#'   and fixed-ratio infinite-SD references when applicable. \code{tolerance} is
#'   the original comparison tolerance; \code{infinite_sd_path} describes the
#'   original unrestricted-model direction; \code{interpretation} states the
#'   identification limitation. Added fields \code{component_only = TRUE},
#'   \code{target = "marginal_likelihood_limit"},
#'   \code{finite_scale_estimate = FALSE}, and
#'   \code{uncertainty_available = FALSE} distinguish this result from a fit.
#' @details Endpoint nuisance parameters retain all three original starts,
#'   polishing, native/automatic-differentiation comparisons, and curvature
#'   diagnostics. A finite feasible value can exclude a candidate even if its
#'   optimizer or Hessian does not qualify as an accepted fitted reference;
#'   those diagnostics remain available. Conversely, completion and absence of
#'   a witness do not establish global identification. This operation supplies
#'   no intervals, posterior inference, new parameter count, or finite scale.
#'   Priors remain excluded. Original fitting and constructor guards remain.
#' @keywords internal
.msr_factorized_limit_component_audit <- function(spec, parameters, candidate_nll) {
  .msr_factorized_limit_component_validate(spec)
  if (!is.numeric(candidate_nll) || is.complex(candidate_nll) ||
      !is.null(dim(candidate_nll)) || length(candidate_nll) != 1L || !is.finite(candidate_nll))
    stop("Supply one finite numeric candidate negative log likelihood.", call. = FALSE)
  kernel_factory <- function(...) tryCatch(.msr_factorized_kernel_component(...),
    error = function(e)e)
  out <- .msr_factorized_limit_audit_engine(spec, parameters, candidate_nll,
    kernel_factory, .msr_factorized_path_component, .msr_factorized_audit_density,
    function(...) .msr_factorized_audit_density(..., coordinate = "physical"))
  c(out, list(component_only = TRUE, target = "marginal_likelihood_limit",
    finite_scale_estimate = FALSE, uncertainty_available = FALSE))
}

# Conditional profile targets share the numerical audit, while the fixed
# block, zero reduction and structurally irrelevant limits remain explicit.
.msr_fixed_sd_profile_problem <- function(spec, sd, parameters, component,
    allow_irrelevant = FALSE, cold = FALSE) {
  component <- match.arg(component, c("response", "scale"))
  block <- paste0("log_", component, "_sd")
  if (allow_irrelevant && !(component == "scale" && !length(spec$scales[[1L]]$random) &&
      ncol(spec$scales[[1L]]$X) == 0L && !length(spec$parameters$scale_beta)))
    stop("Only an exact kernel limit can omit a positive scale SD.", call. = FALSE)
  if (sd == 0) spec <- .msr_marginal_reduce(spec, component)
  parameters$scale_z <- parameters$response_z <- numeric()
  parameters[[block]] <- if (sd == 0 || allow_irrelevant) numeric() else log(sd)
  object <- if (component == "scale" && cold)
    .msr_factorized_adfun(spec, parameters, rel.tol = 1e-10, step = .5) else
      .msr_factorized_optimizer(spec, parameters, step = 1)
  fixed <- which(names(object$par) == block)
  free <- setdiff(seq_along(object$par), fixed)
  if (length(fixed) != as.integer(sd > 0 && !allow_irrelevant))
    stop(paste0("The ", component, "-SD profile requires exactly one fixed SD coordinate."), call. = FALSE)
  full <- function(x) {
    if (length(x) != length(free) || any(!is.finite(x)))
      stop("Invalid nuisance parameter vector.", call. = FALSE)
    p <- object$par
    p[free] <- as.numeric(x)
    if (length(fixed)) p[fixed] <- log(sd)
    p
  }
  list(spec = spec, object = object, full = full, free = free, fixed = fixed,
    par = object$par[free], fn = function(x) object$fn(full(x)),
    gr = function(x) object$gr(full(x))[free],
    he = function(x) object$he(full(x))[free, free, drop = FALSE],
    parameters = function(x) {
      p <- full(x)
      object$fn(p) # Restore accepted state after any rejected trial.
      object$env$parList(x = p, par = object$env$last.par)
    })
}

.msr_fixed_sd_profile_attempt <- function(spec, sd, initial, iter.max, component,
    allow_irrelevant = FALSE) {
  component <- match.arg(component, c("response", "scale"))
  block <- paste0("log_", component, "_sd")
  build <- function(s, p, cold = FALSE) {
    if (component == "response") .msr_response_profile_problem(s, sd, p) else
      .msr_fixed_sd_profile_problem(s, sd, p, "scale", allow_irrelevant, cold)
  }
  warnings <- character(); candidate <- Inf; parameters <- NULL; par <- NULL
  imposed_initial <- initial
  imposed_initial$scale_z <- imposed_initial$response_z <- numeric()
  imposed_initial[[block]] <- if (sd == 0 || allow_irrelevant) numeric() else log(sd)
  stage <- "objective construction"
  answer <- tryCatch(withCallingHandlers({
    problem <- build(spec, imposed_initial)
    stage <- "nuisance optimization"
    fit <- stats::nlminb(problem$par, problem$fn, problem$gr,
      control = list(iter.max = iter.max, eval.max = 3L * iter.max, rel.tol = 1e-11))
    candidate <- fit$objective; par <- fit$par
    parameters <- problem$parameters(par)
    stage <- "Newton certification"
    fit$refinement <- .msr_marginal_polish(par, problem$fn, problem$gr, problem$he, maxit = 6L)
    if (isTRUE(fit$refinement$passed) && fit$refinement$objective <= fit$objective + 1e-10) {
      fit$optimizer_par <- fit$par; fit$optimizer_objective <- fit$objective
      fit$par <- fit$refinement$par; fit$objective <- fit$refinement$objective
    }
    par <- fit$par; value <- as.numeric(problem$fn(par)); candidate <- value
    gradient <- as.numeric(problem$gr(par)); H <- problem$he(par)
    geometry <- .msr_marginal_curvature(H)
    parameters <- problem$parameters(par)
    full_par <- problem$full(par)
    full_gradient <- as.numeric(problem$object$gr(full_par))
    names(full_gradient) <- names(full_par)
    full_hessian <- problem$object$he(full_par)
    stage <- "cold nuisance audit"
    cold <- build(spec, parameters, cold = TRUE)
    cold_value <- as.numeric(cold$fn(par)); cold_gradient <- as.numeric(cold$gr(par))
    ordinary <- .msr_factorized_objective(cold$spec, FALSE, rel.tol = 1e-10, step = .5)(parameters)
    coordinate <- if (component == "scale" && length(parameters$log_scale_sd) &&
        exp(parameters$log_scale_sd) >= 8) .msr_factorized_scale_objective(cold$spec, FALSE)(parameters) else ordinary
    error <- max(abs(c(value - cold_value, value - ordinary)))
    if (component == "scale") error <- max(error, abs(value - coordinate))
    gradient_error <- max(abs(gradient - cold_gradient))
    audited <- all(is.finite(c(value, cold_value, ordinary, gradient, cold_gradient, full_gradient, coordinate))) && error <= 1e-6
    stationary <- audited && (fit$convergence == 0L || isTRUE(fit$refinement$passed)) &&
      max(abs(gradient)) <= 1e-4 && max(abs(cold_gradient)) <= 1e-4 && gradient_error <= 1e-4 &&
      all(is.finite(full_hessian))
    list(eligible_nuisance = stationary && geometry$passed, stationary = stationary,
      audited_value = audited, objective = as.numeric(ordinary), par = par, parameters = parameters,
      optim = fit, nuisance_hessian = H, full_hessian = full_hessian, full_gradient = full_gradient,
      full_curvature = .msr_marginal_curvature(full_hessian),
      parameter_blocks = names(full_par),
      free_indices = problem$free, fixed_indices = problem$fixed, full_par = full_par,
      nuisance_gradient_max = max(abs(gradient)), condition = geometry$condition,
      minimum_eigenvalue = geometry$minimum, cold_error = abs(value - cold_value),
      reference_error = abs(value - ordinary), gradient_error = gradient_error,
      fixed_sd = sd, zero_sd = sd == 0, global_fit_eligible = FALSE)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) list(eligible_nuisance = FALSE, stationary = FALSE,
    audited_value = FALSE, objective = candidate, parameters = parameters, par = par,
    fixed_sd = sd, zero_sd = sd == 0, global_fit_eligible = FALSE,
    error_stage = stage, error = conditionMessage(e)))
  names(answer)[names(answer) == "zero_sd"] <- paste0("zero_sd_", component)
  if (component == "scale") {
    answer$sd_coordinate <- if (allow_irrelevant) "structurally_omitted_kernel_limit" else
      if (sd == 0) "exact_zero_reduction" else "fixed_log_scale_sd"
    if (exists("coordinate", inherits = FALSE)) answer$coordinate_error <- abs(candidate - coordinate)
  }
  answer$initial <- imposed_initial; answer$supplied_initial <- initial
  answer$warnings <- unique(warnings)
  answer
}


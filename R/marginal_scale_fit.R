# Numerical gates are fixed in this prototype; they are not interval-calibration
# claims. No likelihood curvature is regularized to pass these checks.
.msr_marginal_curvature <- function(H) {
  if (!length(H)) return(list(passed = TRUE, condition = 1, minimum = Inf))
  if (any(!is.finite(H))) return(list(passed = FALSE, condition = Inf, minimum = NA_real_))
  values <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  condition <- if (min(values) > 0) max(values) / min(values) else Inf
  list(passed = is.finite(condition) && condition <= 1e12,
       condition = condition, minimum = min(values))
}

# An optimizer's singular-convergence message is not a stationarity certificate.
# Independently refine with an unmodified positive-definite Hessian. Record the
# original optimizer status and require a tighter gradient for this certificate.
.msr_marginal_polish <- function(par, fn, gr, he, maxit = 12L) {
  value <- as.numeric(fn(par))
  gradient <- as.numeric(gr(par))
  iterations <- 0L
  for (iteration in seq_len(maxit)) {
    H <- he(par)
    if (any(!is.finite(c(value, gradient))) || !.msr_marginal_curvature(H)$passed) break
    if (max(abs(gradient)) <= 1e-7) break
    step <- as.numeric(solve(H, gradient))
    accepted <- FALSE
    for (backtrack in 0:24) {
      candidate <- par - step * 2^(-backtrack)
      next_value <- as.numeric(fn(candidate))
      next_gradient <- as.numeric(gr(candidate))
      roundoff <- 10 * .Machine$double.eps * (1 + abs(value))
      if (all(is.finite(c(next_value, next_gradient))) &&
          next_value <= value + roundoff &&
          (next_value <= value - 1e-4 * sum(gradient * step) * 2^(-backtrack) ||
           max(abs(next_gradient)) < max(abs(gradient)))) {
        par <- candidate; value <- next_value; gradient <- next_gradient
        accepted <- TRUE; iterations <- iteration
        break
      }
    }
    if (!accepted) break
  }
  # Re-evaluate at the accepted point, since a rejected trial can change a TMB
  # object's working state. No Hessian repair or boundary floor is applied.
  value <- as.numeric(fn(par)); gradient <- as.numeric(gr(par))
  geometry <- .msr_marginal_curvature(he(par))
  list(par = par, objective = value, gradient = max(abs(gradient)),
    passed = all(is.finite(c(value, gradient))) && max(abs(gradient)) <= 1e-7 && geometry$passed,
    iterations = iterations, condition = geometry$condition)
}

.msr_marginal_inner_audit <- function(spec, parameters, starts = 3L) {
  .msr_scale_surface_ram(spec)
  joint <- RTMB::MakeADFun(.msr_marginal_objective(spec, TRUE), parameters, silent = TRUE)
  full <- unname(unlist(parameters)); random <- which(names(joint$par) %in% c("scale_z", "response_z"))
  if (!length(random)) return(list(passed = TRUE, gradient = 0, condition = 1,
    identity_error = 0, competing = FALSE, lower_mode = FALSE, starts = list()))
  if (length(random) > 256L)
    stop("The private prototype's dense inner-curvature audit supports at most 256 latent effects.", call. = FALSE)
  H <- joint$he(full)[random, random, drop = FALSE]
  geometry <- .msr_marginal_curvature(H)
  gradient <- max(abs(joint$gr(full)[random]))
  current <- joint$fn(full)
  replace_random <- function(z) { p <- full; p[random] <- z; p }
  fn <- function(z) joint$fn(replace_random(z))
  gr <- function(z) joint$gr(replace_random(z))[random]
  he <- function(z) joint$he(replace_random(z))[random, random, drop = FALSE]
  attempts <- lapply(seq_len(starts), function(i) {
    initial <- if (i == 1L) rep(0, length(random)) else
      .5 * sin(seq_along(random) + 2 * i)
    tryCatch({
      fit <- stats::nlminb(initial, fn, gr,
        control = list(iter.max = 300L, eval.max = 600L, rel.tol = 1e-12))
      fit$optimizer_convergence <- fit$convergence
      fit$refinement <- .msr_marginal_polish(fit$par, fn, gr, he)
      if (isTRUE(fit$refinement$passed)) {
        fit$par <- fit$refinement$par; fit$objective <- fit$refinement$objective
      }
      fit$stationary <- fit$convergence == 0L || isTRUE(fit$refinement$passed)
      fit$gradient <- max(abs(gr(fit$par)))
      fit$distance <- if (geometry$passed)
        sqrt(max(0, drop(crossprod(fit$par - full[random], H %*% (fit$par - full[random]))))) else Inf
      fit$competing <- FALSE
      if (fit$stationary && fit$gradient <= 1e-4 && fit$distance > 1) {
        alternative <- joint$he(replace_random(fit$par))[random, random, drop = FALSE]
        if (.msr_marginal_curvature(alternative)$passed && geometry$passed) {
          # Compare approximate modal volumes, not only peak densities.
          log_ratio <- current - fit$objective +
            sum(log(diag(chol(H)))) - sum(log(diag(chol(alternative))))
          fit$competing <- log_ratio >= log(.01)
        }
      }
      fit
    }, error = function(e) list(error = conditionMessage(e), objective = Inf,
      convergence = 1L, gradient = Inf, distance = Inf, competing = FALSE))
  })
  usable <- vapply(attempts, function(x) isTRUE(x$stationary) && x$gradient <= 1e-4, TRUE)
  lower <- any(vapply(attempts[usable], function(x) x$objective < current - 1e-5, TRUE))
  competing <- any(vapply(attempts, function(x) isTRUE(x$competing), TRUE))
  list(passed = geometry$passed && gradient <= 1e-4 && any(usable) && !lower && !competing,
    gradient = gradient, condition = geometry$condition, lower_mode = lower,
    competing = competing, starts = attempts,
    laplace_value = current + sum(log(diag(chol(H)))) - length(random) / 2 * log(2 * pi))
}

.msr_marginal_attempt <- function(spec, initial, control) {
  warnings <- character()
  candidate_value <- Inf
  out <- tryCatch(withCallingHandlers({
    obj <- .msr_marginal_adfun(spec)
    fit <- stats::nlminb(initial, obj$fn, obj$gr,
      control = list(iter.max = control$iter.max, eval.max = 2L * control$iter.max, rel.tol = 1e-11))
    candidate_value <- fit$objective
    if (is.finite(fit$objective) && max(abs(obj$gr(fit$par))) > 1e-4) {
      polish <- stats::optim(fit$par, obj$fn, obj$gr, method = "BFGS",
        control = list(maxit = control$iter.max, reltol = 1e-12))
      if (is.finite(polish$value) && polish$value <= fit$objective) {
        fit$par <- polish$par; fit$objective <- polish$value
        fit$convergence <- polish$convergence; fit$message <- "BFGS polish"
      }
    }
    value <- as.numeric(obj$fn(fit$par))
    candidate_value <- value
    cold <- .msr_marginal_adfun(spec)
    cold_value <- as.numeric(cold$fn(fit$par))
    gradient <- as.numeric(cold$gr(fit$par))
    outer_hessian <- function(theta) {
      if (length(cold$env$random)) stats::optimHess(theta, cold$fn, cold$gr,
        control = list(ndeps = rep(1e-4, length(theta)))) else cold$he(theta)
    }
    # Polish only a cold-reproducible neighborhood. A discontinuity or failed
    # inner search must remain a rejection, not become a convergence override.
    fit$refinement <- NULL
    if (all(is.finite(c(cold_value, gradient))) && abs(value - cold_value) <= 1e-5) {
      fit$refinement <- tryCatch(.msr_marginal_polish(fit$par, cold$fn,
        cold$gr, outer_hessian), error = function(e) list(passed = FALSE,
          error = conditionMessage(e)))
      if (isTRUE(fit$refinement$passed) && fit$refinement$objective <= cold_value + 1e-10) {
        fit$optimizer_par <- fit$par; fit$optimizer_objective <- fit$objective
        fit$par <- fit$refinement$par; fit$objective <- fit$refinement$objective
      }
    }
    # Even an unsuccessful refinement evaluates other outer points. Rebuild
    # before extracting a conditional mode, so it belongs to this exact theta.
    value <- as.numeric(obj$fn(fit$par)); candidate_value <- value
    cold <- .msr_marginal_adfun(spec)
    cold_value <- as.numeric(cold$fn(fit$par)); gradient <- as.numeric(cold$gr(fit$par))
    parameters <- cold$env$parList(x = fit$par, par = cold$env$last.par)
    random <- length(parameters$scale_z) + length(parameters$response_z)
    inner <- .msr_marginal_inner_audit(spec, parameters, control$inner_starts)
    identity_error <- if (random) abs(cold_value - inner$laplace_value) else
      abs(cold_value - .msr_marginal_objective(spec)(parameters))
    H <- outer_hessian(fit$par)
    geometry <- .msr_marginal_curvature(H)
    stationary <- (fit$convergence == 0L || isTRUE(fit$refinement$passed)) &&
      all(is.finite(c(value, gradient))) &&
      max(abs(gradient)) <= 1e-4 && abs(value - cold_value) <= 1e-5 &&
      identity_error <= 1e-6 && isTRUE(inner$passed)
    list(eligible = stationary && geometry$passed, stationary = stationary,
      objective = cold_value, parameters = parameters, par = fit$par,
      optim = fit, gradient_max = max(abs(gradient)), cold_error = abs(value - cold_value),
      identity_error = identity_error, inner = inner, hessian = H,
      hessian_method = if (random) "Numerical derivative of AD marginal gradient" else "AD joint Hessian",
      condition = geometry$condition, minimum_eigenvalue = geometry$minimum)
  }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) list(eligible = FALSE, stationary = FALSE,
      objective = candidate_value, error = conditionMessage(e)))
  out$warnings <- unique(warnings); out$initial <- initial
  out
}

# Private ML fitter. Parameter SDs are outer MLEs. Latent effects alone are
# integrated. Return compact numeric diagnostics, never differentiation tapes.
.msr_marginal_fit <- function(spec, control = list()) {
  .msr_scale_check_capability(spec, "marginal_likelihood", "laplace")
  allowed <- c("starts", "iter.max", "inner_starts")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% allowed))))
    stop("Marginal prototype control supports starts, iter.max and inner_starts.", call. = FALSE)
  ctl <- utils::modifyList(list(starts = 3L, iter.max = 500L, inner_starts = 3L), control)
  for (nm in names(ctl)) validate_scalar_numeric(ctl[[nm]], nm, integerish = TRUE, lower = 1)
  out <- list(specification = spec, control = ctl, attempts = list(), selected = NA_integer_,
    parameters = NULL, annual = NULL,
    uncertainty = list(available = FALSE, reason = "Annual prediction intervals and population confidence intervals are not calibrated."))
  class(out) <- "msr_marginal_prototype"
  if (length(spec$identification)) {
    out$diagnostics <- list(passed = FALSE, reason = spec$identification,
      classification = "known unbounded likelihood")
    out$inference_metadata <- .msr_scale_inference_metadata("marginal_likelihood", "laplace")
    return(out)
  }
  present <- c(scale = length(spec$parameters$log_scale_sd) > 0L,
               response = length(spec$parameters$log_response_sd) > 0L)
  blocks <- names(present)[present]
  faces <- if (!length(blocks)) list(character()) else
    lapply(0:(2^length(blocks) - 1L), function(mask)
      blocks[as.logical(intToBits(mask)[seq_along(blocks)])])
  attempts <- list()
  for (face in seq_along(faces)) {
    reduced <- .msr_marginal_reduce(spec, faces[[face]])
    template <- .msr_marginal_adfun(reduced)$par
    for (i in seq_len(ctl$starts)) {
      initial <- template
      if (i > 1L) {
        initial <- initial + .05 * sin(seq_along(initial) + i)
        initial[names(initial) == "scale_beta"] <- initial[names(initial) == "scale_beta"] +
          rep(c(-.6, .6), length.out = ctl$starts - 1L)[i - 1L]
        for (nm in c("log_scale_sd", "log_response_sd"))
          initial[names(initial) == nm] <- log(rep(c(.05, .6), length.out = ctl$starts - 1L)[i - 1L])
      }
      attempt <- .msr_marginal_attempt(reduced, initial, ctl)
      attempt$face <- face; attempt$zero_sd <- reduced$zero_sd
      attempts[[length(attempts) + 1L]] <- attempt
    }
  }
  out$attempts <- attempts
  eligible <- which(vapply(attempts, function(a) isTRUE(a$eligible), TRUE))
  if (!length(eligible)) {
    out$diagnostics <- list(passed = FALSE, reason = "No start passed the fixed marginal-likelihood numerical gates.")
    out$inference_metadata <- .msr_scale_inference_metadata("marginal_likelihood", "laplace")
    return(out)
  }
  best <- eligible[which.min(vapply(attempts[eligible], function(a) a$objective, 0))]
  selected <- attempts[[best]]
  stationary <- vapply(attempts, function(a) isTRUE(a$stationary), TRUE)
  # A failed audit cannot make a better objective disappear from the search.
  better <- any(vapply(attempts, function(a)
    is.finite(a$objective) && a$objective < selected$objective - 1e-5, TRUE))
  complete <- all(vapply(seq_along(faces), function(f)
    any(stationary & vapply(attempts, function(a) a$face == f, TRUE)), TRUE))
  same_face <- eligible[vapply(attempts[eligible], function(a) a$face == selected$face, TRUE)]
  competing <- any(vapply(attempts[same_face], function(a) {
    difference <- as.numeric(a$par - selected$par)
    distance <- sqrt(max(0, drop(crossprod(difference, selected$hessian %*% difference))))
    a$objective - selected$objective <= 2 && distance > 1
  }, TRUE))
  out$selected <- best; out$parameters <- selected$parameters
  out$log_likelihood <- -selected$objective
  out$zero_sd <- selected$zero_sd
  reduced <- .msr_marginal_reduce(spec, names(selected$zero_sd)[selected$zero_sd])
  out$annual <- .msr_marginal_scale_predictions(reduced, selected$parameters)
  out$group_sd <- c(scale = if (selected$zero_sd[["scale"]]) 0 else
    if (length(selected$parameters$log_scale_sd)) exp(selected$parameters$log_scale_sd) else NA_real_,
    response = if (selected$zero_sd[["response"]]) 0 else
    if (length(selected$parameters$log_response_sd)) exp(selected$parameters$log_response_sd) else NA_real_)
  out$diagnostics <- list(passed = complete && !better && !competing,
    reason = if (!complete) "At least one variance face has no stationary audited fit." else
      if (better) "A numerically ineligible fit has a better objective." else
        if (competing) "Distinct outer optima lie within two log-likelihood units." else
          "Fixed numerical gates passed; interval calibration is unestablished.",
    boundary_search_complete = complete, better_ineligible_fit = better,
    competing_optima = competing,
    gradient_max = selected$gradient_max, inner_gradient_max = selected$inner$gradient,
    outer_condition = selected$condition, cold_error = selected$cold_error,
    identity_error = selected$identity_error, attempted_starts = length(attempts),
    eligible_starts = length(eligible), boundary = any(selected$zero_sd),
    criteria = "Optimizer convergence 0 or independent Newton refinement with gradient <=1e-7 and positive-definite unrepaired curvature; outer/inner gradients <=1e-4; Hessians positive definite with condition <=1e12; cold objective error <=1e-5; Laplace identity error <=1e-6; inner modal audits; all variance faces stationary; no better ineligible fit or discovered competing outer optimum.")
  out$inference_metadata <- .msr_scale_inference_metadata("marginal_likelihood", "laplace")
  out
}

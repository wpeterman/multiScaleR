.msr_factorized_attempt <- function(spec, initial, iter.max, object = NULL,
                                    integration_step = 1, cold_audit = "ad") {
  warnings <- character(); candidate <- Inf
  candidate_parameters <- candidate_par <- candidate_optimizer <- NULL
  stage <- "objective construction"
  obj <- NULL; event_start <- 0L
  result <- tryCatch(withCallingHandlers({
    obj <- if (!is.null(object) && exists("obj", envir = object, inherits = FALSE))
      object$obj else .msr_factorized_optimizer(spec, step = integration_step)
    if (!is.null(object)) object$obj <- obj
    event_start <- length(obj$env$integration_events)
    # Certification and final auditing can request exactly the same curvature.
    # Cache only that numeric matrix at an identical parameter vector.
    hessian_point <- hessian_value <- hessian_generation <- NULL
    generation <- function() if (is.null(obj$env$integration_generation)) 0L else
      obj$env$integration_generation
    he <- function(par) {
      if (is.null(hessian_point) || !identical(as.numeric(par), hessian_point) ||
          !identical(generation(), hessian_generation)) {
        hessian_value <<- obj$he(par)
        hessian_point <<- as.numeric(par)
        hessian_generation <<- generation()
      }
      hessian_value
    }
    stage <- "optimization"
    fit <- stats::nlminb(initial, obj$fn, obj$gr,
      control = list(iter.max = iter.max, eval.max = 3L * iter.max, rel.tol = 1e-11))
    candidate <- fit$objective
    candidate_par <- fit$par; candidate_optimizer <- fit
    candidate_parameters <- obj$env$parList(x = fit$par, par = obj$env$last.par)
    stage <- "Newton certification"
    fit$refinement <- .msr_marginal_polish(fit$par, obj$fn, obj$gr, he, maxit = 6L)
    if (isTRUE(fit$refinement$passed) && fit$refinement$objective <= fit$objective + 1e-10) {
      fit$optimizer_par <- fit$par; fit$optimizer_objective <- fit$objective
      fit$par <- fit$refinement$par; fit$objective <- fit$refinement$objective
    }
    value <- as.numeric(obj$fn(fit$par)); candidate <- value
    gradient <- as.numeric(obj$gr(fit$par))
    parameters <- obj$env$parList(x = fit$par, par = obj$env$last.par)
    candidate_parameters <- parameters; candidate_par <- fit$par; candidate_optimizer <- fit
    stage <- "outer curvature"
    H <- he(fit$par); geometry <- .msr_marginal_curvature(H)
    # New, stricter integration and different subdivision boundaries audit the
    # fitted point without relying on the optimizing object's integration state.
    stage <- "cold integration"
    cold_factory <- if (cold_audit == "score") .msr_factorized_score_adfun else .msr_factorized_adfun
    cold <- cold_factory(spec, parameters, rel.tol = 1e-10, step = .5)
    cold_value <- as.numeric(cold$fn(fit$par)); cold_gradient <- as.numeric(cold$gr(fit$par))
    stage <- "ordinary integration"
    ordinary <- .msr_factorized_objective(spec, FALSE, rel.tol = 1e-10, step = .5)(parameters)
    # A separate integration coordinate resolves response transitions that
    # become extremely narrow in standardized z at a large annual SD.
    stage <- "coordinate audit"
    scale_reference <- if (length(parameters$log_scale_sd) &&
      exp(parameters$log_scale_sd) >= 8)
      .msr_factorized_scale_objective(spec, FALSE)(parameters) else ordinary
    error <- max(abs(c(value - cold_value, value - ordinary, value - scale_reference)))
    gradient_error <- max(abs(gradient - cold_gradient))
    stationary <- (fit$convergence == 0L || isTRUE(fit$refinement$passed)) &&
      all(is.finite(c(value, gradient, cold_value, cold_gradient, ordinary))) &&
      max(abs(gradient)) <= 1e-4 && max(abs(cold_gradient)) <= 1e-4 &&
      error <= 1e-6 && gradient_error <= 1e-4
    list(eligible = stationary && geometry$passed, stationary = stationary,
      objective = value, parameters = parameters, par = fit$par, optim = fit,
      gradient_max = max(abs(gradient)), hessian = H, condition = geometry$condition,
      minimum_eigenvalue = geometry$minimum, cold_error = abs(value - cold_value),
      reference_error = abs(value - ordinary), gradient_error = gradient_error,
      coordinate_error = abs(value - scale_reference),
      hessian_method = "AD Hessian of adaptive factorized integral")
  }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
    error = function(e) list(eligible = FALSE, stationary = FALSE, objective = candidate,
      parameters = candidate_parameters, par = candidate_par, optim = candidate_optimizer,
      error_stage = stage, error = conditionMessage(e)))
  result$warnings <- unique(warnings); result$initial <- initial
  result$cold_audit <- cold_audit
  result$integration_step <- integration_step
  result$integration_events <- if (!is.null(obj) && length(obj$env$integration_events) > event_start)
    obj$env$integration_events[seq.int(event_start + 1L, length(obj$env$integration_events))] else list()
  result$integration_final_step <- if (!is.null(obj) && !is.null(obj$env$integration_current_step))
    obj$env$integration_current_step else integration_step
  result
}

# Certify that no conditional mode lies beyond the numerical search grid.
# A Gaussian density is bounded by its residual-density peak; a discrete
# probability is bounded by one. Multiplication by the standard-normal prior
# then bounds the entire exterior, without truncating its integration measure.
.msr_factorized_modes <- function(log_density, log_likelihood_bound) {
  if (!is.function(log_density) || length(log_likelihood_bound) != 1L ||
      !is.finite(log_likelihood_bound)) stop("A finite likelihood-density bound is required.", call. = FALSE)
  radius <- 12
  repeat {
    grid <- seq(-radius, radius, by = .125)
    values <- vapply(grid, log_density, 0)
    if (anyNA(values) || any(values == Inf) || !any(is.finite(values)))
      stop("Conditional mode search has invalid log densities.", call. = FALSE)
    tail_bound <- log_likelihood_bound - .5 * log(2 * pi) - radius^2 / 2
    if (tail_bound < max(values) - 1) break
    if (radius >= 384)
      stop("Conditional mode search could not certify its exterior density bound.", call. = FALSE)
    radius <- 2 * radius
  }
  peaks <- which(values >= c(-Inf, utils::head(values, -1)) &
    values >= c(utils::tail(values, -1), -Inf))
  candidates <- lapply(peaks, function(j) stats::optimize(function(z)-log_density(z),
    c(grid[max(1L,j-1L)],grid[min(length(grid),j+1L)]), tol = 1e-10))
  best <- candidates[[which.min(vapply(candidates, `[[`, 0, "objective"))]]
  list(candidates = candidates, best = best, radius = radius, exterior_log_bound = tail_bound)
}

.msr_factorized_predictions <- function(spec, parameters) {
  sm <- spec$scales[[1L]]; term <- sm$random
  parameters$scale_z <- numeric()
  if (!length(term)) {
    answer <- .msr_marginal_scale_predictions(spec, parameters)
    answer$conditional_mean_log <- answer$log_estimate
    answer$conditional_mode_z <- 0
    return(answer)
  }
  term <- term[[1L]]; modes <- means <- numeric(length(term$levels))
  for (i in seq_along(modes)) {
    g <- .msr_factorized_group(spec, i)
    log_density <- function(z) {
      p <- parameters; p$scale_z <- z
      -sum(.msr_factorized_response_nll(g, p)) + stats::dnorm(z, log = TRUE)
    }
    log_bound <- if (spec$family == "gaussian")
      -length(g$y) * (parameters$log_dispersion + .5 * log(2 * pi)) else 0
    search <- .msr_factorized_modes(log_density, log_bound)
    candidates <- search$candidates; best <- search$best
    modes[i] <- best$minimum; shift <- -best$objective
    cuts <- sort(unique(c(-Inf, seq(-8,8,1), vapply(candidates, `[[`, 0, "minimum"), Inf)))
    integrate <- function(moment) sum(vapply(seq_len(length(cuts)-1L), function(j)
      stats::integrate(function(z)vapply(z,function(t)
        moment(t) * exp(log_density(t)-shift),0), cuts[j],cuts[j+1L],
        rel.tol=1e-9,abs.tol=1e-12,subdivisions=500L)$value,0))
    means[i] <- integrate(identity) / integrate(function(z)1)
  }
  parameters$scale_z <- modes
  answer <- .msr_marginal_scale_predictions(spec, parameters)
  fixed <- drop(sm$X[sm$representatives,,drop=FALSE] %*% parameters$scale_beta[sm$beta_index]) + log(spec$unit_conv)
  answer$conditional_mean_log <- fixed + exp(parameters$log_scale_sd[term$sd_index]) * means[term$index[sm$representatives]]
  answer$conditional_mode_z <- modes[term$index[sm$representatives]]
  answer
}

.msr_factorized_fit <- function(spec, control = list(), initial_parameters = list(), checkpoint = NULL) {
  .msr_factorized_validate(spec)
  if(!is.null(checkpoint)&&!is.function(checkpoint)) stop("checkpoint must be a function or NULL.",call.=FALSE)
  allowed <- c("starts", "iter.max", "integration_step", "cold_audit")
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyDuplicated(names(control)) || any(!names(control) %in% allowed))))
    stop("Factorized control supports starts, iter.max, integration_step and cold_audit.", call. = FALSE)
  ctl <- utils::modifyList(list(starts=3L,iter.max=250L,integration_step=1L,cold_audit="ad"),control)
  if (!is.character(ctl$cold_audit) || length(ctl$cold_audit) != 1L ||
      is.na(ctl$cold_audit) || !ctl$cold_audit %in% c("ad", "score"))
    stop("cold_audit must be 'ad' or 'score'.", call. = FALSE)
  for(nm in setdiff(names(ctl), "cold_audit")) validate_scalar_numeric(ctl[[nm]],nm,integerish=TRUE,lower=1)
  if (!ctl$integration_step %in% c(1L,2L,4L,8L))
    stop("integration_step must be one of 1, 2, 4 or 8.", call. = FALSE)
  if (length(spec$response$random)) {
    index <- spec$response$random[[1L]]$index
    sizes <- tabulate(index)
    if (max(sizes) <= 1L) spec$identification <- c(spec$identification,
      "Gaussian residual and random-intercept variances cannot be separated with singleton response groups.")
    adjusted <- spec$y - spec$offset
    if (max(sizes) > 1L && all(vapply(split(adjusted, index),
      function(y) length(unique(y)) == 1L, TRUE)))
      spec$identification <- c(spec$identification,
        "Gaussian likelihood is unbounded as residual SD approaches zero for constant responses within every response group.")
  }
  out <- list(specification=spec,control=ctl,attempts=list(),selected=NA_integer_,
    parameters=NULL,annual=NULL,uncertainty=list(available=FALSE,
      reason="Conditional summaries omit population-parameter uncertainty; calibrated intervals are unavailable."))
  class(out) <- "msr_factorized_prototype"
  if(length(spec$identification)) {
    out$diagnostics <- list(passed=FALSE,reason=spec$identification)
    out$inference_metadata <- .msr_scale_inference_metadata("marginal_likelihood", "factorized")
    return(out)
  }
  faces <- if(length(spec$scales[[1L]]$random)) list(character(),"scale") else list(character())
  if (length(spec$response$random)) faces <- c(faces, lapply(faces, function(z)c(z,"response")))
  attempts <- list()
  for(face in seq_along(faces)) {
    reduced <- .msr_marginal_reduce(spec,faces[[face]])
    parameters <- reduced$parameters; parameters$scale_z <- parameters$response_z <- numeric()
    flatten <- function(p) {
      value <- unlist(p,use.names=FALSE)
      names(value) <- rep(names(p),lengths(p))
      value
    }
    template <- flatten(parameters)
    # Reuse the optimizing tape within one variance face. Construct it inside
    # the guarded attempt so failures remain recorded. Every cold audit is new.
    object <- new.env(parent = emptyenv())
    initials <- lapply(seq_len(ctl$starts),function(i) {
      p <- template
      if(i>1L) {
        p <- p + .05*sin(seq_along(p)+i)
        p[names(p)=="scale_beta"] <- p[names(p)=="scale_beta"] +
          rep(c(-.6,.6),length.out=ctl$starts-1L)[i-1L]
        p[names(p)=="log_scale_sd"] <- log(rep(c(.05,.6),length.out=ctl$starts-1L)[i-1L])
        p[names(p)=="log_response_sd"] <- log(rep(c(.05,.6),length.out=ctl$starts-1L)[i-1L])
      }
      p
    })
    for(p in initial_parameters) {
      p$scale_z <- p$response_z <- numeric()
      if(isTRUE(reduced$zero_sd[["scale"]])) p$log_scale_sd <- numeric()
      if(isTRUE(reduced$zero_sd[["response"]])) p$log_response_sd <- numeric()
      p <- flatten(p)
      if(length(p)==length(template) && identical(names(p),names(template)) && all(is.finite(p)))
        initials[[length(initials)+1L]] <- p
    }
    initials <- initials[!duplicated(lapply(initials,unname))]
    for(initial in initials) {
      # A wider optimizing subdivision changes work allocation only. Every
      # attempt receives fresh width-0.5 cold and ordinary-R audits. AD is the
      # default cold method; conditional score expectations are opt-in.
      a <- .msr_factorized_attempt(reduced,initial,ctl$iter.max,object,ctl$integration_step,ctl$cold_audit)
      a$face <- face; a$zero_sd <- reduced$zero_sd
      attempts[[length(attempts)+1L]] <- a
      if(!is.null(checkpoint)) checkpoint(attempts)
      gc(FALSE)
    }
  }
  out$attempts <- attempts
  eligible <- which(vapply(attempts,function(a)isTRUE(a$eligible),TRUE))
  if(!length(eligible)) {
    out$diagnostics <- list(passed=FALSE,reason="No start passed fixed factorized numerical and integration gates.")
    out$inference_metadata <- .msr_scale_inference_metadata("marginal_likelihood", "factorized")
    return(out)
  }
  values <- vapply(attempts[eligible],`[[`,0,"objective")
  tied <- eligible[values<=min(values)+1e-6]
  best <- tied[which.max(vapply(attempts[tied],function(a)sum(a$zero_sd),0))]
  selected <- attempts[[best]]
  stationary <- vapply(attempts,function(a)isTRUE(a$stationary),TRUE)
  complete <- all(vapply(seq_along(faces),function(f)
    any(stationary & vapply(attempts,function(a)a$face==f,TRUE)),TRUE))
  better <- any(vapply(attempts,function(a)is.finite(a$objective)&&a$objective<selected$objective-1e-5,TRUE))
  same <- eligible[vapply(attempts[eligible],function(a)a$face==selected$face,TRUE)]
  competing <- any(vapply(attempts[same],function(a) {
    d <- as.numeric(a$par-selected$par)
    a$objective-selected$objective<=2 && sqrt(max(0,drop(crossprod(d,selected$hessian%*%d))))>1
  },TRUE))
  out$selected <- best; out$parameters <- selected$parameters
  out$zero_sd <- selected$zero_sd; out$log_likelihood <- -selected$objective
  out$group_sd <- if(selected$zero_sd[["scale"]])0 else
    if(length(selected$parameters$log_scale_sd))exp(selected$parameters$log_scale_sd) else NA_real_
  if (length(spec$response$random)) out$response_group_sd <-
    if (selected$zero_sd[["response"]]) 0 else exp(selected$parameters$log_response_sd)
  out$diagnostics <- list(passed=complete&&!better&&!competing,
    reason=if(!complete)"A variance face has no stationary audited fit." else
      if(better)"A better numerically ineligible objective was found." else
        if(competing)"Competing outer optima lie within two log-likelihood units." else
          "Fixed numerical and integration gates passed; interval calibration is unestablished.",
    boundary_search_complete=complete,better_ineligible_fit=better,competing_optima=competing,
    gradient_max=selected$gradient_max,outer_condition=selected$condition,
    cold_error=selected$cold_error,reference_error=selected$reference_error,
    gradient_error=selected$gradient_error,attempted_starts=length(attempts),eligible_starts=length(eligible))
  out$diagnostics$integration_requested_step <- ctl$integration_step
  out$diagnostics$cold_audit <- ctl$cold_audit
  out$diagnostics$integration_errors <- sum(vapply(attempts,
    function(a) length(a$integration_events), 0L))
  out$diagnostics$integration_refinements <- sum(vapply(attempts, function(a)
    sum(vapply(a$integration_events, function(e) !is.na(e$to), TRUE)), 0L))
  out$diagnostics$integration_finer_attempts <- sum(vapply(attempts,
    function(a) a$integration_final_step < ctl$integration_step, TRUE))
  # A nearly stationary, finite-curvature point can still approximate a kernel
  # or Poisson limit. Such a witness rejects identification without adopting
  # an ineligible limiting model as the selected finite-scale estimate.
  out$limits <- tryCatch(.msr_factorized_limit_audit(spec, selected$parameters, selected$objective),
    error = function(e) list(completed = FALSE, witnesses = character(), error = conditionMessage(e)))
  out$diagnostics$limiting_audit_complete <- isTRUE(out$limits$completed)
  out$diagnostics$limiting_witnesses <- out$limits$witnesses
  if (!isTRUE(out$limits$completed) || length(out$limits$witnesses)) {
    out$diagnostics$passed_before_limiting_audit <- out$diagnostics$passed
    out$diagnostics$passed <- FALSE
    out$diagnostics$reason <- paste(out$diagnostics$reason,
      if (!isTRUE(out$limits$completed)) "Required limiting-model evaluation failed." else
        paste("Finite-scale identification is unresolved near", paste(out$limits$witnesses, collapse=", ")))
  }
  reduced <- .msr_marginal_reduce(spec,names(selected$zero_sd)[selected$zero_sd])
  annual <- tryCatch(.msr_factorized_predictions(reduced,selected$parameters), error = function(e)e)
  if (inherits(annual, "error")) {
    # Preserve the completed fit and all attempts when a derived summary fails.
    out$prediction_error <- conditionMessage(annual)
    out$diagnostics$fit_passed_before_prediction <- out$diagnostics$passed
    out$diagnostics$passed <- FALSE
    out$diagnostics$reason <- paste("Annual scale summaries failed:", out$prediction_error)
  } else out$annual <- annual
  out$inference_metadata <- .msr_scale_inference_metadata("marginal_likelihood", "factorized")
  out
}


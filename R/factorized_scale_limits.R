# Fixed-kernel limits give feasible likelihood witnesses. They do not provide
# finite annual-scale estimates, confidence intervals, or a global certificate.
.msr_factorized_limit_spec <- function(spec, kernel, poisson = FALSE) {
  .msr_factorized_validate(spec)
  if (!kernel %in% c("broad", "nearest")) stop("Unknown kernel limit.", call. = FALSE)
  if (poisson && spec$family != "nbinom2") stop("Poisson limit requires NB2.", call. = FALSE)
  # Shared or independently fixed scales have no annual SD block to remove.
  # Their fixed-kernel witnesses still need exactly the same likelihood audit.
  out <- if (length(spec$scales[[1L]]$random))
    .msr_marginal_reduce(spec, "scale") else spec
  keep <- if (kernel == "nearest") spec$delta[[1L]] == 0 else
    matrix(TRUE, nrow(spec$V[[1L]]), ncol(spec$V[[1L]]))
  out$V[[1L]] <- matrix(rowSums(spec$V[[1L]] * keep), ncol = 1L)
  out$C[[1L]] <- matrix(rowSums(spec$C[[1L]] * keep), ncol = 1L)
  out$delta[[1L]] <- matrix(0, length(spec$y), 1L)
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
  if (length(spec$response$random) && !length(parameters$log_response_sd))
    spec <- .msr_marginal_reduce(spec, "response")
  references <- list(); witnesses <- character()
  for (kernel in c("broad", "nearest")) for (poisson in
    if (spec$family == "nbinom2") c(FALSE, TRUE) else FALSE) {
    label <- paste(kernel, if (poisson) "poisson" else spec$family, sep = "-")
    s <- .msr_factorized_limit_spec(spec, kernel, poisson)
    p <- parameters
    p$scale_beta <- p$scale_z <- p$log_scale_sd <- numeric()
    if (poisson) p$log_dispersion <- numeric()
    attempts <- list()
    obj <- tryCatch(.msr_factorized_adfun(s, p), error = function(e)e)
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
        ordinary <- .msr_factorized_objective(s, FALSE)(q)
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
      value <- .msr_factorized_objective(s, FALSE, 1e-10, .5)(q)
      coordinate <- if (length(q$log_scale_sd) && exp(q$log_scale_sd) >= 8)
        .msr_factorized_scale_objective(s, FALSE)(q) else value
      list(completed = all(is.finite(c(value, coordinate))) && abs(value - coordinate) <= 1e-6,
        objective = value, coordinate_error = abs(value - coordinate))
    }, error = function(e) list(completed = FALSE, error = conditionMessage(e)))
    references[["fitted-scales-poisson"]] <- a
    if (isTRUE(a$completed) && a$objective <= candidate_nll + 1e-5)
      witnesses <- c(witnesses, "fitted-scales-poisson")
  }
  list(completed = all(vapply(references, function(a)isTRUE(a$completed), TRUE)),
    witnesses = witnesses, references = references, tolerance = 1e-5,
    interpretation = "Feasible limiting-model exclusion witnesses; absence does not prove global identification.")
}

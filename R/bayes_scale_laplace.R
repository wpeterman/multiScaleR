# Joint Laplace approximation of the existing posterior. Keeping this objective
# identical to NUTS makes approximation error distinguishable from model changes.
.msr_bayes_laplace <- function(out, draws, seed, control) {
  out$inference_metadata <- NULL
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Laplace inference requires the optional package 'RTMB'.", call. = FALSE)
  validate_scalar_numeric(draws, "laplace_draws", integerish = TRUE, lower = 2)
  if (!is.null(seed))
    validate_scalar_numeric(seed, "seed", integerish = TRUE, lower = 1,
                            upper = .Machine$integer.max)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% c("starts", "iter.max")))))
    stop("Laplace control supports only starts and iter.max.", call. = FALSE)
  ctl <- utils::modifyList(list(starts = 6L, iter.max = 1000L), control)
  validate_scalar_numeric(ctl$starts, "starts", integerish = TRUE, lower = 1)
  validate_scalar_numeric(ctl$iter.max, "iter.max", integerish = TRUE, lower = 1)
  spec <- out$specification
  obj <- RTMB::MakeADFun(.msr_bayes_objective(spec, ad = TRUE),
                         spec$parameters, silent = TRUE)
  initial <- lapply(seq_len(ctl$starts), function(i) {
    p <- obj$par
    if (i == 1L) return(p)
    # Vary distance, SD, and latent effects together without consuming RNG.
    p <- p + .15 * sin(seq_along(p) + i)
    p[names(p) == "scale_beta"] <- p[names(p) == "scale_beta"] +
      rep(c(-1, 1, -.5, .5, 0), length.out = ctl$starts - 1L)[i - 1L]
    for (nm in c("log_scale_sd", "log_response_sd"))
      p[names(p) == nm] <- p[names(p) == nm] +
        rep(c(-.7, .7, -.3, .3, 0), length.out = ctl$starts - 1L)[i - 1L]
    p
  })
  attempts <- lapply(initial, function(p) {
    tryCatch({
      a <- stats::nlminb(p, obj$fn, obj$gr,
        control = list(iter.max = ctl$iter.max, eval.max = 2 * ctl$iter.max,
                       rel.tol = 1e-10))
      # Tight stationarity matters more than the optimizer's relative stopping
      # rule, which can stop early for large negative log likelihoods.
      if (max(abs(obj$gr(a$par))) > 1e-4) {
        b <- stats::optim(a$par, obj$fn, obj$gr, method = "BFGS",
          control = list(maxit = ctl$iter.max, reltol = 1e-12))
        if (is.finite(b$value) && b$value <= a$objective) {
          a$par <- b$par
          a$objective <- b$value
          a$convergence <- b$convergence
          a$message <- "BFGS polish"
        }
      }
      # Relative objective tolerances can leave a small gradient even after
      # BFGS. Refine with positive-curvature Newton steps and a line search;
      # never repair indefinite curvature or weaken the stationarity gate.
      a$newton_steps <- 0L
      for (step in seq_len(8L)) {
        gradient <- as.numeric(obj$gr(a$par))
        if (max(abs(gradient)) <= 1e-6) {
          a$convergence <- 0L
          break
        }
        R <- tryCatch(chol(obj$he(a$par)), error = function(e) NULL)
        if (is.null(R)) break
        direction <- as.numeric(backsolve(R, forwardsolve(t(R), gradient)))
        accepted <- FALSE
        for (half in 0:12) {
          candidate <- a$par - direction * 2^-half
          value <- obj$fn(candidate)
          # Allow roundoff in a nearly identical objective only when the
          # gradient improves. This is not an objective-increasing restart.
          roundoff <- 1e-12 * max(1, abs(a$objective))
          if (is.finite(value) && value <= a$objective + roundoff &&
              max(abs(obj$gr(candidate))) < max(abs(gradient))) {
            a$par <- candidate
            a$objective <- value
            a$newton_steps <- a$newton_steps + 1L
            accepted <- TRUE
            break
          }
        }
        if (!accepted) break
      }
      if (max(abs(obj$gr(a$par))) <= 1e-6) a$convergence <- 0L
      a$gradient_max <- max(abs(obj$gr(a$par)))
      a$usable <- is.finite(a$objective) && a$convergence == 0L &&
        is.finite(a$gradient_max) && a$gradient_max <= 1e-4
      a
    }, error = function(e) list(usable = FALSE, objective = Inf,
                                message = conditionMessage(e)))
  })
  good <- which(vapply(attempts, function(a) isTRUE(a$usable), logical(1)))
  if (!length(good))
    .msr_bayes_laplace_failure("Laplace optimization found no converged stationary mode.", attempts)
  best <- good[which.min(vapply(attempts[good], function(a) a$objective, 0))]
  mode <- attempts[[best]]$par
  hessian <- obj$he(mode)
  hessian <- (hessian + t(hessian)) / 2
  geometry <- tryCatch(.msr_bayes_laplace_geometry(hessian),
    error = function(e) .msr_bayes_laplace_failure(conditionMessage(e), attempts))
  covariance <- geometry$covariance
  # Distances are measured in posterior standard deviations. This screens for
  # competing modes; deterministic starts cannot prove a mode is global.
  distances <- vapply(attempts[good], function(a) {
    d <- a$par - mode
    sqrt(max(0, drop(crossprod(d, hessian %*% d))))
  }, 0)
  objectives <- vapply(attempts[good], function(a) a$objective, 0)
  competing <- any(objectives - attempts[[best]]$objective <= 2 & distances > 1)
  samples <- .msr_bayes_laplace_draws(mode, geometry$root, draws, seed)
  ends <- cumsum(lengths(spec$parameters))
  begins <- ends - lengths(spec$parameters)
  out$draws <- lapply(seq_along(ends), function(j)
    samples[, begins[j] + seq_len(ends[j] - begins[j]), drop = FALSE])
  names(out$draws) <- names(spec$parameters)
  out$inference <- "joint Laplace Gaussian approximation"
  out$laplace <- list(mode = mode, covariance = covariance, hessian = hessian,
    starts = attempts, selected_start = best, draws = draws, seed = seed,
    coordinates = "Unconstrained, noncentered; log SDs include prior Jacobians")
  out$diagnostics <- list(passed = !competing, type = "Laplace numerical checks",
    gradient_max = attempts[[best]]$gradient_max,
    min_eigenvalue = geometry$min_eigenvalue, condition_number = geometry$condition_number,
    competing_modes = competing, successful_starts = length(good), starts = ctl$starts,
    criteria = paste("Optimizer convergence 0; max absolute gradient <= 1e-4;",
      "positive-definite Hessian; condition number <= 1e12;",
      "no distinct discovered mode within 2 log-posterior units.",
      "These checks do not establish Gaussian approximation accuracy."))
  if (competing)
    warning("Laplace approximation found competing posterior modes; inspect diagnostics and compare with full Bayesian sampling.",
            call. = FALSE)
  out$inference_metadata <- .msr_scale_inference_metadata("posterior", "joint_laplace")
  out
}

.msr_bayes_laplace_failure <- function(message, attempts) {
  # Preserve failed starts for reproducible benchmark and user diagnostics.
  stop(structure(list(message = message, call = NULL, starts = attempts),
    class = c("msr_laplace_error", "error", "condition")))
}

.msr_bayes_laplace_geometry <- function(hessian) {
  if (any(!is.finite(hessian)))
    stop("Laplace Hessian contains nonfinite values.", call. = FALSE)
  ev <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
  condition <- max(ev) / min(ev)
  if (min(ev) <= 0 || !is.finite(condition) || condition > 1e12)
    stop("Laplace Hessian is not positive definite or is too ill-conditioned; no covariance was regularized.",
         call. = FALSE)
  # If H = R'R, columns of R^-1 transform independent normals to H^-1.
  root <- backsolve(chol(hessian), diag(nrow(hessian)))
  list(root = root, covariance = tcrossprod(root),
       min_eigenvalue = min(ev), condition_number = condition)
}

.msr_bayes_laplace_draws <- function(mode, root, draws, seed) {
  if (!is.null(seed)) {
    exists_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (exists_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (exists_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  sweep(matrix(stats::rnorm(draws * length(mode)), draws) %*% t(root),
        2L, mode, "+")
}

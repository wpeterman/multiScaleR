# Likelihood profiles of kernel sigma for the joint_mixed backend.

# Cold native refit shared by the joint_mixed fit (every start and its
# terminal audit) and by its likelihood profiles. Each call builds a fresh
# glmmTMB fit from the original formula, family, zero-inflation and
# dispersion formulas, and contrasts, with the kernel covariates evaluated at
# log scales `kp` (engine log-distance units). Nothing carries over between
# calls, so the result depends only on `kp`.
.msr_mixed_refit_fn <- function(model, engine, covs, data) {
  function(kp) {
    nd <- data
    kd <- .msr_eval_kernel_engine(engine, unname(kp))
    for (nm in covs) nd[[nm]] <- kd$value[, nm]
    glmmTMB::glmmTMB(formula = stats::formula(model), data = nd,
      family = stats::family(model), ziformula = model$modelInfo$allForm$ziformula,
      dispformula = model$modelInfo$allForm$dispformula,
      contrasts = model$modelInfo$contrasts, na.action = stats::na.fail)
  }
}


# Profile evaluator for a joint_mixed fit, in the form that
# .msr_profile_likelihood_intervals() and profile_sigma() expect: the fitted
# log scales `par`, their search bounds, and `at(j, z, reoptimize)`, which
# returns the negative log marginal likelihood (Laplace approximation, as in
# the fit) with log scale j held at z.
#
# * One scale, or `reoptimize = FALSE`: every other scale stays at its
#   estimate, and the value is the cold native refit at those scales. All
#   coefficients, the random-intercept SD, and the NB2 dispersion are
#   re-estimated. This is the refit the fit audits each start with.
# * Several scales with `reoptimize = TRUE`: the other scales are also
#   re-estimated, on the joint RTMB objective with scale j held at z
#   (equal L-BFGS-B bounds), from the optimum of every eligible start. A cold
#   native refit at the resulting scales audits the value, and the lower of
#   the two is returned.
#
# A failed evaluation stops with an error rather than returning a value that
# could pass for a likelihood. With `reoptimize = TRUE`, a value below the
# stored optimum by more than `tolerance` also stops, because the fit missed a
# better optimum. The tolerance exceeds the fit's 1e-5 native audit, so the
# stored optimum itself never trips it. Values are cached by scale and z.
.msr_mixed_profiler <- function(object, tolerance = 1e-4) {
  for (pkg in c("RTMB", "TMB", "glmmTMB"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Profiling a joint_mixed fit requires the optional package ", pkg, ".", call. = FALSE)
  mf <- object$mixed_fit
  model <- object$fitted_mod_original
  p <- nrow(object$scale_est)
  if (!is.list(mf) || is.null(mf$engine) || !is.data.frame(mf$prediction_data) ||
      !inherits(model, "glmmTMB") || length(mf$log_scale) != p ||
      length(mf$lower) != p || length(mf$upper) != p ||
      !is.finite(object$optim_results$value))
    stop("This joint_mixed fit lacks the stored data needed to profile sigma. ",
         "Refit it with the current version of multiScaleR.", call. = FALSE)
  covs <- mf$engine$covariates
  par <- unname(mf$log_scale)
  lower <- unname(mf$lower)
  upper <- unname(mf$upper)
  optimum <- object$optim_results$value
  unit_conv <- object$kernel_inputs$unit_conv
  refit <- .msr_mixed_refit_fn(model, mf$engine, covs, mf$prediction_data)

  # Negative log-likelihood of the cold native refit, or NA when the refit
  # errors or nlminb does not report convergence. Warnings about the native
  # Hessian are ignored: the profile needs only the maximized likelihood.
  native_value <- function(kp) {
    m <- tryCatch(suppressWarnings(refit(kp)), error = function(e) NULL)
    if (!inherits(m, "glmmTMB") || !identical(m$fit$convergence, 0L)) return(NA_real_)
    value <- -as.numeric(stats::logLik(m))
    if (is.finite(value)) value else NA_real_
  }

  # One RTMB tape, built on first use from the native fit at the estimate.
  # The design does not depend on the scales (the fit audits this), so the
  # same tape serves every profile point.
  tape <- NULL
  joint_tape <- function() {
    if (is.null(tape)) {
      initial <- suppressWarnings(refit(par))
      args <- .msr_native_mixed_design(initial, mf$engine, covs, par,
                                       rownames(mf$prediction_data))
      obj <- do.call(.msr_mixed_rtmb, args)
      ki <- which(names(obj$par) == "kernel_par")
      if (length(ki) != p) stop("Joint scale parameter mapping failed.", call. = FALSE)
      # Starts: the selected start's optimum first, then every other eligible
      # start's optimum, keeping one start per distinct set of log scales.
      starts <- list()
      for (i in unique(c(mf$selected, seq_along(mf$start_results)))) {
        o <- mf$start_results[[i]]
        s <- o$optim$par
        if (isTRUE(o$eligible) && identical(names(s), names(obj$par)) && all(is.finite(s)))
          starts[[length(starts) + 1L]] <- unname(s)
      }
      if (!length(starts)) starts <- list(unname(obj$par))
      keys <- vapply(starts, function(s) paste(round(s[ki], 3), collapse = ","), "")
      tape <<- list(obj = obj, ki = ki, starts = starts[!duplicated(keys)])
    }
    tape
  }

  joint_value <- function(j, z) {
    tp <- joint_tape()
    lo <- rep(-Inf, length(tp$obj$par))
    hi <- rep(Inf, length(tp$obj$par))
    lo[tp$ki] <- lower
    hi[tp$ki] <- upper
    lo[tp$ki[j]] <- hi[tp$ki[j]] <- z
    best <- NULL
    for (s in tp$starts) {
      s[tp$ki] <- pmin(pmax(s[tp$ki], lo[tp$ki]), hi[tp$ki])
      opt <- tryCatch(suppressWarnings(stats::optim(s, tp$obj$fn, tp$obj$gr,
        method = "L-BFGS-B", lower = lo, upper = hi,
        control = list(maxit = 300, factr = 1e5, pgtol = 1e-7))),
        error = function(e) NULL)
      if (!is.null(opt) && identical(opt$convergence, 0L) && is.finite(opt$value) &&
          (is.null(best) || opt$value < best$value)) best <- opt
    }
    if (is.null(best)) return(NA_real_)
    kp <- unname(best$par[tp$ki])
    kp[j] <- z
    min(best$value, native_value(kp), na.rm = TRUE)
  }

  cache <- new.env(parent = emptyenv())
  at <- function(j, z, reoptimize = TRUE) {
    joint <- reoptimize && p > 1L
    key <- paste(j, joint, sprintf("%.17g", z))
    if (exists(key, envir = cache, inherits = FALSE))
      return(get(key, envir = cache, inherits = FALSE))
    kp <- par
    kp[j] <- z
    value <- if (joint) joint_value(j, z) else native_value(kp)
    if (!is.finite(value))
      stop("The native refit failed at sigma = ",
           format(signif(exp(z) * unit_conv, 6)), " for ", covs[j], ".", call. = FALSE)
    if (reoptimize && value < optimum - tolerance)
      stop("Profiling found a better optimum; refit with additional starts before inference.",
           call. = FALSE)
    assign(key, value, envir = cache)
    value
  }
  list(par = par, lower = lower, upper = upper, at = at)
}

# Fix selected log parameters and reoptimize the remainder. Retain independent
# starts so profiling does not silently inherit one nuisance-parameter basin.
.msr_profile_at <- function(object, index, value, reoptimize = TRUE) {
  pf <- object$profile_fit
  obj <- pf$contract$objective
  fixed <- pf$par
  fixed[index] <- value
  free <- if (reoptimize) setdiff(seq_along(fixed), index) else
    if (obj$family == "nbinom2") length(fixed) else integer()
  if (!length(free)) return(obj$fn(fixed))
  expand <- function(par) { z <- fixed; z[free] <- par; z }
  reduced <- list(fn = function(par) obj$fn(expand(par)),
                  gr = function(par) obj$gr(expand(par))[free])
  if (is.function(obj$reset_start)) reduced$reset_start <- obj$reset_start
  if (is.function(obj$cold_evaluate)) reduced$cold_evaluate <- function(par) {
    z <- obj$cold_evaluate(expand(par)); z$gradient <- z$gradient[free]; z
  }
  starts <- unique(rbind(pf$par, pf$starts)[, free, drop = FALSE])
  fits <- .msr_multistart(reduced, starts, pf$lower[free], pf$upper[free],
    gradient_tolerance = pf$control$gradient_tolerance,
    control = list(iter.max = pf$control$iter.max, eval.max = pf$control$eval.max,
                   rel.tol = 1e-11, x.tol = 1e-9))
  if (is.null(fits$best)) stop("No nuisance fit passed the profile convergence checks.", call. = FALSE)
  if (fits$best$objective < object$optim_results$value - 1e-5)
    stop("Profiling found a better optimum; refit with additional starts before inference.", call. = FALSE)
  fits$best$objective
}

.msr_profile_likelihood_intervals <- function(object, level = .95) {
  validate_scalar_numeric(level, "level", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)
  pf <- object$profile_fit
  p <- nrow(object$scale_est)
  target <- object$optim_results$value + stats::qchisq(level, 1) / 2
  limits <- matrix(NA_real_, p, 2L)
  bounded <- matrix(FALSE, p, 2L)
  for (j in seq_len(p)) {
    fn <- function(z) .msr_profile_at(object, j, z) - target
    for (side in 1:2) {
      bound <- if (side == 1L) pf$lower[j] else pf$upper[j]
      grid <- seq(pf$par[j], bound, length.out = 17L)
      # Stop at the first crossing of the connected region. Evaluating every
      # remaining grid point does not refine this endpoint and is expensive.
      cross <- NA_integer_
      for (k in seq_along(grid)) {
        if (fn(grid[k]) >= 0) {
          cross <- k
          break
        }
      }
      if (is.na(cross)) {
        limits[j, side] <- bound
        bounded[j, side] <- TRUE
      } else if (cross == 1L) {
        stop("Stored optimum failed the likelihood-profile check.", call. = FALSE)
      } else {
        limits[j, side] <- stats::uniroot(fn, sort(grid[c(cross - 1L, cross)]),
                                        tol = 1e-5)$root
      }
    }
  }
  limits <- exp(limits) * object$kernel_inputs$unit_conv
  rownames(limits) <- rownames(object$scale_est)
  colnames(limits) <- c("2.5%", "97.5%")
  dimnames(bounded) <- dimnames(limits)
  list(intervals = as.data.frame(limits),
       method = stats::setNames(ifelse(rowSums(bounded) > 0, "profile-bound-limited", "profile"),
                                rownames(limits)),
       bound_limited = bounded, level = level)
}

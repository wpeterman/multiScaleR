# Experimental matrix interface for independent observations. Supported links:
# Gaussian identity, Poisson/NB2 log, and binomial logit. X_fixed includes any
# intercept. Weights multiply log likelihoods (frequency/case weights), not
# Gaussian precision weights. The fitting sample must already be fixed.
#
# Beta is solved by IRLS at each scale evaluation. NB2 log(theta) is optimized
# alongside log(kernel scales), with an analytic score; it is not held fixed.
# Gaussian residual SD is profiled at its ML estimate. This is distinct from
# the spatial scale, which the public API calls sigma.
.msr_profile_glm_objective <- function(engine, y, family = "gaussian",
                                       X_fixed = matrix(1, engine$n, 1L),
                                       offset = rep(0, engine$n),
                                       weights = rep(1, engine$n),
                                       trials = rep(1, engine$n)) {
  family <- match.arg(family, c("gaussian", "poisson", "binomial", "nbinom2"))
  n <- engine$n
  valid_vector <- function(x) is.numeric(x) && length(x) == n && all(is.finite(x))
  if (!all(vapply(list(y, offset, weights, trials), valid_vector, logical(1))) ||
      any(weights <= 0) || any(trials < 1 | trials != floor(trials)))
    stop("Responses, offsets, positive weights, and integer trials must match the fitting rows.",
         call. = FALSE)
  if (family != "gaussian" && any(y < 0 | y != floor(y)))
    stop("Count responses must be nonnegative integers.", call. = FALSE)
  if (family == "binomial" && any(y > trials))
    stop("Binomial successes cannot exceed trials.", call. = FALSE)
  if (!is.matrix(X_fixed) || !is.numeric(X_fixed) || nrow(X_fixed) != n ||
      any(!is.finite(X_fixed)))
    stop("X_fixed must be a finite numeric matrix with one row per observation.",
         call. = FALSE)
  if (family == "nbinom2" && !requireNamespace("MASS", quietly = TRUE))
    stop("NB2 profiling requires the suggested package MASS.", call. = FALSE)
  kernel_npar <- length(engine$parameter_names)
  parameter_names <- c(engine$parameter_names,
                       if (family == "nbinom2") "log_theta")
  cache <- new.env(parent = emptyenv())
  cache$par <- NULL
  cache$kernel_par <- NULL
  cache$evaluations <- 0L
  evaluate <- function(par) {
    if (identical(par, cache$par)) return(cache$value)
    if (!is.numeric(par) || length(par) != length(parameter_names) ||
        any(!is.finite(par)))
      stop("Invalid profile parameter vector.", call. = FALSE)
    if (!is.null(names(par)) && !identical(names(par), parameter_names))
      stop("Profile parameter names are out of order.", call. = FALSE)
    kernel_par <- unname(par[seq_len(kernel_npar)])
    # Dispersion-only profile steps reuse the same landscape summaries.
    if (!identical(kernel_par, cache$kernel_par)) {
      cache$kernel <- .msr_eval_kernel_engine(engine, kernel_par)
      cache$kernel_par <- kernel_par
    }
    kd <- cache$kernel
    X <- cbind(X_fixed, kd$value)
    theta <- if (family == "nbinom2") exp(par[kernel_npar + 1L]) else NA_real_
    if (family == "nbinom2" && (!is.finite(theta) || theta <= 0))
      stop("NB2 theta must be positive and finite.", call. = FALSE)
    fam <- switch(family, gaussian = stats::gaussian(),
                  poisson = stats::poisson(), binomial = stats::binomial(),
                  nbinom2 = MASS::negative.binomial(theta))
    fit <- stats::glm.fit(
      X, if (family == "binomial") y / trials else y,
      weights = weights * if (family == "binomial") trials else 1,
      offset = offset, family = fam,
      control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
    if (family == "nbinom2" && !isTRUE(fit$converged) && fit$rank == ncol(X))
      fit <- .msr_nb2_refine_glm(fit, X, y, theta, offset, weights)
    if (!isTRUE(fit$converged) || fit$rank != ncol(X) ||
        any(!is.finite(fit$coefficients)))
      stop("Profile GLM did not converge to a full-rank finite fit.", call. = FALSE)
    mu <- fit$fitted.values
    residual_sd <- NA_real_
    if (family == "gaussian") {
      residual_sd <- sqrt(sum(weights * (y - mu)^2) / sum(weights))
      ll <- stats::dnorm(y, mu, residual_sd, log = TRUE)
      eta_score <- weights * (mu - y) / residual_sd^2
    } else if (family == "poisson") {
      ll <- stats::dpois(y, mu, log = TRUE)
      eta_score <- weights * (mu - y)
    } else if (family == "binomial") {
      ll <- stats::dbinom(y, trials, mu, log = TRUE)
      eta_score <- weights * (trials * mu - y)
    } else {
      ll <- stats::dnbinom(y, size = theta, mu = mu, log = TRUE)
      eta_score <- weights * (mu - y) / (1 + mu / theta)
    }
    beta_kernel <- utils::tail(fit$coefficients, engine$p)
    gradient <- colSums(kd$d_scale * eta_score) * beta_kernel
    if (engine$kernel == "expow")
      gradient <- c(gradient, colSums(kd$d_shape * eta_score) * beta_kernel)
    if (family == "nbinom2") {
      # Exact NB2 score for log(theta), including the normalizing constants.
      theta_score <- .msr_nb2_theta_score(y, mu, theta, weights)
      gradient <- c(gradient, theta_score)
    }
    nll <- -sum(weights * ll)
    if (!is.finite(nll) || any(!is.finite(gradient)))
      stop("Nonfinite profile likelihood or score.", call. = FALSE)
    value <- list(value = nll, gradient = unname(gradient),
                  beta = fit$coefficients, theta = theta,
                  residual_sd = residual_sd, fit = fit, design = X)
    cache$evaluations <- cache$evaluations + 1L
    cache$par <- par
    cache$value <- value
    value
  }
  list(fn = function(par) evaluate(par)$value,
       gr = function(par) evaluate(par)$gradient, evaluate = evaluate,
       evaluations = function() cache$evaluations,
       parameter_names = parameter_names, engine = engine, family = family,
       y = y, X_fixed = X_fixed, offset = offset, weights = weights, trials = trials)
}

# NB2 log-link Fisher scoring can oscillate even when the conditional
# likelihood is strictly concave. Use its observed curvature only as a
# fallback, then rebuild the ordinary GLM result at the recovered solution.
.msr_nb2_refine_glm <- function(fit, X, y, theta, offset, weights) {
  start <- fit$coefficients
  if (any(!is.finite(start))) start <- rep(0, ncol(X))
  eta <- function(beta) drop(X %*% beta) + offset
  fn <- function(beta) {
    linear <- eta(beta)
    if (any(linear > log(.Machine$double.xmax) - 1)) return(Inf)
    -sum(weights * stats::dnbinom(y, mu = exp(linear), size = theta, log = TRUE))
  }
  # Logistic forms avoid overflow in the score and observed information.
  gr <- function(beta) {
    linear <- eta(beta) - log(theta)
    drop(crossprod(X, weights * (theta * stats::plogis(linear) -
                                  y * stats::plogis(-linear))))
  }
  he <- function(beta) {
    linear <- eta(beta) - log(theta)
    w <- weights * (theta + y) * stats::plogis(linear) * stats::plogis(-linear)
    crossprod(X, X * w)
  }
  refined <- stats::nlminb(start, fn, gr, he,
    control = list(iter.max = 500L, eval.max = 1000L, rel.tol = 1e-12, x.tol = 1e-10))
  if (!is.finite(refined$objective) || max(abs(gr(refined$par))) > 1e-6)
    stop("Conditional NB2 refinement did not pass the coefficient score check.", call. = FALSE)
  stats::glm.fit(X, y, weights = weights, offset = offset,
    family = MASS::negative.binomial(theta), start = refined$par,
    control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
}

# Avoid subtracting log(theta) terms near the Poisson limit. For moderate
# counts the rising-factorial identity also avoids cancellation in digamma.
.msr_nb2_theta_score <- function(y, mu, theta, weights) {
  if (theta > 1000 && max(y) <= 10000) {
    sums <- c(0, cumsum(1 / (1 + seq.int(0, max(y)) / theta)))
    delta <- sums[y + 1L]
  } else {
    delta <- theta * (digamma(y + theta) - digamma(theta))
  }
  -sum(weights * (delta - theta * log1p(mu / theta) +
                    (mu - y) / (1 + mu / theta)))
}

# Use the same convergence screen for both experimental backends. At a bound,
# only a gradient pointing toward an allowed decrease violates stationarity.
.msr_projected_gradient <- function(par, gradient, lower, upper, tolerance = 1e-7) {
  g <- gradient
  g[par <= lower + tolerance & g > 0] <- 0
  g[par >= upper - tolerance & g < 0] <- 0
  g
}

# No Hessian during the starts: uncertainty is computed only for the selected
# fit by the caller. A failed start is retained and cannot win the comparison.
.msr_multistart <- function(objective, starts, lower, upper,
                            gradient_tolerance = 1e-3,
                            polish = TRUE,
                            control = list(iter.max = 500L, eval.max = 1000L,
                                           rel.tol = 1e-11, x.tol = 1e-9)) {
  if (!is.matrix(starts) || !is.numeric(starts) || !nrow(starts) ||
      ncol(starts) != length(lower) || length(lower) != length(upper) ||
      anyNA(lower) || anyNA(upper) || any(lower >= upper) ||
      any(!is.finite(starts)) ||
      any(sweep(starts, 2L, lower, "<") | sweep(starts, 2L, upper, ">")) ||
      length(gradient_tolerance) != 1L || !is.finite(gradient_tolerance) ||
      gradient_tolerance <= 0)
    stop("Invalid starts, bounds, or gradient tolerance.", call. = FALSE)
  fits <- vector("list", nrow(starts))
  records <- vector("list", nrow(starts))
  for (i in seq_len(nrow(starts))) {
    t0 <- proc.time()[["elapsed"]]
    result <- tryCatch({
      if (is.function(objective$reset_start)) objective$reset_start()
      fit <- stats::nlminb(unname(starts[i, ]), objective$fn, objective$gr,
                           lower = lower, upper = upper, control = control)
      g <- .msr_projected_gradient(fit$par, objective$gr(fit$par), lower, upper)
      primary <- fit
      primary_gradient <- max(abs(g))
      polished <- FALSE
      polish_message <- ""
      # A separate bounded optimizer can verify a near-stationary nlminb
      # termination. Keep its original diagnostics rather than relabeling it.
      if (polish && is.finite(fit$objective) &&
          (fit$convergence != 0L || primary_gradient > gradient_tolerance)) {
        attempt <- tryCatch(stats::optim(
          fit$par, objective$fn, objective$gr, method = "L-BFGS-B",
          lower = lower, upper = upper,
          control = list(maxit = 500L, factr = 1e7, pgtol = 1e-6)),
          error = function(e) list(convergence = 99L, message = conditionMessage(e)))
        polish_message <- paste(attempt$convergence, attempt$message)
        if (attempt$convergence == 0L && is.finite(attempt$value)) {
          pg <- .msr_projected_gradient(attempt$par, objective$gr(attempt$par),
                                        lower, upper)
          if (max(abs(pg)) <= gradient_tolerance &&
              attempt$value <= fit$objective + 1e-8) {
            fit <- list(par = attempt$par, objective = attempt$value,
                        convergence = attempt$convergence, message = attempt$message,
                        primary = primary)
            g <- pg
            polished <- TRUE
          }
        }
      }
      if (is.function(objective$cold_evaluate)) {
        cold <- objective$cold_evaluate(fit$par)
        cold_gap <- abs(cold$value - fit$objective)
        if (!is.finite(cold_gap) || cold_gap > 1e-5)
          stop("Independent cold coefficient audit disagrees with the warm likelihood.", call. = FALSE)
        fit$cold_gap <- cold_gap
        g <- .msr_projected_gradient(fit$par, cold$gradient, lower, upper)
      }
      list(fit = fit, gradient = max(abs(g)), error = "",
           primary = primary, primary_gradient = primary_gradient,
           polished = polished, polish_message = polish_message)
    }, error = function(e) list(fit = NULL, gradient = Inf,
                               error = conditionMessage(e)))
    fit <- result$fit
    fits[i] <- list(fit)
    eligible <- !is.null(fit) && fit$convergence == 0L &&
      is.finite(fit$objective) && result$gradient <= gradient_tolerance
    records[[i]] <- data.frame(
      start = i, value = if (is.null(fit)) NA_real_ else fit$objective,
      convergence = if (is.null(fit)) NA_integer_ else fit$convergence,
      projected_gradient = result$gradient, eligible = eligible,
      primary_value = if (is.null(result$primary)) NA_real_ else result$primary$objective,
      primary_convergence = if (is.null(result$primary)) NA_integer_ else result$primary$convergence,
      primary_gradient = if (is.null(result$primary_gradient)) Inf else result$primary_gradient,
      primary_message = if (is.null(result$primary)) result$error else result$primary$message,
      polished = isTRUE(result$polished),
      polish_message = if (is.null(result$polish_message)) "" else result$polish_message,
      seconds = proc.time()[["elapsed"]] - t0,
      message = if (is.null(fit)) result$error else fit$message)
    if (is.function(objective$cold_evaluate))
      records[[i]]$cold_gap <- if (is.null(fit$cold_gap)) NA_real_ else fit$cold_gap
  }
  records <- do.call(rbind, records)
  eligible <- which(records$eligible)
  winner <- if (length(eligible))
    eligible[which.min(records$value[eligible])] else NA_integer_
  list(best = if (is.na(winner)) NULL else fits[[winner]],
       winner = winner, records = records, fits = fits, starts = starts)
}

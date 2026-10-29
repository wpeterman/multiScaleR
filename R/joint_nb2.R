# The joint backend uses rho = log(1 + alpha), alpha = 1/theta. Its lower
# bound is the exact Poisson model. No large finite theta stands in for it.
.msr_analytic_backend <- function(object) {
  isTRUE(object$backend %in% c("profile_glm", "joint_nb2"))
}

# Cache the count-only part of the likelihood. Moderate counts use an exact
# rising factorial; large counts use a stable log-gamma ratio without allocating
# a vector whose length is the largest response.
.msr_nb2_normalizer <- function(y, weights, cache = FALSE) {
  small <- y <= 10000
  k <- if (any(small) && max(y[small]) > 0) seq_len(max(y[small])) - 1 else numeric()
  frequency <- numeric(length(k))
  if (length(k)) {
    counts <- numeric(length(k) + 1L)
    totals <- rowsum(weights[small], y[small], reorder = FALSE)
    counts[as.integer(rownames(totals)) + 1L] <- totals[, 1L]
    frequency <- rev(cumsum(rev(counts[-1L])))
  }
  yl <- y[!small]; wl <- weights[!small]
  evaluate <- function(alpha) {
    if (alpha == 0) return(c(value = 0, derivative = sum(weights * y * (y - 1) / 2)))
    value <- sum(frequency * log1p(k * alpha))
    derivative <- sum(frequency * k / (1 + k * alpha))
    if (length(yl)) {
      theta <- 1 / alpha
      if (theta < 100) {
        r <- lgamma(theta + yl) - lgamma(theta) - yl * log(theta)
        dr <- (digamma(theta) - digamma(theta + yl) + yl / theta) * theta^2
      } else {
        q <- alpha * yl; lq <- log1p(q)
        g <- ((1 + q) * lq - q) / q
        h <- (q - lq) / q^2
        s <- q < 1e-3; u <- q[s]
        g[s] <- u/2-u^2/6+u^3/12-u^4/20+u^5/30-u^6/42+u^7/56
        h[s] <- .5-u/3+u^2/4-u^3/5+u^4/6-u^5/7+u^6/8
        r <- yl * g - .5 * lq
        dr <- yl^2 * h - .5 * yl / (1 + q)
        # Stirling remainder, accurate well beyond double precision at theta 100.
        m <- c(1, 3, 5, 7, 9); cc <- c(1/12, -1/360, 1/1260, -1/1680, 1/1188)
        for (j in seq_along(m)) {
          r <- r + cc[j] * alpha^m[j] * expm1(-m[j] * lq)
          dr <- dr + cc[j] * m[j] * alpha^(m[j]-1) * expm1(-(m[j]+1) * lq)
        }
      }
      value <- value + sum(wl * r); derivative <- derivative + sum(wl * dr)
    }
    c(value = value, derivative = derivative)
  }
  if (!cache) return(evaluate)
  # One dispersion value only: bounded state, owned by this response/objective.
  # Conditional Newton iterations reuse it exactly without changing their path.
  last_alpha <- NULL; last_value <- NULL
  cached <- function(alpha) {
    if (!identical(alpha, last_alpha)) {
      value <- evaluate(alpha)
      last_alpha <<- alpha; last_value <<- value
    }
    last_value
  }
  attr(cached, "log_factorial") <- lgamma(y + 1)
  cached
}

.msr_nb2_terms <- function(y, eta, alpha, weights, normalizer) {
  if (!is.finite(alpha) || alpha < 0 || any(!is.finite(eta)))
    stop("Nonfinite predictor or negative NB2 dispersion.", call. = FALSE)
  if (alpha == 0) {
    mu <- exp(eta)
    value <- -sum(weights * stats::dpois(y, mu, log = TRUE))
    # exp(eta) can underflow even for a finite log mean. Keep the log-mean
    # likelihood for positive counts in that tail.
    tail <- mu == 0 & y > 0
    if (any(tail)) value <- -sum(weights[!tail] * stats::dpois(y[!tail], mu[!tail], log = TRUE)) +
      sum(weights[tail] * (lgamma(y[tail] + 1) - y[tail] * eta[tail]))
    score <- weights * (mu - y)
    da <- sum(weights * (y - (y - mu)^2)) / 2
    curvature <- weights * mu
  } else {
    lt <- log(alpha) + eta
    lp <- pmax(lt, 0) + log1p(exp(-abs(lt)))
    frac <- stats::plogis(lt); inv <- stats::plogis(-lt)
    eminus <- eta - lp
    high <- lt > 0
    eminus[high] <- -log(alpha) - log1p(exp(-lt[high]))
    # Division by alpha is safe for ordinary trials. A series avoids cancellation
    # in the dispersion score as alpha approaches zero.
    da_rows <- (y * frac + (frac - lp) / alpha) / alpha
    low <- lt < log(1e-4)
    u <- exp(lt[low]); mu_low <- exp(eta[low])
    da_rows[low] <- y[low] * mu_low / (1 + u) + mu_low^2 *
      (-.5 + 2*u/3 - 3*u^2/4 + 4*u^3/5 - 5*u^4/6 + 6*u^5/7)
    norm <- normalizer(alpha)
    log_factorial <- attr(normalizer, "log_factorial", exact = TRUE)
    if (is.null(log_factorial)) log_factorial <- lgamma(y + 1)
    value <- sum(weights * (log_factorial - y * eminus + lp / alpha)) - norm[[1L]]
    score <- weights * (frac / alpha - y * inv)
    da <- sum(weights * da_rows) - norm[[2L]]
    curvature <- weights * frac * inv * (1 / alpha + y)
    mu <- exp(eta)
  }
  list(value = unname(value), score = score, alpha_score = unname(da),
       curvature = curvature, mu = mu)
}

# Conditional coefficients are checked independently of the joint optimizer.
# Poisson IRLS is only an initializer; Newton steps solve the requested NB2
# likelihood, including small positive alpha that IRLS can round to Poisson.
.msr_nb2_coefficients <- function(X, y, offset, weights, alpha, normalizer) {
  family <- if (alpha <= 1e-8) stats::poisson() else MASS::negative.binomial(1 / alpha)
  fit <- stats::glm.fit(X, y, weights = weights, offset = offset, family = family,
    control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
  beta <- unname(fit$coefficients)
  if (any(!is.finite(beta)) || fit$rank != ncol(X))
    stop("Nonfinite or rank-deficient conditional NB2 fit.", call. = FALSE)
  for (iteration in seq_len(100L)) {
    a <- .msr_nb2_terms(y, drop(X %*% beta) + offset, alpha, weights, normalizer)
    score <- drop(crossprod(X, a$score))
    if (max(c(0, abs(score))) <= 1e-7) break
    step <- solve(crossprod(X, X * a$curvature), score)
    rate <- 1
    repeat {
      candidate <- beta - rate * step
      b <- .msr_nb2_terms(y, drop(X %*% candidate) + offset, alpha, weights, normalizer)
      if (is.finite(b$value) && b$value <= a$value - 1e-4 * rate * sum(score * step) + 1e-10) break
      rate <- rate / 2
      if (rate < 2^-30) stop("Conditional NB2 line search failed.", call. = FALSE)
    }
    beta <- candidate
  }
  a <- .msr_nb2_terms(y, drop(X %*% beta) + offset, alpha, weights, normalizer)
  if (max(c(0, abs(crossprod(X, a$score)))) > 1e-6)
    stop("Conditional NB2 coefficient score did not converge.", call. = FALSE)
  names(beta) <- colnames(X)
  fit$coefficients <- beta; fit$fitted.values <- a$mu
  fit$linear.predictors <- drop(X %*% beta) + offset
  fit$converged <- TRUE
  list(beta = beta, fit = fit, terms = a)
}

.msr_joint_nb2_objective <- function(engine, y, family = "nbinom2",
                                    X_fixed, offset, weights, trials = rep(1, length(y)),
                                    coefficient_start = "cold") {
  coefficient_start <- match.arg(coefficient_start, c("cold", "warm"))
  # Reuse validation, but never evaluate the older finite-theta objective here.
  obj <- .msr_profile_glm_objective(engine, y, family, X_fixed, offset, weights, trials)
  if (family != "nbinom2") stop("joint_nb2 requires an NB2 log-link model.", call. = FALSE)
  nk <- length(engine$parameter_names); nb <- ncol(X_fixed) + engine$p
  ib <- nk + seq_len(nb); outer <- c(seq_len(nk), nk + nb + 1L)
  normalizer <- .msr_nb2_normalizer(y, weights,
    cache = identical(engine$implementation, "stream"))
  worker <- if (coefficient_start == "warm") .msr_nb2_warm_coefficients() else NULL
  coefficients <- if (is.null(worker)) .msr_nb2_coefficients else worker$solve
  # Keep one solver closure, rather than serializing it again inside worker.
  if (!is.null(worker)) worker$solve <- NULL
  cache <- new.env(parent = emptyenv())
  cache$kernel <- cache$joint <- cache$profile <- NULL
  design <- function(kp) {
    if (!identical(kp, cache$kernel)) {
      cache$kd <- .msr_eval_kernel_engine(engine, kp)
      cache$X <- cbind(X_fixed, cache$kd$value); cache$kernel <- kp
    }
    list(X = cache$X, kd = cache$kd)
  }
  joint_eval <- function(par) {
    par <- unname(par)
    if (identical(par, cache$joint)) return(cache$joint_result)
    rho <- utils::tail(par, 1L)
    # Permit only floating-point roundoff outside the lower bound.
    if (!is.finite(rho) || rho < -1e-14) stop("Invalid NB2 rho.", call. = FALSE)
    rho <- max(0, rho); alpha <- expm1(rho)
    z <- design(par[seq_len(nk)]); beta <- par[ib]
    a <- .msr_nb2_terms(y, drop(z$X %*% beta) + offset, alpha, weights, normalizer)
    bk <- utils::tail(beta, engine$p)
    gk <- colSums(z$kd$d_scale * a$score) * bk
    if (engine$kernel == "expow") gk <- c(gk, colSums(z$kd$d_shape * a$score) * bk)
    result <- list(value = a$value,
      gradient = c(gk, drop(crossprod(z$X, a$score)), a$alpha_score * (1 + alpha)),
      beta = stats::setNames(beta, colnames(z$X)), design = z$X,
      theta = if (alpha == 0) Inf else 1 / alpha, alpha = alpha,
      fit = list(family = if (alpha == 0) stats::poisson() else MASS::negative.binomial(1/alpha),
                 fitted.values = a$mu), curvature = a$curvature)
    cache$joint <- par; cache$joint_result <- result
    result
  }
  profile_eval <- function(par) {
    par <- unname(par)
    if (identical(par, cache$profile)) return(cache$profile_result)
    rho <- utils::tail(par, 1L)
    if (!is.finite(rho) || rho < -1e-14) stop("Invalid NB2 rho.", call. = FALSE)
    z <- design(par[seq_len(nk)])
    fit <- coefficients(z$X, y, offset, weights, expm1(max(0, rho)), normalizer)
    jp <- c(par[seq_len(nk)], unname(fit$beta), max(0, rho))
    result <- joint_eval(jp)
    result$joint_par <- jp; result$gradient <- result$gradient[outer]
    cache$profile <- par; cache$profile_result <- result
    result
  }
  obj$parameter_names <- c(engine$parameter_names, "rho")
  obj$evaluate <- profile_eval; obj$fn <- function(par) profile_eval(par)$value
  obj$gr <- function(par) profile_eval(par)$gradient
  obj$joint <- list(evaluate = joint_eval, fn = function(par) joint_eval(par)$value,
    gr = function(par) joint_eval(par)$gradient, beta_index = ib, profile_index = outer)
  obj$boundary_nb2 <- TRUE
  if (!is.null(worker)) {
    reset_start <- function() {
      worker$reset()
      cache$profile <- NULL
    }
    obj$reset_start <- reset_start
    obj$cold_evaluate <- function(par) {
      reset_start()
      profile_eval(par)
    }
    obj$coefficient_audit <- function() list(counts = worker$state$counts,
      failures = worker$state$failures, cold_recoveries = worker$state$cold_recoveries)
    # The returned functions share this frame. The constructor's output list
    # is not needed by them and would retain a second copy of their bodies.
    on.exit(rm("obj", envir = environment()), add = TRUE)
  }
  obj$coefficient_start <- coefficient_start
  obj
}

.msr_joint_nb2_bounds <- function(contract, inputs, par, control, model) {
  b <- .msr_profile_starts(contract, inputs, par, control, model)
  nk <- length(b$lower) - 1L
  scales <- b$starts[, seq_len(nk), drop = FALSE]
  # Search each dispersed scale from both the exact boundary and an interior
  # dispersion. Retain the current six scale starts, giving twelve joint starts.
  b$starts <- rbind(cbind(scales, 0), cbind(scales, log1p(exp(-b$starts[, nk + 1L]))))
  b$lower[nk + 1L] <- 0; b$upper[nk + 1L] <- log1p(1 / control$theta_bounds[1L])
  colnames(b$starts) <- contract$objective$parameter_names
  b
}

.msr_joint_nb2_multistart <- function(obj, bounds, control) {
  if (any(!is.finite(bounds$starts)) ||
      any(sweep(bounds$starts, 2L, bounds$lower, "<") |
          sweep(bounds$starts, 2L, bounds$upper, ">")))
    stop("Invalid starts: kernel parameters must lie within the fitting bounds.", call. = FALSE)
  nk <- length(obj$engine$parameter_names); ib <- obj$joint$beta_index
  idx <- obj$joint$profile_index
  lower <- c(bounds$lower[seq_len(nk)], rep(-Inf, length(ib)), utils::tail(bounds$lower, 1L))
  upper <- c(bounds$upper[seq_len(nk)], rep(Inf, length(ib)), utils::tail(bounds$upper, 1L))
  fits <- vector("list", nrow(bounds$starts)); records <- vector("list", length(fits))
  for (i in seq_along(fits)) {
    t0 <- proc.time()[["elapsed"]]; warnings <- character()
    ans <- tryCatch(withCallingHandlers({
      if (is.function(obj$reset_start)) obj$reset_start()
      start <- obj$evaluate(bounds$starts[i, ])$joint_par
      # L-BFGS-B needs finite trial values. Reject unrepresentable trials with a
      # finite barrier; a barrier point cannot pass final independent checks.
      safe <- function(p) {
        z <- tryCatch(obj$joint$evaluate(p), error = function(e) NULL)
        if (is.null(z) || !is.finite(z$value) || any(!is.finite(z$gradient)))
          list(value = 1e100, gradient = rep(0, length(p))) else z
      }
      fit <- stats::optim(start, function(p) safe(p)$value, function(p) safe(p)$gradient,
        method = "L-BFGS-B", lower = lower, upper = upper,
        control = list(maxit = control$iter.max, factr = 1e4, pgtol = 1e-7))
      exact <- obj$joint$evaluate(fit$par)
      pg <- .msr_projected_gradient(fit$par, exact$gradient, lower, upper)
      independent <- if (is.function(obj$cold_evaluate)) obj$cold_evaluate(fit$par[idx]) else
        obj$evaluate(fit$par[idx])
      gap <- abs(exact$value - independent$value)
      list(par = unname(fit$par[idx]), joint_par = fit$par,
        objective = exact$value, convergence = fit$convergence,
        message = if (is.null(fit$message)) "" else fit$message,
        gradient = max(abs(pg)), conditional_gap = gap,
        eligible = fit$convergence == 0L && is.finite(exact$value) &&
          all(is.finite(pg)) && max(abs(pg)) <= control$gradient_tolerance && gap <= 1e-5)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) list(par = bounds$starts[i, ], objective = NA_real_,
      convergence = NA_integer_, message = conditionMessage(e), gradient = Inf,
      conditional_gap = NA_real_, eligible = FALSE))
    fits[[i]] <- ans
    records[[i]] <- data.frame(start = i, value = ans$objective, convergence = ans$convergence,
      projected_gradient = ans$gradient, conditional_gap = ans$conditional_gap,
      eligible = ans$eligible, seconds = proc.time()[["elapsed"]] - t0,
      message = ans$message, warnings = paste(unique(warnings), collapse = "; "))
  }
  records <- do.call(rbind, records); good <- which(records$eligible)
  winner <- if (length(good)) good[which.min(records$value[good])] else NA_integer_
  list(best = if (is.na(winner)) NULL else fits[[winner]], winner = winner,
       fits = fits, records = records, starts = bounds$starts)
}

# A boundary maximum has no regular two-sided dispersion covariance. Keep
# joint Wald inference unavailable, while exposing coefficient uncertainty
# conditional on both the selected kernel and the exact Poisson face.
.msr_poisson_boundary_information <- function(obj, par, fits) {
  fitted <- obj$evaluate(par)
  nk <- length(obj$engine$parameter_names); nb <- length(fitted$beta)
  jp <- fitted$joint_par
  names(jp) <- c(obj$engine$parameter_names, paste0("beta:", names(fitted$beta)), "rho")
  H <- matrix(NA_real_, length(jp), length(jp), dimnames = list(names(jp), names(jp)))
  ib <- nk + seq_len(nb)
  conditional <- .msr_information_inverse(crossprod(fitted$design,
    fitted$design * fitted$curvature))$vcov
  good <- which(fits$records$eligible & fits$records$value <=
    fits$best$objective + stats::qchisq(.95, 1) / 2)
  multiple <- any(vapply(good, function(i)
    max(abs(fits$fits[[i]]$par[seq_len(nk)] - par[seq_len(nk)])) > .5, logical(1)))
  list(available = FALSE, vcov = H, hessian = H, par = jp,
    beta_index = ib, profile_index = c(seq_len(nk), length(jp)),
    nuisance_vcov = NULL, conditional_vcov = conditional,
    boundary = TRUE, weak = TRUE, multiple_modes = multiple,
    poisson_limit = TRUE, poisson_compatible = TRUE, poisson_deviance = 0,
    condition_number = Inf, log_scale_se = stats::setNames(rep(NA_real_, nk), obj$engine$parameter_names),
    reasons = c("exact Poisson dispersion boundary",
      "Poisson limit within the diagnostic likelihood threshold",
      if (multiple) "competing scale modes within the 95% likelihood threshold"))
}

.msr_joint_nb2_zero_effect <- function(setup) {
  obj <- setup$obj
  normalizer <- .msr_nb2_normalizer(obj$y, obj$weights)
  evaluate <- function(rho) {
    alpha <- expm1(max(0, rho))
    a <- .msr_nb2_coefficients(obj$X_fixed, obj$y, obj$offset, obj$weights,
                              alpha, normalizer)$terms
    list(value = a$value, gradient = a$alpha_score * (1 + alpha))
  }
  objective <- list(fn = function(rho) evaluate(rho)$value,
                    gr = function(rho) evaluate(rho)$gradient)
  fits <- .msr_multistart(objective,
    unique(setup$starts[, ncol(setup$starts), drop = FALSE]), 0, utils::tail(setup$upper, 1L),
    gradient_tolerance = setup$control$gradient_tolerance,
    control = list(iter.max = setup$control$iter.max, eval.max = setup$control$eval.max))
  if (is.null(fits$best)) stop("No zero-effect fit passed convergence checks.", call. = FALSE)
  fits$best$objective
}

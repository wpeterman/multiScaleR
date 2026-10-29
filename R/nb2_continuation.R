# Opt-in coefficient continuation, ported from the validated development runner.
# Each objective owns its state; every independent start and cold audit resets it.
.msr_nb2_warm_coefficients <- function(maxit = 20L) {
  cold <- .msr_nb2_cold_fallback
  terms <- .msr_nb2_terms
  state <- new.env(parent = emptyenv())
  state$beta <- NULL
  state$counts <- c(calls = 0L, cold = 0L, warm = 0L, fallback = 0L,
    iterations = 0L, resets = 0L)
  state$failures <- integer()
  state$cold_recoveries <- list()
  reset <- function() {
    state$beta <- NULL
    state$counts['resets'] <- state$counts['resets'] + 1L
  }
  solve <- function(X, y, offset, weights, alpha, normalizer) {
    state$counts['calls'] <- state$counts['calls'] + 1L
    finish_cold <- function() {
      state$counts['cold'] <- state$counts['cold'] + 1L
      # Clear stale coefficients even if the cold fallback itself fails.
      state$beta <- NULL
      z <- cold(X, y, offset, weights, alpha, normalizer)
      if (!is.null(z$initialization_audit))
        state$cold_recoveries[[length(state$cold_recoveries) + 1L]] <- z$initialization_audit
      state$beta <- unname(z$beta)
      z
    }
    if (is.null(state$beta)) return(finish_cold())
    attempt <- tryCatch({
      beta <- state$beta
      if (length(beta) != ncol(X) || any(!is.finite(beta)))
        stop('invalid warm coefficients')
      a <- terms(y, drop(X %*% beta) + offset, alpha, weights, normalizer)
      initial <- a$value
      for (iteration in seq_len(maxit)) {
        score <- drop(crossprod(X, a$score))
        if (any(!is.finite(score)) || !is.finite(a$value)) stop('nonfinite warm score')
        if (max(c(0, abs(score))) <= 1e-7) break
        state$counts['iterations'] <- state$counts['iterations'] + 1L
        step <- base::solve(crossprod(X, X * a$curvature), score)
        descent <- sum(score * step)
        if (!is.finite(descent) || descent <= 0) stop('invalid Newton direction')
        rate <- 1
        repeat {
          candidate <- beta - rate * step
          b <- terms(y, drop(X %*% candidate) + offset, alpha, weights, normalizer)
          if (is.finite(b$value) && b$value <= a$value - 1e-4 * rate * descent + 1e-10) break
          rate <- rate / 2
          if (rate < 2^-30) stop('warm line search failed')
        }
        beta <- candidate; a <- b
      }
      score <- drop(crossprod(X, a$score))
      if (any(!is.finite(score)) || max(c(0, abs(score))) > 1e-7)
        stop('warm score did not converge')
      if (!is.finite(a$value) || a$value > initial + 1e-8 || any(!is.finite(a$mu)))
        stop('invalid warm likelihood')
      # Match glm.fit's QR tolerance for epsilon=1e-11. Check its expected
      # information weights, including its Poisson initializer for tiny alpha.
      eta <- drop(X %*% beta) + offset
      logw <- eta
      if (alpha > 1e-8) {
        t <- log(alpha) + eta
        logw <- eta - (pmax(t, 0) + log1p(exp(-abs(t))))
      }
      wx <- X * sqrt(weights * exp(logw))
      if (any(!is.finite(wx)) || qr(wx, tol = 1e-14, LAPACK = FALSE)$rank != ncol(X))
        stop('rank-deficient warm design')
      names(beta) <- colnames(X)
      list(beta = beta, terms = a)
    }, error = function(e) e)
    if (inherits(attempt, 'error')) {
      reason <- conditionMessage(attempt)
      if (!reason %in% names(state$failures)) state$failures[reason] <- 0L
      state$failures[reason] <- state$failures[reason] + 1L
      state$counts['fallback'] <- state$counts['fallback'] + 1L
      return(finish_cold())
    }
    state$counts['warm'] <- state$counts['warm'] + 1L
    state$beta <- unname(attempt$beta)
    attempt
  }
  # The returned solver does not recurse. Its enclosing frame need not keep
  # a second reference to its own body when fitted objectives are serialized.
  on.exit(rm("solve", envir = environment()), add = TRUE)
  list(solve = solve, reset = reset, state = state)
}


.msr_nb2_cold_fallback <- function(X, y, offset, weights, alpha, normalizer) {
  original <- tryCatch(.msr_nb2_coefficients(
    X, y, offset, weights, alpha, normalizer), error = identity)
  if (!inherits(original, 'error')) return(original)
  terms <- .msr_nb2_terms
  if (any(!is.finite(X)) || any(!is.finite(offset)) ||
      any(!is.finite(weights)) || any(weights < 0) ||
      qr(X * sqrt(weights), tol = 1e-14)$rank != ncol(X))
    stop('Cold fallback: invalid or rank-deficient design')
  # The projection accommodates offsets and models without an intercept.
  beta <- unname(stats::lm.wfit(X, log(y + .5) - offset, weights)$coefficients)
  a <- terms(y, drop(X %*% beta) + offset, alpha, weights, normalizer)
  initial <- a$value
  for (iteration in seq_len(200L)) {
    score <- drop(crossprod(X, a$score))
    if (any(!is.finite(score)) || !is.finite(a$value)) stop('Cold fallback: nonfinite score')
    if (max(c(0, abs(score))) <= 1e-7) break
    step <- solve(crossprod(X, X * a$curvature), score)
    descent <- sum(score * step)
    if (!is.finite(descent) || descent <= 0) stop('Cold fallback: invalid direction')
    rate <- 1
    repeat {
      candidate <- beta - rate * step
      b <- terms(y, drop(X %*% candidate) + offset, alpha, weights, normalizer)
      if (is.finite(b$value) && b$value <= a$value - 1e-4 * rate * descent + 1e-10) break
      rate <- rate / 2
      if (rate < 2^-30) stop('Cold fallback: failed line search')
    }
    beta <- candidate; a <- b
  }
  score <- max(c(0, abs(crossprod(X, a$score))))
  eta <- drop(X %*% beta) + offset
  wx <- X * sqrt(weights * if (alpha == 0) a$mu else stats::plogis(log(alpha) + eta) / alpha)
  direct <- -sum(weights * if (alpha == 0) stats::dpois(y, a$mu, log = TRUE) else
    stats::dnbinom(y, mu = a$mu, size = 1/alpha, log = TRUE))
  if (!is.finite(score) || score > 1e-7 || any(!is.finite(a$mu)) ||
      any(!is.finite(wx)) || qr(wx, tol = 1e-14)$rank != ncol(X) ||
      !is.finite(direct) || abs(direct - a$value) > 1e-5 || a$value > initial + 1e-8)
    stop('Cold fallback: terminal audit failed')
  names(beta) <- colnames(X)
  list(beta = beta, terms = a, fit = list(coefficients = beta,
    fitted.values = a$mu, linear.predictors = eta, converged = TRUE,
    rank = ncol(X), family = if (alpha == 0) stats::poisson() else MASS::negative.binomial(1/alpha)),
    initialization_audit = list(reference_error = conditionMessage(original),
      iterations = iteration, score = score, direct_gap = abs(direct - a$value)))
}


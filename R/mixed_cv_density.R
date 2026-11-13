# Evaluate the response log density and its first two derivatives on the link
# scale. The binomial and NB2 expressions retain informative tails when a
# response-scale probability has rounded to 0 or 1.
.msr_mixed_cv_response <- function(family, dispersion, observed, trials, eta) {
  softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
  if (family == "binomial") {
    p <- stats::plogis(eta)
    q <- stats::plogis(-eta)
    ll <- lchoose(trials, observed) - observed * softplus(-eta) -
      (trials - observed) * softplus(eta)
    return(c(log_density = sum(ll), gradient = sum(observed * q - (trials - observed) * p),
      information = sum(trials * p * q)))
  }
  if (family == "nbinom2") {
    theta <- dispersion
    delta <- eta - log(theta)
    p <- stats::plogis(delta)
    q <- stats::plogis(-delta)
    ll <- -lbeta(theta, observed + 1) - log(theta + observed) -
      observed * softplus(-delta) - theta * softplus(delta)
    return(c(log_density = sum(ll), gradient = sum(observed * q - theta * p),
      information = sum((observed + theta) * p * q)))
  }
  if (family == "Gamma") {
    shape <- 1 / dispersion^2
    scaled <- exp(log(observed) - eta)
    ll <- shape * log(shape) - lgamma(shape) + (shape - 1) * log(observed) -
      shape * (eta + scaled)
    return(c(log_density = sum(ll), gradient = shape * sum(scaled - 1),
      information = shape * sum(scaled)))
  }
  mu <- exp(eta)
  c(log_density = sum(observed * eta - mu - lgamma(observed + 1)),
    gradient = sum(observed - mu), information = sum(mu))
}

.msr_mixed_cv_log_integral <- function(eta, y, n, family, dispersion, sd, tolerance) {
  if (!length(y)) return(0)
  if (sd == 0) return(.msr_mixed_cv_response(family, dispersion, y, n, eta)[["log_density"]])
  state <- function(z) .msr_mixed_cv_response(family, dispersion, y, n, eta + sd * z)
  score <- function(z) {
    value <- sd * state(z)[["gradient"]] - z
    max(-.Machine$double.xmax, min(.Machine$double.xmax, value))
  }
  # Each supported response is log-concave in the single shared effect.
  # Bracket its unique mode, then center and scale adaptive integration there.
  left <- right <- 0
  if (score(0) < 0) {
    left <- -1
    while (score(left) < 0 && is.finite(2 * left)) left <- 2 * left
  } else {
    right <- 1
    while (score(right) > 0 && is.finite(2 * right)) right <- 2 * right
  }
  mode <- stats::uniroot(score, c(left, right), tol = 1e-10)$root
  at_mode <- state(mode)
  width <- 1 / sqrt(1 + sd^2 * at_mode[["information"]])
  height <- stats::dnorm(mode, log = TRUE) + at_mode[["log_density"]]
  if (!is.finite(height) || !is.finite(width) || width <= 0)
    stop("A mixed-model group has a nonrepresentable conditional density.", call. = FALSE)
  integrand <- function(x) vapply(mode + width * x, function(z) {
    exp(stats::dnorm(z, log = TRUE) + state(z)[["log_density"]] - height)
  }, numeric(1))
  integral <- stats::integrate(integrand, -Inf, Inf, rel.tol = tolerance,
    abs.tol = 0, subdivisions = 200L)
  answer <- height + log(width) + log(integral$value)
  if (!is.finite(answer))
    stop("Mixed-model group integration returned a nonfinite log density.", call. = FALSE)
  answer
}

.msr_mixed_cv_glmm_group <- function(train, test, family, dispersion, sd, tolerance) {
  integral <- function(x, shift = 0) .msr_mixed_cv_log_integral(
    x$eta + shift, x$y, x$n, family, dispersion, sd, tolerance)
  train_density <- integral(train)
  joint <- list(eta = c(train$eta, test$eta), y = c(train$y, test$y), n = c(train$n, test$n))
  score <- train_density - integral(joint)
  prediction <- if (family == "binomial") {
    vapply(test$eta, function(eta) {
      probe <- function(y) integral(list(eta = c(train$eta, eta),
        y = c(train$y, y), n = c(train$n, 1)))
      stats::plogis(probe(1) - probe(0))
    }, numeric(1))
  } else {
    # Exponential tilting of the Normal prior evaluates E[exp(b)|training].
    exp(test$eta + sd^2 / 2 + integral(train, sd^2) - train_density)
  }
  list(prediction = prediction, neg_log_score = score)
}

.msr_mixed_cv_gaussian_group <- function(train, test, covariance, sd) {
  # Whiten the random coefficients. Only coefficient-sized matrices are
  # factored, even when a group contains many held-out observations.
  eigen <- eigen(covariance, symmetric = TRUE)
  tolerance <- 64 * .Machine$double.eps * ncol(covariance) * max(abs(eigen$values))
  if (any(eigen$values < -tolerance))
    stop("Native random-effect covariance is not positive semidefinite.", call. = FALSE)
  L <- sweep(eigen$vectors, 2L, sqrt(pmax(eigen$values, 0)), "*")
  A <- train$Z %*% L / sd
  B <- test$Z %*% L
  posterior <- chol(diag(ncol(L)) + crossprod(A))
  mean <- backsolve(posterior, forwardsolve(t(posterior), crossprod(A, (train$y - train$eta) / sd)))
  prediction <- test$eta + drop(B %*% mean)
  loading <- t(forwardsolve(t(posterior), t(B))) / sd
  innovation <- (test$y - prediction) / sd
  updated <- chol(diag(ncol(L)) + crossprod(loading))
  mode <- backsolve(updated, forwardsolve(t(updated), crossprod(loading, innovation)))
  quadratic <- sum((innovation - drop(loading %*% mode))^2) + sum(mode^2)
  logdet <- 2 * (length(innovation) * log(sd) + sum(log(diag(updated))))
  list(prediction = prediction,
    neg_log_score = (length(innovation) * log(2 * pi) + logdet + quadratic) / 2)
}

wald_plan <- function(replicates) {
  plan <- expand.grid(family = c("gaussian", "poisson", "binomial", "nbinom2"),
    scenario = c("interior", "weak", "near_upper", "exact_upper", "exact_lower", "poisson_face"),
    n = c(75L, 240L), replicate = seq_len(replicates),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  plan <- plan[plan$scenario != "poisson_face" | plan$family == "nbinom2", ]
  plan$seed <- 1621000L + (plan$replicate - 1L) * 1000L +
    match(plan$scenario, c("interior", "weak", "near_upper", "exact_upper", "exact_lower", "poisson_face")) * 100L +
    match(plan$n, c(75L, 240L)) * 10L + match(plan$family, c("gaussian", "poisson", "binomial", "nbinom2"))
  plan$id <- with(plan, sprintf("%s-%s-%03d-%03d", family, scenario, n, replicate))
  rownames(plan) <- NULL
  plan
}

wald_case <- function(case) {
  set.seed(case$seed)
  n <- case$n
  d <- seq(.015, 1, length.out = 23L)
  count <- matrix(sample(1:5, n * length(d), replace = TRUE), n)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(12 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .25), n)
  value <- radial * count
  habitat <- function(log_scale) {
    weight <- exp(-(d^2 - min(d)^2) / (2 * exp(2 * log_scale)))
    predictor <- drop(value %*% weight) / drop(count %*% weight)
    (predictor - mean(predictor)) / stats::sd(predictor)
  }
  sigma <- switch(case$scenario, near_upper = .73, exact_upper = .8, exact_lower = .05, .24)
  effect <- if (case$scenario == "weak") .12 else .95
  theta <- if (case$scenario == "poisson_face") Inf else 2
  data <- data.frame(hab = habitat(log(.23)), z = rnorm(n), effort = runif(n, .8, 1.2))
  offset <- log(data$effort)
  eta <- .4 + effect * habitat(log(sigma)) + .3 * data$z + offset
  trials <- rep(5, n)
  data$y <- switch(case$family, gaussian = rnorm(n, eta, .65), poisson = rpois(n, exp(eta)),
    binomial = rbinom(n, trials, plogis(eta)),
    nbinom2 = if (is.finite(theta)) rnbinom(n, mu = exp(eta), size = theta) else rpois(n, exp(eta)))
  data$trials <- trials
  inputs <- list(raw_cov = NULL, d_list = NULL, min_D = 50, max_D = 800,
    unit_conv = 1000, kernel = "gaussian",
    binned = list(bin_rep = d, vsum = list(hab = value), csum = list(hab = count),
      point_ids = rownames(data), n_points = n))
  formula <- if (case$family == "binomial") cbind(y, trials - y) ~ hab + z + offset(log(effort)) else
    y ~ hab + z + offset(log(effort))
  list(data = data, inputs = inputs, formula = formula, habitat = habitat,
    offset = offset, truth = list(log_scale = log(sigma), effect = effect, theta = theta))
}

wald_capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = unique(messages))
}

wald_glm <- function(X, y, offset, family, trials) {
  family_object <- switch(family, gaussian = stats::gaussian(), poisson = stats::poisson(),
    binomial = stats::binomial())
  weights <- if (family == "binomial") trials else rep(1, length(y))
  fit <- stats::glm.fit(X, y / weights, weights = weights, offset = offset,
    family = family_object, control = stats::glm.control(epsilon = 1e-11, maxit = 100L))
  if (!isTRUE(fit$converged) || any(!is.finite(fit$coefficients)))
    stop("Reference GLM did not converge with finite coefficients.")
  mu <- fit$fitted.values
  density <- switch(family, gaussian = stats::dnorm(y, mu, sqrt(mean((y - mu)^2)), log = TRUE),
    poisson = stats::dpois(y, mu, log = TRUE), binomial = stats::dbinom(y, trials, mu, log = TRUE))
  if (any(!is.finite(density))) stop("Reference GLM has a nonfinite density.")
  list(nll = -sum(density), beta = fit$coefficients)
}

# Profile the native normalized NB2 law in alpha = 1/theta. Use a shifted
# log coordinate that includes alpha = 0 literally and allows optimization
# close to either package's dispersion boundary. Coefficients use independent
# Fisher scoring and density-based backtracking, initialized from Poisson.
wald_nb_reference <- function(X, y, offset, poisson, joint) {
  fixed_alpha <- function(alpha) {
    if (alpha == 0) return(poisson$nll)
    beta <- poisson$beta
    eta <- drop(X %*% beta) + offset
    size <- 1 / alpha
    objective <- function(link) -sum(stats::dnbinom(y, mu = exp(link), size = size, log = TRUE))
    old <- objective(eta)
    for (iteration in seq_len(80L)) {
      mu <- exp(eta)
      inv <- 1 / (1 + alpha * mu)
      score <- drop(crossprod(X, (y - mu) * inv))
      if (max(abs(score)) <= 5e-8) return(old)
      direction <- drop(solve(crossprod(X, X * (mu * inv)), score))
      adjustment <- drop(X %*% direction)
      step <- 1
      repeat {
        new <- objective(eta + step * adjustment)
        if (is.finite(new) && new <= old + 64 * .Machine$double.eps * max(1, abs(old))) break
        step <- step / 2
        if (step < 2^-24) stop("Independent NB2 coefficient line search failed.")
      }
      change <- max(abs(step * adjustment))
      eta <- eta + step * adjustment
      old <- new
      if (change <= 5e-10) return(old)
    }
    stop("Independent NB2 coefficient score did not converge.")
  }
  shift <- .01
  alpha_min <- if (joint) 0 else 1e-5
  alpha_max <- 1e4
  nodes <- log(c(alpha_min, .001, .02, .15, .7, 3, 20, 300, alpha_max) + shift)
  cached <- new.env(parent = emptyenv())
  evaluate <- function(x) {
    key <- sprintf("%.17g", x)
    if (!exists(key, cached, inherits = FALSE)) {
      alpha <- if (x == nodes[1L]) alpha_min else if (x == tail(nodes, 1L)) alpha_max else exp(x) - shift
      assign(key, fixed_alpha(alpha), cached)
    }
    get(key, cached, inherits = FALSE)
  }
  values <- vapply(nodes, evaluate, 0)
  k <- length(nodes)
  middle <- seq_len(k - 2L) + 1L
  local <- middle[values[middle] <= values[middle - 1L] & values[middle] <= values[middle + 1L]]
  # Examine both end intervals because a narrow interior minimum can occur
  # beyond the first or last dispersion node.
  ranges <- unique(rbind(c(1L, 2L), c(k - 1L, k), cbind(local - 1L, local + 1L)))
  fits <- apply(ranges, 1L, function(i)
    stats::optimize(evaluate, nodes[i], tol = 1e-8))
  fitted <- c(values, vapply(fits, function(x) x$objective, 0))
  parameters <- c(nodes, vapply(fits, function(x) x$minimum, 0))
  best <- which.min(fitted)
  list(nll = fitted[best], alpha = if (best == 1L) alpha_min else exp(parameters[best]) - shift,
    poisson_nll = poisson$nll)
}

wald_reference <- function(case, family, backend) {
  cache <- new.env(parent = emptyenv())
  y <- case$data$y
  trials <- case$data$trials
  offset <- case$offset
  poisson_fit <- family == "nbinom2"
  evaluate <- function(log_scale) {
    key <- sprintf("%.17g", log_scale)
    if (!exists(key, cache, inherits = FALSE)) {
      evaluated <- wald_capture({
        X <- cbind(1, case$habitat(log_scale), case$data$z)
        ordinary <- wald_glm(X, y, offset, if (poisson_fit) "poisson" else family, trials)
        if (poisson_fit) wald_nb_reference(X, y, offset, ordinary, backend == "joint_nb2") else ordinary
      })
      assign(key, list(log_scale = log_scale, result = evaluated), cache)
    }
    a <- get(key, cache, inherits = FALSE)$result
    if (inherits(a$value, "error")) stop(conditionMessage(a$value), call. = FALSE)
    a$value$nll
  }
  list(fn = evaluate, audit = function() {
    points <- as.list(cache, all.names = TRUE)
    lapply(points, function(x) data.frame(log_scale = x$log_scale,
      nll = if (inherits(x$result$value, "error")) NA_real_ else x$result$value$nll,
      alpha = if (!inherits(x$result$value, "error") && length(x$result$value$alpha)) x$result$value$alpha else NA_real_,
      warnings = paste(x$result$warnings, collapse = " | "),
      error = if (inherits(x$result$value, "error")) conditionMessage(x$result$value) else ""))
  })
}

wald_global <- function(fn, grid, nll) {
  k <- length(grid)
  middle <- seq_len(k - 2L) + 1L
  local <- middle[nll[middle] < nll[middle - 1L] & nll[middle] < nll[middle + 1L]]
  # Refine boundary-adjacent intervals in addition to all detected minima.
  limits <- unique(rbind(c(1L, 2L), c(k - 1L, k), cbind(local - 1L, local + 1L)))
  refined <- apply(limits, 1L, function(i) stats::optimize(fn, grid[i], tol = 1e-9))
  candidates <- c(nll, vapply(refined, function(x) x$objective, 0))
  points <- c(grid, vapply(refined, function(x) x$minimum, 0))
  winner <- which.min(candidates)
  list(best = candidates[winner], optimum = points[winner], minima = points[-seq_along(grid)])
}

wald_profile_set <- function(fn, grid, nll, best, level) {
  cutoff <- best + stats::qchisq(level, 1) / 2
  accepted <- nll <= cutoff
  edge <- diff(c(FALSE, accepted, FALSE))
  first <- which(edge == 1L)
  last <- which(edge == -1L) - 1L
  if (!length(first)) stop("Profile set has no accepted component.")
  root <- function(i, j) stats::uniroot(function(x) fn(x) - cutoff, grid[c(i, j)], tol = 1e-8)$root
  lower <- vapply(first, function(i) if (i == 1L) grid[i] else root(i - 1L, i), 0)
  upper <- vapply(last, function(i) if (i == length(grid)) grid[i] else root(i, i + 1L), 0)
  data.frame(lower = exp(lower) * 1000, upper = exp(upper) * 1000,
    lower_bound = first == 1L, upper_bound = last == length(grid))
}

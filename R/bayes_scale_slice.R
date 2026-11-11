# A sweep uses exact Gaussian updates where available and slice updates for
# the remaining full conditional distributions. The likelihood is unchanged.
.msr_bayes_slice <- function(out, chains, iter, warmup, seed, cores = 1L) {
  spec <- out$specification
  sizes <- lengths(spec$parameters)
  variables <- unlist(Map(function(nm, n)
    if (n) paste0(nm, "[", seq_len(n), "]") else character(),
    names(sizes), sizes), use.names = FALSE)
  saved <- iter - warmup
  cores <- min(cores, chains)
  kinds <- RNGkind()
  # Seed by chain ID before dispatch so worker scheduling cannot change draws.
  seeds <- .msr_bayes_prediction_seed(seed, sample.int(.Machine$integer.max, chains))
  pool <- .msr_bayes_pool(cores, .msr_bayes_slice_context, spec)
  if (is.null(pool)) {
    matrices <- lapply(seeds, .msr_bayes_slice_run, .msr_bayes_slice_context(spec),
      iter, warmup, kinds)
  } else {
    on.exit(.msr_bayes_pool_stop(pool), add = TRUE)
    matrices <- .msr_bayes_pool_map(pool, as.list(seeds), .msr_bayes_slice_run,
      iter, warmup, kinds)
  }
  a <- array(NA_real_, c(saved, chains, sum(sizes)), dimnames = list(NULL, NULL, variables))
  for (chain in seq_len(chains)) a[, chain, ] <- matrices[[chain]]
  out$chains <- a
  out$draws <- .msr_bayes_parameter_blocks(matrix(a, nrow = saved * chains), spec$parameters)
  out$inference <- "full Bayesian slice sampling"
  out$sampling <- list(chains = chains, iter = iter, warmup = warmup, seed = seed, cores = cores)
  out$inference_metadata <- .msr_scale_inference_metadata("posterior", "slice")
  status <- .msr_bayes_package_status("posterior", "1.6.0")
  if (status$passed) {
    tab <- .msr_bayes_chain_summary(a)
    passed <- .msr_bayes_chain_checks(tab, chains)
    out$diagnostics <- list(passed = passed, parameters = tab,
      criteria = .msr_bayes_chain_criteria())
    if (!passed)
      warning("Bayesian sampling diagnostics need attention; inspect summary(fit)$diagnostics before interpreting posterior summaries.",
        call. = FALSE)
  } else {
    out$diagnostics <- list(passed = NA, parameters = NULL,
      reason = "Install posterior (>= 1.6.0) to compute chain diagnostics.")
  }
  out
}

.msr_bayes_slice_context <- function(spec) {
  list(spec = spec, shared = .msr_bayes_slice_shared_bins(spec))
}

.msr_bayes_slice_run <- function(seed, context, iter, warmup, kinds) {
  .msr_bayes_prediction_seed(seed, {
    set.seed(seed, kind = kinds[1L], normal.kind = kinds[2L], sample.kind = kinds[3L])
    .msr_bayes_slice_chain(context$spec, iter, warmup, context$shared)
  })
}

# With a shared scale, observations whose nearest populated bin is the same
# have the same shifted weights. Compute each such vector once per proposal.
.msr_bayes_slice_shared_bins <- function(spec) {
  lapply(spec$covariates, function(nm) {
    sm <- spec$scales[[nm]]
    if (is.null(spec$C[[nm]]) || length(sm$random) || ncol(sm$X) != 1L ||
        any(sm$X[, 1L] != sm$X[1L, 1L])) return(NULL)
    groups <- split(seq_along(spec$y), max.col(spec$C[[nm]] > 0, ties.method = "first"))
    lapply(groups, function(rows) list(rows = rows,
      V = spec$V[[nm]][rows, , drop = FALSE], C = spec$C[[nm]][rows, , drop = FALSE],
      delta = spec$delta[[nm]][rows[1L], ]))
  })
}

.msr_bayes_slice_kernel <- function(spec, par, j, shared) {
  nm <- spec$covariates[j]
  log_s <- .msr_bayes_log_scale(spec$scales[[nm]], par)
  if (anyNA(log_s)) return(rep(NA_real_, length(spec$y)))
  if (is.null(shared[[j]])) return(.msr_bayes_kernel_value(spec, nm, log_s))
  multiplier <- exp(-(if (spec$kernel == "gaussian") 2 else 1) * log_s[1L])
  if (spec$kernel == "gaussian") multiplier <- multiplier / 2
  raw <- numeric(length(spec$y))
  for (group in shared[[j]]) {
    w <- if (is.infinite(multiplier)) as.numeric(group$delta == 0) else
      exp(-group$delta * multiplier)
    raw[group$rows] <- drop(group$V %*% w) / drop(group$C %*% w)
  }
  (raw - spec$center[[nm]]) / spec$spread[[nm]]
}

.msr_bayes_slice_chain <- function(spec, iter, warmup, shared) {
  priors <- spec$priors
  loglik <- function(eta, log_dispersion = par$log_dispersion) {
    if (any(!is.finite(eta))) return(-Inf)
    if (length(log_dispersion) && (!is.finite(exp(log_dispersion)) ||
        exp(log_dispersion) == 0)) return(-Inf)
    if (spec$family %in% c("poisson", "nbinom2") && any(!is.finite(exp(eta)))) return(-Inf)
    value <- sum(.msr_scale_response_loglik(spec$family, spec$y, eta, log_dispersion, spec$trials))
    if (is.finite(value)) value else -Inf
  }
  response_effect <- function(term, z, sd) {
    if (is.null(term$gaussian_block)) exp(sd) * z[term$index] else
      .msr_gaussian_effect(term$gaussian_block, z, sd)
  }
  initialize <- function(p) {
    k <- stats::setNames(lapply(seq_along(spec$covariates),
      function(j) .msr_bayes_slice_kernel(spec, p, j, shared)), spec$covariates)
    X <- .msr_bayes_slice_design(spec, k)
    offset <- spec$offset
    for (term in spec$response$random)
      offset <- offset + response_effect(term, p$response_z[term$z_index],
        p$log_response_sd[term$sd_index])
    eta <- if (is.null(X)) NULL else offset + drop(X %*% p$beta)
    list(kernel = k, X = X, offset = offset, eta = eta,
      ll = if (is.null(eta)) -Inf else loglik(eta, p$log_dispersion))
  }
  # Modest random perturbations give chains distinct starts. If a perturbation
  # leaves the model's numerical domain, use the validated reference point.
  par <- lapply(spec$parameters, function(x) x + stats::rnorm(length(x), sd = 0.1))
  state <- initialize(par)
  if (!is.finite(state$ll)) {
    par <- spec$parameters
    state <- initialize(par)
  }
  if (!is.finite(state$ll))
    stop("Slice sampling needs a finite likelihood at the reference scales and starting priors.",
      call. = FALSE)

  update_scale <- function(j, block, index, center = NULL, spread = 1) {
    evaluate <- function(value) {
      candidate <- par
      candidate[[block]][index] <- value
      k <- state$kernel
      k[[j]] <- .msr_bayes_slice_kernel(spec, candidate, j, shared)
      X <- .msr_bayes_slice_design(spec, k)
      if (is.null(X)) return(list(ll = -Inf))
      eta <- state$offset + drop(X %*% par$beta)
      list(ll = loglik(eta), kernel = k, X = X, eta = eta)
    }
    current <- par[[block]][index]
    accepted <- if (block == "log_scale_sd")
      .msr_bayes_slice_scalar(current, state$ll, evaluate, function(x)
        .msr_bayes_log_sd_prior(x, priors$scale_group_sd)) else if (length(current) == 1L)
      .msr_bayes_slice_scalar(current, state$ll, evaluate, function(x)
        stats::dnorm(x, center, spread, log = TRUE), width = spread) else
      .msr_bayes_slice_ellipse(current, center, spread, state$ll, evaluate)
    par[[block]][index] <<- accepted$value
    state[names(accepted$evaluated)] <<- accepted$evaluated
  }
  update_response_group <- function(term, sd_update) {
    zi <- term$z_index
    si <- term$sd_index
    z <- par$response_z[zi]
    sd <- par$log_response_sd[si]
    old <- response_effect(term, z, sd)
    evaluate <- function(value) {
      effect <- if (sd_update) response_effect(term, z, value) else
        response_effect(term, value, sd)
      change <- effect - old
      eta <- state$eta + change
      list(ll = loglik(eta), offset = state$offset + change, eta = eta)
    }
    if (sd_update) {
      accepted <- .msr_bayes_slice_scalar(sd, state$ll, evaluate, function(x)
        .msr_bayes_log_sd_prior(x, priors$response_group_sd))
      par$log_response_sd[si] <<- accepted$value
    } else {
      accepted <- .msr_bayes_slice_ellipse(z, 0, 1, state$ll, evaluate)
      par$response_z[zi] <<- accepted$value
    }
    state[names(accepted$evaluated)] <<- accepted$evaluated
  }

  # A fixed anchor keeps the reference independent of the current beta block.
  beta_start <- par$beta
  if (spec$family != "gaussian") {
    reference <- .msr_bayes_beta_geometry(state$X, spec$y, spec$trials, state$offset,
      spec$family, par$log_dispersion, priors$beta_sd, spec$parameters$beta)
    if (!reference$fallback) beta_start <- reference$mode
  }
  draws <- matrix(NA_real_, iter - warmup, sum(lengths(par)))
  for (i in seq_len(iter)) {
    if (spec$family == "gaussian") {
      par$beta <- .msr_bayes_slice_gaussian_beta(state$X, spec$y - state$offset,
        par$log_dispersion, priors$beta_sd)
      state$eta <- state$offset + drop(state$X %*% par$beta)
      state$ll <- loglik(state$eta)
    } else {
      geometry <- .msr_bayes_beta_geometry(state$X, spec$y, spec$trials, state$offset,
        spec$family, par$log_dispersion, priors$beta_sd, beta_start)
      accepted <- .msr_bayes_slice_beta(par$beta, state, geometry, priors$beta_sd, loglik)
      par$beta <- accepted$beta
      state$eta <- accepted$eta
      state$ll <- accepted$likelihood
    }
    for (j in seq_along(spec$scales)) {
      sm <- spec$scales[[j]]
      update_scale(j, "scale_beta", sm$beta_index, sm$prior_mean, sm$prior_sd)
      for (term in sm$random) {
        update_scale(j, "scale_z", term$z_index, center = 0)
        update_scale(j, "log_scale_sd", term$sd_index)
      }
    }
    for (term in spec$response$random) {
      update_response_group(term, FALSE)
      update_response_group(term, TRUE)
    }
    if (length(par$log_dispersion)) {
      accepted <- .msr_bayes_slice_scalar(par$log_dispersion, state$ll,
        function(x) list(ll = loglik(state$eta, x)), function(x)
          stats::dnorm(x, priors$dispersion_log_mean, priors$dispersion_log_sd, log = TRUE),
        width = priors$dispersion_log_sd)
      par$log_dispersion <- accepted$value
      state$ll <- accepted$evaluated$ll
    }
    if (i > warmup) draws[i - warmup, ] <- unlist(par, use.names = FALSE)
  }
  draws
}

.msr_bayes_slice_design <- function(spec, kernel) {
  if (any(!is.finite(unlist(kernel, use.names = FALSE)))) return(NULL)
  if (is.null(spec$response_design)) return(.msr_scale_response_matrix(spec, kernel))
  tryCatch(suppressWarnings(.msr_scale_response_matrix(spec, kernel)), error = function(e) {
    # A candidate can leave the domain of a kernel expression, such as log(x).
    if (identical(conditionMessage(e), "Dynamic response design is nonfinite.")) return(NULL)
    stop(e)
  })
}

.msr_bayes_slice_gaussian_beta <- function(X, y, log_sd, prior_sd) {
  # Factor the precision of the exact Gaussian conditional, then draw without
  # forming the inverse covariance.
  inverse_sd <- exp(-log_sd)
  scaled <- X * inverse_sd
  precision <- crossprod(scaled)
  diag(precision) <- diag(precision) + 1 / prior_sd^2
  R <- chol(precision)
  rhs <- crossprod(scaled, y * inverse_sd)
  backsolve(R, forwardsolve(t(R), rhs) + stats::rnorm(ncol(X)))[, 1L]
}

.msr_bayes_slice_beta <- function(beta, state, geometry, prior_sd, loglik) {
  # A Gaussian reference close to the conditional posterior allows larger
  # moves than a prior ellipse. Subtract its log density from the slice target.
  # Its center and precision depend only on the other parameter blocks.
  if (geometry$fallback) {
    evaluate_prior <- function(b) {
      eta <- state$offset + drop(state$X %*% b)
      value <- loglik(eta)
      list(ll = value, likelihood = value, eta = eta, beta = b)
    }
    return(.msr_bayes_slice_ellipse(beta, 0, prior_sd, state$ll, evaluate_prior)$evaluated)
  }
  correction <- function(b, z) (sum(z^2) - sum((b / prior_sd)^2)) / 2
  evaluate <- function(z) {
    b <- geometry$mode + drop(backsolve(geometry$R, z))
    eta <- state$offset + drop(state$X %*% b)
    value <- loglik(eta)
    list(ll = value + correction(b, z), likelihood = value, eta = eta, beta = b)
  }
  z <- drop(geometry$R %*% (beta - geometry$mode))
  .msr_bayes_slice_ellipse(z, 0, 1, state$ll + correction(beta, z), evaluate)$evaluated
}

.msr_bayes_log_sd_prior <- function(x, sd) {
  stats::dnorm(exp(x), 0, sd, log = TRUE) + log(2) + x
}

.msr_bayes_slice_scalar <- function(current, ll, evaluate, prior, width = 1) {
  height <- ll + prior(current) - stats::rexp(1L)
  left <- current - width * stats::runif(1L)
  right <- left + width
  # Randomly allocate a finite step-out budget to keep the update reversible.
  steps_left <- sample.int(20L, 1L) - 1L
  steps_right <- 19L - steps_left
  above <- function(x) isTRUE(evaluate(x)$ll + prior(x) > height)
  while (steps_left > 0L && above(left)) {
    left <- left - width
    steps_left <- steps_left - 1L
  }
  while (steps_right > 0L && above(right)) {
    right <- right + width
    steps_right <- steps_right - 1L
  }
  for (i in seq_len(1000L)) {
    candidate <- stats::runif(1L, left, right)
    evaluated <- evaluate(candidate)
    if (isTRUE(evaluated$ll + prior(candidate) >= height))
      return(list(value = candidate, evaluated = evaluated))
    if (candidate < current) left <- candidate else right <- candidate
  }
  stop("A scalar slice update did not find a finite accepted value.", call. = FALSE)
}

.msr_bayes_slice_ellipse <- function(current, center, spread, ll, evaluate) {
  height <- ll - stats::rexp(1L)
  direction <- stats::rnorm(length(current)) * spread
  angle <- stats::runif(1L, 0, 2 * pi)
  left <- angle - 2 * pi
  right <- angle
  centered <- current - center
  for (i in seq_len(1000L)) {
    candidate <- center + cos(angle) * centered + sin(angle) * direction
    evaluated <- evaluate(candidate)
    if (isTRUE(evaluated$ll >= height))
      return(list(value = candidate, evaluated = evaluated))
    if (angle < 0) left <- angle else right <- angle
    angle <- stats::runif(1L, left, right)
  }
  stop("An elliptical slice update did not find a finite accepted value.", call. = FALSE)
}

# Expand a nested response-effect approximation into a joint proposal. The
# cross-covariance follows the derivative of the conditional random-effect mode.
.msr_bayes_nested_expand <- function(node, conditional, joint, hyper) {
  conditional$fn(node$mode)
  q <- conditional$env$last.par
  stopifnot(identical(names(q), names(joint$par)))
  random <- which(names(q) == "response_z")
  fixed <- setdiff(seq_along(q), c(random, hyper))
  nuisance <- setdiff(seq_along(q), hyper)
  H <- joint$he(q)
  Hrr <- H[random, random, drop = FALSE]
  geometry <- .msr_bayes_laplace_geometry(Hrr)
  B <- -geometry$covariance %*% H[random, fixed, drop = FALSE]
  V <- matrix(0, length(nuisance), length(nuisance))
  f <- match(fixed, nuisance); r <- match(random, nuisance)
  V[f, f] <- node$covariance
  V[r, f] <- B %*% node$covariance
  V[f, r] <- t(V[r, f, drop = FALSE])
  V[r, r] <- geometry$covariance + B %*% node$covariance %*% t(B)
  node$nested <- list(fixed_mode = node$mode, fixed_hessian = node$hessian,
    random_hessian = Hrr, random_slope = B, fixed_indices = fixed, random_indices = random)
  node$mode <- q
  node$covariance <- (V + t(V)) / 2
  node$hessian <- chol2inv(chol(node$covariance))
  expanded <- .msr_bayes_laplace_geometry(node$hessian)
  node$root <- expanded$root
  node$condition_number <- expanded$condition_number
  node
}

.msr_bayes_logadd <- function(a, b) {
  m <- pmax(a, b)
  ans <- m + log(exp(a - m) + exp(b - m))
  ans[is.infinite(m) & m < 0] <- -Inf
  ans
}

# A heavy-tailed conditional proposal retains support beyond the local Gaussian.
.msr_bayes_importance_density <- function(u, node, df = 7, inflation = 1.2) {
  d <- ncol(u)
  z <- sweep(u, 2L, node$mode, "-") %*% t(chol(node$hessian)) / inflation
  radius <- rowSums(z * z)
  log_det <- sum(log(diag(chol(node$hessian)))) - d * log(inflation)
  shared <- lgamma((df + d) / 2) - lgamma(df / 2) - d / 2 * log(df * pi) -
    (df + d) / 2 * log1p(radius / df)
  separate <- rowSums(stats::dt(z, df = df, log = TRUE))
  log_det + .msr_bayes_logadd(log(.5) + shared, log(.5) + separate)
}

.msr_bayes_importance_tdraw <- function(node, hyper, n, df = 7) {
  nuisance <- setdiff(seq_along(node$mode), hyper)
  z <- matrix(stats::rnorm(n * length(nuisance)), n)
  shared <- stats::runif(n) < .5
  scale <- matrix(sqrt(stats::rchisq(length(z), df) / df), n)
  if (any(shared)) scale[shared, ] <- sqrt(stats::rchisq(sum(shared), df) / df)
  z <- (1.2 * z / scale) %*% t(node$root)
  sweep(z, 2L, node$mode[nuisance], "+")
}

.msr_bayes_importance_candidates <- function(grid, hyper, n, adapted = NULL) {
  if (identical(grid$proposal_type, "latent_prior"))
    return(.msr_bayes_latent_prior_draws(grid, n))
  if (!is.null(adapted)) {
    if (!is.null(adapted$components)) {
      ids <- sample.int(length(adapted$components), n, replace = TRUE, prob = adapted$weights)
      ans <- matrix(NA_real_, n, length(grid$nodes[[1L]]$mode))
      for (j in unique(ids)) {
        rows <- which(ids == j)
        ans[rows, ] <- .msr_bayes_importance_candidates(adapted$components[[j]], hyper, length(rows))
      }
      return(ans)
    }
    use_adapted <- stats::runif(n) < .8
    ans <- matrix(NA_real_, n, length(grid$nodes[[1L]]$mode))
    for (flag in c(FALSE, TRUE)) {
      rows <- which(use_adapted == flag)
      if (length(rows)) ans[rows, ] <- .msr_bayes_importance_candidates(
        if (flag) adapted else grid, hyper, length(rows))
    }
    return(ans)
  }
  ids <- sample.int(length(grid$x), n, replace = TRUE, prob = grid$weights)
  df <- if (is.null(grid$df)) 7 else grid$df
  d <- length(grid$nodes[[1L]]$mode)
  ans <- matrix(NA_real_, n, d)
  nuisance <- setdiff(seq_len(d), hyper)
  for (i in unique(ids)) {
    rows <- which(ids == i)
    left <- if (i > 1L) grid$x[i] - grid$x[i - 1L] else 0
    right <- if (i < length(grid$x)) grid$x[i + 1L] - grid$x[i] else 0
    go_left <- stats::runif(length(rows)) < left / (left + right)
    fraction <- sqrt(stats::runif(length(rows)))
    x <- ifelse(go_left, grid$x[i] - left + left * fraction,
                         grid$x[i] + right - right * fraction)
    x <- pmin(1 - .Machine$double.eps, pmax(.Machine$double.eps, x))
    ans[rows, hyper] <- log(grid$prior_sd * stats::qnorm((1 + x) / 2))
    ans[rows, nuisance] <- .msr_bayes_importance_tdraw(grid$nodes[[i]], hyper, length(rows), df)
  }
  if (!is.null(grid$tail_fraction) && grid$tail_fraction > 0) {
    rows <- which(stats::runif(n) < grid$tail_fraction)
    if (length(rows)) {
      ans[rows, hyper] <- stats::rnorm(length(rows), log(2 * grid$prior_sd), .6)
      x <- 2 * stats::pnorm(exp(ans[rows, hyper]) / grid$prior_sd) - 1
      i <- findInterval(x, grid$x, all.inside = TRUE)
      fraction <- (x - grid$x[i]) / (grid$x[i + 1L] - grid$x[i])
      left <- (1 - fraction) * grid$density[i]
      right <- fraction * grid$density[i + 1L]
      selected <- i + as.integer(stats::runif(length(rows)) >= left / (left + right))
      for (j in unique(selected)) {
        take <- rows[selected == j]
        ans[take, nuisance] <- .msr_bayes_importance_tdraw(grid$nodes[[j]], hyper, length(take), df)
      }
    }
  }
  ans
}

.msr_bayes_importance_logq <- function(samples, grid, hyper, adapted = NULL) {
  if (identical(grid$proposal_type, "latent_prior"))
    return(.msr_bayes_latent_prior_logq(samples, grid))
  if (!is.null(adapted$components)) {
    value <- rep(-Inf, nrow(samples))
    for (j in seq_along(adapted$components)) value <- .msr_bayes_logadd(value,
      log(adapted$weights[j]) + .msr_bayes_importance_logq(samples, adapted$components[[j]], hyper))
    return(value)
  }
  if (!is.null(adapted)) return(.msr_bayes_logadd(
    log(.2) + .msr_bayes_importance_logq(samples, grid, hyper),
    log(.8) + .msr_bayes_importance_logq(samples, adapted, hyper)))
  tau <- exp(samples[, hyper])
  df <- if (is.null(grid$df)) 7 else grid$df
  x <- 2 * stats::pnorm(tau / grid$prior_sd) - 1
  i <- findInterval(x, grid$x, all.inside = TRUE)
  fraction <- (x - grid$x[i]) / (grid$x[i + 1L] - grid$x[i])
  logq <- numeric(nrow(samples))
  for (j in unique(i)) {
    rows <- which(i == j)
    a <- grid$nodes[[j]]; a$mode <- a$mode[-hyper]
    b <- grid$nodes[[j + 1L]]; b$mode <- b$mode[-hyper]
    logq[rows] <- .msr_bayes_logadd(
      log(grid$density[j]) + log1p(-fraction[rows]) +
        .msr_bayes_importance_density(samples[rows, -hyper, drop = FALSE], a, df),
      log(grid$density[j + 1L]) + log(fraction[rows]) +
        .msr_bayes_importance_density(samples[rows, -hyper, drop = FALSE], b, df))
  }
  # Convert the proposal in prior-CDF coordinates back to log-SD coordinates.
  prior_h <- stats::dnorm(tau, 0, grid$prior_sd, log = TRUE) + log(2) + samples[, hyper]
  if (is.null(grid$tail_fraction) || grid$tail_fraction == 0) return(logq + prior_h)
  density_x <- (1 - fraction) * grid$density[i] + fraction * grid$density[i + 1L]
  # A log-SD tail component extends the proposal beyond the numerical upper
  # limit of the prior-CDF inverse. The conditional proposal still uses the
  # nearest endpoint components, and the full posterior corrects their error.
  logq - log(density_x) + .msr_bayes_logadd(
    log1p(-grid$tail_fraction) + log(density_x) + prior_h,
    log(grid$tail_fraction) + stats::dnorm(samples[, hyper], log(2 * grid$prior_sd), .6, log = TRUE))
}

# Pilot weights train a proposal only. Tempering and covariance shrinkage here
# stabilize proposal fitting; production weights still use the exact posterior.
# Independent production samples avoid recycling adaptive training observations.
.msr_bayes_importance_adapt <- function(samples, log_weights, grid, hyper, power = .75) {
  x <- 2 * stats::pnorm(exp(samples[, hyper]) / grid$prior_sd) - 1
  u <- samples[, -hyper, drop = FALSE]
  weight <- exp(power * (log_weights - max(log_weights)))
  adapted <- grid[intersect(c("x", "density", "weights", "prior_sd", "tail_fraction", "df"), names(grid))]
  adapted$nodes <- lapply(grid$nodes, function(node)
    node[c("mode", "covariance", "hessian", "root")])
  adapted$skipped_nodes <- 0L
  for (j in seq_along(grid$x)) {
    local <- weight * exp(-.5 * ((x - grid$x[j]) / .12)^2)
    local <- local / sum(local)
    ess <- 1 / sum(local^2)
    if (!is.finite(ess) || ess < ncol(u) + 5) {
      adapted$skipped_nodes <- adapted$skipped_nodes + 1L; next
    }
    offset <- x - grid$x[j]
    center <- sum(local * offset)
    spread <- sum(local * (offset - center)^2)
    if (spread <= 1e-12) { adapted$skipped_nodes <- adapted$skipped_nodes + 1L; next }
    mean_u <- colSums(u * local)
    slope <- colSums(sweep(u, 2L, mean_u, "-") * (local * (offset - center))) / spread
    mode <- mean_u - slope * center
    residual <- sweep(u, 2L, mode, "-") - offset %o% slope
    covariance <- crossprod(residual * sqrt(local))
    # This is a deliberately conservative proposal covariance, not a repaired
    # likelihood Hessian or a modified posterior uncertainty estimate.
    covariance <- .8 * covariance + .2 * grid$nodes[[j]]$covariance
    geometry <- tryCatch(.msr_bayes_laplace_geometry(chol2inv(chol(covariance))),
      error = function(e) NULL)
    if (is.null(geometry)) { adapted$skipped_nodes <- adapted$skipped_nodes + 1L; next }
    adapted$nodes[[j]]$mode[-hyper] <- mode
    adapted$nodes[[j]]$covariance <- geometry$covariance
    adapted$nodes[[j]]$root <- geometry$root
    adapted$nodes[[j]]$hessian <- chol2inv(chol(geometry$covariance))
  }
  adapted
}

.msr_bayes_weighted_quantile <- function(x, weights, probs) {
  o <- order(x)
  cumulative <- cumsum(weights[o]) / sum(weights)
  # Inverse CDF uses the first value whose cumulative mass is at least p.
  # left.open prevents an exact cumulative probability from jumping to the
  # next value, which matters for small discrete weighted populations.
  x[o][pmin(length(x), findInterval(probs, cumulative, left.open = TRUE) + 1L)]
}

.msr_bayes_importance_targets <- function(samples, spec) {
  sizes <- lengths(spec$parameters)
  ends <- cumsum(sizes); begins <- ends - sizes
  blocks <- lapply(seq_along(sizes), function(i)
    samples[, begins[i] + seq_len(sizes[i]), drop = FALSE])
  names(blocks) <- names(sizes)
  annual <- lapply(spec$scales, function(sm) {
    rows <- sm$representatives
    value <- blocks$scale_beta[, sm$beta_index, drop = FALSE] %*% t(sm$X[rows, , drop = FALSE])
    for (term in sm$random)
      value <- value + blocks$scale_z[, term$z_index[term$index[rows]], drop = FALSE] *
        exp(blocks$log_scale_sd[, term$sd_index])
    value
  })
  cbind(samples, do.call(cbind, annual), exp(blocks$log_scale_sd), exp(blocks$log_response_sd))
}

.msr_bayes_importance_checks <- function(samples, log_weights, spec) {
  if (anyNA(log_weights) || any(log_weights == Inf) || !any(is.finite(log_weights)))
    stop("Invalid posterior importance ratios.", call. = FALSE)
  weights <- exp(log_weights - max(log_weights)); weights <- weights / sum(weights)
  # Zero posterior density is permitted; replace -Inf only for Pareto fitting.
  diagnostic_weights <- pmax(weights, .Machine$double.xmin)
  k <- if (diff(range(diagnostic_weights)) <= 1e-12 * max(diagnostic_weights)) 0 else
    posterior::pareto_khat(diagnostic_weights, tail = "right", r_eff = 1)
  targets <- .msr_bayes_importance_targets(samples, spec)
  min_tail <- Inf; max_split <- 0
  half <- seq_len(nrow(samples)) <= nrow(samples) / 2
  for (j in seq_len(ncol(targets))) {
    x <- targets[, j]
    limits <- .msr_bayes_weighted_quantile(x, weights, c(.025, .5, .975))
    for (q in limits) {
      indicator <- as.numeric(x <= q)
      p <- min(1, max(0, sum(weights * indicator)))
      variance <- sum(weights^2 * (indicator - p)^2)
      if (variance > 0 && p > 0 && p < 1) min_tail <- min(min_tail, p * (1 - p) / variance)
      pa <- sum(weights[half] * indicator[half]) / sum(weights[half])
      pb <- sum(weights[!half] * indicator[!half]) / sum(weights[!half])
      max_split <- max(max_split, abs(pa - pb))
    }
  }
  ess <- 1 / sum(weights^2)
  passed <- is.finite(k) && k <= .7 && ess >= 1000 && min_tail >= 400 &&
    max(weights) <= .01 && max_split <= .05
  list(passed = passed, pareto_k = as.numeric(k), ess = ess,
    min_tail_ess = min_tail, max_weight = max(weights), split_cdf_error = max_split,
    weights = weights)
}

# Importance proposals can visit extreme but valid finite parameters where an
# AD density evaluates to NaN although the ordinary full target is defined.
# Recheck only those values. Negative infinity is valid zero posterior mass;
# unresolved targets and invalid proposal densities still fail the same checks.
.msr_bayes_importance_target <- function(samples, obj, spec) {
  value <- vapply(seq_len(nrow(samples)), function(i) -obj$fn(samples[i, ]), 0)
  invalid <- which(is.na(value) | value == Inf)
  counts <- c(evaluations = length(value), ad_invalid = length(invalid),
    recovered_finite = 0L, recovered_zero = 0L, unresolved = 0L)
  if (length(invalid)) {
    ordinary <- .msr_bayes_objective(spec, ad = FALSE)
    blocks <- .msr_bayes_parameter_blocks(samples[invalid, , drop = FALSE], spec$parameters)
    for (j in seq_along(invalid)) {
      parameters <- lapply(blocks, function(x) unname(x[j, ]))
      checked <- tryCatch(-ordinary(parameters), error = function(e) NA_real_)
      if (length(checked) == 1L && is.finite(checked)) {
        value[invalid[j]] <- checked
        counts["recovered_finite"] <- counts["recovered_finite"] + 1L
      } else if (isTRUE(checked == -Inf)) {
        value[invalid[j]] <- -Inf
        counts["recovered_zero"] <- counts["recovered_zero"] + 1L
      }
    }
    counts["unresolved"] <- sum(is.na(value[invalid]) | value[invalid] == Inf)
  }
  list(log_density = value, diagnostics = counts)
}

.msr_bayes_importance <- function(out, obj, draws, seed, control) {
  out$inference_metadata <- NULL
  minimum <- if (is.null(control$importance_min)) 16384L else control$importance_min
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  grid <- out$marginal; hyper <- grid$hyper_index
  indices <- seq_along(obj$par)
  conditional_beta <- isTRUE(control$conditional_beta)
  conditional_response <- !conditional_beta && control$nested_response &&
    length(out$specification$response$random) == 1L
  if (conditional_beta) {
    projected <- .msr_bayes_beta_project(grid, obj)
    grid <- projected$grid; hyper <- projected$hyper; indices <- projected$indices
    beta_start <- projected$start
  }
  if (conditional_response) {
    if (!is.null(grid$parameter_indices) && length(grid$parameter_indices) < length(obj$par)) {
      indices <- grid$parameter_indices
    } else {
      projected <- .msr_bayes_response_project(grid, obj)
      grid <- projected$grid; hyper <- projected$hyper; indices <- projected$indices
    }
  }
  # A uniform component in prior-CDF coordinates guarantees positive proposal
  # support even when a very small approximate node density underflows.
  nodes <- lapply(grid$nodes, function(node) node[c("mode", "root", "covariance", "hessian")])
  normalized <- .msr_bayes_marginal_weights(grid$x, log(.99 * grid$density + .01))
  # Degrees of freedom affect only the proposal. Every mixture component keeps
  # its own setting so the exact proposal density follows sampling and training.
  df <- if (is.null(control$importance_df)) 7 else control$importance_df
  grid <- c(normalized, list(prior_sd = grid$prior_sd, nodes = nodes, tail_fraction = .01, df = df))
  response_fallbacks <- 0L; response_max_step <- 0
  beta_fallbacks <- 0L; beta_max_step <- 0
  augment <- function(outer) {
    if (conditional_beta) {
      value <- .msr_bayes_beta_candidates(outer, indices, out$specification, beta_start)
      beta_fallbacks <<- beta_fallbacks + value$fallbacks
      beta_max_step <<- max(beta_max_step, value$max_step)
      return(value)
    }
    if (!conditional_response) return(list(samples = outer, logq = rep(0, nrow(outer))))
    value <- .msr_bayes_response_candidates(outer, indices, out$specification)
    response_fallbacks <<- response_fallbacks + value$fallbacks
    response_max_step <<- max(response_max_step, value$max_step)
    value
  }
  pilot_proposal <- NULL
  if (conditional_beta) {
    prior_grid <- .msr_bayes_latent_prior_proposal(out, indices)
    pilot_proposal <- list(components = list(prior_grid, grid), weights = c(.2, .8))
  }
  pilot_outer <- .msr_bayes_importance_candidates(grid, hyper, 16384L, adapted = pilot_proposal)
  pilot <- augment(pilot_outer)
  pilot_logq <- .msr_bayes_importance_logq(pilot_outer, grid, hyper, adapted = pilot_proposal) + pilot$logq
  pilot <- pilot$samples
  pilot_target <- .msr_bayes_importance_target(pilot, obj, out$specification)
  pilot_logp <- pilot_target$log_density
  target_evaluation <- pilot_target$diagnostics
  learned <- .msr_bayes_importance_adapt(pilot_outer, pilot_logp - pilot_logq, grid, hyper)
  adapted <- list(components = list(grid, learned), weights = c(.2, .8))
  if (conditional_beta)
    adapted <- list(components = list(prior_grid, grid, learned), weights = c(.2, .2, .6))
  samples <- matrix(numeric(), 0L, length(obj$par)); log_weights <- numeric()
  history <- list(); total_production <- 0L
  for (round in seq.int(0L, control$importance_refinements)) {
   samples <- matrix(numeric(), 0L, length(obj$par)); log_weights <- numeric()
   for (n in c(16384L, 32768L, 65536L, 131072L)) {
    if (n > control$importance_max) break
    # A successful small batch may grow to an explicitly requested precision
    # budget without changing its proposal. Failed batches can still train the
    # next proposal; only its fresh production population is retained.
    if (round < control$importance_refinements && n > 16384L && !check$passed) break
    outer <- .msr_bayes_importance_candidates(grid, hyper, n - nrow(samples), adapted)
    extra <- augment(outer)
    logq <- .msr_bayes_importance_logq(outer, grid, hyper, adapted) + extra$logq
    extra <- extra$samples
    total_production <- total_production + nrow(extra)
    production_target <- .msr_bayes_importance_target(extra, obj, out$specification)
    logp <- production_target$log_density
    target_evaluation <- target_evaluation + production_target$diagnostics
    samples <- rbind(samples, extra)
    log_weights <- c(log_weights, logp - logq)
    check <- .msr_bayes_importance_checks(samples, log_weights, out$specification)
    history[[length(history) + 1L]] <- data.frame(refinement = round, candidates = n,
      pareto_k = check$pareto_k, ess = check$ess, min_tail_ess = check$min_tail_ess,
      max_weight = check$max_weight, split_cdf_error = check$split_cdf_error,
      passed = check$passed, budget_met = n >= minimum)
    if (check$passed && n >= minimum) break
   }
   if ((check$passed && nrow(samples) >= minimum) || round == control$importance_refinements) break
   learned <- .msr_bayes_importance_adapt(samples[, indices, drop = FALSE], log_weights,
     adapted$components[[length(adapted$components)]], hyper, power = .9)
   adapted$components <- c(adapted$components, list(learned))
   adapted$weights <- c(.4 * adapted$weights, .6)
   if (conditional_beta) {
     # Keep defensive latent support through every refinement, rather than
     # letting repeated adaptation reduce its mixture probability.
     adapted$weights[-1L] <- .8 * adapted$weights[-1L] / sum(adapted$weights[-1L])
     adapted$weights[1L] <- .2
   }
  }
  ids <- sample.int(nrow(samples), draws, replace = TRUE, prob = check$weights)
  sizes <- lengths(out$parameters); ends <- cumsum(sizes); begins <- ends - sizes
  out$draws <- lapply(seq_along(sizes), function(j)
    samples[ids, begins[j] + seq_len(sizes[j]), drop = FALSE])
  names(out$draws) <- names(out$parameters)
  out$importance <- list(candidates = samples, weights = check$weights,
    log_weights = log_weights, resample_ids = ids, history = do.call(rbind, history),
    seed = seed, proposal = paste0("Equal multivariate/coordinatewise Student t, df ", df,
      ", scale 1.2; .01 uniform prior-CDF support and .01 log-SD normal tail"),
    adapted_proposal = adapted, pilot_candidates = nrow(pilot),
    proposal_indices = indices, conditional_response = conditional_response,
    production_candidates = total_production,
    response_proposal = list(fallbacks = response_fallbacks, max_step = response_max_step),
    proposal_diagnostics = out$diagnostics, target_evaluation = target_evaluation)
  if (conditional_beta) {
    out$importance$conditional_beta <- TRUE
    out$importance$beta_proposal <- list(fallbacks = beta_fallbacks,
      max_step = beta_max_step, start = beta_start, latent_prior_fraction = .2)
    out$importance$proposal <- paste(out$importance$proposal,
      "; conditional Gaussian/t7 response coefficients and .2 latent-prior defensive mixture")
  }
  out$inference <- "importance-corrected marginal Laplace"
  out$diagnostics <- c(check[setdiff(names(check), "weights")],
    list(type = "Marginal importance diagnostics", candidates = nrow(samples),
      min_candidates = minimum,
      proposal_df = df,
      production_candidates = total_production,
      proposal_passed = out$importance$proposal_diagnostics$passed,
      criteria = paste("Pareto k <= .7; raw importance ESS >= 1000; minimum indicator-tail ESS >= 400;",
        "maximum weight <= .01; split-half CDF difference <= .05.",
        "Proposal flags are retained separately. These finite-sample checks do not prove accuracy.")))
  if (!check$passed)
    warning("Marginal importance diagnostics need attention; use full NUTS or increase importance_max.", call. = FALSE)
  out$inference_metadata <- .msr_scale_inference_metadata("posterior", "importance_corrected")
  out
}

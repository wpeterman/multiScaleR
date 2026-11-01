# Marginal covariance for an outer proposal excluding fixed response coefficients.
# This changes the proposal dimension, not the posterior or its prior factors.
.msr_bayes_beta_project <- function(grid, obj) {
  keep <- which(names(obj$par) != "beta")
  nuisance <- setdiff(seq_along(obj$par), grid$hyper_index)
  small <- match(setdiff(keep, grid$hyper_index), nuisance)
  beta <- which(names(obj$par) == "beta")
  starts <- do.call(rbind, lapply(grid$nodes, function(node) node$mode[beta]))
  ans <- grid[c("x", "density", "weights", "prior_sd")]
  ans$nodes <- lapply(grid$nodes, function(node) {
    covariance <- node$covariance[small, small, drop = FALSE]
    hessian <- chol2inv(chol(covariance))
    geometry <- .msr_bayes_laplace_geometry(hessian)
    list(mode = node$mode[keep], covariance = geometry$covariance,
         root = geometry$root, hessian = hessian)
  })
  list(grid = ans, indices = keep, hyper = match(grid$hyper_index, keep),
       start = colSums(starts * grid$weights))
}

# Normal priors give positive conditional curvature for all supported GLMs.
# An unfinished center is still a normalized proposal; its exact density is
# evaluated below. Nonfinite geometry uses the normalized coefficient prior.
.msr_bayes_beta_geometry <- function(X, y, trials, offset, family, dispersion,
                                     prior_sd, start) {
  p <- ncol(X)
  precision <- 1 / prior_sd^2
  evaluate <- function(b) {
    eta <- drop(X %*% b) + offset
    if (family == "gaussian") {
      variance <- exp(2 * dispersion)
      score <- (y - eta) / variance
      curvature <- rep(1 / variance, length(eta))
    } else if (family == "poisson") {
      mu <- exp(eta)
      score <- y - mu; curvature <- mu
    } else if (family == "binomial") {
      prob <- stats::plogis(eta)
      score <- y - trials * prob
      curvature <- trials * prob * (1 - prob)
    } else {
      theta <- exp(dispersion)
      prob <- stats::plogis(eta - dispersion)
      score <- y - (y + theta) * prob
      curvature <- (y + theta) * prob * (1 - prob)
    }
    list(gradient = precision * b - drop(crossprod(X, score)),
         hessian = diag(precision, p) + crossprod(X * sqrt(curvature)))
  }
  b <- start
  bad <- FALSE
  for (step in seq_len(20L)) {
    geometry <- evaluate(b)
    R <- if (any(!is.finite(geometry$hessian))) NULL else
      tryCatch(chol(geometry$hessian), error = function(e) NULL)
    if (is.null(R) || any(!is.finite(R)) || any(!is.finite(geometry$gradient))) { bad <- TRUE; break }
    direction <- backsolve(R, forwardsolve(t(R), geometry$gradient))
    if (any(!is.finite(direction))) { bad <- TRUE; break }
    if (max(abs(direction)) < 1e-8) break
    b <- b - if (family == "gaussian") direction else pmax(-1, pmin(1, direction))
  }
  geometry <- evaluate(b)
  R <- if (any(!is.finite(geometry$hessian))) NULL else
    tryCatch(chol(geometry$hessian), error = function(e) NULL)
  fallback <- bad || is.null(R) || any(!is.finite(R)) || any(!is.finite(geometry$gradient))
  list(mode = b, hessian = geometry$hessian, R = R, fallback = fallback,
       max_step = if (fallback) 0 else
         max(abs(backsolve(R, forwardsolve(t(R), geometry$gradient)))))
}

.msr_bayes_beta_candidates <- function(outer, indices, spec, start, given = NULL) {
  samples <- matrix(0, nrow(outer), sum(lengths(spec$parameters)))
  samples[, indices] <- outer
  beta <- setdiff(seq_len(ncol(samples)), indices)
  if (!is.null(given) && !identical(dim(given), dim(samples)))
    stop("Given coefficient evaluations must preserve the full parameter layout.", call. = FALSE)
  logq <- numeric(nrow(samples)); fallbacks <- 0L; max_step <- 0
  p <- length(beta); fixed <- ncol(spec$X_fixed)
  for (begin in seq.int(1L, nrow(samples), by = 256L)) {
    rows <- seq.int(begin, min(begin + 255L, nrow(samples)))
    blocks <- .msr_bayes_parameter_blocks(samples[rows, , drop = FALSE], spec$parameters)
    kernels <- if (is.null(spec$response_design))
      lapply(seq_along(spec$covariates), function(j) {
        b <- blocks; b$beta[, fixed + j] <- 1
        sweep(.msr_bayes_batch_eta(spec, b), 2L, spec$offset, "-")
      }) else NULL
    for (i in seq_along(rows)) {
      X <- if (is.null(spec$response_design))
        cbind(spec$X_fixed, do.call(cbind, lapply(kernels, function(k) k[i, ]))) else {
          # One kernel evaluation per candidate, rather than probing a response
          # coefficient for each expanded column or repeating every kernel.
          candidate <- lapply(blocks, function(x) unname(x[i, ]))
          k <- .msr_bayes_components(spec, candidate, include_response = FALSE)$kernel
          .msr_scale_response_matrix(spec, k)
        }
      offset <- spec$offset
      for (term in spec$response$random)
        offset <- offset + blocks$response_z[i, term$z_index][term$index] *
          exp(blocks$log_response_sd[i, term$sd_index])
      dispersion <- if (ncol(blocks$log_dispersion)) blocks$log_dispersion[i, 1L] else 0
      geometry <- .msr_bayes_beta_geometry(X, spec$y, spec$trials, offset,
        spec$family, dispersion, spec$priors$beta_sd, start)
      if (geometry$fallback) {
        fallbacks <- fallbacks + 1L
        value <- if (is.null(given)) stats::rnorm(p, sd = spec$priors$beta_sd) else given[rows[i], beta]
        samples[rows[i], beta] <- value
        logq[rows[i]] <- sum(stats::dnorm(value, 0, spec$priors$beta_sd, log = TRUE))
        next
      }
      max_step <- max(max_step, geometry$max_step)
      if (!is.null(given)) {
        value <- given[rows[i], beta]
        if (spec$family == "gaussian") {
          z <- drop(geometry$R %*% (value - geometry$mode))
          logq[rows[i]] <- sum(log(diag(geometry$R))) - p / 2 * log(2 * pi) - sum(z^2) / 2
        } else logq[rows[i]] <- .msr_bayes_importance_density(matrix(value, 1L), geometry)
      } else if (spec$family == "gaussian") {
        z <- stats::rnorm(p)
        value <- geometry$mode + drop(backsolve(geometry$R, z))
        logq[rows[i]] <- sum(log(diag(geometry$R))) - p / 2 * log(2 * pi) - sum(z^2) / 2
      } else {
        z <- stats::rnorm(p)
        scale <- if (stats::runif(1) < .5) rep(sqrt(stats::rchisq(1, 7) / 7), p) else
          sqrt(stats::rchisq(p, 7) / 7)
        value <- geometry$mode + drop(backsolve(geometry$R, 1.2 * z / scale))
        logq[rows[i]] <- .msr_bayes_importance_density(matrix(value, 1L), geometry)
      }
      samples[rows[i], beta] <- value
    }
  }
  list(samples = samples, logq = logq, fallbacks = fallbacks, max_step = max_step)
}

# A fixed defensive component covers latent-scale regions missed by local
# covariance estimates. Its priors are proposal factors from the existing
# model, not modified model priors. Dispersion uses a separate t proposal.
.msr_bayes_latent_prior_proposal <- function(out, indices) {
  spec <- out$specification
  ans <- list(proposal_type = "latent_prior", parameters = out$parameters,
    indices = indices, priors = spec$priors,
    scales = lapply(spec$scales, function(sm) sm[c("beta_index", "prior_mean", "prior_sd")]))
  if (length(out$parameters$log_dispersion)) {
    pos <- which(rep(names(out$parameters), lengths(out$parameters)) == "log_dispersion")
    nuisance <- setdiff(seq_len(sum(lengths(out$parameters))), out$marginal$hyper_index)
    means <- vapply(out$marginal$nodes, function(node) node$mode[pos], 0)
    variances <- vapply(out$marginal$nodes, function(node)
      node$covariance[match(pos, nuisance), match(pos, nuisance)], 0)
    center <- sum(out$marginal$weights * means)
    spread <- 1.2 * sqrt(sum(out$marginal$weights * (variances + (means - center)^2)))
    if (!is.finite(spread) || spread <= 0)
      stop("Invalid dispersion proposal variance.", call. = FALSE)
    ans$dispersion <- list(mean = center, scale = spread)
  }
  ans
}

.msr_bayes_latent_prior_draws <- function(grid, n) {
  blocks <- lapply(grid$parameters, function(p) matrix(0, n, length(p)))
  for (sm in grid$scales) for (j in seq_along(sm$beta_index))
    blocks$scale_beta[, sm$beta_index[j]] <- stats::rnorm(n, sm$prior_mean[j], sm$prior_sd[j])
  for (nm in c("scale_z", "response_z")) if (length(blocks[[nm]]))
    blocks[[nm]][] <- stats::rnorm(length(blocks[[nm]]))
  if (length(blocks$log_scale_sd)) blocks$log_scale_sd[] <-
    log(abs(stats::rnorm(length(blocks$log_scale_sd), sd = grid$priors$scale_group_sd)))
  if (length(blocks$log_response_sd)) blocks$log_response_sd[] <-
    log(abs(stats::rnorm(length(blocks$log_response_sd), sd = grid$priors$response_group_sd)))
  if (length(blocks$log_dispersion)) blocks$log_dispersion[] <-
    grid$dispersion$mean + grid$dispersion$scale * stats::rt(n, 7)
  do.call(cbind, blocks)[, grid$indices, drop = FALSE]
}

.msr_bayes_latent_prior_logq <- function(samples, grid) {
  full <- matrix(0, nrow(samples), sum(lengths(grid$parameters)))
  full[, grid$indices] <- samples
  blocks <- .msr_bayes_parameter_blocks(full, grid$parameters)
  value <- numeric(nrow(samples))
  for (sm in grid$scales) for (j in seq_along(sm$beta_index))
    value <- value + stats::dnorm(blocks$scale_beta[, sm$beta_index[j]],
      sm$prior_mean[j], sm$prior_sd[j], log = TRUE)
  for (nm in c("scale_z", "response_z")) if (length(blocks[[nm]]))
    value <- value + rowSums(stats::dnorm(blocks[[nm]], log = TRUE))
  for (nm in c("log_scale_sd", "log_response_sd")) if (length(blocks[[nm]])) {
    sd <- grid$priors[[if (nm == "log_scale_sd") "scale_group_sd" else "response_group_sd"]]
    value <- value + rowSums(stats::dnorm(exp(blocks[[nm]]), 0, sd, log = TRUE) + log(2) + blocks[[nm]])
  }
  if (length(blocks$log_dispersion)) value <- value +
    stats::dt((blocks$log_dispersion[, 1L] - grid$dispersion$mean) / grid$dispersion$scale,
      7, log = TRUE) - log(grid$dispersion$scale)
  value
}

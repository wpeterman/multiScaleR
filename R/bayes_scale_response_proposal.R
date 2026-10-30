# One random-intercept term has independent site effects conditional on the
# other parameters. Update its proposal at each sampled parameter vector,
# rather than extending the mode's linear approximation across the SD funnel.
.msr_bayes_response_project <- function(grid, obj) {
  keep <- which(names(obj$par) != "response_z")
  old_nuisance <- setdiff(seq_along(obj$par), grid$hyper_index)
  small_nuisance <- match(setdiff(keep, grid$hyper_index), old_nuisance)
  ans <- grid[c("x", "density", "weights", "prior_sd")]
  ans$nodes <- lapply(grid$nodes, function(node) {
    covariance <- node$covariance[small_nuisance, small_nuisance, drop = FALSE]
    hessian <- chol2inv(chol(covariance))
    geometry <- .msr_bayes_laplace_geometry(hessian)
    list(mode = node$mode[keep], covariance = geometry$covariance,
         root = geometry$root, hessian = hessian)
  })
  list(grid = ans, indices = keep, hyper = match(grid$hyper_index, keep))
}

.msr_bayes_parameter_blocks <- function(samples, parameters) {
  sizes <- lengths(parameters); ends <- cumsum(sizes); begins <- ends - sizes
  ans <- lapply(seq_along(sizes), function(i)
    samples[, begins[i] + seq_len(sizes[i]), drop = FALSE])
  names(ans) <- names(parameters)
  ans
}

# Calculate fixed-response predictors for a modest batch of candidates. The
# distance-bin loop avoids a candidates x observations x bins array.
.msr_bayes_batch_eta <- function(spec, blocks) {
  if (length(spec$surface$buffers)) {
    # SQ uses raw cells, not binned sums. Iterate over this bounded candidate
    # batch without allocating a candidates x observations x cells array.
    eta <- matrix(0, nrow(blocks$beta), length(spec$y))
    for (i in seq_len(nrow(eta))) {
      p <- lapply(blocks, function(x) unname(x[i, ]))
      eta[i, ] <- .msr_bayes_components(spec, p, include_response = FALSE)$eta
    }
    return(eta)
  }
  eta <- blocks$beta[, seq_len(ncol(spec$X_fixed)), drop = FALSE] %*% t(spec$X_fixed)
  eta <- sweep(eta, 2L, spec$offset, "+")
  for (j in seq_along(spec$covariates)) {
    nm <- spec$covariates[j]; sm <- spec$scales[[nm]]
    log_s <- blocks$scale_beta[, sm$beta_index, drop = FALSE] %*% t(sm$X)
    for (term in sm$random)
      log_s <- log_s + blocks$scale_z[, term$z_index[term$index], drop = FALSE] *
        exp(blocks$log_scale_sd[, term$sd_index])
    power <- if (spec$kernel == "gaussian") 2 else 1
    precision <- exp(-power * log_s) / if (power == 2) 2 else 1
    numerator <- denominator <- matrix(0, nrow(eta), ncol(eta))
    for (bin in seq_len(ncol(spec$delta[[nm]]))) {
      w <- exp(-sweep(precision, 2L, spec$delta[[nm]][, bin], "*"))
      numerator <- numerator + sweep(w, 2L, spec$V[[nm]][, bin], "*")
      denominator <- denominator + sweep(w, 2L, spec$C[[nm]][, bin], "*")
    }
    kernel <- (numerator / denominator - spec$center[[nm]]) / spec$spread[[nm]]
    eta <- eta + kernel * blocks$beta[, ncol(spec$X_fixed) + j]
  }
  eta
}

.msr_bayes_response_geometry <- function(spec, blocks) {
  term <- spec$response$random[[1L]]
  sigma <- exp(blocks$log_response_sd[, term$sd_index])
  base <- .msr_bayes_batch_eta(spec, blocks)
  n <- nrow(base); groups <- length(term$levels)
  y <- matrix(rep(spec$y, each = n), n)
  trials <- matrix(rep(spec$trials, each = n), n)
  aggregate_sites <- function(value) t(rowsum(t(value), term$index, reorder = TRUE))
  z <- matrix(0, n, groups)
  fallback <- rep(FALSE, n)
  evaluate <- function(z) {
    eta <- base + z[, term$index, drop = FALSE] * sigma
    if (spec$family == "gaussian") {
      variance <- exp(2 * blocks$log_dispersion[, 1L])
      score <- (y - eta) / variance
      curvature <- matrix(1, n, ncol(base)) / variance
    } else if (spec$family == "poisson") {
      mu <- exp(eta)
      score <- y - mu; curvature <- mu
    } else if (spec$family == "binomial") {
      p <- stats::plogis(eta)
      score <- y - trials * p; curvature <- trials * p * (1 - p)
    } else {
      theta <- exp(blocks$log_dispersion[, 1L])
      p <- stats::plogis(eta - blocks$log_dispersion[, 1L])
      score <- y - (y + theta) * p
      curvature <- (y + theta) * p * (1 - p)
    }
    list(gradient = z - aggregate_sites(score) * sigma,
         hessian = 1 + aggregate_sites(curvature) * sigma^2)
  }
  for (step in seq_len(20L)) {
    geometry <- evaluate(z)
    bad <- rowSums(!is.finite(geometry$gradient) | !is.finite(geometry$hessian) |
      geometry$hessian <= 0) > 0
    fallback <- fallback | bad
    if (any(fallback)) {
      z[fallback, ] <- 0
      geometry$gradient[fallback, ] <- 0
      geometry$hessian[fallback, ] <- 1
    }
    direction <- geometry$gradient / geometry$hessian
    if (max(abs(direction)) <= 1e-8) break
    # Strictly concave conditional likelihoods give positive curvature. Bound
    # each Newton update to prevent an extreme initial log-link prediction.
    z <- z - pmax(-1, pmin(1, direction))
  }
  geometry <- evaluate(z)
  fallback <- fallback | rowSums(!is.finite(geometry$gradient) |
    !is.finite(geometry$hessian) | geometry$hessian <= 0) > 0
  if (any(fallback)) {
    z[fallback, ] <- 0
    geometry$gradient[fallback, ] <- 0
    geometry$hessian[fallback, ] <- 1
  }
  # Stationarity is useful for proposal efficiency, but is not required to
  # define this normalized proposal. Exact importance ratios still correct it.
  list(mode = z, hessian = geometry$hessian, fallbacks = sum(fallback),
       fallback_rows = fallback,
       max_step = max(abs(geometry$gradient / geometry$hessian)))
}

.msr_bayes_response_candidates <- function(outer, indices, spec) {
  samples <- matrix(0, nrow(outer), sum(lengths(spec$parameters)))
  samples[, indices] <- outer
  random <- setdiff(seq_len(ncol(samples)), indices)
  logq <- numeric(nrow(samples))
  fallbacks <- 0L; max_step <- 0
  for (begin in seq.int(1L, nrow(samples), by = 256L)) {
    rows <- seq.int(begin, min(begin + 255L, nrow(samples)))
    blocks <- .msr_bayes_parameter_blocks(samples[rows, , drop = FALSE], spec$parameters)
    geometry <- .msr_bayes_response_geometry(spec, blocks)
    fallbacks <- fallbacks + geometry$fallbacks
    max_step <- max(max_step, geometry$max_step)
    standard <- matrix(stats::rt(length(rows) * length(random), df = 20), length(rows))
    if (any(geometry$fallback_rows)) standard[geometry$fallback_rows, ] <-
      stats::rnorm(sum(geometry$fallback_rows) * length(random))
    samples[rows, random] <- geometry$mode + standard / sqrt(geometry$hessian)
    logq[rows] <- rowSums(stats::dt(standard, df = 20, log = TRUE) + .5 * log(geometry$hessian))
    if (any(geometry$fallback_rows)) logq[rows[geometry$fallback_rows]] <- rowSums(
      stats::dnorm(standard[geometry$fallback_rows, , drop = FALSE], log = TRUE))
  }
  list(samples = samples, logq = logq, fallbacks = fallbacks, max_step = max_step)
}

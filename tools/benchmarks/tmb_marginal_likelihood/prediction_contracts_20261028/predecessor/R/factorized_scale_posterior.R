# Private conditional reconstruction for independent annual log-scale effects.
# Population parameters are supplied, not estimated or sampled by these helpers.
.msr_factorized_distribution <- function(log_density, log_likelihood_bound,
                                          knots = numeric(), rel.tol = 1e-10,
                                          step = .5) {
  if (!is.numeric(knots) || any(!is.finite(knots)) ||
      length(rel.tol) != 1L || !is.finite(rel.tol) || rel.tol <= 0 || rel.tol >= 1 ||
      length(step) != 1L || !is.finite(step) || step <= 0)
    stop("Finite knots and positive integration controls are required.", call. = FALSE)
  search <- .msr_factorized_modes(log_density, log_likelihood_bound)
  # Physical log-distance knots help resolve response transitions when a large
  # population SD compresses them into a narrow standardized-effect interval.
  knots <- knots[abs(knots) < search$radius]
  grid <- sort(unique(c(seq(-search$radius, search$radius, by = step), knots)))
  values <- vapply(grid, log_density, 0)
  if (anyNA(values) || any(values == Inf))
    stop("Invalid conditional log density.", call. = FALSE)
  peaks <- which(values >= c(-Inf, utils::head(values, -1L)) &
    values >= c(utils::tail(values, -1L), -Inf))
  extra <- lapply(peaks, function(j) stats::optimize(function(z) -log_density(z),
    c(grid[max(1L, j - 1L)], grid[min(length(grid), j + 1L)]), tol = 1e-10))
  candidates <- c(search$candidates, extra)
  shift <- -min(vapply(candidates, `[[`, 0, "objective"))
  cuts <- sort(unique(c(-Inf, grid,
    vapply(candidates, `[[`, 0, "minimum"), Inf)))
  density <- function(z) vapply(z, function(x) exp(log_density(x) - shift), 0)
  integral <- function(lo, hi) {
    if (lo == hi) return(list(value = 0, abs.error = 0))
    stats::integrate(density, lo, hi, rel.tol = rel.tol, abs.tol = 0,
      subdivisions = 1000L, stop.on.error = TRUE)
  }
  piece <- function(lo, hi) integral(lo, hi)$value
  pieces <- lapply(seq_len(length(cuts) - 1L), function(j) integral(cuts[j], cuts[j + 1L]))
  mass <- vapply(pieces, `[[`, 0, "value")
  errors <- vapply(pieces, `[[`, 0, "abs.error")
  total <- sum(mass)
  if (!all(is.finite(c(mass, errors))) || any(mass < 0) ||
      !is.finite(total) || total <= 0 || sum(errors) > rel.tol * total)
    stop("Invalid conditional normalization.", call. = FALSE)
  left <- c(0, cumsum(mass))
  right <- c(rev(cumsum(rev(mass))), 0)
  lower_probability <- left / total
  upper_probability <- right / total
  lower_probability[c(1L, length(lower_probability))] <- c(0, 1)
  upper_probability[c(1L, length(upper_probability))] <- c(1, 0)
  cdf <- function(q, lower.tail = TRUE) {
    if (!is.numeric(q) || anyNA(q) || !is.logical(lower.tail) ||
        length(lower.tail) != 1L || is.na(lower.tail))
      stop("Numeric nonmissing quantiles and one logical tail flag are required.", call. = FALSE)
    vapply(q, function(x) {
      if (x == -Inf) return(if (lower.tail) 0 else 1)
      if (x == Inf) return(if (lower.tail) 1 else 0)
      j <- findInterval(x, cuts, all.inside = TRUE)
      value <- if (lower.tail) (left[j] + piece(cuts[j], x)) / total else
        (right[j + 1L] + piece(x, cuts[j + 1L])) / total
      if (!is.finite(value) || value < -1e-12 || value > 1 + 1e-12)
        stop("Invalid conditional probability.", call. = FALSE)
      # Only probability roundoff is clipped; tails are integrated to infinity.
      max(0, min(1, value))
    }, 0)
  }
  quantile <- function(p) {
    if (!is.numeric(p) || anyNA(p) || any(p < 0 | p > 1))
      stop("Probabilities must be nonmissing numbers between zero and one.", call. = FALSE)
    vapply(p, function(prob) {
      if (prob == 0) return(-Inf)
      if (prob == 1) return(Inf)
      boundary <- if (prob <= .5) which(lower_probability >= prob)[1L] else
        which(upper_probability <= 1 - prob)[1L]
      j <- max(1L, min(length(mass), boundary - 1L))
      lo <- cuts[j]; hi <- cuts[j + 1L]
      if (!is.finite(lo)) {
        lo <- -search$radius
        while (cdf(lo) > prob) {
          lo <- 2 * lo
          if (!is.finite(lo)) stop("Conditional lower quantile cannot be bracketed.", call. = FALSE)
        }
      }
      if (!is.finite(hi)) {
        hi <- search$radius
        while (cdf(hi, FALSE) > 1 - prob) {
          hi <- 2 * hi
          if (!is.finite(hi)) stop("Conditional upper quantile cannot be bracketed.", call. = FALSE)
        }
      }
      target <- if (prob <= .5) function(z) cdf(z) - prob else
        function(z) cdf(z, FALSE) - (1 - prob)
      stats::uniroot(target, c(lo, hi), tol = 1e-9, maxiter = 1000L)$root
    }, 0)
  }
  list(cdf = cdf, quantile = quantile, log_normalizer = log(total) + shift,
    cuts = cuts, interval_probability = mass / total,
    relative_error = sum(errors) / total, infinite_support = TRUE, rel.tol = rel.tol)
}

.msr_factorized_annual_conditional <- function(spec, parameters, group,
                                               rel.tol = 1e-10, cache = NULL) {
  .msr_factorized_validate(spec)
  term <- spec$scales[[1L]]$random
  if (length(term) != 1L || length(group) != 1L || !is.numeric(group) ||
      !is.finite(group) || group != as.integer(group) ||
      group < 1L || group > length(term[[1L]]$levels))
    stop("Select one existing annual random-effect group.", call. = FALSE)
  tau <- exp(parameters$log_scale_sd[term[[1L]]$sd_index])
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0 || length(parameters$scale_z))
    stop("Supply population parameters with positive finite annual SD and no latent annual effects.", call. = FALSE)
  g <- .msr_factorized_group(spec, group)
  log_density <- function(z) {
    p <- parameters; p$scale_z <- z
    -sum(.msr_factorized_response_nll(g, p)) + stats::dnorm(z, log = TRUE)
  }
  if (!is.null(cache)) {
    if (!inherits(cache, "msr_factorized_conditional_cache") ||
        !identical(cache$specification, spec))
      stop("Conditional density cache must match the prepared specification.", call. = FALSE)
    log_density <- cache$log_density(parameters, group)
  }
  bound <- if (spec$family == "gaussian")
    -length(g$y) * (parameters$log_dispersion + .5 * log(2 * pi)) else 0
  sm <- g$scales[[1L]]
  locations <- unique(drop(sm$X %*% parameters$scale_beta[sm$beta_index]))
  knots <- as.vector(outer(seq(-30, 30, by = 1), locations, "-") / tau)
  .msr_factorized_distribution(log_density, bound, knots, rel.tol)
}

# A separately shifted ordinary-double reference. This does not alter the
# optimizing objective, its derivatives, or the production auditing default.
.msr_factorized_shifted_reference <- function(spec, parameters, rel.tol = 1e-10) {
  .msr_factorized_validate(spec)
  parameters$scale_z <- numeric()
  term <- spec$scales[[1L]]$random
  if (!length(term)) return(list(value = sum(.msr_factorized_response_nll(spec, parameters)),
    log_group_normalizers = numeric(), relative_error = numeric(), infinite_support = TRUE))
  groups <- lapply(seq_along(term[[1L]]$levels), function(i)
    .msr_factorized_annual_conditional(spec, parameters, i, rel.tol))
  normalizers <- vapply(groups, `[[`, 0, "log_normalizer")
  list(value = -sum(normalizers), log_group_normalizers = normalizers,
    relative_error = vapply(groups, `[[`, 0, "relative_error"), infinite_support = TRUE)
}

# Inverse-CDF reconstruction is deterministic for supplied uniforms. This is a
# conditional distribution at one population vector, not a posterior sampler.
# Integrating population uncertainty requires valid draws from that posterior.
.msr_factorized_reconstruct <- function(spec, parameters, uniforms = NULL,
                                         probs = c(.025, .5, .975), rel.tol = 1e-10,
                                         cache = NULL) {
  .msr_factorized_validate(spec)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) || any(probs < 0 | probs > 1))
    stop("Quantile probabilities must lie between zero and one.", call. = FALSE)
  if (length(parameters$scale_z))
    stop("Omit integrated latent annual effects from population parameters.", call. = FALSE)
  sm <- spec$scales[[1L]]; term <- sm$random; rows <- sm$representatives
  fixed <- drop(sm$X[rows, , drop = FALSE] %*% parameters$scale_beta[sm$beta_index]) + log(spec$unit_conv)
  groups <- if (length(term)) length(term[[1L]]$levels) else 0L
  if (!is.null(uniforms) && (!is.numeric(uniforms) || length(uniforms) != groups ||
      anyNA(uniforms) || any(uniforms <= 0 | uniforms >= 1)))
    stop("Supply one open-unit-interval uniform per annual random-effect group.", call. = FALSE)
  output <- list(covariate = spec$covariates, labels = sm$labels, rows = rows, probs = probs,
    log_quantiles = matrix(rep(fixed, each = length(probs)), nrow = length(probs)),
    log_draw = if (is.null(uniforms)) NULL else fixed,
    latent_draw = if (is.null(uniforms)) NULL else numeric(groups),
    conditional_only = TRUE)
  normalizers <- numeric(groups)
  if (groups) for (i in seq_len(groups)) {
    conditional <- .msr_factorized_annual_conditional(spec, parameters, i, rel.tol, cache)
    selected <- which(term[[1L]]$index[rows] == i)
    tau <- exp(parameters$log_scale_sd[term[[1L]]$sd_index])
    output$log_quantiles[, selected] <- outer(conditional$quantile(probs) * tau,
      fixed[selected], "+")
    normalizers[i] <- conditional$log_normalizer
    if (!is.null(uniforms)) {
      output$latent_draw[i] <- conditional$quantile(uniforms[i])
      output$log_draw[selected] <- fixed[selected] + tau * output$latent_draw[i]
    }
  }
  output$scale_quantiles <- exp(output$log_quantiles)
  output$scale_draw <- if (is.null(output$log_draw)) NULL else exp(output$log_draw)
  output$log_group_normalizers <- normalizers
  output
}

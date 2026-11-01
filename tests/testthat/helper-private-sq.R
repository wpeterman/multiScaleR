# Neutral, tiny software fixtures. No SQ fit admission or scientific cohort.
private_sq_fixture <- function(kernel = "gaussian", family = "gaussian", response = FALSE,
                               annual = TRUE) {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3, annual = FALSE, kernel = kernel)
  d <- f$inputs$binned$bin_rep
  ids <- f$inputs$binned$point_ids
  f$inputs$d_list <- rep(list(d), length(ids))
  f$inputs$raw_cov <- lapply(seq_along(ids), function(i)
    data.frame(elev = sin(d * (i + 1)) + i * d^2))
  f$inputs$raw_cov[[1]]$elev[] <- 2
  for (i in 2:4) {
    f$inputs$d_list[[i]][2] <- d[1]
    f$inputs$raw_cov[[i]]$elev[1:2] <- c(0, i)
  }
  f$inputs$raw_cov[[2]]$elev[3] <- NA_real_
  f$inputs$kernel_dat <- data.frame(row.names = ids)
  f$inputs$scale_vars <- msr_vars(rough = surface_var("elev", "sq", weighted = TRUE))
  f$inputs$binned <- NULL
  if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
  f$data$batch <- interaction(f$data$year, as.integer(f$data$site) <= 2)
  formula <- if (response) count ~ rough + temp + (1 | batch) + offset(log(effort)) else
    count ~ rough + temp + offset(log(effort))
  scales <- list(rough = if (annual) ~1 + (1 | year) else ~1)
  spec <- .msr_scale_prepare(formula, f$inputs, f$data, scales, family, "key",
    reference_scale = 250)
  list(spec = spec, inputs = f$inputs, data = f$data, formula = formula, scales = scales)
}


# Independent base-R reference, using retained cell values and physical units.
private_sq_reference <- function(s, p, response = TRUE) {
  sm <- s$scales[[1L]]
  log_scale <- drop(sm$X %*% p$scale_beta[sm$beta_index])
  for (term in sm$random) {
    stopifnot(is.null(term$gaussian_block))
    log_scale <- log_scale + exp(p$log_scale_sd[term$sd_index]) *
      p$scale_z[term$z_index][term$index]
  }
  raw <- vapply(seq_along(s$y), function(i) {
    b <- s$surface$buffers$rough[[i]]
    sigma <- exp(log_scale[i])
    w <- if (s$kernel == "gaussian")
      exp(-(b$d^2 - min(b$d)^2) / (2 * sigma^2)) else
      exp(-(b$d - min(b$d)) / sigma)
    v <- b$v - b$v[1L]
    sqrt(sum(w * (v - sum(w * v) / sum(w))^2) / sum(w))
  }, 0)
  k <- (raw - s$center[[1L]]) / s$spread[[1L]]
  eta <- drop(s$X_fixed %*% p$beta[seq_len(ncol(s$X_fixed))]) + s$offset +
    p$beta[ncol(s$X_fixed) + 1L] * k
  if (response) for (term in s$response$random) {
    stopifnot(is.null(term$gaussian_block))
    eta <- eta + exp(p$log_response_sd[term$sd_index]) *
      p$response_z[term$z_index][term$index]
  }
  list(eta = eta, raw = raw, kernel = k, physical = exp(log_scale) * s$unit_conv)
}

private_sq_parameters <- function(s) {
  p <- s$parameters; p$beta <- c(.4, -.2, .7)
  p$scale_beta[] <- log(.35)
  if (length(p$scale_z)) p$scale_z[] <- c(-.4, .1, .5)
  p$log_scale_sd[] <- log(.3)
  p$log_response_sd[] <- log(.25)
  if (length(p$response_z)) p$response_z[] <- seq(-.2, .2, length.out = length(p$response_z))
  if (length(p$log_dispersion)) p$log_dispersion[] <- if (s$family == "gaussian") log(.8) else log(3)
  p
}

private_sq_unpack <- function(x, template) {
  n <- lengths(template); at <- cumsum(c(0L, n))
  for (i in seq_along(template)) template[[i]] <- if (n[i]) x[at[i] + seq_len(n[i])] else numeric()
  template
}

private_sq_joint_reference <- function(s, p) {
  eta <- private_sq_reference(s, p)$eta
  ll <- switch(s$family,
    gaussian = dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE),
    poisson = dpois(s$y, exp(eta), log = TRUE),
    binomial = dbinom(s$y, s$trials, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(s$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE))
  -sum(ll) - sum(dnorm(p$scale_z, log = TRUE)) - sum(dnorm(p$response_z, log = TRUE))
}

private_sq_fd_gradient <- function(fn, x) vapply(seq_along(x), function(j) {
  a <- b <- x; a[j] <- a[j] + 1e-5; b[j] <- b[j] - 1e-5
  (fn(a) - fn(b)) / 2e-5
}, 0)

# Richardson extrapolation limits truncation error without using AD scores.
private_sq_fd_hessian <- function(fn, x) {
  central <- function(step) {
    out <- matrix(0, length(x), length(x)); base <- fn(x)
    h <- step * (1 + abs(x))
    for (j in seq_along(x)) for (k in seq_len(j)) {
      a <- b <- c <- d <- x
      if (j == k) {
        a[j] <- a[j] + h[j]; b[j] <- b[j] - h[j]
        v <- (fn(a) - 2 * base + fn(b)) / h[j]^2
      } else {
        a[j] <- a[j] + h[j]; a[k] <- a[k] + h[k]
        b[j] <- b[j] + h[j]; b[k] <- b[k] - h[k]
        c[j] <- c[j] - h[j]; c[k] <- c[k] + h[k]
        d[j] <- d[j] - h[j]; d[k] <- d[k] - h[k]
        v <- (fn(a) - fn(b) - fn(c) + fn(d)) / (4 * h[j] * h[k])
      }
      out[j, k] <- out[k, j] <- v
    }
    out
  }
  (4 * central(5e-4) - central(1e-3)) / 3
}

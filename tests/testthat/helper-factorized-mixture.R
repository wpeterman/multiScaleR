# Limiting group-density arithmetic only. No fitted finite-scale candidate.
mixture_fixture <- function(kernel, family, surface, response = FALSE) {
  if (surface) return(private_sq_fixture(kernel, family, response = response))
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3, annual = FALSE, kernel = kernel)
  if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
  f$data$batch <- interaction(f$data$year, as.integer(f$data$site) <= 2)
  form <- if (response) count ~ habitat + temp + (1 | batch) + offset(log(effort)) else
    count ~ habitat + temp + offset(log(effort))
  spec <- .msr_scale_prepare(form, f$inputs, f$data, list(habitat = ~1 + (1 | year)),
    family, "key", reference_scale = 250)
  list(spec = spec)
}

# Native density leaves, independent of the limiting/objective implementation.
mixture_native_nll <- function(s, p, kernel) {
  if (length(s$surface$buffers)) {
    raw <- vapply(s$surface$buffers[[1L]], function(b) {
      v <- if (kernel == "nearest") b$v[b$d == min(b$d)] else b$v
      sqrt(mean((v - mean(v))^2))
    }, 0)
  } else {
    keep <- if (kernel == "nearest") s$delta[[1L]] == 0 else matrix(TRUE, nrow(s$V[[1L]]), ncol(s$V[[1L]]))
    raw <- rowSums(s$V[[1L]] * keep) / rowSums(s$C[[1L]] * keep)
  }
  eta <- drop(s$X_fixed %*% p$beta[seq_len(ncol(s$X_fixed))]) + s$offset +
    p$beta[ncol(s$X_fixed) + 1L] * (raw - s$center[[1L]]) / s$spread[[1L]]
  index <- s$scales[[1L]]$random[[1L]]$index
  groups <- split(seq_along(s$y), index)
  if (length(s$response$random)) {
    r <- s$response$random[[1L]]$index
    return(vapply(groups, function(rows) {
      cov <- diag(exp(2 * p$log_dispersion), length(rows)) +
        outer(r[rows], r[rows], "==") * exp(2 * p$log_response_sd)
      ch <- chol(cov)
      .5 * (length(rows) * log(2 * pi) + 2 * sum(log(diag(ch))) +
        sum(backsolve(ch, s$y[rows] - eta[rows], transpose = TRUE)^2))
    }, 0))
  }
  leaf <- -switch(s$family,
    gaussian = dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE),
    poisson = dpois(s$y, exp(eta), log = TRUE),
    binomial = dbinom(s$y, s$trials, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(s$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE))
  vapply(groups, function(rows) sum(leaf[rows]), 0)
}

mixture_native_total <- function(a, b, weight) {
  if (weight == 0) return(sum(b))
  if (weight == 1) return(sum(a))
  -sum(vapply(seq_along(a), function(i) {
    values <- c(log(weight) - a[i], log1p(-weight) - b[i])
    shift <- max(values)
    shift + log(sum(exp(values - shift)))
  }, 0))
}


# Higher-accuracy marginal integration for independent annual scale effects.
# Population coefficients and SDs remain likelihood parameters, without priors.
.msr_factorized_validate <- function(spec) {
  .msr_scale_check_capability(spec, "marginal_likelihood", "factorized")
}

.msr_factorized_group <- function(spec, group) {
  term <- spec$scales[[1L]]$random[[1L]]
  rows <- which(term$index == group)
  for (nm in c("X_fixed", "y", "offset", "trials"))
    spec[[nm]] <- if (is.matrix(spec[[nm]])) spec[[nm]][rows, , drop = FALSE] else spec[[nm]][rows]
  for (nm in c("V", "C", "delta")) spec[[nm]] <- lapply(spec[[nm]], function(x) x[rows, , drop = FALSE])
  if (length(spec$surface$buffers))
    spec$surface$buffers <- lapply(spec$surface$buffers, function(x) x[rows])
  spec$scales[[1L]]$X <- spec$scales[[1L]]$X[rows, , drop = FALSE]
  term$index <- rep(1L, length(rows)); term$z_index <- 1L
  term$levels <- term$levels[group]
  spec$scales[[1L]]$random[[1L]] <- term
  for (j in seq_along(spec$response$random)) {
    response <- spec$response$random[[j]]
    observed <- sort(unique(response$index[rows]))
    response$index <- match(response$index[rows], observed)
    response$levels <- response$levels[observed]
    response$z_index <- seq_along(observed)
    spec$response$random[[j]] <- response
  }
  spec
}

# Integrate Gaussian response intercepts exactly at fixed kernel scales.
# Centering residuals avoids subtracting two large quadratic forms. The
# determinant includes every normalizing constant and no outer-SD Jacobian.
.msr_factorized_gaussian_response_nll <- function(spec, eta, parameters, ad = FALSE) {
  term <- spec$response$random[[1L]]
  log_v <- 2 * parameters$log_dispersion
  log_w <- 2 * parameters$log_response_sd[term$sd_index]
  logadd <- if (ad) RTMB::logspace_add else function(a, b)
    max(a, b) + log1p(exp(-abs(a - b)))
  answer <- 0
  for (group in unique(term$index)) {
    rows <- which(term$index == group); n <- length(rows)
    residual <- spec$y[rows] - eta[rows]
    average <- sum(residual) / n
    log_total <- logadd(log_v, log(n) + log_w)
    quadratic <- sum((residual - average)^2) * exp(-log_v) +
      n * average^2 * exp(-log_total)
    answer <- answer + .5 * (n * log(2 * pi) + (n - 1) * log_v + log_total + quadratic)
  }
  answer
}

# The nearest supported bin always has weight one. Saturation below the
# floating-point weight range avoids overflow in first and second derivatives.
.msr_factorized_response_nll <- function(spec, parameters, ad = FALSE) {
  if (length(spec$surface$buffers)) {
    eta <- .msr_bayes_components(spec, parameters, include_response = FALSE)$eta
    if (length(spec$response$random)) {
      if (spec$family != "gaussian") stop("Gaussian response integration required.", call. = FALSE)
      return(.msr_factorized_gaussian_response_nll(spec, eta, parameters, ad))
    }
    return(-.msr_scale_response_loglik(spec$family, spec$y, eta,
      parameters$log_dispersion, spec$trials, ad))
  }
  sm <- spec$scales[[1L]]
  log_s <- drop(sm$X %*% parameters$scale_beta[sm$beta_index])
  for (term in sm$random)
    log_s <- log_s + exp(parameters$log_scale_sd[term$sd_index]) *
      parameters$scale_z[term$z_index][term$index]
  power <- if (spec$kernel == "gaussian") 2 else 1
  log_exponent <- log(spec$delta[[1L]]) - power * log_s - if (power == 2) log(2) else 0
  capped <- if (ad) RTMB::Vectorize(function(x) {
    "if" <- RTMB::ADoverload("if")
    if (x > log(745)) RTMB::advector(log(745)) else x
  })(log_exponent) else pmin(log_exponent, log(745))
  dim(capped) <- dim(log_exponent)
  weights <- exp(-exp(capped))
  sum_rows <- if (ad) RTMB::rowSums else base::rowSums
  kernel <- sum_rows(weights * spec$V[[1L]]) / sum_rows(weights * spec$C[[1L]])
  kernel <- (kernel - spec$center[[1L]]) / spec$spread[[1L]]
  eta <- drop(spec$X_fixed %*% parameters$beta[seq_len(ncol(spec$X_fixed))]) +
    spec$offset + parameters$beta[ncol(spec$X_fixed) + 1L] * kernel
  if (length(spec$response$random)) {
    if (spec$family != "gaussian") stop("Gaussian response integration required.", call. = FALSE)
    return(.msr_factorized_gaussian_response_nll(spec, eta, parameters, ad))
  }
  -.msr_scale_response_loglik(spec$family, spec$y, eta,
    parameters$log_dispersion, spec$trials, ad)
}

.msr_factorized_objective <- function(spec, ad = FALSE, rel.tol = 1e-8, step = 1) {
  .msr_factorized_validate(spec)
  term <- spec$scales[[1L]]$random
  if (!length(term)) return(function(parameters) sum(.msr_factorized_response_nll(spec, parameters, ad)))
  groups <- lapply(seq_along(term[[1L]]$levels), function(i) .msr_factorized_group(spec, i))
  cuts <- c(-Inf, seq(-8, 8, by = step), Inf)
  normal <- if (ad) RTMB::dnorm else stats::dnorm
  integrate <- if (ad) RTMB::integrate else stats::integrate
  vectorize <- if (ad) RTMB::Vectorize else base::Vectorize
  function(parameters) {
    answer <- 0
    for (g in groups) {
      # This upper bound keeps the integrand <=1, including Gaussian densities
      # above one. Ordinary integration retains a purely relative tolerance.
      # The AD integrator needs a tiny absolute tolerance to avoid undefined
      # derivatives when a whole subdivision underflows. Cold ordinary checks
      # still require the original relative accuracy and objective agreement.
      # All normal-prior tails are integrated, not truncated.
      shift <- -.5 * log(2 * pi)
      if (spec$family == "gaussian") shift <- shift -
        length(g$y) * (parameters$log_dispersion + .5 * log(2 * pi))
      integrand <- vectorize(function(z) {
        p <- parameters; p$scale_z <- z
        exp(-sum(.msr_factorized_response_nll(g, p, ad)) + normal(z, log = TRUE) - shift)
      })
      value <- 0
      for (j in seq_len(length(cuts) - 1L)) value <- value +
        integrate(integrand, cuts[j], cuts[j + 1L], rel.tol = rel.tol,
          abs.tol = if (ad) 1e-100 else 0, subdivisions = 500L, stop.on.error = TRUE)$value
      answer <- answer - log(value) - shift
    }
    answer
  }
}

.msr_factorized_adfun <- function(spec, parameters = spec$parameters, rel.tol = 1e-8, step = 1) {
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Factorized likelihood requires the optional package 'RTMB'.", call. = FALSE)
  .msr_factorized_validate(spec)
  parameters$scale_z <- numeric()
  parameters$response_z <- numeric()
  RTMB::MakeADFun(.msr_factorized_objective(spec, TRUE, rel.tol, step), parameters, silent = TRUE)
}

# An independent coordinate system for very large annual SDs. The middle
# intervals resolve the response transition in log-scale units; shifted normal
# coordinates resolve both infinite tails. The Jacobian is included exactly.
# At very small SDs the standard-normal coordinate system is preferable.
.msr_factorized_scale_objective <- function(spec, ad = FALSE, rel.tol = 1e-10, step = 2) {
  .msr_factorized_validate(spec)
  term <- spec$scales[[1L]]$random
  if (!length(term)) return(.msr_factorized_objective(spec, ad, rel.tol))
  groups <- lapply(seq_along(term[[1L]]$levels), function(i) .msr_factorized_group(spec, i))
  normal <- if (ad) RTMB::dnorm else stats::dnorm
  integrate <- if (ad) RTMB::integrate else stats::integrate
  vectorize <- if (ad) RTMB::Vectorize else base::Vectorize
  function(parameters) {
    tau <- exp(parameters$log_scale_sd)
    answer <- 0
    for (g in groups) {
      mu <- sum(drop(g$scales[[1L]]$X %*% parameters$scale_beta)) / length(g$y)
      shift <- -.5 * log(2 * pi)
      if (spec$family == "gaussian") shift <- shift -
        length(g$y) * (parameters$log_dispersion + .5 * log(2 * pi))
      density <- function(z) {
        p <- parameters; p$scale_z <- z
        exp(-sum(.msr_factorized_response_nll(g, p, ad)) + normal(z, log = TRUE) - shift)
      }
      cuts <- seq(-30, 30, by = step)
      central <- vectorize(function(v) density((v - mu) / tau) / tau)
      piece <- function(f, lo, hi) integrate(f, lo, hi, rel.tol = rel.tol,
        abs.tol = if (ad) 1e-100 else 0, subdivisions = 1000L, stop.on.error = TRUE)$value
      value <- 0
      for (j in seq_len(length(cuts) - 1L)) value <- value + piece(central, cuts[j], cuts[j + 1L])
      left <- vectorize(function(t) density((-30 - mu) / tau + t))
      right <- vectorize(function(t) density((30 - mu) / tau + t))
      value <- value + piece(left, -Inf, 0) + piece(right, 0, Inf)
      answer <- answer - log(value) - shift
    }
    answer
  }
}

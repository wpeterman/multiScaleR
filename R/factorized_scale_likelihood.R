# Higher-accuracy marginal integration for independent annual scale effects.
# Population coefficients and SDs remain likelihood parameters, without priors.
.msr_factorized_validate <- function(spec) {
  .msr_scale_check_capability(spec, "marginal_likelihood", "factorized")
}

.msr_factorized_group <- function(spec, group) {
  term <- spec$scales[[1L]]$random[[1L]]
  rows <- which(term$index == group)
  # Keep observation provenance aligned with the actual group payload. `rows`
  # in the specification are source-point indices, not local observation IDs.
  for (nm in c("rows", "observation_ids"))
    if (length(spec[[nm]])) spec[[nm]] <- spec[[nm]][rows]
  if (!is.null(spec$data)) spec$data <- spec$data[rows, , drop = FALSE]
  for (nm in c("X_fixed", "y", "offset", "trials"))
    spec[[nm]] <- if (is.matrix(spec[[nm]])) spec[[nm]][rows, , drop = FALSE] else spec[[nm]][rows]
  for (nm in c("V", "C", "delta")) spec[[nm]] <- lapply(spec[[nm]], function(x) x[rows, , drop = FALSE])
  if (!is.null(spec$response_design))
    spec$response_design <- .msr_response_design_rows(spec$response_design, rows)
  if (length(spec$surface$buffers)) {
    spec$surface$buffers <- lapply(spec$surface$buffers, function(x) x[rows])
    spec$surface$memory <- .msr_scale_surface_size(spec$surface)
  }
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
  eta <- if (is.null(spec$response_design))
    drop(spec$X_fixed %*% parameters$beta[seq_len(ncol(spec$X_fixed))]) +
      spec$offset + parameters$beta[ncol(spec$X_fixed) + 1L] * kernel else
    .msr_response_design_eta(spec$response_design, list(kernel), parameters$beta) + spec$offset
  if (length(spec$response$random)) {
    if (spec$family != "gaussian") stop("Gaussian response integration required.", call. = FALSE)
    return(.msr_factorized_gaussian_response_nll(spec, eta, parameters, ad))
  }
  -.msr_scale_response_loglik(spec$family, spec$y, eta,
    parameters$log_dispersion, spec$trials, ad)
}

# One-group likelihood arithmetic, not a fitting/admission interface. The
# original normalized integration algorithm is shared by the guarded wrapper.
.msr_factorized_group_check <- function(group, ad = FALSE) {
  if (length(group$covariates) != 1L || length(group$scales) != 1L ||
      length(group$scales[[1L]]$random) != 1L)
    stop("Group integration requires one kernel and one scale grouping term.", call. = FALSE)
  term <- group$scales[[1L]]$random[[1L]]
  if (!length(group$y) || length(term$levels) != 1L ||
      length(term$index) != length(group$y) || anyNA(term$index) ||
      any(term$index != 1L) || !identical(as.integer(term$z_index), 1L) ||
      !is.null(term$gaussian_block))
    stop("Supply one recoded iid scale group, not a complete multi-group model.", call. = FALSE)
  if (!group$family %in% c("gaussian", "poisson", "binomial", "nbinom2") ||
      length(group$response$random) > 1L ||
      (length(group$response$random) && (group$family != "gaussian" ||
        !is.null(group$response$random[[1L]]$gaussian_block))))
    stop("Unsupported group response structure for normalized quadrature.", call. = FALSE)
  if ((length(group$rows) && length(group$rows) != length(group$y)) ||
      (length(group$observation_ids) && length(group$observation_ids) != length(group$y)) ||
      (!is.null(group$data) && nrow(group$data) != length(group$y)))
    stop("Group observation metadata does not match its density payload.", call. = FALSE)
  if (ad) .msr_scale_surface_ram(group)
  invisible(group)
}

.msr_factorized_group_integral <- function(group, ad = FALSE, rel.tol = 1e-8, step = 1) {
  .msr_factorized_group_check(group, ad)
  force(group)
  cuts <- c(-Inf, seq(-8, 8, by = step), Inf)
  normal <- if (ad) RTMB::dnorm else stats::dnorm
  integrate <- if (ad) RTMB::integrate else stats::integrate
  vectorize <- if (ad) RTMB::Vectorize else base::Vectorize
  function(parameters) {
    # The upper bound keeps even Gaussian observation integrands at most one.
    # Cold ordinary integration stays purely relative. AD's tiny absolute
    # tolerance handles underflowed pieces; both normal tails remain unbounded.
    shift <- -.5 * log(2 * pi)
    if (group$family == "gaussian") shift <- shift -
      length(group$y) * (parameters$log_dispersion + .5 * log(2 * pi))
    integrand <- vectorize(function(z) {
      p <- parameters; p$scale_z <- z
      exp(-sum(.msr_factorized_response_nll(group, p, ad)) + normal(z, log = TRUE) - shift)
    })
    value <- 0
    for (j in seq_len(length(cuts) - 1L)) value <- value +
      integrate(integrand, cuts[j], cuts[j + 1L], rel.tol = rel.tol,
        abs.tol = if (ad) 1e-100 else 0, subdivisions = 500L, stop.on.error = TRUE)$value
    list(value = value, shift = shift)
  }
}

.msr_factorized_objective <- function(spec, ad = FALSE, rel.tol = 1e-8, step = 1) {
  .msr_factorized_validate(spec)
  if (ad) .msr_scale_surface_ram(spec)
  term <- spec$scales[[1L]]$random
  if (!length(term)) return(function(parameters) sum(.msr_factorized_response_nll(spec, parameters, ad)))
  groups <- lapply(seq_along(term[[1L]]$levels), function(i) .msr_factorized_group(spec, i))
  integrals <- lapply(groups, function(g) .msr_factorized_group_integral(g, ad, rel.tol, step))
  function(parameters) {
    answer <- 0
    for (integral in integrals) {
      part <- integral(parameters)
      answer <- answer - log(part$value) - part$shift
    }
    answer
  }
}

.msr_factorized_adfun <- function(spec, parameters = spec$parameters, rel.tol = 1e-8, step = 1) {
  if (!requireNamespace("RTMB", quietly = TRUE))
    stop("Factorized likelihood requires the optional package 'RTMB'.", call. = FALSE)
  .msr_factorized_validate(spec)
  .msr_scale_surface_ram(spec)
  parameters$scale_z <- numeric()
  parameters$response_z <- numeric()
  RTMB::MakeADFun(.msr_factorized_objective(spec, TRUE, rel.tol, step), parameters, silent = TRUE)
}

# An independent coordinate system for very large annual SDs. The middle
# intervals resolve the response transition in log-scale units; shifted normal
# coordinates resolve both infinite tails. The Jacobian is included exactly.
# At very small SDs the standard-normal coordinate system is preferable.
# Same normalized group density in physical log-scale coordinates. This
# component supplies no fitted model or finite-identification certificate.
.msr_factorized_group_scale_integral <- function(group, ad = FALSE, rel.tol = 1e-10, step = 2) {
  .msr_factorized_group_check(group, ad)
  force(group)
  normal <- if (ad) RTMB::dnorm else stats::dnorm
  integrate <- if (ad) RTMB::integrate else stats::integrate
  vectorize <- if (ad) RTMB::Vectorize else base::Vectorize
  function(parameters) {
    tau <- exp(parameters$log_scale_sd)
    mu <- sum(drop(group$scales[[1L]]$X %*% parameters$scale_beta)) / length(group$y)
    shift <- -.5 * log(2 * pi)
    if (group$family == "gaussian") shift <- shift -
      length(group$y) * (parameters$log_dispersion + .5 * log(2 * pi))
    density <- function(z) {
      p <- parameters; p$scale_z <- z
      exp(-sum(.msr_factorized_response_nll(group, p, ad)) + normal(z, log = TRUE) - shift)
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
    list(value = value, shift = shift)
  }
}

.msr_factorized_scale_objective <- function(spec, ad = FALSE, rel.tol = 1e-10, step = 2) {
  .msr_factorized_validate(spec)
  if (ad) .msr_scale_surface_ram(spec)
  term <- spec$scales[[1L]]$random
  if (!length(term)) return(.msr_factorized_objective(spec, ad, rel.tol))
  groups <- lapply(seq_along(term[[1L]]$levels), function(i) .msr_factorized_group(spec, i))
  integrals <- lapply(groups, function(g) .msr_factorized_group_scale_integral(g, ad, rel.tol, step))
  function(parameters) {
    answer <- 0
    for (integral in integrals) {
      part <- integral(parameters)
      answer <- answer - log(part$value) - part$shift
    }
    answer
  }
}

# Ordinary point evaluation only. Reuse the strict neutral full-model policy
# and original group measures; fitting/AD wrappers keep their own admission.
# A function value at supplied parameters is not an estimate or an interval.
.msr_factorized_density_component <- function(spec, coordinate = c("standard", "physical"),
    rel.tol = NULL, step = NULL) {
  .msr_factorized_limit_component_validate(spec)
  coordinate <- match.arg(coordinate)
  if (is.null(rel.tol)) rel.tol <- if (coordinate == "standard") 1e-8 else 1e-10
  if (is.null(step)) step <- if (coordinate == "standard") 1 else 2
  if (is.complex(rel.tol) || is.complex(step) || !is.null(dim(rel.tol)) || !is.null(dim(step)))
    stop("rel.tol and step must be finite positive numeric scalars.", call. = FALSE)
  validate_scalar_numeric(rel.tol, "rel.tol", positive = TRUE)
  validate_scalar_numeric(step, "step", positive = TRUE)
  sm <- spec$scales[[1L]]
  if (length(sm$random) && (is.null(.msr_factorized_location_direction(sm, spec$parameters)) ||
      nrow(sm$X) != length(spec$y)))
    stop("Density components require one canonical annual SD and a supported prepared scale design.", call. = FALSE)
  sizes <- lengths(spec$parameters)
  sizes[c("scale_z", "response_z")] <- 0L
  integrals <- if (length(sm$random)) lapply(seq_along(sm$random[[1L]]$levels), function(i) {
    group <- .msr_factorized_group(spec, i)
    if (coordinate == "standard") .msr_factorized_group_integral(group, FALSE, rel.tol, step) else
      .msr_factorized_group_scale_integral(group, FALSE, rel.tol, step)
  }) else list()
  force(spec)
  answer <- function(parameters) {
    if (!is.list(parameters) || is.null(names(parameters)) || anyDuplicated(names(parameters)) ||
        !identical(lengths(parameters), sizes) ||
        !all(vapply(parameters, function(x) is.numeric(x) && !is.complex(x) &&
          !is.object(x) && is.null(dim(x)), TRUE)) ||
        any(!is.finite(unlist(parameters, use.names = FALSE))))
      stop("Supply matching finite numeric population blocks without realized latent effects.", call. = FALSE)
    if (!length(integrals)) return(sum(.msr_factorized_response_nll(spec, parameters, FALSE)))
    value <- 0
    for (integral in integrals) {
      part <- integral(parameters)
      value <- value - log(part$value) - part$shift
    }
    value
  }
  attr(answer, "component") <- list(target = "normalized_full_marginal_likelihood",
    coordinate = coordinate, rel.tol = rel.tol, step = step, population_layout = sizes,
    finite_scale_estimate = FALSE, uncertainty_available = FALSE)
  answer
}

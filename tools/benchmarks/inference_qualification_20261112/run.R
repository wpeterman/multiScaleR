# Independent references for the public analytic scale backends.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, args[1L] %in% c("glm", "mixed", "bootstrap"))
phase <- args[1L]; output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(compile = FALSE, quiet = TRUE)
for (p in c("MASS", if (phase == "mixed") c("glmmTMB", "RTMB", "TMB")))
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing ", p)

families <- switch(phase, glm = c("gaussian", "poisson", "binomial", "nbinom2"),
  mixed = c("poisson", "nbinom2", "zero_inflated_poisson"), bootstrap = "nbinom2")
plan <- expand.grid(family = families, scenario = c("interior", "weak", "near_bound"),
  replicate = seq_len(switch(phase, glm = 8L, mixed = 6L, bootstrap = 4L)),
  stringsAsFactors = FALSE)
plan$seed <- switch(phase, glm = 930000L, mixed = 940000L, bootstrap = 950000L) +
  seq_len(nrow(plan)) * 100L
plan$id <- sprintf("%s-%s-%02d", plan$family, plan$scenario, plan$replicate)
grid <- exp(seq(log(.05), log(.8), length.out = 19L))
warning_log <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = unique(messages))
}
atomic <- function(x, file) {
  tmp <- tempfile(tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, file)) stop("Cannot publish ", file)
}
sources <- c(list.files("R", pattern = "[.]R$", full.names = TRUE),
  list.files("src", pattern = "[.](c|cpp|h)$", full.names = TRUE),
  "DESCRIPTION", "NAMESPACE", "tools/benchmarks/inference_qualification_20261112/run.R")
config <- list(phase = phase, plan = plan, R = R.version.string,
  sources = stats::setNames(vapply(sources, digest::digest, "", file = TRUE, algo = "sha256"), sources))
metadata <- file.path(output, "configuration.rds")
if (file.exists(metadata)) stopifnot(identical(readRDS(metadata), config)) else atomic(config, metadata)

simulate_case <- function(case, mixed = FALSE) {
  set.seed(case$seed)
  n <- if (mixed) 100L else 110L
  sigma <- if (case$scenario == "near_bound") .73 else .24
  effect <- if (case$scenario == "weak") .12 else .95
  size <- if (case$scenario == "near_bound") Inf else 2
  d <- seq(.015, 1, length.out = 21L)
  count <- matrix(sample(1:5, n * length(d), replace = TRUE), n)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(12 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .25), n)
  V <- radial * count
  # Evaluate the Gaussian kernel directly in ordinary R.
  habitat <- function(scale) {
    w <- exp(-(d^2 - min(d)^2) / (2 * scale^2))
    as.numeric(base::scale(drop(V %*% w) / drop(count %*% w)))
  }
  data <- data.frame(hab = habitat(.23), z = rnorm(n), effort = runif(n, .8, 1.2),
    site = factor(rep(seq_len(20L), each = 5L, length.out = n)))
  eta <- .4 + effect * habitat(sigma) + .3 * data$z + log(data$effort)
  if (mixed) {
    sd <- if (case$scenario == "weak") .4 else .8
    eta <- eta + rnorm(nlevels(data$site), sd = sd)[data$site]
  }
  y <- switch(case$family, gaussian = rnorm(n, eta, .65),
    poisson = rpois(n, exp(eta)), binomial = rbinom(n, 5L, plogis(eta)),
    nbinom2 = if (is.infinite(size)) rpois(n, exp(eta)) else rnbinom(n, size = size, mu = exp(eta)),
    zero_inflated_poisson = rpois(n, exp(eta)) * rbinom(n, 1L, .70))
  data$y <- y
  data$trials <- rep(5L, n)
  ids <- rownames(data)
  inputs <- list(raw_cov = NULL, d_list = NULL, min_D = 50, max_D = 800,
    unit_conv = 1000, kernel = "gaussian",
    binned = list(bin_rep = d, vsum = list(hab = V), csum = list(hab = count),
      point_ids = ids, n_points = n))
  formula <- if (case$family == "binomial")
    cbind(y, trials - y) ~ hab + z + offset(log(effort)) else
    if (mixed) y ~ hab + z + offset(log(effort)) + (1 | site) else
      y ~ hab + z + offset(log(effort))
  list(data = data, inputs = inputs, formula = formula, habitat = habitat,
    truth = list(sigma = sigma * 1000, effect = effect, size = size))
}

native_at <- function(data, formula, family, mixed, independent_nb2 = TRUE) {
  if (mixed) {
    fam <- if (family == "nbinom2") glmmTMB::nbinom2() else poisson()
    zi <- if (family == "zero_inflated_poisson") ~1 else ~0
    fit <- glmmTMB::glmmTMB(formula, data = data, family = fam, ziformula = zi)
    return(list(model = fit, nll = -as.numeric(logLik(fit)),
      converged = isTRUE(fit$sdr$pdHess) && isTRUE(fit$fit$convergence == 0)))
  }
  fam <- switch(family, gaussian = gaussian(), poisson = poisson(),
    binomial = binomial(), nbinom2 = MASS::negative.binomial(2))
  if (family != "nbinom2" || !independent_nb2) {
    fit <- glm(formula, data = data, family = fam)
    return(list(model = fit, nll = -as.numeric(logLik(fit)), converged = fit$converged))
  }
  # Separate base-R optimization of the normalized NB2 likelihood with a
  # literal Poisson face. This reference never evaluates package objectives.
  start <- glm(formula, data = data, family = poisson())
  X <- model.matrix(start); offset <- model.offset(model.frame(start)); y <- data$y
  fn <- function(p) {
    eta <- drop(X %*% p[-length(p)]) + offset
    alpha <- p[length(p)]
    if (any(!is.finite(eta) | eta > 600)) return(1e60)
    ll <- if (alpha == 0) dpois(y, exp(eta), log = TRUE) else
      dnbinom(y, size = 1 / alpha, mu = exp(eta), log = TRUE)
    if (any(!is.finite(ll))) return(1e60)
    -sum(ll)
  }
  k <- ncol(X)
  opts <- lapply(c(.05, .5, 2), function(a) optim(c(coef(start), a), fn,
    method = "L-BFGS-B", lower = c(rep(-Inf, k), 0), upper = c(rep(Inf, k), 1e4),
    control = list(factr = 1e5, pgtol = 1e-7, maxit = 300)))
  valid <- vapply(opts, function(x) x$convergence == 0L && is.finite(x$value), logical(1))
  candidates <- c(-as.numeric(logLik(start)), vapply(opts[valid], "[[", 0, "value"))
  list(nll = min(candidates), converged = isTRUE(start$converged),
    poisson_nll = candidates[1L], interior_attempts = length(opts), interior_converged = sum(valid))
}

run_case <- function(case) {
  mixed <- phase == "mixed"
  f <- simulate_case(case, mixed)
  initial <- warning_log(native_at(f$data, f$formula, case$family, mixed, independent_nb2 = FALSE))
  methods <- if (mixed) "joint_mixed" else if (phase == "bootstrap") "joint_nb2" else
    if (case$family == "nbinom2") c("profile_glm", "joint_nb2") else "profile_glm"
  score <- function(scale) {
    data <- f$data; data$hab <- f$habitat(scale)
    ref <- warning_log(native_at(data, f$formula, case$family, mixed))
    list(scale = scale, nll = if (inherits(ref$value, "error") || !isTRUE(ref$value$converged))
      NA_real_ else ref$value$nll, warnings = ref$warnings,
      error = if (inherits(ref$value, "error")) conditionMessage(ref$value))
  }
  profile <- if (phase != "bootstrap") lapply(grid, score) else NULL
  reference <- vapply(profile, function(x) x$nll, 0)
  attempts <- lapply(methods, function(backend) {
    record <- list(method = backend, status = "failed", reason = NA_character_,
      sigma = NA_real_, nll = NA_real_, native_gap = NA_real_, grid_excess = NA_real_,
      wald_90_covers = NA, alpha_zero = NA, profile_accepts_truth = NA,
      bootstrap_accepts_truth = NA, bootstrap_p = NA_real_, bootstrap_ambiguous = NA)
    if (inherits(initial$value, "error")) {
      record$reason <- conditionMessage(initial$value)
      return(record)
    }
    fit <- warning_log(multiScale_optim(initial$value$model, f$inputs,
      backend = backend, n_cores = 1, verbose = FALSE,
      profile_control = list(starts = matrix(c(.075, .20, .43, .71), ncol = 1L))))
    record$warnings <- fit$warnings
    if (inherits(fit$value, "error")) {
      record$reason <- conditionMessage(fit$value)
      record$starts <- fit$value$start_results
      return(record)
    }
    x <- fit$value
    record$status <- "fitted"
    record$sigma <- x$scale_est$Mean[1L]; record$nll <- -as.numeric(logLik(x$opt_mod))
    record$alpha_zero <- if (backend == "joint_nb2") x$dispersion_est$exact_poisson else NA
    if (phase != "bootstrap") {
      at <- score(record$sigma / 1000)
      record$native_gap <- record$nll - at$nll
      if (any(is.finite(reference))) record$grid_excess <- record$nll - min(reference, na.rm = TRUE)
      record$independent <- at
      s <- x$scale_est$SE[1L]
      if (is.finite(s)) record$wald_90_covers <- abs(log(record$sigma / f$truth$sigma)) <= qnorm(.95) * s / record$sigma
    } else {
      a <- warning_log(scale_inference(x, sigma_values = f$truth$sigma,
        level = .90, bootstrap = 19, bootstrap_values = f$truth$sigma,
        seed = case$seed + 1L, n_starts = 2))
      if (inherits(a$value, "error")) {
        record$status <- "diagnostic_failed"
        record$reason <- conditionMessage(a$value)
      } else {
        record$inference <- a$value[c("grid", "diagnostics", "bootstrap_records")]
        row <- a$value$grid[which.min(abs(a$value$grid$sigma - f$truth$sigma)), , drop = FALSE]
        record$profile_accepts_truth <- row$profile_accept
        record$bootstrap_accepts_truth <- row$bootstrap_accept
        record$bootstrap_p <- row$bootstrap_p
        record$bootstrap_ambiguous <- row$mc_ambiguous
        if (is.na(record$bootstrap_p)) record$status <- "diagnostic_failed"
      }
    }
    record
  })
  list(case = case, truth = f$truth, reference = profile, attempts = attempts)
}
for (i in seq_len(nrow(plan))) {
  path <- file.path(output, paste0(plan$id[i], ".rds"))
  if (!file.exists(path)) {
    cat(plan$id[i], "\n")
    atomic(run_case(plan[i, ]), path)
  }
}
keys <- c("method", "status", "reason", "sigma", "nll", "native_gap", "grid_excess",
  "wald_90_covers", "alpha_zero", "profile_accepts_truth", "bootstrap_accepts_truth",
  "bootstrap_p", "bootstrap_ambiguous")
results <- lapply(file.path(output, paste0(plan$id, ".rds")), readRDS)
table <- do.call(rbind, lapply(results, function(x) do.call(rbind, lapply(x$attempts, function(y) {
  data.frame(x$case, true_sigma = x$truth$sigma, as.data.frame(y[keys]))
}))))
utils::write.csv(table, file.path(output, "summary.csv"), row.names = FALSE)
print(with(table, table(method, status)))

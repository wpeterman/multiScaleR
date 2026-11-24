# Design and per-dataset evaluation for the joint_mixed Wald coverage study.
# Sourced by run.R in the driver and in every worker.

jm_plan <- function(replicates) {
  sigmas <- c(.12, .24, .45)
  plan <- expand.grid(family = c("poisson", "nbinom2"), sigma = sigmas,
    replicate = seq_len(replicates), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  plan$seed <- 2611000L + (plan$replicate - 1L) * 100L +
    match(plan$sigma, sigmas) * 10L + match(plan$family, c("poisson", "nbinom2"))
  plan$id <- sprintf("%s-s%03d-%03d", plan$family, round(plan$sigma * 1000), plan$replicate)
  plan
}

# Binned radial landscape, as in tools/benchmarks/analytic_wald_intervals,
# with points grouped into sites that share a normal random intercept.
# Distances are in internal units; unit_conv = 1000 gives map units.
jm_case <- function(case, n = 200L, groups = 40L, effect = .8, re_sd = .5, size = 2) {
  set.seed(case$seed)
  d <- seq(.015, 1, length.out = 23L)
  count <- matrix(sample(1:5, n * length(d), replace = TRUE), n)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(12 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .25), n)
  value <- radial * count
  # Same standardized Gaussian kernel mean that the package computes.
  habitat <- function(log_scale) {
    weight <- exp(-(d^2 - min(d)^2) / (2 * exp(2 * log_scale)))
    predictor <- drop(value %*% weight) / drop(count %*% weight)
    (predictor - mean(predictor)) / stats::sd(predictor)
  }
  site <- factor(sample(rep(seq_len(groups), length.out = n)))
  b <- rnorm(groups, 0, re_sd)
  data <- data.frame(hab = habitat(log(.23)), z = rnorm(n), effort = runif(n, .8, 1.2),
                     site = site)
  eta <- .4 + effect * habitat(log(case$sigma)) + .3 * data$z + log(data$effort) +
    b[as.integer(site)]
  data$y <- if (case$family == "poisson") rpois(n, exp(eta)) else
    rnbinom(n, mu = exp(eta), size = size)
  inputs <- list(raw_cov = NULL, d_list = NULL, min_D = 50, max_D = 800,
    unit_conv = 1000, kernel = "gaussian",
    binned = list(bin_rep = d, vsum = list(hab = value), csum = list(hab = count),
      point_ids = rownames(data), n_points = n))
  list(data = data, inputs = inputs, habitat = habitat)
}

jm_capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = unique(messages))
}

jm_native <- function(data, family) {
  fam <- if (family == "poisson") stats::poisson() else glmmTMB::nbinom2()
  glmmTMB::glmmTMB(y ~ hab + z + (1 | site) + offset(log(effort)), data = data, family = fam)
}

jm_run_case <- function(case) {
  sample <- jm_case(case)
  started <- proc.time()[["elapsed"]]
  start <- jm_capture(jm_native(sample$data, case$family))
  fit <- if (inherits(start$value, "error")) start else jm_capture(multiScale_optim(
    start$value, sample$inputs, backend = "joint_mixed", n_cores = 1, verbose = FALSE))
  elapsed <- proc.time()[["elapsed"]] - started
  ok <- !inherits(fit$value, "error")
  truth <- case$sigma * 1000
  z95 <- stats::qnorm(.975); z90 <- stats::qnorm(.95)
  out <- data.frame(case[c("id", "family", "sigma", "replicate", "seed")],
    true_sigma = truth, fit_ok = ok,
    fit_error = if (ok) "" else conditionMessage(fit$value),
    warnings = paste(unique(c(start$warnings, fit$warnings)), collapse = " | "),
    seconds = elapsed, stringsAsFactors = FALSE)
  if (!ok) return(out)
  m <- fit$value
  jm <- m$diagnostics$joint_mixed
  selected <- m$mixed_fit$start_results[[m$mixed_fit$selected]]
  V <- selected$joint_sdreport$cov.fixed
  k <- which(rownames(V) == "kernel_par")
  # Raw log-scale SE straight from the sdreport, before any withholding.
  raw_log_se <- sqrt(V[k, k])
  estimate <- m$scale_est$Mean
  public <- summary(m)$opt_scale
  reported_log_se <- unname(jm$log_scale_se)
  wald <- function(se, z) estimate * exp(c(-1, 1) * z * se)
  raw95 <- wald(raw_log_se, z95); raw90 <- wald(raw_log_se, z90)
  rep95 <- if (is.finite(reported_log_se)) wald(reported_log_se, z95) else c(NA_real_, NA_real_)
  # Likelihood-ratio reference at the generating scale: the native glmmTMB
  # fit at the true sigma re-optimizes every other parameter.
  at_truth <- sample$data
  at_truth$hab <- sample$habitat(log(case$sigma))
  native_truth <- jm_capture(jm_native(at_truth, case$family))
  nll_truth <- if (inherits(native_truth$value, "error")) NA_real_ else
    -as.numeric(stats::logLik(native_truth$value))
  lr <- 2 * (nll_truth - m$optim_results$value)
  # Spread of eligible starts inside the 95% likelihood threshold.
  outcomes <- m$mixed_fit$start_results
  eligible <- vapply(outcomes, function(x) isTRUE(x$eligible), TRUE)
  objective <- vapply(outcomes, function(x) if (isTRUE(x$eligible)) x$objective else NA_real_, 0)
  near <- which(eligible & objective <= min(objective, na.rm = TRUE) + stats::qchisq(.95, 1) / 2)
  spread <- max(abs(vapply(outcomes[near], function(x) x$log_scale, 0) - m$mixed_fit$log_scale))
  coef <- glmmTMB::fixef(m$opt_mod)$cond
  cbind(out, data.frame(
    eligible_starts = jm$eligible_starts, attempted_starts = jm$attempted_starts,
    estimate = estimate, log_error = log(estimate / truth),
    at_lower = m$mixed_fit$log_scale - m$mixed_fit$lower < 1e-5,
    at_upper = m$mixed_fit$upper - m$mixed_fit$log_scale < 1e-5,
    status = jm$scale_se_status, reason = unname(jm$scale_se_reasons),
    pd_hessian = isTRUE(jm$pd_hessian), raw_log_se = raw_log_se,
    reported_log_se = reported_log_se, start_spread = spread,
    public_lower = public[1L, "2.5%"], public_upper = public[1L, "97.5%"],
    public_matches = if (is.finite(reported_log_se))
      isTRUE(all.equal(unname(unlist(public[1L, c("2.5%", "97.5%")])), rep95, tolerance = 1e-10)) else
      all(is.na(public[1L, c("2.5%", "97.5%")])),
    raw_lower95 = raw95[1L], raw_upper95 = raw95[2L],
    raw_cover95 = truth >= raw95[1L] & truth <= raw95[2L],
    raw_cover90 = truth >= raw90[1L] & truth <= raw90[2L],
    reported_cover95 = if (is.finite(reported_log_se)) truth >= rep95[1L] & truth <= rep95[2L] else NA,
    raw_above = truth > raw95[2L], raw_below = truth < raw95[1L],
    nll_hat = m$optim_results$value, nll_truth = nll_truth, lr_truth = lr,
    lr_cover95 = lr <= stats::qchisq(.95, 1), lr_cover90 = lr <= stats::qchisq(.90, 1),
    beta_hab = unname(coef["hab"]),
    re_sd = unname(attr(glmmTMB::VarCorr(m$opt_mod)$cond$site, "stddev")),
    theta = if (case$family == "nbinom2") stats::sigma(m$opt_mod) else NA_real_,
    stringsAsFactors = FALSE))
}

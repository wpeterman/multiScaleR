# Run from the source root. Shards resume per dataset and use one R thread.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1L] %in% c("1", "2"))
  stop("Supply shard 1 or 2 and one output directory.")
shard <- as.integer(args[1L])
output <- args[2L]
study <- "tools/benchmarks/analytic_wald_intervals"
source(file.path(study, "reference.R"))
pkgload::load_all(compile = FALSE, quiet = TRUE)
stopifnot(requireNamespace("MASS", quietly = TRUE), !is.null(getLoadedDLLs()$multiScaleR))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
records <- file.path(output, paste0("shard", shard))
dir.create(records, showWarnings = FALSE)
atomic <- function(x, file) {
  tmp <- tempfile(tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, file)) stop("Cannot save ", file)
}
plan <- wald_plan(48L)
files <- c(list.files("R", full.names = TRUE), list.files("src", pattern = "[.](c|cpp|h)$", full.names = TRUE),
  "DESCRIPTION", "NAMESPACE", file.path(study, c("run.R", "reference.R")))
config <- list(plan = plan, R = R.version.string, locale = Sys.getlocale(),
  MASS = as.character(packageVersion("MASS")),
  sources = stats::setNames(vapply(files, digest::digest, "", file = TRUE, algo = "sha256"), files))
config_path <- file.path(records, "config.rds")
if (file.exists(config_path)) stopifnot(identical(config, readRDS(config_path))) else atomic(config, config_path)

run_backend <- function(case, sample, initial, backend) {
  fit <- if (inherits(initial$value, "error")) initial else wald_capture(
    multiScale_optim(initial$value, sample$inputs, backend = backend, n_cores = 1, verbose = FALSE))
  reference <- wald_reference(sample, case$family, backend)
  fn <- reference$fn
  truth <- sample$truth$log_scale
  bounds <- log(c(.05, .8))
  fit_ok <- !inherits(fit$value, "error")
  model <- if (fit_ok) fit$value else NULL
  estimate <- if (fit_ok) model$scale_est$Mean[1L] else NA_real_
  package_nll <- if (fit_ok) model$optim_results$value else NA_real_
  levels <- c(.90, .95)
  local_wald <- lapply(levels, function(level) {
    available <- fit_ok && isTRUE(model$uncertainty$available)
    logse <- if (available) sqrt(model$uncertainty$vcov[1L, 1L]) else NA_real_
    if (available) {
      endpoints <- estimate * exp(c(-1, 1) * stats::qnorm((1 + level) / 2) * logse)
      if (level == .95) {
        public <- summary(model)$opt_scale[1L, c("2.5%", "97.5%")]
        stopifnot(isTRUE(all.equal(as.numeric(public), endpoints, tolerance = 1e-12)))
      }
    } else endpoints <- c(NA_real_, NA_real_)
    data.frame(level = level, wald_available = available, wald_lower = endpoints[1L], wald_upper = endpoints[2L],
      wald_cover = if (available) exp(truth) * 1000 >= endpoints[1L] & exp(truth) * 1000 <= endpoints[2L] else NA,
      wald_outside_search = if (available) endpoints[1L] < 50 || endpoints[2L] > 800 else NA,
      wald_log_se = logse)
  })
  wald <- do.call(rbind, local_wald)
  profile <- wald_capture({
    selected <- if (fit_ok) log(estimate / 1000) else numeric()
    base_grid <- seq(bounds[1L], bounds[2L], length.out = 41L)
    fine_grid <- seq(bounds[1L], bounds[2L], length.out = 81L)
    primary <- sort(unique(c(base_grid, truth, selected)))
    primary_nll <- vapply(primary, fn, 0)
    optimum <- wald_global(fn, primary, primary_nll)
    refined <- sort(unique(c(fine_grid, primary, optimum$minima)))
    check <- wald_global(fn, refined, vapply(refined, fn, 0))
    refined <- sort(unique(c(refined, check$minima)))
    final <- vapply(refined, fn, 0)
    if (check$best > min(final) + 1e-7) stop("Scale refinement missed a better likelihood.")
    profile_sets <- lapply(levels, function(level) wald_profile_set(fn, refined, final, check$best, level))
    base_sets <- lapply(levels, function(level) {
      grid <- sort(unique(c(primary, optimum$minima)))
      wald_profile_set(fn, grid, vapply(grid, fn, 0), optimum$best, level)
    })
    null <- wald_glm(cbind(1, sample$data$z), sample$data$y, sample$offset,
      if (case$family == "nbinom2") "poisson" else case$family, sample$data$trials)
    null_nll <- if (case$family == "nbinom2") wald_nb_reference(cbind(1, sample$data$z),
      sample$data$y, sample$offset, null, backend == "joint_nb2")$nll else null$nll
    list(nll = check$best, optimum = check$optimum, sets = profile_sets, base = base_sets,
      same_scale = if (fit_ok) fn(selected) else NA_real_, nll_truth = fn(truth),
      null_lr = 2 * (null_nll - check$best), refine_gap = optimum$best - check$best)
  })
  ref_ok <- !inherits(profile$value, "error")
  p <- if (ref_ok) profile$value else NULL
  same <- if (ref_ok) p$same_scale - package_nll else NA_real_
  gap <- if (ref_ok) package_nll - p$nll else NA_real_
  summary <- data.frame(case, backend = backend, fit_available = fit_ok,
    fit_message = if (fit_ok) "" else conditionMessage(fit$value),
    fit_warnings = paste(unique(c(initial$warnings, fit$warnings)), collapse = " | "),
    estimated_scale = estimate, true_scale = exp(truth) * 1000, effect = sample$truth$effect,
    true_theta = sample$truth$theta, nll = package_nll,
    fit_boundary = if (fit_ok) isTRUE(model$uncertainty$boundary) else NA,
    fit_scale_boundary = if (fit_ok) min(abs(log(estimate / 1000) - bounds)) < 1e-5 else NA,
    fit_poisson_face = if (fit_ok && backend == "joint_nb2") isTRUE(model$dispersion_est$exact_poisson) else NA,
    fit_theta = if (fit_ok && case$family == "nbinom2") model$dispersion_est$theta else NA_real_,
    fit_weak = if (fit_ok) isTRUE(model$uncertainty$weak) else NA,
    fit_multiple = if (fit_ok) isTRUE(model$uncertainty$multiple_modes) else NA,
    fit_wald_reasons = if (fit_ok) paste(model$uncertainty$reasons, collapse = " | ") else "",
    profile_available = ref_ok, profile_message = if (ref_ok) "" else conditionMessage(profile$value),
    profile_warnings = paste(profile$warnings, collapse = " | "),
    profile_estimate = if (ref_ok) exp(p$optimum) * 1000 else NA_real_,
    profile_scale_boundary = if (ref_ok) min(abs(p$optimum - bounds)) < 1e-5 else NA,
    same_scale_nll_gap = same, global_nll_gap = gap,
    profile_minimum_grid_refinement = if (ref_ok) p$refine_gap else NA_real_,
    profile_truth_lr = if (ref_ok) 2 * (p$nll_truth - p$nll) else NA_real_,
    null_effect_lr = if (ref_ok) p$null_lr else NA_real_)
  metrics <- do.call(rbind, lapply(seq_along(levels), function(j) {
    d <- if (ref_ok) p$sets[[j]] else NULL
    base <- if (ref_ok) p$base[[j]] else NULL
    profile_cover <- if (ref_ok) summary$profile_truth_lr <= stats::qchisq(levels[j], 1) else NA
    profile_set_cover <- if (ref_ok)
      any(summary$true_scale >= d$lower - 1e-6 & summary$true_scale <= d$upper + 1e-6) else NA
    data.frame(profile_cover = profile_cover, profile_set_cover = profile_set_cover,
      profile_components = if (ref_ok) nrow(d) else NA_integer_,
      profile_lower_bound = if (ref_ok) any(d$lower_bound) else NA,
      profile_upper_bound = if (ref_ok) any(d$upper_bound) else NA,
      profile_full = if (ref_ok) nrow(d) == 1L && all(c(d$lower_bound, d$upper_bound)) else NA,
      profile_width_log = if (ref_ok) sum(log(d$upper / d$lower)) else NA_real_,
      profile_refined_count_change = if (ref_ok) nrow(d) != nrow(base) else NA,
      profile_refined_endpoint_max = if (ref_ok && nrow(d) == nrow(base))
        max(abs(log(c(d$lower, d$upper) / c(base$lower, base$upper)))) else NA_real_)
  }))
  package_profile <- if (fit_ok && case$replicate == 1L) wald_capture(summary(model, profile = TRUE)$opt_scale) else NULL
  if (!is.null(package_profile)) {
    ok <- !inherits(package_profile$value, "error")
    d <- if (ok) package_profile$value else NULL
    h <- if (ref_ok) p$sets[[2L]] else NULL
    connected <- if (ref_ok) h[estimate >= h$lower - 1e-6 & estimate <= h$upper + 1e-6, ] else NULL
    package_profile_record <- data.frame(case, backend = backend,
      returned = ok, reason = if (ok) "" else conditionMessage(package_profile$value),
      warnings = paste(package_profile$warnings, collapse = " | "),
      lower = if (ok) d[1L, "2.5%"] else NA_real_, upper = if (ok) d[1L, "97.5%"] else NA_real_,
      connected_lower = if (!is.null(connected) && nrow(connected) == 1L) connected$lower else NA_real_,
      connected_upper = if (!is.null(connected) && nrow(connected) == 1L) connected$upper else NA_real_)
  } else package_profile_record <- NULL
  # Keep complete summaries and every observed optimization/reference failure.
  audit <- do.call(rbind, reference$audit())
  failed <- if (!is.null(audit)) audit[nzchar(audit$error) | nzchar(audit$warnings), ] else NULL
  starts <- if (fit_ok) model$profile_fit$start_results else fit$value$start_results
  list(metrics = cbind(summary[rep(1L, length(levels)), ], wald, metrics),
    sets = if (ref_ok) do.call(rbind, lapply(seq_along(levels), function(j)
      cbind(id = case$id, backend = backend, level = levels[j], p$sets[[j]]))) else NULL,
    failed_reference_points = failed, starts = starts, package_profile = package_profile_record)
}

run_case <- function(case) {
  data <- wald_case(case)
  fam <- switch(case$family, gaussian = gaussian(), poisson = poisson(),
    binomial = binomial(), nbinom2 = MASS::negative.binomial(2))
  initial <- wald_capture(glm(data$formula, data = data$data, family = fam))
  backends <- c("profile_glm", if (case$family == "nbinom2") "joint_nb2")
  lapply(backends, function(backend) run_backend(case, data, initial, backend))
}

# Alternate shards within every cell, including both sample sizes.
todo <- plan[plan$replicate %% 2L == shard %% 2L, ]
for (i in seq_len(nrow(todo))) {
  case <- todo[i, , drop = FALSE]
  destination <- file.path(records, paste0(case$id, ".rds"))
  if (!file.exists(destination)) atomic(run_case(case), destination)
  if (i %% 12L == 0L || i == nrow(todo)) cat("Shard", shard, "finished", i, "of", nrow(todo), "datasets\n")
}

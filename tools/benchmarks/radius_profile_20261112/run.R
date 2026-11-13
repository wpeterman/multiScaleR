# Run from the package root: Rscript tools/benchmarks/radius_profile_20261112/run.R
# Set an optional first argument to write an independent copy of the results.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1L]] else "tools/benchmarks/radius_profile_20261112"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(compile = FALSE, quiet = TRUE)
stopifnot(requireNamespace("survival", quietly = TRUE))
strata <- survival::strata

scenarios <- c("gaussian_interior", "gaussian_lower_bound",
  "gaussian_flat", "choice_expow")
seeds <- c(1401L, 1402L, 1403L)
radius_grid <- list(diversity = c(40, 100, 150), rough = c(100, 180))
combinations <- expand.grid(rough = radius_grid$rough,
  diversity = radius_grid$diversity, KEEP.OUT.ATTRS = FALSE)
bounds <- log(c(.04, .24))
shape_bounds <- log(c(.75, 50))
objective_tolerance <- 1e-5

capture <- function(expr) {
  warnings <- character()
  result <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(result = result, warnings = unique(warnings))
}

fixture <- function(scenario, seed) {
  set.seed(seed)
  n <- 60L
  distance <- c(.025, rep(c(.075, .125, .175, .225), each = 5L))
  cells <- length(distance)
  cover <- t(replicate(n, sample(1:3, cells, replace = TRUE,
    prob = stats::rexp(3L))))
  terrain <- matrix(stats::rnorm(n * cells), nrow = n) * stats::runif(n, .5, 2)
  habitat <- if (scenario == "gaussian_flat")
    matrix(rep(stats::rnorm(n), cells), nrow = n) else
      matrix(stats::rnorm(n * cells), nrow = n)
  kernel <- if (scenario == "choice_expow") "expow" else "gaussian"
  raw <- function(radius) {
    keep <- distance <= radius / 1000
    vapply(seq_len(n), function(i) {
      freq <- prop.table(table(cover[i, keep]))
      -sum(freq * log(freq))
    }, numeric(1))
  }
  rough <- function(radius)
    apply(terrain[, distance <= radius / 1000, drop = FALSE], 1L, stats::sd)
  smooth <- function(sigma, shape = 2) {
    weights <- if (kernel == "gaussian") exp(-.5 * (distance / sigma)^2) else
      exp(-(distance / sigma)^shape)
    drop(habitat %*% (weights / sum(weights)))
  }
  true_sigma <- if (scenario == "gaussian_lower_bound") .015 else .104
  signal <- drop(.8 * scale(raw(150)) - .65 * scale(rough(180)) +
    1.25 * scale(smooth(true_sigma, 2.5)))
  ids <- sprintf("row_%03d", seq_len(n))
  step <- rep(seq_len(n / 3L), each = 3L)
  used <- integer(n)
  for (i in split(seq_len(n), step))
    used[sample(i, 1L, prob = exp(signal[i]))] <- 1L
  initial_values <- cbind(rough = rough(100), diversity = raw(100), habitat = smooth(.12))
  prepared <- as.data.frame(scale(initial_values))
  rownames(prepared) <- ids
  data <- data.frame(y = signal + stats::rnorm(n, sd = .3), used = used,
    time = 1, step = step, prepared, row.names = ids)
  raw_cov <- stats::setNames(lapply(seq_len(n), function(i)
    cbind(cover = cover[i, ], elevation = terrain[i, ], hab = habitat[i, ])), ids)
  d_list <- stats::setNames(rep(list(distance), n), ids)
  inputs <- structure(list(kernel_dat = prepared, d_list = d_list, raw_cov = raw_cov,
    kernel = kernel, sigma = c(.1, .1, .12),
    shape = if (kernel == "expow") rep(2, 3L) else NULL,
    min_D = 40, max_D = 240, unit_conv = 1000, resolution = 20, n_cols = 50L,
    cell_size = c(20, 20), n_covs = 3L, cell_data_stored = TRUE,
    scale_vars = msr_vars(rough = surface_var("elevation", "sq"),
      diversity = landscape_var("cover", "shdi"), habitat = kernel_var("hab")),
    binned = multiScaleR:::.msr_build_kernel_bins(d_list, raw_cov,
      "habitat", "hab", nbins = 16L, point_ids = ids),
    scl_params = list(mean = colMeans(initial_values), sd = apply(initial_values, 2L, stats::sd))),
    class = "multiScaleR_data")
  list(data = data, inputs = inputs, raw = raw, rough = rough, smooth = smooth,
    true_sigma = 1000 * true_sigma, true_shape = if (kernel == "expow") 2.5 else NA_real_)
}

# Use native unweighted least squares and its closed-form Gaussian likelihood.
# For choice sets use Breslow's native C fit. One event per stratum makes it
# equivalent to the exact partial likelihood, checked against coxph() below.
native_target <- function(f, hard) {
  n <- nrow(f$data)
  h1 <- f$raw(hard$diversity)
  h2 <- f$rough(hard$rough)
  if (any(!is.finite(c(h1, h2))) || stats::sd(h1) == 0 || stats::sd(h2) == 0)
    stop("Undefined or constant hard-radius predictor.")
  hard_values <- scale(cbind(diversity = h1, rough = h2))
  if (f$inputs$kernel == "gaussian") return(function(theta) {
    X <- cbind(1, drop(scale(f$smooth(exp(theta[1L])))), hard_values)
    fit <- stats::lm.fit(X, f$data$y)
    if (fit$rank != ncol(X)) stop("Rank-deficient Gaussian fit.")
    -n / 2 * (log(2 * pi) + 1 + log(sum(fit$residuals^2) / n))
  })

  y <- survival::Surv(f$data$time, f$data$used)
  control <- survival::coxph.control(iter.max = 30L, eps = 1e-10)
  function(theta) {
    X <- cbind(habitat = drop(scale(f$smooth(exp(theta[1L]), exp(theta[2L])))),
      hard_values)
    fit <- survival::coxph.fit(X, y, strata = f$data$step, offset = rep(0, n),
      init = NULL, control = control, weights = rep(1, n), method = "breslow",
      rownames = rownames(f$data), resid = FALSE)
    if (any(!is.finite(fit$coefficients))) stop("Invalid native choice coefficients.")
    fit$loglik[2L]
  }
}

reference <- function(f, hard) {
  target <- native_target(f, hard)
  sigma <- seq(bounds[1L], bounds[2L], length.out = 121L)
  shape <- if (f$inputs$kernel == "gaussian") NULL else
    seq(shape_bounds[1L], shape_bounds[2L], length.out = 41L)
  points <- if (is.null(shape)) matrix(sigma, ncol = 1L) else
    as.matrix(expand.grid(sigma = sigma[seq(1L, 121L, by = 3L)], shape = shape))
  safe <- function(p) {
    out <- capture(target(p))
    if (inherits(out$result, "error") || length(out$warnings) ||
        length(out$result) != 1L || !is.finite(out$result)) return(NA_real_)
    out$result
  }
  scores <- apply(points, 1L, safe)
  dense <- which.max(scores)
  if (!length(dense)) stop("The independent dense search had no finite objective.")
  best <- list(par = points[dense, ], ll = scores[dense])
  finite <- scores[is.finite(scores)]
  spread <- diff(range(finite))
  if (is.null(shape) && spread > 1e-10) {
    local <- which(vapply(seq_along(sigma), function(i)
      i > 1L && i < length(sigma) && is.finite(scores[i]) &&
        all(is.finite(scores[c(i - 1L, i + 1L)])) &&
        scores[i] >= max(scores[c(i - 1L, i + 1L)]), logical(1)))
    for (i in local) {
      fit <- stats::optimize(function(x) {
        ll <- safe(x)
        if (is.finite(ll)) ll else -1e50
      }, sigma[c(i - 1L, i + 1L)], maximum = TRUE, tol = 1e-10)
      if (fit$objective > best$ll) best <- list(par = fit$maximum, ll = fit$objective)
    }
  } else if (!is.null(shape)) {
    # Refine the dense leaders and the best point from every shape band.
    bands <- cut(points[, 2L], breaks = 5L, include.lowest = TRUE)
    band_best <- vapply(split(seq_len(nrow(points)), bands), function(i) {
      rank <- order(scores[i], decreasing = TRUE, na.last = NA)
      if (length(rank)) i[rank[1L]] else NA_integer_
    }, integer(1))
    leaders <- order(scores, decreasing = TRUE, na.last = NA)
    leaders <- unique(c(head(leaders, 6L), band_best[!is.na(band_best)]))
    for (i in leaders) {
      fit <- stats::optim(points[i, ], function(x) {
        ll <- safe(x)
        if (is.finite(ll)) -ll else 1e50
      }, method = "L-BFGS-B", lower = c(bounds[1L], shape_bounds[1L]),
      upper = c(bounds[2L], shape_bounds[2L]),
      control = list(factr = 1e4, pgtol = 1e-7, ndeps = c(1e-5, 1e-5)))
      if (fit$convergence == 0L && -fit$value > best$ll)
        best <- list(par = fit$par, ll = -fit$value)
    }
  }
  list(sigma = 1000 * exp(best$par[1L]),
    shape = if (length(best$par) == 2L) exp(best$par[2L]) else NA_real_,
    log_likelihood = best$ll, dense_log_likelihood = scores[dense],
    dense_range = spread, dense_attempted = length(scores),
    dense_failed = sum(!is.finite(scores)), target = target)
}

all_candidates <- all_starts <- all_replicates <- list()
for (scenario in scenarios) for (seed in seeds) {
  key <- paste(scenario, seed, sep = "_")
  cat("Running", key, "\n")
  flush.console()
  f <- fixture(scenario, seed)
  choice <- scenario == "choice_expow"
  initial <- if (choice)
    survival::coxph(survival::Surv(time, used) ~ habitat + diversity + rough + strata(step),
      data = f$data, method = "exact", model = TRUE) else
    stats::lm(y ~ habitat + diversity + rough, data = f$data)
  base_target <- native_target(f, list(diversity = 100, rough = 100))
  agreement <- abs(base_target(log(c(.12, 2))) - as.numeric(stats::logLik(initial)))
  stopifnot(agreement <= 1e-8)
  timer <- proc.time()[["elapsed"]]
  fitted <- capture(multiScale_optim(initial, f$inputs, radius_grid = radius_grid,
    backend = "refit", verbose = FALSE))
  package_seconds <- proc.time()[["elapsed"]] - timer
  succeeded <- inherits(fitted$result, "multiScaleR")
  tab <- if (succeeded) diagnostics(fitted$result)$radius_profile else fitted$result$radius_profile
  if (is.null(tab)) {
    tab <- data.frame("radius:diversity" = combinations$diversity,
      "radius:rough" = combinations$rough, "sigma:habitat" = NA_real_,
      log_likelihood = NA_real_, status = "fit_error",
      message = conditionMessage(fitted$result), selected = FALSE, check.names = FALSE)
    if (choice) tab[["shape:habitat"]] <- NA_real_
  }
  stopifnot(nrow(tab) == nrow(combinations),
    identical(tab[["radius:diversity"]], combinations$diversity),
    identical(tab[["radius:rough"]], combinations$rough))
  for (i in seq_len(nrow(tab))) {
    row <- tab[i, ]
    hard <- combinations[i, ]
    ref <- capture(reference(f, hard))
    ok <- !inherits(ref$result, "error")
    r <- ref$result
    p_sigma <- row[["sigma:habitat"]]
    p_shape <- if (choice) row[["shape:habitat"]] else NA_real_
    audit <- if (ok && is.finite(p_sigma))
      r$target(log(c(p_sigma / 1000, if (choice) p_shape))) else NA_real_
    all_candidates[[paste(key, i)]] <- data.frame(scenario, seed, candidate = i,
      diversity = hard$diversity, rough = hard$rough,
      package_status = row$status, package_message = row$message,
      package_selected = row$selected, package_sigma = p_sigma, package_shape = p_shape,
      package_loglik = row$log_likelihood, independent_at_package = audit,
      reference_status = if (ok) "scored" else "failed",
      reference_message = if (ok) paste(ref$warnings, collapse = " | ") else conditionMessage(r),
      dense_attempted = if (ok) r$dense_attempted else 0L,
      dense_failed = if (ok) r$dense_failed else 0L,
      dense_range = if (ok) r$dense_range else NA_real_,
      dense_loglik = if (ok) r$dense_log_likelihood else NA_real_,
      reference_sigma = if (ok) r$sigma else NA_real_,
      reference_shape = if (ok) r$shape else NA_real_,
      reference_loglik = if (ok) r$log_likelihood else NA_real_,
      package_minus_reference = if (ok) row$log_likelihood - r$log_likelihood else NA_real_,
      package_minus_independent_at_package = row$log_likelihood - audit)
    starts <- row$smooth_starts[[1L]]
    if (!is.null(starts)) {
      if (!"initial:shape:habitat" %in% names(starts)) starts[["initial:shape:habitat"]] <- NA_real_
      all_starts[[paste(key, i)]] <- cbind(scenario, seed, candidate = i, starts)
    }
  }
  these <- do.call(rbind, all_candidates[paste(key, seq_len(nrow(tab)))])
  reference_eligible <- which(these$reference_status == "scored")
  reference_best <- reference_eligible[which.max(these$reference_loglik[reference_eligible])]
  selected <- which(these$package_selected)
  opt <- if (succeeded) fitted$result else NULL
  h <- if (succeeded) diagnostics(opt)$refit_hessian else NULL
  all_replicates[[key]] <- data.frame(scenario, seed,
    true_sigma = f$true_sigma, true_shape = f$true_shape, package_seconds,
    initial_native_loglik_gap = agreement,
    package_status = if (succeeded) "selected" else "failed",
    package_warnings = paste(fitted$warnings, collapse = " | "),
    package_selected_candidate = if (length(selected)) selected else NA_integer_,
    reference_selected_candidate = reference_best,
    selected_loglik_gap = if (length(selected))
      these$package_loglik[selected] - these$reference_loglik[reference_best] else NA_real_,
    selected_sigma = if (length(selected)) these$package_sigma[selected] else NA_real_,
    selected_shape = if (length(selected)) these$package_shape[selected] else NA_real_,
    selected_sigma_se = if (succeeded) opt$scale_est["habitat", "SE"] else NA_real_,
    selected_shape_se = if (succeeded && choice) opt$shape_est["habitat", "SE"] else NA_real_,
    curvature_status = if (is.null(h)) NA_character_ else h$status,
    active_boundary = if (is.null(h)) NA else h$bounds$active_boundary)
}
candidates <- do.call(rbind, all_candidates)
starts <- do.call(rbind, all_starts)
replicates <- do.call(rbind, all_replicates)
utils::write.csv(candidates, file.path(output, "candidates.csv"), row.names = FALSE)
utils::write.csv(starts, file.path(output, "starts.csv"), row.names = FALSE)
utils::write.csv(replicates, file.path(output, "replicates.csv"), row.names = FALSE)
writeLines(capture.output({
  print(sessionInfo())
  print(tools::md5sum(c("R/radius_grid.R", "R/multiScale_optim.R",
    "R/kernel_scale_fn.R", "R/smooth_glm.R",
    "tools/benchmarks/radius_profile_20261112/run.R")))
  cat("Objective agreement threshold:", objective_tolerance, "\n")
}), file.path(output, "session.txt"))
print(replicates[, c("scenario", "seed", "package_seconds", "selected_loglik_gap",
  "selected_sigma", "selected_shape", "curvature_status")], row.names = FALSE)
print(with(candidates, table(scenario, package_status, reference_status)))
cat("Max absolute package/reference gap:",
  max(abs(candidates$package_minus_reference), na.rm = TRUE), "\n")

# Fail the script after writing evidence if the declared comparison fails.
matched <- candidates$package_status == "scored"
stopifnot(nrow(candidates) == length(scenarios) * length(seeds) * nrow(combinations),
  identical(matched, candidates$reference_status == "scored"),
  sum(matched) == length(scenarios) * length(seeds) * 4L,
  all(abs(candidates$package_minus_reference[matched]) <= objective_tolerance),
  all(abs(candidates$package_minus_independent_at_package[matched]) <= objective_tolerance),
  all(replicates$package_selected_candidate == replicates$reference_selected_candidate),
  all(is.na(replicates$selected_sigma_se[replicates$scenario %in%
    c("gaussian_lower_bound", "gaussian_flat")])))
cat("All declared point-optimization and curvature availability checks passed.\n")

# Compare the optional prefix cache with the scalar metric in public fits.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
output <- args[[1L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(compile = FALSE, quiet = TRUE)
ns <- asNamespace("multiScaleR")
ids <- sprintf("site%03d", seq_len(48L))
radii <- seq(30, 210, by = 15)
d <- seq(.015, .215, length.out = 81L)
d[seq(2L, length(d), by = 3L)] <- d[seq(1L, length(d) - 1L, by = 3L)]
distances <- stats::setNames(rep(list(d), length(ids)), ids)

fixture <- function(seed, scenario) {
  set.seed(seed)
  values <- t(replicate(length(ids), stats::rnorm(length(d)) *
    exp(stats::runif(1L, -.9, .9)) + sin(45 * d + stats::runif(1L, 0, 6))))
  if (scenario == "missing")
    values[matrix(stats::runif(length(values)) < .15, nrow(values))] <- NA_real_
  if (scenario == "large_datum") values <- values + 1e10
  raw <- function(radius) apply(values[, d <= radius / 1000, drop = FALSE],
    1L, stats::sd, na.rm = TRUE)
  z <- raw(105)
  y <- drop(scale(raw(165))) + stats::rnorm(length(ids), sd = .35)
  kernel_dat <- data.frame(rough = drop(scale(z)), row.names = ids)
  data <- cbind(data.frame(y, row.names = ids), kernel_dat)
  buffers <- stats::setNames(lapply(seq_along(ids), function(i)
    cbind(elevation = values[i, ])), ids)
  inputs <- structure(list(kernel_dat = kernel_dat, d_list = distances,
    raw_cov = buffers, kernel = "gaussian", sigma = .105, min_D = 30, max_D = 210,
    unit_conv = 1000, resolution = 10, n_cols = 25L, cell_size = c(10, 10),
    n_covs = 1L, cell_data_stored = TRUE,
    scale_vars = msr_vars(rough = surface_var("elevation", "sq")),
    scl_params = list(mean = c(rough = mean(z)), sd = c(rough = sd(z)))),
    class = "multiScaleR_data")
  list(data = data, inputs = inputs, raw = raw)
}

fit_with_cache <- function(fixture, cached, grid) {
  old <- options(multiScaleR.experimental_sq_cache = cached)
  on.exit(options(old))
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(withCallingHandlers({
    data <- fixture$data
    model <- stats::lm(y ~ rough, data = data)
    multiScale_optim(model, fixture$inputs, radius_grid = if (grid) radii else NULL,
      verbose = FALSE)
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(fit = fit, seconds = proc.time()[["elapsed"]] - start,
    warnings = paste(unique(warnings), collapse = " | "))
}

scenarios <- c("ordinary", "missing", "large_datum")
seeds <- 1801:1808
plan <- expand.grid(scenario = scenarios, seed = seeds, stringsAsFactors = FALSE)
result <- vector("list", nrow(plan) * 2L)
point <- vector("list", nrow(plan))
for (i in seq_len(nrow(plan))) {
  f <- fixture(plan$seed[i], plan$scenario[i])
  prepared <- ns$.msr_prepare_surface_sq(distances, f$inputs$raw_cov, f$inputs$scale_vars)
  reference <- vapply(radii, f$raw, numeric(length(ids)))
  cached <- vapply(radii, function(radius) vapply(prepared, function(p)
    ns$.msr_eval_surface_neighbors(p$rough, radius / 1000), numeric(1)), numeric(length(ids)))
  point[[i]] <- cbind(plan[i, ], missing_mask_equal = identical(is.na(reference), is.na(cached)),
    max_abs_sq_error = max(abs(reference - cached), na.rm = TRUE))
  for (grid in c(FALSE, TRUE)) {
    order <- c(FALSE, TRUE)
    if ((i + grid) %% 2L) order <- rev(order)
    runs <- lapply(order, function(enabled) fit_with_cache(f, enabled, grid))
    runs <- runs[match(c(FALSE, TRUE), order)]
    ok <- vapply(runs, function(x) !inherits(x$fit, "error"), logical(1))
    row <- cbind(plan[i, ], search = if (grid) "declared_grid" else "continuous",
      status_scalar = if (ok[1L]) "ok" else conditionMessage(runs[[1L]]$fit),
      status_cache = if (ok[2L]) "ok" else conditionMessage(runs[[2L]]$fit),
      warnings_scalar = runs[[1L]]$warnings, warnings_cache = runs[[2L]]$warnings,
      seconds_scalar = runs[[1L]]$seconds, seconds_cache = runs[[2L]]$seconds,
      radius_scalar = NA_real_, radius_cache = NA_real_, ll_difference = NA_real_,
      max_fitted_difference = NA_real_, grid_status_equal = NA, max_grid_ll_difference = NA_real_)
    if (all(ok)) {
      a <- runs[[1L]]$fit; b <- runs[[2L]]$fit
      row$radius_scalar <- a$scale_est["rough", "Mean"]
      row$radius_cache <- b$scale_est["rough", "Mean"]
      row$ll_difference <- as.numeric(logLik(b$opt_mod) - logLik(a$opt_mod))
      row$max_fitted_difference <- max(abs(fitted(b$opt_mod) - fitted(a$opt_mod)))
      if (grid) {
        a_grid <- diagnostics(a)$radius_profile
        b_grid <- diagnostics(b)$radius_profile
        row$grid_status_equal <- identical(a_grid$status, b_grid$status)
        scored <- a_grid$status == "scored" & b_grid$status == "scored"
        row$max_grid_ll_difference <- max(abs(a_grid$log_likelihood[scored] -
          b_grid$log_likelihood[scored]))
      }
    }
    result[[2L * (i - 1L) + as.integer(grid) + 1L]] <- row
  }
}
write.csv(do.call(rbind, point), file.path(output, "metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, result), file.path(output, "fits.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))

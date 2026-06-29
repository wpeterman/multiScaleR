# Benchmark the held-covariate reuse path in profile_sigma().
#
# This script compares the current implementation with a legacy-equivalent
# profile loop that recomputes every covariate at every grid point. It uses
# small landscape and surface metric objects so the benchmark can be re-run on a
# laptop while still exercising the expensive covariate families that motivated
# the optimization.
#
# Run from the package root:
#   Rscript tools/benchmarks/profile_sigma_covariate_reuse_benchmark.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package `pkgload` is required to run this benchmark.", call. = FALSE)
}
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package `terra` is required to run this benchmark.", call. = FALSE)
}

pkgload::load_all(".", helpers = FALSE, quiet = TRUE, export_all = TRUE)

set.seed(20260629)

bench_dir <- file.path("tools", "benchmarks")
results_csv <- file.path(bench_dir, "profile_sigma_covariate_reuse_results.csv")
summary_md <- file.path(bench_dir, "profile_sigma_covariate_reuse_summary.md")

make_points <- function(raster, n_x = 6, n_y = 4) {
  ext <- terra::ext(raster)
  x <- seq(ext$xmin + 90, ext$xmax - 90, length.out = n_x)
  y <- seq(ext$ymin + 90, ext$ymax - 90, length.out = n_y)
  pts <- terra::vect(cbind(rep(x, times = n_y), rep(y, each = n_x)),
                     crs = terra::crs(raster))
  pts
}

make_landscape_case <- function() {
  landcover <- terra::rast(nrows = 60, ncols = 60,
                           xmin = 0, xmax = 600,
                           ymin = 0, ymax = 600,
                           crs = "EPSG:3857")
  values <- matrix(sample(1:4, terra::ncell(landcover), replace = TRUE,
                          prob = c(0.35, 0.30, 0.20, 0.15)),
                   nrow = 60, ncol = 60)
  values[seq(1, 60, by = 5), ] <- 1L
  values[, seq(3, 60, by = 7)] <- 3L
  terra::values(landcover) <- as.vector(values)
  names(landcover) <- "landcover"

  pts <- make_points(landcover)
  scale_vars <- msr_vars(
    lc_ed = landscape_var("landcover", metric = "ed"),
    lc_contag = landscape_var("landcover", metric = "contag"),
    lc_ent = landscape_var("landcover", metric = "ent"),
    lc_ai = landscape_var("landcover", metric = "ai")
  )

  kernel_inputs <- kernel_prep(pts = pts,
                               raster_stack = landcover,
                               max_D = 170,
                               scale_vars = scale_vars,
                               verbose = FALSE)
  dat <- data.frame(kernel_inputs$kernel_dat)
  dat$y <- 0.2 * scale(dat$lc_ed)[, 1] - 0.1 * scale(dat$lc_contag)[, 1] +
    0.3 * scale(dat$lc_ent)[, 1] + 0.15 * scale(dat$lc_ai)[, 1] +
    stats::rnorm(nrow(dat), 0, 0.5)
  mod <- lm(y ~ lc_ed + lc_contag + lc_ent + lc_ai, data = dat)
  make_profile_object(
    mod,
    kernel_inputs,
    sigma = c(lc_ed = 60, lc_contag = 90, lc_ent = 120, lc_ai = 150)
  )
}

make_surface_case <- function() {
  elevation <- terra::rast(nrows = 60, ncols = 60,
                           xmin = 0, xmax = 600,
                           ymin = 0, ymax = 600,
                           crs = "EPSG:3857")
  xy <- terra::crds(elevation)
  z <- 40 + 0.04 * xy[, 1] - 0.03 * xy[, 2] +
    8 * sin(xy[, 1] / 35) + 5 * cos(xy[, 2] / 45) +
    stats::rnorm(nrow(xy), 0, 1.5)
  terra::values(elevation) <- z
  names(elevation) <- "elevation"

  pts <- make_points(elevation)
  scale_vars <- msr_vars(
    elev_sq = surface_var("elevation", metric = "sq", weighted = TRUE),
    elev_sa = surface_var("elevation", metric = "sa", weighted = TRUE),
    elev_ssk = surface_var("elevation", metric = "ssk", weighted = TRUE),
    elev_sku = surface_var("elevation", metric = "sku", weighted = TRUE)
  )

  kernel_inputs <- kernel_prep(pts = pts,
                               raster_stack = elevation,
                               max_D = 170,
                               scale_vars = scale_vars,
                               verbose = FALSE)
  dat <- data.frame(kernel_inputs$kernel_dat)
  dat$y <- 0.2 * scale(dat$elev_sq)[, 1] - 0.1 * scale(dat$elev_sa)[, 1] +
    0.3 * scale(dat$elev_ssk)[, 1] + 0.15 * scale(dat$elev_sku)[, 1] +
    stats::rnorm(nrow(dat), 0, 0.5)
  mod <- lm(y ~ elev_sq + elev_sa + elev_ssk + elev_sku, data = dat)
  make_profile_object(
    mod,
    kernel_inputs,
    sigma = c(elev_sq = 60, elev_sa = 90, elev_ssk = 120, elev_sku = 150)
  )
}

make_profile_object <- function(mod, kernel_inputs, sigma) {
  opt_context <- build_opt_context(
    fitted_mod = mod,
    cov_df = kernel_inputs$raw_cov,
    scale_vars = kernel_inputs$scale_vars,
    unit_conv = kernel_inputs$unit_conv,
    resolution = kernel_inputs$resolution,
    n_cols = kernel_inputs$n_cols,
    binned = kernel_inputs$binned
  )

  scale_est <- data.frame(
    Mean = unname(sigma),
    SE = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    row.names = names(sigma),
    check.names = FALSE
  )

  out <- list(
    opt_mod = mod,
    kernel_inputs = kernel_inputs,
    opt_context = opt_context,
    join_by = NULL,
    scale_est = scale_est,
    shape_est = NULL,
    min_D = kernel_inputs$min_D,
    max_D = kernel_inputs$max_D
  )
  class(out) <- "multiScaleR"
  out
}

profile_sigma_legacy <- function(x,
                                 n_pts = 5,
                                 spacing = "linear",
                                 sigma_range = c(30, 130)) {
  kernel_inputs <- x$kernel_inputs
  opt_context <- x$opt_context
  kernel <- kernel_inputs$kernel
  unit_conv <- kernel_inputs$unit_conv
  scale_est <- x$scale_est
  covs <- row.names(scale_est)
  n_covs <- length(covs)
  opt_par <- scale_est$Mean / unit_conv
  sigma_seq <- seq(sigma_range[1], sigma_range[2], length.out = n_pts)

  mod <- x$opt_mod
  k_base <- .msr_parameter_count(mod)
  n <- .msr_model_nobs(mod)
  k_total <- k_base + n_covs

  rows <- vector("list", n_covs)
  for (j in seq_len(n_covs)) {
    ll_vec <- numeric(length(sigma_seq))
    aicc_vec <- numeric(length(sigma_seq))
    for (i in seq_along(sigma_seq)) {
      par_i <- opt_par
      par_i[j] <- sigma_seq[i] / unit_conv
      neg_ll <- kernel_scale_fn(
        par = par_i,
        d_list = kernel_inputs$d_list,
        cov_df = kernel_inputs$raw_cov,
        kernel = kernel,
        fitted_mod = opt_context$fitted_mod,
        join_by = x$join_by,
        mod_return = NULL,
        opt_context = opt_context
      )
      ll_val <- -neg_ll
      aic_val <- -2 * ll_val + 2 * k_total
      aicc_vec[i] <- if ((n - k_total - 1) > 0) {
        aic_val + (2 * k_total * (k_total + 1)) / (n - k_total - 1)
      } else {
        NA_real_
      }
      ll_vec[i] <- ll_val
    }
    rows[[j]] <- data.frame(variable = covs[j],
                            sigma = sigma_seq,
                            LL = ll_vec,
                            AICc = aicc_vec,
                            stringsAsFactors = FALSE)
  }

  out <- list(
    profiles = do.call(rbind, rows),
    opt_sigma = stats::setNames(scale_est$Mean, covs),
    metric = "AICc",
    spacing = spacing,
    sigma_grid = sigma_seq
  )
  class(out) <- "sigma_profile"
  row.names(out$profiles) <- NULL
  out
}

time_one <- function(expr) {
  gc()
  unname(system.time(force(expr))[["elapsed"]])
}

benchmark_case <- function(label, x, reps = 3, n_pts = 5) {
  # Warm both paths and confirm the profile values agree before timing.
  legacy <- profile_sigma_legacy(x, n_pts = n_pts)
  current <- profile_sigma(x, n_pts = n_pts, spacing = "linear",
                           sigma_range = c(30, 130), verbose = FALSE)
  stopifnot(isTRUE(all.equal(legacy$profiles, current$profiles,
                             tolerance = 1e-10, check.attributes = FALSE)))

  rows <- vector("list", reps * 2)
  idx <- 1L
  for (replicate in seq_len(reps)) {
    rows[[idx]] <- data.frame(
      case = label,
      implementation = "legacy_full_recompute",
      replicate = replicate,
      elapsed_sec = time_one(profile_sigma_legacy(x, n_pts = n_pts)),
      n_profile_points = n_pts,
      n_covariates = length(x$scale_est$Mean)
    )
    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      case = label,
      implementation = "current_reuse_held_covariates",
      replicate = replicate,
      elapsed_sec = time_one(profile_sigma(x, n_pts = n_pts,
                                           spacing = "linear",
                                           sigma_range = c(30, 130),
                                           verbose = FALSE)),
      n_profile_points = n_pts,
      n_covariates = length(x$scale_est$Mean)
    )
    idx <- idx + 1L
  }
  do.call(rbind, rows)
}

results <- rbind(
  benchmark_case("landscape_metrics_four_covariates", make_landscape_case()),
  benchmark_case("weighted_surface_metrics_four_covariates", make_surface_case())
)

write.csv(results, results_csv, row.names = FALSE)

summary_rows <- do.call(
  rbind,
  lapply(split(results, results$case), function(dat) {
    legacy <- stats::median(dat$elapsed_sec[dat$implementation == "legacy_full_recompute"])
    current <- stats::median(dat$elapsed_sec[dat$implementation == "current_reuse_held_covariates"])
    data.frame(
      case = dat$case[[1]],
      legacy_median_sec = legacy,
      current_median_sec = current,
      speedup = legacy / current,
      n_profile_points = dat$n_profile_points[[1]],
      n_covariates = dat$n_covariates[[1]]
    )
  })
)

fmt <- function(x) format(round(x, 3), nsmall = 3)
md <- c(
  "# profile_sigma() Held-Covariate Reuse Benchmark",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "This benchmark compares the current `profile_sigma()` implementation with a",
  "legacy-equivalent loop that recomputes every covariate at every profile grid",
  "point. The profile values are checked for equality before timing.",
  "",
  "| Case | Covariates | Points per profile | Legacy median (s) | Current median (s) | Speedup |",
  "|------|------------|--------------------|-------------------|--------------------|---------|",
  apply(summary_rows, 1, function(row) {
    sprintf("| `%s` | %s | %s | %s | %s | %sx |",
            row[["case"]],
            row[["n_covariates"]],
            row[["n_profile_points"]],
            fmt(as.numeric(row[["legacy_median_sec"]])),
            fmt(as.numeric(row[["current_median_sec"]])),
            fmt(as.numeric(row[["speedup"]])))
  }),
  "",
  sprintf("Raw timings: `%s`.", results_csv)
)

writeLines(md, summary_md)
print(summary_rows)

# Internal helpers for continuous surface texture metric support.
#
# Surface texture metrics summarize the within-neighborhood heterogeneity of a
# *continuous* raster (elevation, NDVI, canopy height, temperature, moisture,
# and similar surfaces). Unlike the landscape metrics in `landscape_metrics.R`,
# they require continuous values rather than integer-like class codes, so they
# are exposed through `surface_var()` rather than `landscape_var()` to keep the
# two raster contracts separate.
#
# Supported metrics:
#   `sa`  average roughness     = mean absolute deviation from the neighborhood
#                                 mean.
#   `sq`  RMS roughness         = sample standard deviation of the neighborhood
#                                 values (matches `geodiv::sq()`, N - 1).
#   `ssk` skewness              = adjusted skewness of the neighborhood values
#                                 (matches `geodiv::ssk(adj = TRUE)`).
#   `sku` kurtosis              = excess kurtosis of the neighborhood values
#                                 (matches `geodiv::sku(excess = TRUE)`).
#   `sdq` root mean square slope = RMS of the local forward-difference gradient
#                                 (matches `geodiv::sdq()` up to a unit factor;
#                                 see `.surface_slope_metric()`).
#
# `sa`, `sq`, `ssk`, and `sku` are pointwise functions of the neighborhood value
# distribution. `sdq` additionally needs each cell's neighbors, so (like the
# landscape edge/adjacency metrics) it requires cell IDs. All definitions
# intentionally match the `geodiv` package so that `geodiv` can be used as a
# validation oracle in the tests.


# Validate that buffered raster values are usable as a continuous surface.
#
# This is the surface-metric counterpart to
# `.landscape_validate_categorical_values()`, but with the opposite contract:
# any finite continuous value is acceptable. The check only rejects non-finite
# values (e.g., Inf) that would silently propagate through the moment formulas.
# NA values are allowed and dropped during computation.
.surface_validate_continuous_values <- function(values,
                                                metric,
                                                context = "Surface texture metric") {
  finite_values <- as.vector(values)
  finite_values <- finite_values[!is.na(finite_values)]

  if (length(finite_values) == 0) {
    return(invisible(NULL))
  }

  if (any(!is.finite(finite_values))) {
    stop(
      sprintf(
        "%s `%s` requires finite continuous raster values; found non-finite values such as Inf.",
        context,
        metric
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


# Adjusted skewness, matching geodiv::ssk(adj = TRUE): the third standardized
# moment (using the sample standard deviation, N - 1) times the bias-adjustment
# factor sqrt(N (N - 1)) / (N - 2). Returns NA when fewer than three finite
# values are available or the surface is flat (zero spread).
.surface_skewness <- function(z, n = length(z)) {
  if (n < 3) {
    return(NA_real_)
  }
  s <- stats::sd(z)
  if (!is.finite(s) || s == 0) {
    return(NA_real_)
  }
  zbar <- mean(z)
  val_unadj <- (sum((z - zbar)^3) / n) / s^3
  (sqrt(n * (n - 1)) / (n - 2)) * val_unadj
}


# Excess kurtosis, matching geodiv::sku(excess = TRUE): the fourth moment over N
# divided by the sample standard deviation (N - 1) to the fourth power, minus 3.
# Returns NA when fewer than two finite values are available or the surface is
# flat.
.surface_kurtosis <- function(z, n = length(z)) {
  if (n < 2) {
    return(NA_real_)
  }
  s <- stats::sd(z)
  if (!is.finite(s) || s == 0) {
    return(NA_real_)
  }
  zbar <- mean(z)
  (sum((z - zbar)^4) / n) / s^4 - 3
}


# Compute a single pointwise surface texture metric from a vector of continuous
# values. NA values are dropped before computation. Returns NA_real_ when there
# are too few finite values to define the metric.
.surface_metric <- function(values, metric) {
  metric <- match.arg(metric, c("sa", "sq", "ssk", "sku"))

  z <- values[!is.na(values)]
  n <- length(z)
  if (n == 0) {
    return(NA_real_)
  }

  switch(
    metric,
    sa = {
      zbar <- mean(z)
      sum(abs(z - zbar)) / n
    },
    sq = if (n < 2) NA_real_ else stats::sd(z),
    ssk = .surface_skewness(z, n),
    sku = .surface_kurtosis(z, n)
  )
}


# Root mean square slope from a reconstructed local grid (a matrix with NA
# outside the buffer). Uses forward differences to the right and down neighbors;
# a base cell contributes its squared gradient magnitude only if both neighbors
# are present (finite). The result is the *true* slope (rise over run), i.e.,
# the differences are divided by `resolution`. geodiv::sdq() hard-codes a cell
# spacing of 1 (slope in value units per pixel), so multiplying this value by
# `resolution` reproduces geodiv on the same grid. Returns NA when no base cell
# has both neighbors present.
.surface_slope_metric <- function(values, resolution) {
  validate_scalar_numeric(resolution, "resolution", positive = TRUE)
  if (!is.matrix(values)) {
    values <- as.matrix(values)
  }
  nr <- nrow(values)
  nc <- ncol(values)
  if (nr < 2 || nc < 2) {
    return(NA_real_)
  }

  # Forward differences aligned to the (nr - 1) x (nc - 1) base-cell grid.
  right_diff <- values[-nr, -nc, drop = FALSE] - values[-nr, -1, drop = FALSE]
  down_diff <- values[-nr, -nc, drop = FALSE] - values[-1, -nc, drop = FALSE]

  slope_sq <- (right_diff^2 + down_diff^2) / resolution^2
  valid <- is.finite(slope_sq)
  if (!any(valid)) {
    return(NA_real_)
  }

  sqrt(mean(slope_sq[valid]))
}


# Point-extraction path for the pointwise metrics (sa, sq, ssk, sku): compute
# the metric for every column of buffered raster values within `radius`.
# Mirrors `.landscape_composition_by_buffer()` so the dispatch in
# `.msr_eval_scale_vars()` is uniform across metric families. `validate` is set
# FALSE by the optimizer to skip the per-evaluation re-scan once the source
# raster has been checked at `kernel_prep()` time.
.surface_metric_by_buffer <- function(d,
                                      r_stack.df,
                                      radius,
                                      metric,
                                      validate = TRUE) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)

  if (length(d) == 0) {
    stop("`d` must contain at least one distance.", call. = FALSE)
  }

  values <- as.matrix(r_stack.df)
  storage.mode(values) <- "double"
  if (nrow(values) != length(d)) {
    stop("`r_stack.df` must have one row for each distance in `d`.", call. = FALSE)
  }
  if (isTRUE(validate)) {
    .surface_validate_continuous_values(
      values,
      metric = metric,
      context = "Surface texture metric"
    )
  }

  keep <- is.finite(d) & d <= radius
  out <- vapply(
    seq_len(ncol(values)),
    function(j) .surface_metric(values[keep, j], metric),
    numeric(1)
  )

  stats::setNames(out, colnames(values))
}


# Point-extraction path for the neighbor-based metric (sdq): reconstruct the
# local grid from the buffered cell IDs and compute the slope metric. Like the
# landscape edge/adjacency metrics, this needs `cells` and `n_cols`. `validate`
# behaves as in `.surface_metric_by_buffer()`.
.surface_slope_by_buffer <- function(d,
                                     r_stack.df,
                                     cells,
                                     radius,
                                     resolution,
                                     n_cols,
                                     metric = "sdq",
                                     validate = TRUE) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_numeric(resolution, "resolution", positive = TRUE)
  validate_scalar_numeric(n_cols, "n_cols", integerish = TRUE, positive = TRUE)

  if (length(d) == 0) {
    stop("`d` must contain at least one distance.", call. = FALSE)
  }
  if (length(cells) != length(d)) {
    stop("`cells` must have the same length as `d`.", call. = FALSE)
  }

  values <- as.matrix(r_stack.df)
  storage.mode(values) <- "double"
  if (nrow(values) != length(d)) {
    stop("`r_stack.df` must have one row for each distance in `d`.", call. = FALSE)
  }
  if (isTRUE(validate)) {
    .surface_validate_continuous_values(
      values,
      metric = metric,
      context = "Surface texture metric"
    )
  }

  keep <- is.finite(d) & d <= radius
  out <- vapply(
    seq_len(ncol(values)),
    function(j) {
      local_matrix <- .landscape_cells_to_matrix(
        values = values[keep, j],
        cells = cells[keep],
        n_cols = n_cols
      )
      .surface_slope_metric(local_matrix, resolution)
    },
    numeric(1)
  )

  stats::setNames(out, colnames(values))
}


# Raster-projection path for `sq`.
#
# Sq is the sample standard deviation within each circular neighborhood, which
# decomposes into neighborhood sums:
#   var = (sum(z^2) - sum(z)^2 / n) / (n - 1)
#   sq  = sqrt(var)
# Each of sum(z), sum(z^2), and n (the count of finite cells) is a circular
# convolution, so the whole surface is produced with three FFT passes. This
# matches the FFT idiom used by the kernel-smoothing and landscape-metric
# projections and is exact (verified to floating-point tolerance against the
# per-point computation). Cells with fewer than two finite neighbors return NA.
.surface_sq_raster_fft <- function(raster, radius, na.rm = TRUE) {
  resolution <- .landscape_validate_single_layer_raster(raster, "Surface")
  values <- terra::as.matrix(raster, wide = TRUE)
  .surface_validate_continuous_values(
    values,
    metric = "sq",
    context = "Surface texture projection"
  )

  kernel <- .landscape_circle_kernel(radius = radius, resolution = resolution)

  valid <- matrix(as.numeric(!is.na(values)),
                  nrow = nrow(values),
                  ncol = ncol(values))
  z0 <- values
  z0[is.na(z0)] <- 0

  n <- .landscape_clean_count_matrix(
    fft_convolution(valid, kernel, fun = "sum", na.rm = na.rm)
  )
  s1 <- fft_convolution(z0, kernel, fun = "sum", na.rm = na.rm)
  s2 <- fft_convolution(z0 * z0, kernel, fun = "sum", na.rm = na.rm)

  variance <- (s2 - (s1 * s1) / n) / (n - 1)
  # Guard against tiny negative variances from floating-point cancellation.
  variance[!is.finite(variance) | variance < 0] <- 0
  result <- sqrt(variance)
  result[!is.finite(n) | n < 2] <- NA_real_

  out <- terra::rast(raster)
  terra::values(out) <- as.vector(result)
  names(out) <- paste0(names(raster), "_sq")
  out
}


# Raster-projection path for the higher standardized moments `ssk` and `sku`.
#
# Both are functions of the neighborhood central moments, which decompose into
# convolved power sums (sum(z), sum(z^2), sum(z^3), sum(z^4), and the count).
# To control the catastrophic cancellation that raw power sums of large-mean
# surfaces would otherwise produce, the values are first centered on the global
# raster mean; skewness and kurtosis are shift-invariant, so this does not
# change the result but keeps the convolved sums small. Cells with too few
# finite neighbors (or a flat neighborhood) return NA.
.surface_moment_raster_fft <- function(raster, radius, metric, na.rm = TRUE) {
  metric <- match.arg(metric, c("ssk", "sku"))
  resolution <- .landscape_validate_single_layer_raster(raster, "Surface")
  values <- terra::as.matrix(raster, wide = TRUE)
  .surface_validate_continuous_values(
    values,
    metric = metric,
    context = "Surface texture projection"
  )

  kernel <- .landscape_circle_kernel(radius = radius, resolution = resolution)

  global_mean <- mean(values, na.rm = TRUE)
  if (!is.finite(global_mean)) {
    global_mean <- 0
  }
  centered <- values - global_mean

  valid <- matrix(as.numeric(!is.na(values)),
                  nrow = nrow(values),
                  ncol = ncol(values))
  z0 <- centered
  z0[is.na(z0)] <- 0

  n <- .landscape_clean_count_matrix(
    fft_convolution(valid, kernel, fun = "sum", na.rm = na.rm)
  )
  s1 <- fft_convolution(z0, kernel, fun = "sum", na.rm = na.rm)
  s2 <- fft_convolution(z0^2, kernel, fun = "sum", na.rm = na.rm)
  s3 <- fft_convolution(z0^3, kernel, fun = "sum", na.rm = na.rm)

  mu <- s1 / n
  variance <- (s2 - (s1 * s1) / n) / (n - 1)
  variance[!is.finite(variance) | variance < 0] <- 0
  s <- sqrt(variance)

  # Central third moment from the (centered) power sums.
  m3 <- s3 - 3 * mu * s2 + 2 * n * mu^3

  if (metric == "ssk") {
    result <- (sqrt(n * (n - 1)) / (n - 2)) * (m3 / n) / s^3
    result[!is.finite(n) | n < 3 | s == 0] <- NA_real_
  } else {
    s4 <- fft_convolution(z0^4, kernel, fun = "sum", na.rm = na.rm)
    m4 <- s4 - 4 * mu * s3 + 6 * mu^2 * s2 - 3 * n * mu^4
    result <- (m4 / n) / s^4 - 3
    result[!is.finite(n) | n < 2 | s == 0] <- NA_real_
  }
  result[!is.finite(result)] <- NA_real_

  out <- terra::rast(raster)
  terra::values(out) <- as.vector(result)
  names(out) <- paste0(names(raster), "_", metric)
  out
}


# Raster-projection path for `sdq` (RMS slope).
#
# Each cell's squared gradient magnitude (forward differences to the right and
# down neighbors, divided by resolution to give true slope) is precomputed on
# the full raster, then the circular neighborhood mean is taken via FFT and
# square-rooted. Because the per-cell slope uses each cell's actual raster
# neighbors, a window-edge base cell may use a neighbor just outside the
# circular window; the point-extraction path instead restricts to in-buffer
# neighbors, so the two paths agree closely but not bit-for-bit at boundaries
# (the same boundary behavior as the landscape edge metrics). Cells with no
# valid base cell in the window return NA.
.surface_sdq_raster_fft <- function(raster, radius, na.rm = TRUE) {
  resolution <- .landscape_validate_single_layer_raster(raster, "Surface")
  values <- terra::as.matrix(raster, wide = TRUE)
  .surface_validate_continuous_values(
    values,
    metric = "sdq",
    context = "Surface texture projection"
  )

  nr <- nrow(values)
  nc <- ncol(values)

  right_diff <- matrix(NA_real_, nrow = nr, ncol = nc)
  down_diff <- matrix(NA_real_, nrow = nr, ncol = nc)
  if (nc > 1) {
    right_diff[, -nc] <- values[, -nc, drop = FALSE] - values[, -1, drop = FALSE]
  }
  if (nr > 1) {
    down_diff[-nr, ] <- values[-nr, , drop = FALSE] - values[-1, , drop = FALSE]
  }

  # A base cell needs both forward neighbors, so slope_sq is NA in the last row
  # and column and wherever a neighbor is NA.
  slope_sq <- (right_diff^2 + down_diff^2) / resolution^2
  valid <- matrix(as.numeric(is.finite(slope_sq)), nrow = nr, ncol = nc)
  slope_sq[!is.finite(slope_sq)] <- 0

  kernel <- .landscape_circle_kernel(radius = radius, resolution = resolution)
  num <- fft_convolution(slope_sq, kernel, fun = "sum", na.rm = na.rm)
  den <- .landscape_clean_count_matrix(
    fft_convolution(valid, kernel, fun = "sum", na.rm = na.rm)
  )

  mean_slope_sq <- num / den
  mean_slope_sq[!is.finite(mean_slope_sq) | mean_slope_sq < 0] <- 0
  result <- sqrt(mean_slope_sq)
  result[!is.finite(den) | den < 1] <- NA_real_

  out <- terra::rast(raster)
  terra::values(out) <- as.vector(result)
  names(out) <- paste0(names(raster), "_sdq")
  out
}


# Raster-projection path for `sa`.
#
# Sa is the mean absolute deviation from the *neighborhood* mean. Because that
# mean changes with the output location, Sa is not a convolution and has no
# exact FFT decomposition (the absolute value couples every cell to the output
# location through that window's own mean). It is therefore projected with an
# exact masked focal pass, implemented in compiled code (`surface_sa_focal_cpp`)
# rather than via an R closure through `terra::focal()`, which is far faster
# while remaining exact. The circular kernel from `.landscape_circle_kernel()`
# supplies the in-window mask; NA and non-finite neighbors are dropped, and a
# window with no finite cells returns NA. `terra::as.matrix(wide = TRUE)` and
# `as.vector(t(...))` round-trip the raster so cell order is preserved.
.surface_sa_raster_cpp <- function(raster, radius, na.rm = TRUE) {
  resolution <- .landscape_validate_single_layer_raster(raster, "Surface")
  values <- terra::as.matrix(raster, wide = TRUE)
  .surface_validate_continuous_values(
    values,
    metric = "sa",
    context = "Surface texture projection"
  )

  kernel <- .landscape_circle_kernel(radius = radius, resolution = resolution)
  result <- surface_sa_focal_cpp(values, kernel)

  out <- terra::rast(raster)
  terra::values(out) <- as.vector(t(result))
  names(out) <- paste0(names(raster), "_sa")
  out
}


# --- Kernel-weighted surface moments -----------------------------------------
#
# The metrics above summarize a hard-edged circular neighborhood: every cell
# within `radius` counts equally, and the optimized parameter is that radius.
# The weighted variants below instead weight each cell by the model's distance
# kernel (gaussian, exp, expow) and optimize the kernel scale (sigma), exactly
# like a `kernel_var()` mean. This estimates the scale of effect of
# heterogeneity itself with a likelihood that is smooth in sigma, because the
# weights vary smoothly with sigma. Weighted skewness and kurtosis use the plain
# standardized central moments (the small-sample bias adjustment that
# `geodiv::ssk()` applies is a discrete-count correction that does not transfer
# to continuous weights). Weighted Sq uses the weighted root-mean-square
# deviation (population-style, normalized by the sum of weights).

# Compute a single weighted surface moment from aligned value and weight
# vectors. Cells that are NA, or carry a non-positive/non-finite weight, are
# dropped. Returns NA when no positive weight remains or the weighted spread is
# zero (for the standardized higher moments).
.surface_weighted_metric <- function(values, weights, metric) {
  fin <- is.finite(values) & is.finite(weights) & weights > 0
  v <- values[fin]
  w <- weights[fin]
  s0 <- sum(w)
  if (length(v) == 0 || !is.finite(s0) || s0 <= 0) {
    return(NA_real_)
  }

  mu <- sum(w * v) / s0
  dev <- v - mu
  m2 <- sum(w * dev^2) / s0

  switch(
    metric,
    sa = sum(w * abs(dev)) / s0,
    sq = if (m2 <= 0) 0 else sqrt(m2),
    ssk = if (m2 <= 0) NA_real_ else (sum(w * dev^3) / s0) / m2^1.5,
    sku = if (m2 <= 0) NA_real_ else (sum(w * dev^4) / s0) / m2^2 - 3
  )
}


# Point-extraction path for weighted moments: weight every buffered cell by the
# model kernel at the current sigma (and shape, for expow) and reduce to the
# weighted moment. `sigma`/`shape` are in the same internal (scaled) units as
# `d`, exactly as the kernel-mean path receives them.
.surface_weighted_by_buffer <- function(d,
                                        r_stack.df,
                                        kernel,
                                        sigma,
                                        shape,
                                        metric,
                                        validate = TRUE) {
  validate_scalar_numeric(sigma, "sigma", positive = TRUE)

  if (length(d) == 0) {
    stop("`d` must contain at least one distance.", call. = FALSE)
  }

  values <- as.matrix(r_stack.df)
  storage.mode(values) <- "double"
  if (nrow(values) != length(d)) {
    stop("`r_stack.df` must have one row for each distance in `d`.", call. = FALSE)
  }
  if (isTRUE(validate)) {
    .surface_validate_continuous_values(
      values,
      metric = metric,
      context = "Surface texture metric"
    )
  }

  weights <- .msr_kernel_shape_weight(
    d = d,
    sigma = sigma,
    shape = if (identical(kernel, "expow")) shape else NULL,
    kernel = kernel
  )

  out <- vapply(
    seq_len(ncol(values)),
    function(j) .surface_weighted_metric(values[, j], weights, metric),
    numeric(1)
  )

  stats::setNames(out, colnames(values))
}


# Build the kernel weight matrix used to project a weighted surface metric. This
# mirrors the focal-window construction in `.msr_kernel_raster_one()`: a window
# wide enough to hold `pct_wt` of the kernel mass, filled with the distance-
# decay weights for the given kernel, sigma, and shape.
.surface_kernel_weight_matrix <- function(raster, sigma, shape, kernel, pct_wt) {
  mx <- kernel_dist(kernel = kernel, sigma = sigma, beta = shape, prob = pct_wt)
  r_res <- terra::res(raster)[1]
  focal_d <- ceiling(mx / r_res) * 2
  if ((focal_d %% 2) == 0) {
    focal_d <- focal_d + 1
  }

  r_wt <- terra::rast(matrix(0, focal_d, focal_d))
  terra::crs(r_wt) <- terra::crs(raster)
  cntr_crd <- terra::xyFromCell(r_wt, ceiling(focal_d^2 / 2))
  cell_crds <- terra::crds(r_wt)
  d_vec <- fields::rdist(cntr_crd, cell_crds)[1, ] * r_res
  terra::values(r_wt) <- scale_type_r(d = d_vec,
                                      kernel = kernel,
                                      sigma = sigma,
                                      shape = shape,
                                      output = "wts")
  terra::as.matrix(r_wt, wide = TRUE)
}


# Raster-projection path for the weighted moments. `sq`, `ssk`, and `sku`
# decompose into weight-convolved power sums (S0 = sum of weights, S1..S4 =
# weighted power sums), centered on the global raster mean for numerical
# stability. `sa` again has no convolution form, so it reuses the compiled
# weighted focal pass with the kernel weight matrix.
.surface_weighted_raster_fft <- function(raster,
                                         kernel,
                                         sigma,
                                         shape,
                                         metric,
                                         pct_wt = 0.975,
                                         na.rm = TRUE) {
  metric <- match.arg(metric, c("sa", "sq", "ssk", "sku"))
  .landscape_validate_single_layer_raster(raster, "Surface")
  values <- terra::as.matrix(raster, wide = TRUE)
  .surface_validate_continuous_values(
    values,
    metric = metric,
    context = "Surface texture projection"
  )

  wt_mat <- .surface_kernel_weight_matrix(
    raster = raster,
    sigma = sigma,
    shape = shape,
    kernel = kernel,
    pct_wt = pct_wt
  )

  out <- terra::rast(raster)

  if (metric == "sa") {
    result <- surface_sa_focal_cpp(values, wt_mat)
    terra::values(out) <- as.vector(t(result))
    names(out) <- paste0(names(raster), "_sa")
    return(out)
  }

  global_mean <- mean(values, na.rm = TRUE)
  if (!is.finite(global_mean)) {
    global_mean <- 0
  }
  centered <- values - global_mean
  valid <- matrix(as.numeric(!is.na(values)),
                  nrow = nrow(values),
                  ncol = ncol(values))
  z0 <- centered
  z0[is.na(z0)] <- 0

  s0 <- fft_convolution(valid, wt_mat, fun = "sum", na.rm = na.rm)
  s1 <- fft_convolution(z0, wt_mat, fun = "sum", na.rm = na.rm)
  s2 <- fft_convolution(z0^2, wt_mat, fun = "sum", na.rm = na.rm)

  mu <- s1 / s0
  m2 <- s2 / s0 - mu^2
  m2[!is.finite(m2) | m2 < 0] <- 0

  if (metric == "sq") {
    result <- sqrt(m2)
  } else if (metric == "ssk") {
    s3 <- fft_convolution(z0^3, wt_mat, fun = "sum", na.rm = na.rm)
    m3 <- s3 / s0 - 3 * mu * (s2 / s0) + 2 * mu^3
    result <- m3 / m2^1.5
  } else {
    s3 <- fft_convolution(z0^3, wt_mat, fun = "sum", na.rm = na.rm)
    s4 <- fft_convolution(z0^4, wt_mat, fun = "sum", na.rm = na.rm)
    m4 <- s4 / s0 - 4 * mu * (s3 / s0) + 6 * mu^2 * (s2 / s0) - 3 * mu^4
    result <- m4 / m2^2 - 3
  }

  result[!is.finite(s0) | s0 <= 0 | m2 <= 0] <- NA_real_
  result[!is.finite(result)] <- NA_real_

  terra::values(out) <- as.vector(result)
  names(out) <- paste0(names(raster), "_", metric)
  out
}


# Projection dispatcher used by `.msr_scale_vars_raster()`. `sq` and the
# higher moments (`ssk`, `sku`) use FFT moment convolutions; `sdq` uses an FFT
# slope convolution; `sa` uses the compiled masked focal path.
.surface_metric_raster_fft <- function(raster, radius, metric, na.rm = TRUE) {
  metric <- match.arg(metric, c("sa", "sq", "ssk", "sku", "sdq"))
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_logical(na.rm, "na.rm")

  switch(
    metric,
    sq = .surface_sq_raster_fft(raster = raster, radius = radius, na.rm = na.rm),
    sa = .surface_sa_raster_cpp(raster = raster, radius = radius, na.rm = na.rm),
    ssk = .surface_moment_raster_fft(raster = raster, radius = radius,
                                     metric = "ssk", na.rm = na.rm),
    sku = .surface_moment_raster_fft(raster = raster, radius = radius,
                                     metric = "sku", na.rm = na.rm),
    sdq = .surface_sdq_raster_fft(raster = raster, radius = radius, na.rm = na.rm)
  )
}

# Internal helpers for continuous surface texture metric support.
#
# Surface texture metrics summarize the within-neighborhood heterogeneity of a
# *continuous* raster (elevation, NDVI, canopy height, temperature, moisture,
# and similar surfaces). Unlike the landscape metrics in `landscape_metrics.R`,
# they require continuous values rather than integer-like class codes, so they
# are exposed through `surface_var()` rather than `landscape_var()` to keep the
# two raster contracts separate.
#
# Tier 1 implements two metrics:
#   `sa` average roughness      = mean absolute deviation from the neighborhood
#                                 mean.
#   `sq` root mean square (RMS) = sample standard deviation of the neighborhood
#        roughness                values (matches `geodiv::sq()`, which calls
#                                 `stats::sd()` with an N - 1 denominator).
#
# Both definitions intentionally match the `geodiv` package so that `geodiv`
# can be used as a validation oracle in the tests.


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


# Compute a single surface texture metric from a vector of continuous values.
# NA values are dropped before computation. Returns NA_real_ when there are too
# few finite values to define the metric (0 values for `sa`; fewer than 2 for
# `sq`, since a sample standard deviation is undefined for a single value).
.surface_metric <- function(values, metric) {
  metric <- match.arg(metric, c("sa", "sq"))

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
    sq = if (n < 2) NA_real_ else stats::sd(z)
  )
}


# Point-extraction path: compute a surface metric for every column of buffered
# raster values within `radius`. Mirrors `.landscape_composition_by_buffer()`
# so that the dispatch in `.msr_eval_scale_vars()` is uniform across metric
# families. `validate` is set FALSE by the optimizer to skip the per-evaluation
# re-scan once the source raster has been checked at `kernel_prep()` time.
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


# Projection dispatcher used by `.msr_scale_vars_raster()`. `sq` uses the FFT
# moment path; `sa` uses the masked focal path.
.surface_metric_raster_fft <- function(raster, radius, metric, na.rm = TRUE) {
  metric <- match.arg(metric, c("sa", "sq"))
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_logical(na.rm, "na.rm")

  if (metric == "sq") {
    .surface_sq_raster_fft(raster = raster, radius = radius, na.rm = na.rm)
  } else {
    .surface_sa_raster_cpp(raster = raster, radius = radius, na.rm = na.rm)
  }
}

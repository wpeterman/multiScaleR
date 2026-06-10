# Tests for continuous surface texture metrics (`surface_var()`): `sa` (average
# roughness) and `sq` (RMS roughness = sample standard deviation).

# A continuous, textured test surface on a square 10 m grid. Point coordinates
# are placed on cell centers (values ending in 5) so that the FFT/focal
# projection neighborhood and the exact point-buffer neighborhood select the
# same cells, making the two paths directly comparable.
surface_test_raster <- function(seed, name = "elevation") {
  set.seed(seed)
  r <- terra::rast(
    nrows = 80,
    ncols = 80,
    xmin = 0,
    xmax = 800,
    ymin = 0,
    ymax = 800,
    crs = "EPSG:3857"
  )
  n <- terra::ncell(r)
  trend <- seq(0, 20, length.out = n)
  terra::values(r) <- 50 + trend + stats::rnorm(n, 0, 5)
  names(r) <- name
  r
}


surface_test_points <- function(raster) {
  terra::vect(
    cbind(c(205, 305, 405, 505, 605), c(205, 585, 425, 245, 615)),
    crs = terra::crs(raster)
  )
}


# Reference metric computed directly from the raster: select all cells whose
# centers fall within `radius` of each point, then apply the metric definition.
surface_cached_metric <- function(raster, pts, radius, metric) {
  cc <- terra::crds(raster)
  vals <- terra::values(raster)[, 1]
  xy <- terra::crds(pts)

  vapply(
    seq_len(nrow(xy)),
    function(i) {
      d <- sqrt((cc[, 1] - xy[i, 1])^2 + (cc[, 2] - xy[i, 2])^2)
      .surface_metric(vals[d <= radius], metric)
    },
    numeric(1)
  )
}


test_that("surface metric definitions match hand-computed values", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(.surface_metric(x, "sa"), 1.2)
  expect_equal(.surface_metric(x, "sq"), stats::sd(x))
  expect_equal(.surface_metric(x, "sq"), sqrt(2.5))

  # NA values are dropped before computation.
  xn <- c(1, 2, NA, 4, 5)
  expect_equal(.surface_metric(xn, "sa"), 1.5)
  expect_equal(.surface_metric(xn, "sq"), stats::sd(c(1, 2, 4, 5)))

  # A single value has zero average roughness and an undefined sample SD.
  expect_equal(.surface_metric(7, "sa"), 0)
  expect_true(is.na(.surface_metric(7, "sq")))

  # An empty (all-NA) neighborhood returns NA for both metrics.
  expect_true(is.na(.surface_metric(c(NA, NA), "sa")))
  expect_true(is.na(.surface_metric(numeric(0), "sq")))
})


test_that("surface metrics match geodiv point sampling", {
  skip_if_not_installed("geodiv")

  elevation <- surface_test_raster(101)
  pts <- surface_test_points(elevation)
  radius <- 120
  cc <- terra::crds(elevation)
  vals <- terra::values(elevation)[, 1]
  xy <- terra::crds(pts)

  for (metric in c("sa", "sq")) {
    geodiv_fun <- if (metric == "sa") geodiv::sa else geodiv::sq

    reference <- vapply(
      seq_len(nrow(xy)),
      function(i) {
        d <- sqrt((cc[, 1] - xy[i, 1])^2 + (cc[, 2] - xy[i, 2])^2)
        masked <- elevation
        mv <- vals
        mv[d > radius] <- NA
        terra::values(masked) <- mv
        geodiv_fun(masked)
      },
      numeric(1)
    )

    cached <- surface_cached_metric(elevation, pts, radius, metric)
    expect_equal(cached, reference, tolerance = 1e-6)
  }
})


test_that("FFT/focal surface projections agree with cached point metrics", {
  elevation <- surface_test_raster(202)
  pts <- surface_test_points(elevation)
  radius <- 95

  for (metric in c("sa", "sq")) {
    cached <- surface_cached_metric(elevation, pts, radius, metric)

    metric_raster <- .surface_metric_raster_fft(
      elevation,
      radius = radius,
      metric = metric
    )
    projected <- terra::extract(metric_raster, pts)[, 2]

    expect_s4_class(metric_raster, "SpatRaster")
    expect_equal(names(metric_raster), paste0("elevation_", metric))
    expect_equal(projected, cached, tolerance = 1e-6)
  }
})


test_that("compiled sa projection matches an independent focal reference", {
  # Use a non-square grid (nrow != ncol) so a row/column transpose bug in the
  # compiled focal pass would surface as a large disagreement rather than pass
  # silently.
  elevation <- terra::rast(
    nrows = 50,
    ncols = 70,
    xmin = 0,
    xmax = 700,
    ymin = 0,
    ymax = 500,
    crs = "EPSG:3857"
  )
  set.seed(717)
  terra::values(elevation) <- 40 + stats::rnorm(terra::ncell(elevation), 0, 6)
  names(elevation) <- "elevation"
  radius <- 65

  compiled <- .surface_metric_raster_fft(elevation, radius = radius, metric = "sa")

  # Independent reference: terra::focal with a masked window and an R closure.
  mask <- .landscape_circle_kernel(radius = radius,
                                   resolution = terra::res(elevation)[[1]])
  mask[mask == 0] <- NA_real_
  reference <- terra::focal(
    elevation,
    w = mask,
    na.policy = "all",
    fillvalue = NA,
    fun = function(x, ...) {
      z <- x[!is.na(x)]
      if (length(z) == 0) NA_real_ else sum(abs(z - mean(z))) / length(z)
    }
  )

  expect_equal(names(compiled), "elevation_sa")
  expect_equal(terra::values(compiled)[, 1],
               terra::values(reference)[, 1],
               tolerance = 1e-8)
})


test_that("compiled by-buffer surface metrics match the metric definition", {
  elevation <- surface_test_raster(303)
  pts <- surface_test_points(elevation)
  pts_sf <- sf::st_as_sf(pts)
  point_xy <- sf::st_coordinates(pts_sf)
  radius <- 110

  extracted <- exactextractr::exact_extract(
    elevation,
    sf::st_buffer(pts_sf, dist = radius),
    include_xy = TRUE,
    progress = FALSE
  )

  for (i in seq_along(extracted)) {
    d <- fields::rdist(
      matrix(point_xy[i, ], nrow = 1),
      as.matrix(extracted[[i]][, c("x", "y")])
    )[1, ]
    value_col <- setdiff(names(extracted[[i]]),
                         c("x", "y", "cell", "coverage_fraction"))[[1]]
    values <- extracted[[i]][, value_col, drop = FALSE]
    keep <- d <= radius

    for (metric in c("sa", "sq")) {
      expect_equal(
        .surface_metric_by_buffer(d, values, radius, metric)[[1]],
        .surface_metric(values[[1]][keep], metric)
      )
    }
  }
})


test_that("surface by-buffer helper validates inputs and handles empty buffers", {
  expect_error(
    .surface_metric_by_buffer(
      d = c(1, 2),
      r_stack.df = matrix(1:3, ncol = 1),
      radius = 1,
      metric = "sq"
    ),
    "one row for each distance"
  )

  # No cells within the radius -> NA, with the source column name preserved.
  out <- .surface_metric_by_buffer(
    d = c(10, 20),
    r_stack.df = data.frame(elev = c(1.5, 2.5)),
    radius = 1,
    metric = "sa"
  )
  expect_named(out, "elev")
  expect_true(is.na(out[["elev"]]))

  # Non-finite raster values are rejected.
  expect_error(
    .surface_metric_by_buffer(
      d = c(1, 2),
      r_stack.df = data.frame(elev = c(1, Inf)),
      radius = 2,
      metric = "sq"
    ),
    "finite continuous"
  )
})


test_that("scale variable specs support surface metrics from a continuous raster", {
  elevation <- surface_test_raster(404)
  pts <- surface_test_points(elevation)
  radius <- 95

  vars <- msr_vars(
    elev_sa = surface_var("elevation", metric = "sa", radius = radius),
    elev_sq = surface_var("elevation", metric = "sq", radius = radius)
  )

  expect_s3_class(vars, "multiScaleR_vars")
  expect_equal(vars$type, c("surface", "surface"))
  expect_error(surface_var("elevation", metric = "bogus"))

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = elevation,
    max_D = radius,
    scale_vars = vars,
    verbose = FALSE
  )

  expect_named(kernel_inputs$kernel_dat, c("elev_sa", "elev_sq"))
  # Fixed-radius covariates are not optimized.
  expect_equal(kernel_inputs$n_covs, 0)

  for (metric in c("sa", "sq")) {
    cov <- paste0("elev_", metric)
    cached <- surface_cached_metric(elevation, pts, radius, metric)
    unscaled <- kernel_inputs$kernel_dat[[cov]] *
      kernel_inputs$scl_params$sd[[cov]] +
      kernel_inputs$scl_params$mean[[cov]]
    expect_equal(as.numeric(unscaled), cached, tolerance = 1e-6)
  }

  projected <- kernel_scale.raster(
    raster_stack = elevation,
    scale_vars = vars,
    verbose = FALSE
  )
  expect_named(projected, c("elev_sa", "elev_sq"))
  for (metric in c("sa", "sq")) {
    cov <- paste0("elev_", metric)
    cached <- surface_cached_metric(elevation, pts, radius, metric)
    expect_equal(terra::extract(projected[[cov]], pts)[, 2], cached,
                 tolerance = 1e-6)
  }
})


test_that("optimized surface specs evaluate during model refits", {
  elevation <- surface_test_raster(505)
  pts <- surface_test_points(elevation)
  vars <- msr_vars(
    elev_sq = surface_var("elevation", metric = "sq")
  )

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = elevation,
    max_D = 120,
    scale_vars = vars,
    verbose = FALSE
  )
  expect_equal(kernel_inputs$n_covs, 1)
  expect_true(kernel_inputs$cell_data_stored)

  dat <- data.frame(
    y = c(1, 0, 1, 0, 1),
    kernel_inputs$kernel_dat
  )
  mod <- glm(y ~ elev_sq, family = binomial(), data = dat)
  opt_context <- build_opt_context(
    fitted_mod = mod,
    cov_df = kernel_inputs$raw_cov,
    scale_vars = kernel_inputs$scale_vars,
    unit_conv = kernel_inputs$unit_conv,
    resolution = kernel_inputs$resolution,
    n_cols = kernel_inputs$n_cols
  )

  neg_ll <- kernel_scale_fn(
    par = 0.5,
    d_list = kernel_inputs$d_list,
    cov_df = kernel_inputs$raw_cov,
    kernel = kernel_inputs$kernel,
    fitted_mod = mod,
    opt_context = opt_context
  )
  final <- kernel_scale_fn(
    par = 0.5,
    d_list = kernel_inputs$d_list,
    cov_df = kernel_inputs$raw_cov,
    kernel = kernel_inputs$kernel,
    fitted_mod = mod,
    opt_context = opt_context,
    mod_return = TRUE
  )

  expect_true(is.finite(neg_ll))
  expect_s3_class(final$mod, "glm")
  expect_named(final$scl_params$mean, "elev_sq")
})


test_that("optimized surface specs are rejected with the expow kernel", {
  elevation <- surface_test_raster(606)
  vars <- msr_vars(
    elev_sq = surface_var("elevation", metric = "sq")
  )
  expect_error(
    .msr_validate_scale_vars(vars, elevation, kernel = "expow"),
    "expow"
  )
})


test_that("kernel preparation rejects flat (zero-variance) surfaces", {
  elevation <- terra::rast(
    nrows = 20,
    ncols = 20,
    xmin = 0,
    xmax = 200,
    ymin = 0,
    ymax = 200,
    crs = "EPSG:3857"
  )
  terra::values(elevation) <- 5
  names(elevation) <- "elevation"
  pts <- terra::vect(cbind(c(55, 85, 115), c(55, 85, 115)),
                     crs = terra::crs(elevation))
  vars <- msr_vars(
    elev_sq = surface_var("elevation", metric = "sq", radius = 40)
  )

  expect_error(
    kernel_prep(
      pts = pts,
      raster_stack = elevation,
      max_D = 40,
      scale_vars = vars,
      verbose = FALSE
    ),
    "zero variance"
  )
})


test_that("surface projections require square raster cells", {
  elevation <- terra::rast(
    nrows = 20,
    ncols = 40,
    xmin = 0,
    xmax = 800,
    ymin = 0,
    ymax = 200,
    crs = "EPSG:3857"
  )
  terra::values(elevation) <- stats::rnorm(terra::ncell(elevation))
  names(elevation) <- "elevation"
  pts <- terra::vect(cbind(c(105, 205), c(95, 105)),
                     crs = terra::crs(elevation))
  vars <- msr_vars(
    elev_sq = surface_var("elevation", metric = "sq", radius = 60)
  )

  expect_error(
    kernel_prep(
      pts = pts,
      raster_stack = elevation,
      max_D = 60,
      scale_vars = vars,
      verbose = FALSE
    ),
    "square raster cells"
  )
})


# --- Tier 2 metrics: ssk (skewness), sku (kurtosis), sdq (RMS slope) ---------

# Reconstruct the local grid for each point's circular buffer and apply the
# slope metric directly, as a reference for sdq.
surface_cached_sdq <- function(raster, pts, radius) {
  cc <- terra::crds(raster)
  vals <- terra::values(raster)[, 1]
  xy <- terra::crds(pts)
  resolution <- terra::res(raster)[[1]]
  n_cols <- terra::ncol(raster)

  vapply(
    seq_len(nrow(xy)),
    function(i) {
      d <- sqrt((cc[, 1] - xy[i, 1])^2 + (cc[, 2] - xy[i, 2])^2)
      keep <- d <= radius
      local_matrix <- .landscape_cells_to_matrix(
        values = vals[keep],
        cells = terra::cellFromXY(raster, cc[keep, , drop = FALSE]),
        n_cols = n_cols
      )
      .surface_slope_metric(local_matrix, resolution)
    },
    numeric(1)
  )
}


test_that("ssk and sku match geodiv point sampling", {
  skip_if_not_installed("geodiv")

  elevation <- surface_test_raster(111)
  pts <- surface_test_points(elevation)
  radius <- 120
  cc <- terra::crds(elevation)
  vals <- terra::values(elevation)[, 1]
  xy <- terra::crds(pts)

  for (metric in c("ssk", "sku")) {
    geodiv_fun <- if (metric == "ssk") geodiv::ssk else geodiv::sku

    reference <- vapply(
      seq_len(nrow(xy)),
      function(i) {
        d <- sqrt((cc[, 1] - xy[i, 1])^2 + (cc[, 2] - xy[i, 2])^2)
        masked <- elevation
        mv <- vals
        mv[d > radius] <- NA
        terra::values(masked) <- mv
        geodiv_fun(masked)
      },
      numeric(1)
    )

    cached <- surface_cached_metric(elevation, pts, radius, metric)
    expect_equal(cached, reference, tolerance = 1e-6)
  }
})


test_that("sdq matches geodiv on a full window, scaled by resolution", {
  skip_if_not_installed("geodiv")

  set.seed(321)
  block <- terra::rast(
    nrows = 25,
    ncols = 30,
    xmin = 0,
    xmax = 300,
    ymin = 0,
    ymax = 250,
    crs = "EPSG:3857"
  )
  terra::values(block) <- 40 + stats::rnorm(terra::ncell(block), 0, 6)
  resolution <- terra::res(block)[[1]]

  # geodiv::sdq uses a cell spacing of 1; .surface_slope_metric returns true
  # slope (divided by resolution), so multiplying by resolution recovers geodiv.
  mine <- .surface_slope_metric(terra::as.matrix(block, wide = TRUE), resolution)
  expect_equal(mine * resolution, geodiv::sdq(block), tolerance = 1e-6)
})


test_that("sdq equals the known slope of a tilted plane", {
  resolution <- 10
  slope <- 0.7
  # z increases by slope * resolution per column step, so the true gradient is
  # `slope` everywhere and sdq must equal it.
  z <- outer(seq_len(30), seq_len(30), function(i, j) slope * j * resolution)
  expect_equal(.surface_slope_metric(z, resolution), slope, tolerance = 1e-8)

  # A flat surface has zero slope.
  expect_equal(.surface_slope_metric(matrix(5, 8, 8), resolution), 0)
})


test_that("ssk and sku FFT projections agree with cached point metrics", {
  # Large-mean surface to confirm the moment projection is numerically stable.
  elevation <- surface_test_raster(212)
  terra::values(elevation) <- terra::values(elevation) + 500
  pts <- surface_test_points(elevation)
  radius <- 95

  for (metric in c("ssk", "sku")) {
    cached <- surface_cached_metric(elevation, pts, radius, metric)
    metric_raster <- .surface_metric_raster_fft(elevation, radius = radius,
                                                metric = metric)
    projected <- terra::extract(metric_raster, pts)[, 2]

    expect_equal(names(metric_raster), paste0("elevation_", metric))
    expect_equal(projected, cached, tolerance = 1e-5)
  }
})


test_that("sdq FFT projection agrees with the point metric within boundary tolerance", {
  elevation <- surface_test_raster(213)
  pts <- surface_test_points(elevation)
  radius <- 105

  cached <- surface_cached_sdq(elevation, pts, radius)
  metric_raster <- .surface_metric_raster_fft(elevation, radius = radius,
                                              metric = "sdq")
  projected <- terra::extract(metric_raster, pts)[, 2]

  expect_equal(names(metric_raster), "elevation_sdq")
  # The point path restricts to in-buffer neighbors while the projection uses
  # each cell's actual raster neighbors, so the two agree to a few percent
  # (boundary behavior shared with the landscape edge metrics), not bit-for-bit.
  expect_equal(projected, cached, tolerance = 0.05)
})


test_that("Tier 2 surface specs integrate through kernel_prep and projection", {
  elevation <- surface_test_raster(214)
  pts <- surface_test_points(elevation)
  radius <- 95

  vars <- msr_vars(
    elev_ssk = surface_var("elevation", metric = "ssk", radius = radius),
    elev_sku = surface_var("elevation", metric = "sku", radius = radius),
    elev_sdq = surface_var("elevation", metric = "sdq", radius = radius)
  )

  # sdq needs cell IDs; ssk/sku do not. kernel_prep must extract them.
  expect_true(.msr_needs_cells(vars))

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = elevation,
    max_D = radius,
    scale_vars = vars,
    verbose = FALSE
  )
  expect_named(kernel_inputs$kernel_dat, c("elev_ssk", "elev_sku", "elev_sdq"))
  expect_equal(kernel_inputs$n_covs, 0)

  # Point-extracted ssk/sku reproduce the direct cached values exactly.
  for (metric in c("ssk", "sku")) {
    cov <- paste0("elev_", metric)
    cached <- surface_cached_metric(elevation, pts, radius, metric)
    unscaled <- kernel_inputs$kernel_dat[[cov]] *
      kernel_inputs$scl_params$sd[[cov]] +
      kernel_inputs$scl_params$mean[[cov]]
    expect_equal(as.numeric(unscaled), cached, tolerance = 1e-6)
  }

  # sdq extracted via cell reconstruction reproduces the cached slope.
  cached_sdq <- surface_cached_sdq(elevation, pts, radius)
  unscaled_sdq <- kernel_inputs$kernel_dat$elev_sdq *
    kernel_inputs$scl_params$sd[["elev_sdq"]] +
    kernel_inputs$scl_params$mean[["elev_sdq"]]
  expect_equal(as.numeric(unscaled_sdq), cached_sdq, tolerance = 1e-6)

  projected <- kernel_scale.raster(
    raster_stack = elevation,
    scale_vars = vars,
    verbose = FALSE
  )
  expect_named(projected, c("elev_ssk", "elev_sku", "elev_sdq"))
})


test_that("higher-moment and slope helpers handle small or flat windows", {
  # Skewness needs at least three values; kurtosis at least two.
  expect_true(is.na(.surface_metric(c(1, 2), "ssk")))
  expect_true(is.na(.surface_metric(1, "sku")))
  # Flat windows have zero spread, so standardized moments are undefined.
  expect_true(is.na(.surface_metric(rep(3, 10), "ssk")))
  expect_true(is.na(.surface_metric(rep(3, 10), "sku")))

  # A perfectly symmetric distribution has exactly zero skewness.
  expect_equal(.surface_metric(c(-3, -1, 0, 1, 3), "ssk"), 0)

  # The slope helper needs at least a 2x2 grid.
  expect_true(is.na(.surface_slope_metric(matrix(1:3, nrow = 1), 10)))
})


# --- Tier 3a: kernel-weighted moments (weighted = TRUE) ----------------------

test_that("surface_var validates the weighted flag", {
  spec <- surface_var("elevation", metric = "sq", weighted = TRUE)
  expect_true(spec$weighted)
  expect_true(spec$optimize)

  vars <- msr_vars(w = surface_var("elevation", metric = "sq", weighted = TRUE))
  expect_true(vars$weighted)

  # Weighting applies to the distribution metrics, not the slope metric.
  expect_error(surface_var("elevation", metric = "sdq", weighted = TRUE),
               "not supported for the slope metric")
  # Weighted metrics optimize the kernel scale, so a fixed radius is rejected.
  expect_error(surface_var("elevation", metric = "sq", weighted = TRUE,
                           radius = 300),
               "must be NULL")
})


test_that("weighted moment computation matches a direct weighted moment", {
  set.seed(414)
  z <- stats::rnorm(400, 20, 4)
  d <- runif(400, 0, 300)
  sigma <- 90
  w <- exp(-(d^2) / (2 * sigma^2))

  s0 <- sum(w)
  mu <- sum(w * z) / s0
  m2 <- sum(w * (z - mu)^2) / s0
  manual <- c(
    sa = sum(w * abs(z - mu)) / s0,
    sq = sqrt(m2),
    ssk = (sum(w * (z - mu)^3) / s0) / m2^1.5,
    sku = (sum(w * (z - mu)^4) / s0) / m2^2 - 3
  )

  for (metric in c("sa", "sq", "ssk", "sku")) {
    expect_equal(.surface_weighted_metric(z, w, metric), manual[[metric]],
                 tolerance = 1e-10)
    # The by-buffer wrapper builds the same Gaussian weights from distances.
    expect_equal(
      .surface_weighted_by_buffer(d, matrix(z, ncol = 1), "gaussian",
                                  sigma, NULL, metric)[[1]],
      manual[[metric]],
      tolerance = 1e-10
    )
  }
})


test_that("weighted moment projections agree with point extraction", {
  elevation <- surface_test_raster(415)
  pts <- surface_test_points(elevation)
  sigma <- 90
  cc <- terra::crds(elevation)
  vals <- terra::values(elevation)[, 1]
  xy <- terra::crds(pts)

  for (metric in c("sa", "sq", "ssk", "sku")) {
    point_vals <- vapply(
      seq_len(nrow(xy)),
      function(i) {
        d <- sqrt((cc[, 1] - xy[i, 1])^2 + (cc[, 2] - xy[i, 2])^2)
        keep <- d <= 300
        .surface_weighted_by_buffer(d[keep], matrix(vals[keep], ncol = 1),
                                    "gaussian", sigma, NULL, metric)[[1]]
      },
      numeric(1)
    )

    metric_raster <- .surface_weighted_raster_fft(
      elevation, kernel = "gaussian", sigma = sigma, shape = NULL,
      metric = metric, pct_wt = 0.99
    )
    projected <- terra::extract(metric_raster, pts)[, 2]

    expect_equal(names(metric_raster), paste0("elevation_", metric))
    # The projection truncates the kernel tail at pct_wt while the point path
    # uses the full buffer, so they agree to a small tail tolerance (the same
    # behavior as kernel-weighted means). The dispersion metrics are compared
    # relatively; the standardized higher moments sit near zero, so they are
    # compared on an absolute scale.
    if (metric %in% c("sa", "sq")) {
      expect_equal(projected, point_vals, tolerance = 0.03)
    } else {
      expect_lt(max(abs(projected - point_vals)), 0.1)
    }
  }
})


test_that("weighted surface metrics are allowed with the expow kernel", {
  elevation <- surface_test_raster(416)
  weighted_vars <- msr_vars(
    elev_sqw = surface_var("elevation", metric = "sq", weighted = TRUE)
  )
  # Weighted metrics use the kernel (including its shape), so expow is allowed.
  expect_silent(.msr_validate_scale_vars(weighted_vars, elevation,
                                         kernel = "expow"))
  # Unweighted optimized surface metrics remain unsupported with expow.
  unweighted_vars <- msr_vars(
    elev_sq = surface_var("elevation", metric = "sq")
  )
  expect_error(.msr_validate_scale_vars(unweighted_vars, elevation,
                                        kernel = "expow"),
               "expow")
})


test_that("weighted surface metrics optimize and project end to end", {
  elevation <- surface_test_raster(417)
  pts <- surface_test_points(elevation)

  vars <- msr_vars(
    elev_sqw = surface_var("elevation", metric = "sq", weighted = TRUE)
  )
  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = elevation,
    max_D = 300,
    kernel = "gaussian",
    scale_vars = vars,
    verbose = FALSE
  )
  expect_equal(kernel_inputs$n_covs, 1)
  expect_false(.msr_needs_cells(kernel_inputs$scale_vars))
  expect_true(all(is.finite(kernel_inputs$kernel_dat$elev_sqw)))

  dat <- data.frame(y = c(1, 0, 1, 0, 1), kernel_inputs$kernel_dat)
  mod <- glm(y ~ elev_sqw, family = binomial(), data = dat)
  opt_context <- build_opt_context(
    fitted_mod = mod,
    cov_df = kernel_inputs$raw_cov,
    scale_vars = kernel_inputs$scale_vars,
    unit_conv = kernel_inputs$unit_conv,
    resolution = kernel_inputs$resolution,
    n_cols = kernel_inputs$n_cols
  )

  # The objective is finite and varies smoothly with the scale parameter.
  obj <- vapply(
    c(0.2, 0.3, 0.4, 0.5),
    function(p) {
      kernel_scale_fn(
        par = p,
        d_list = kernel_inputs$d_list,
        cov_df = kernel_inputs$raw_cov,
        kernel = kernel_inputs$kernel,
        fitted_mod = mod,
        opt_context = opt_context
      )
    },
    numeric(1)
  )
  expect_true(all(is.finite(obj)))

  projected <- kernel_scale.raster(
    raster_stack = elevation,
    sigma = 90,
    kernel = "gaussian",
    scale_vars = vars,
    verbose = FALSE
  )
  expect_named(projected, "elev_sqw")
  expect_true(any(is.finite(terra::values(projected)[, 1])))
})


# --- Tier 3b: sdr (surface area ratio, physical units) -----------------------

# Reference sdr for each point's circular buffer via grid reconstruction.
surface_cached_sdr <- function(raster, pts, radius) {
  cc <- terra::crds(raster)
  vals <- terra::values(raster)[, 1]
  xy <- terra::crds(pts)
  resolution <- terra::res(raster)[[1]]
  n_cols <- terra::ncol(raster)

  vapply(
    seq_len(nrow(xy)),
    function(i) {
      d <- sqrt((cc[, 1] - xy[i, 1])^2 + (cc[, 2] - xy[i, 2])^2)
      keep <- d <= radius
      local_matrix <- .landscape_cells_to_matrix(
        values = vals[keep],
        cells = terra::cellFromXY(raster, cc[keep, , drop = FALSE]),
        n_cols = n_cols
      )
      .surface_area_metric(local_matrix, resolution)
    },
    numeric(1)
  )
}


test_that("sdr equals the analytic surface area ratio of a tilted plane", {
  resolution <- 10
  for (slope in c(0.3, 0.7, 1.5)) {
    z <- outer(seq_len(30), seq_len(30),
               function(i, j) slope * j * resolution)
    expected <- (sqrt(1 + slope^2) - 1) * 100
    expect_equal(.surface_area_metric(z, resolution), expected, tolerance = 1e-8)
  }
  # A flat surface has zero excess surface area.
  expect_equal(.surface_area_metric(matrix(5, 8, 8), resolution), 0)
  # The helper needs at least a 2x2 grid.
  expect_true(is.na(.surface_area_metric(matrix(1:3, nrow = 1), resolution)))
})


test_that("sdr quad-area triangulation matches geodiv::surface_area", {
  skip_if_not_installed("geodiv")

  set.seed(733)
  block <- matrix(stats::rnorm(20 * 24, 50, 8), 20, 24)
  # geodiv::surface_area() rescales z to [0, 1] and uses unit cell spacing;
  # replicate those inputs so the triangulation formula can be compared directly.
  zs <- (block - min(block)) / (max(block) - min(block))
  nr <- nrow(zs)
  nc <- ncol(zs)
  mine <- sum(.surface_akl(
    zs[-nr, -nc], zs[-nr, -1], zs[-1, -nc], zs[-1, -1], 1, 1
  ))
  expect_equal(mine, geodiv::surface_area(terra::rast(block)), tolerance = 1e-6)
})


test_that("sdr FFT projection agrees with the point metric within boundary tolerance", {
  elevation <- surface_test_raster(734)
  pts <- surface_test_points(elevation)
  radius <- 105

  cached <- surface_cached_sdr(elevation, pts, radius)
  metric_raster <- .surface_metric_raster_fft(elevation, radius = radius,
                                              metric = "sdr")
  projected <- terra::extract(metric_raster, pts)[, 2]

  expect_equal(names(metric_raster), "elevation_sdr")
  expect_equal(projected, cached, tolerance = 0.05)
})


test_that("sdr integrates through kernel_prep and projection", {
  elevation <- surface_test_raster(735)
  pts <- surface_test_points(elevation)
  radius <- 95

  vars <- msr_vars(elev_sdr = surface_var("elevation", metric = "sdr",
                                          radius = radius))
  expect_true(.msr_needs_cells(vars))

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = elevation,
    max_D = radius,
    scale_vars = vars,
    verbose = FALSE
  )
  expect_named(kernel_inputs$kernel_dat, "elev_sdr")

  cached_sdr <- surface_cached_sdr(elevation, pts, radius)
  unscaled <- kernel_inputs$kernel_dat$elev_sdr *
    kernel_inputs$scl_params$sd[["elev_sdr"]] +
    kernel_inputs$scl_params$mean[["elev_sdr"]]
  expect_equal(as.numeric(unscaled), cached_sdr, tolerance = 1e-6)

  projected <- kernel_scale.raster(
    raster_stack = elevation,
    scale_vars = vars,
    verbose = FALSE
  )
  expect_named(projected, "elev_sdr")

  # sdr is a slope-family metric, so kernel weighting is not offered.
  expect_error(surface_var("elevation", metric = "sdr", weighted = TRUE),
               "not supported for the slope metric")
})

landscape_test_raster <- function(seed, values, name = "landcover") {
  set.seed(seed)
  landcover <- terra::rast(
    nrows = 80,
    ncols = 80,
    xmin = 0,
    xmax = 800,
    ymin = 0,
    ymax = 800,
    crs = "EPSG:3857"
  )
  terra::values(landcover) <- sample(values,
                                     terra::ncell(landcover),
                                     replace = TRUE)
  names(landcover) <- name
  landcover
}


landscape_test_points <- function(raster) {
  terra::vect(
    cbind(c(165, 265, 425, 615, 705), c(165, 585, 425, 245, 675)),
    crs = terra::crs(raster)
  )
}


landscape_value_column <- function(x) {
  x[, setdiff(names(x), c("x", "y", "cell", "coverage_fraction"))[[1]],
    drop = FALSE]
}


landscape_cached_metric <- function(raster, pts, radius, metric_fun) {
  pts_sf <- sf::st_as_sf(pts)
  point_xy <- sf::st_coordinates(pts_sf)

  extracted <- exactextractr::exact_extract(
    raster,
    sf::st_buffer(pts_sf, dist = radius),
    include_xy = TRUE,
    include_cell = TRUE,
    progress = FALSE
  )

  vapply(
    seq_along(extracted),
    function(i) {
      d <- fields::rdist(
        matrix(point_xy[i, ], nrow = 1),
        as.matrix(extracted[[i]][, c("x", "y")])
      )[1, ]

      metric_fun(
        d = d,
        values = landscape_value_column(extracted[[i]]),
        cells = extracted[[i]]$cell
      )
    },
    numeric(1)
  )
}


landscape_reference_values <- function(reference, metric) {
  rows <- reference[reference$metric == metric, ]
  rows$value[order(rows$plot_id)]
}


test_that("fixed-buffer composition metrics match landscapemetrics point sampling", {
  skip_if_not_installed("landscapemetrics")

  landcover <- landscape_test_raster(101, 1:5)
  pts <- landscape_test_points(landcover)
  pts_sf <- sf::st_as_sf(pts)
  radius <- 80
  metrics <- c("shdi", "shei", "sidi", "siei", "msidi", "msiei",
               "pr", "prd", "rpr", "ta")

  cached <- lapply(
    metrics,
    function(metric) {
      landscape_cached_metric(
        raster = landcover,
        pts = pts,
        radius = radius,
        metric_fun = function(d, values, cells) {
          .landscape_composition_by_buffer(
            d = d,
            r_stack.df = values,
            radius = radius,
            metric = metric,
            resolution = terra::res(landcover)[[1]],
            classes_max = 5
          )[[1]]
        }
      )
    }
  )
  names(cached) <- metrics

  reference <- landscapemetrics::sample_lsm(
    landscape = landcover,
    y = pts_sf,
    shape = "circle",
    size = radius,
    what = paste0("lsm_l_", metrics),
    classes_max = 5,
    verbose = FALSE,
    progress = FALSE
  )

  for (metric in metrics) {
    tolerance <- if (metric %in% c("prd")) {
      5
    } else if (metric %in% c("ta")) {
      0.3
    } else {
      0.05
    }
    expect_equal(
      cached[[metric]],
      landscape_reference_values(reference, metric),
      tolerance = tolerance
    )
  }
})


test_that("FFT composition projections agree with cached point metrics", {
  landcover <- landscape_test_raster(202, 1:5)
  pts <- landscape_test_points(landcover)
  radius <- 80
  metrics <- c("shdi", "shei", "sidi", "siei", "msidi", "msiei",
               "pr", "prd", "rpr", "ta")

  for (metric in metrics) {
    cached <- landscape_cached_metric(
      raster = landcover,
      pts = pts,
      radius = radius,
      metric_fun = function(d, values, cells) {
        .landscape_composition_by_buffer(
          d = d,
          r_stack.df = values,
          radius = radius,
          metric = metric,
          resolution = terra::res(landcover)[[1]],
          classes_max = 5
        )[[1]]
      }
    )

    metric_raster <- .landscape_composition_raster_fft(
      landcover,
      radius = radius,
      metric = metric,
      classes_max = 5
    )
    projected <- terra::extract(metric_raster, pts)[, 2]

    expect_s4_class(metric_raster, "SpatRaster")
    expect_equal(projected, cached, tolerance = 0.01)
  }
})


test_that("fixed-buffer edge metrics match landscapemetrics point sampling", {
  skip_if_not_installed("landscapemetrics")

  landcover <- landscape_test_raster(303, 0:1, name = "forest")
  pts <- landscape_test_points(landcover)
  pts_sf <- sf::st_as_sf(pts)
  radius <- 80
  metrics <- c("ed", "te", "lsi")

  cached <- lapply(
    metrics,
    function(metric) {
      landscape_cached_metric(
        raster = landcover,
        pts = pts,
        radius = radius,
        metric_fun = function(d, values, cells) {
          .landscape_edge_by_buffer(
            d = d,
            r_stack.df = values,
            cells = cells,
            radius = radius,
            resolution = terra::res(landcover)[[1]],
            n_cols = terra::ncol(landcover),
            metric = metric
          )[[1]]
        }
      )
    }
  )
  names(cached) <- metrics

  reference <- landscapemetrics::sample_lsm(
    landscape = landcover,
    y = pts_sf,
    shape = "circle",
    size = radius,
    what = paste0("lsm_l_", metrics),
    verbose = FALSE,
    progress = FALSE
  )

  expect_equal(cached$ed, landscape_reference_values(reference, "ed"),
               tolerance = 25)
  expect_equal(cached$te, landscape_reference_values(reference, "te"),
               tolerance = 50)
  expect_equal(cached$lsi, landscape_reference_values(reference, "lsi"),
               tolerance = 0.15)
})


test_that("FFT edge projections agree with cached point metrics", {
  landcover <- landscape_test_raster(404, 0:1, name = "forest")
  pts <- landscape_test_points(landcover)
  radius <- 80

  for (metric in c("ed", "te", "lsi")) {
    cached <- landscape_cached_metric(
      raster = landcover,
      pts = pts,
      radius = radius,
      metric_fun = function(d, values, cells) {
        .landscape_edge_by_buffer(
          d = d,
          r_stack.df = values,
          cells = cells,
          radius = radius,
          resolution = terra::res(landcover)[[1]],
          n_cols = terra::ncol(landcover),
          metric = metric
        )[[1]]
      }
    )

    metric_raster <- .landscape_edge_raster_fft(
      landcover,
      radius = radius,
      metric = metric
    )
    projected <- terra::extract(metric_raster, pts)[, 2]

    expect_s4_class(metric_raster, "SpatRaster")
    tolerance <- switch(metric, ed = 25, te = 50, lsi = 0.15)
    expect_equal(projected, cached, tolerance = tolerance)
  }
})


test_that("fixed-buffer adjacency metrics match landscapemetrics point sampling", {
  skip_if_not_installed("landscapemetrics")

  landcover <- landscape_test_raster(505, 1:4)
  pts <- landscape_test_points(landcover)
  pts_sf <- sf::st_as_sf(pts)
  radius <- 80
  cached <- lapply(
    c("ai", "pladj", "contag"),
    function(metric) {
      landscape_cached_metric(
        raster = landcover,
        pts = pts,
        radius = radius,
        metric_fun = function(d, values, cells) {
          .landscape_adjacency_by_buffer(
            d = d,
            r_stack.df = values,
            cells = cells,
            radius = radius,
            n_cols = terra::ncol(landcover),
            metric = metric
          )[[1]]
        }
      )
    }
  )
  names(cached) <- c("ai", "pladj", "contag")

  contag_reference <- landscape_cached_metric(
    raster = landcover,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      keep <- d <= radius
      local_matrix <- .landscape_cells_to_matrix(
        values = values[[1]][keep],
        cells = cells[keep],
        n_cols = terra::ncol(landcover)
      )
      landscapemetrics::lsm_l_contag(local_matrix, verbose = FALSE)$value
    }
  )

  reference <- landscapemetrics::sample_lsm(
    landscape = landcover,
    y = pts_sf,
    shape = "circle",
    size = radius,
    what = c("lsm_l_ai", "lsm_l_pladj"),
    verbose = FALSE,
    progress = FALSE
  )

  ai_reference <- landscape_reference_values(reference, "ai")
  expect_lt(max(abs(cached$ai - ai_reference), na.rm = TRUE), 2)
  expect_equal(cached$pladj, landscape_reference_values(reference, "pladj"),
               tolerance = 5)
  expect_equal(cached$contag, contag_reference, tolerance = 0.01)
})


test_that("FFT adjacency projections agree with cached point metrics", {
  landcover <- landscape_test_raster(606, 1:4)
  pts <- landscape_test_points(landcover)
  radius <- 80

  for (metric in c("ai", "pladj", "contag")) {
    cached <- landscape_cached_metric(
      raster = landcover,
      pts = pts,
      radius = radius,
      metric_fun = function(d, values, cells) {
        .landscape_adjacency_by_buffer(
          d = d,
          r_stack.df = values,
          cells = cells,
          radius = radius,
          n_cols = terra::ncol(landcover),
          metric = metric
        )[[1]]
      }
    )

    metric_raster <- .landscape_adjacency_raster_fft(
      landcover,
      radius = radius,
      metric = metric
    )
    projected <- terra::extract(metric_raster, pts)[, 2]

    expect_s4_class(metric_raster, "SpatRaster")
    tolerance <- switch(metric, ai = 1.5, pladj = 1.5, contag = 0.25)
    if (metric == "ai") {
      expect_lt(max(abs(projected - cached), na.rm = TRUE), 2.5)
    } else {
      expect_equal(projected, cached, tolerance = tolerance)
    }
  }
})


test_that("FFT edge and adjacency projections match rare-class point buffers", {
  values <- matrix(1, nrow = 7, ncol = 7)
  values[4, 4] <- 2
  landcover <- terra::rast(
    nrows = 7,
    ncols = 7,
    xmin = 0,
    xmax = 70,
    ymin = 0,
    ymax = 70,
    crs = "EPSG:3857"
  )
  terra::values(landcover) <- as.vector(values)
  names(landcover) <- "landcover"
  pts <- terra::vect(cbind(35, 35), crs = terra::crs(landcover))
  radius <- 10

  cached_ed <- landscape_cached_metric(
    raster = landcover,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        resolution = terra::res(landcover)[[1]],
        n_cols = terra::ncol(landcover),
        metric = "ed"
      )[[1]]
    }
  )
  cached_pladj <- landscape_cached_metric(
    raster = landcover,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        n_cols = terra::ncol(landcover),
        metric = "pladj"
      )[[1]]
    }
  )
  cached_contag <- landscape_cached_metric(
    raster = landcover,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        n_cols = terra::ncol(landcover),
        metric = "contag"
      )[[1]]
    }
  )

  projected_ed <- terra::extract(
    .landscape_edge_raster_fft(landcover, radius = radius, metric = "ed"),
    pts
  )[, 2]
  projected_pladj <- terra::extract(
    .landscape_adjacency_raster_fft(landcover, radius = radius, metric = "pladj"),
    pts
  )[, 2]
  projected_contag <- terra::extract(
    .landscape_adjacency_raster_fft(landcover, radius = radius, metric = "contag"),
    pts
  )[, 2]

  expect_equal(projected_ed, cached_ed, tolerance = 1e-8)
  expect_equal(projected_pladj, cached_pladj, tolerance = 1e-8)
  expect_equal(projected_contag, cached_contag, tolerance = 1e-8)
})


test_that("compiled by-buffer metrics match direct matrix formulas", {
  landcover <- landscape_test_raster(707, 1:4)
  pts <- landscape_test_points(landcover)
  pts_sf <- sf::st_as_sf(pts)
  point_xy <- sf::st_coordinates(pts_sf)
  radius <- 80

  extracted <- exactextractr::exact_extract(
    landcover,
    sf::st_buffer(pts_sf, dist = radius),
    include_xy = TRUE,
    include_cell = TRUE,
    progress = FALSE
  )

  for (i in seq_along(extracted)) {
    d <- fields::rdist(
      matrix(point_xy[i, ], nrow = 1),
      as.matrix(extracted[[i]][, c("x", "y")])
    )[1, ]
    values <- landscape_value_column(extracted[[i]])
    cells <- extracted[[i]]$cell
    keep <- d <= radius
    local_matrix <- .landscape_cells_to_matrix(
      values = values[[1]][keep],
      cells = cells[keep],
      n_cols = terra::ncol(landcover)
    )

    expect_equal(
      .landscape_composition_by_buffer(d, values, radius, "shdi")[[1]],
      .landscape_composition_metric(values[[1]][keep], "shdi")
    )
    expect_equal(
      .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        resolution = terra::res(landcover)[[1]],
        n_cols = terra::ncol(landcover),
        metric = "lsi"
      )[[1]],
      .landscape_edge_metric(
        values = local_matrix,
        resolution = terra::res(landcover)[[1]],
        metric = "lsi"
      )
    )
    expect_equal(
      .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        n_cols = terra::ncol(landcover),
        metric = "contag"
      )[[1]],
      .landscape_adjacency_metric(local_matrix, "contag")
    )
    expect_equal(
      .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        n_cols = terra::ncol(landcover),
        metric = "ai"
      )[[1]],
      .landscape_adjacency_metric(local_matrix, "ai")
    )
  }
})


test_that("scale variable specs support multiple transforms from one source raster", {
  forest <- landscape_test_raster(808, 0:1, name = "forest")
  pts <- landscape_test_points(forest)
  radius <- 80

  vars <- msr_vars(
    forest_kernel = kernel_var("forest"),
    forest_ed = landscape_var("forest", metric = "ed", radius = radius),
    forest_shdi = landscape_var("forest", metric = "shdi", radius = radius)
  )

  expect_s3_class(vars, "multiScaleR_vars")
  expect_error(msr_vars(kernel_var("forest")), "must be named")

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = forest,
    max_D = radius,
    kernel = "gaussian",
    scale_vars = vars,
    verbose = FALSE
  )

  expect_named(kernel_inputs$kernel_dat,
               c("forest_kernel", "forest_ed", "forest_shdi"))
  expect_equal(kernel_inputs$n_covs, 1)
  expect_equal(.msr_optimized_covariates(kernel_inputs$scale_vars,
                                         cov_df = kernel_inputs$raw_cov),
               "forest_kernel")

  cached_ed <- landscape_cached_metric(
    raster = forest,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        resolution = terra::res(forest)[[1]],
        n_cols = terra::ncol(forest),
        metric = "ed"
      )[[1]]
    }
  )

  cached_shdi <- landscape_cached_metric(
    raster = forest,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_composition_by_buffer(
        d = d,
        r_stack.df = values,
        radius = radius,
        metric = "shdi"
      )[[1]]
    }
  )

  ed_unscaled <- kernel_inputs$kernel_dat$forest_ed *
    kernel_inputs$scl_params$sd[["forest_ed"]] +
    kernel_inputs$scl_params$mean[["forest_ed"]]
  shdi_unscaled <- kernel_inputs$kernel_dat$forest_shdi *
    kernel_inputs$scl_params$sd[["forest_shdi"]] +
    kernel_inputs$scl_params$mean[["forest_shdi"]]

  expect_equal(as.numeric(ed_unscaled), cached_ed)
  expect_equal(as.numeric(shdi_unscaled), cached_shdi)

  projected <- kernel_scale.raster(
    raster_stack = forest,
    sigma = radius,
    kernel = "gaussian",
    scale_vars = vars,
    verbose = FALSE
  )

  expect_named(projected, c("forest_kernel", "forest_ed", "forest_shdi"))
  expect_equal(terra::extract(projected[["forest_ed"]], pts)[, 2],
               cached_ed,
               tolerance = 25)
  expect_equal(terra::extract(projected[["forest_shdi"]], pts)[, 2],
               cached_shdi,
               tolerance = 0.01)
})


test_that("optimized landscape specs evaluate during model refits", {
  forest <- landscape_test_raster(909, 0:1, name = "forest")
  pts <- landscape_test_points(forest)
  vars <- msr_vars(
    forest_ed = landscape_var("forest", metric = "ed")
  )

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = forest,
    max_D = 80,
    scale_vars = vars,
    verbose = FALSE
  )
  dat <- data.frame(
    y = c(1, 0, 1, 0, 1),
    kernel_inputs$kernel_dat
  )
  mod <- glm(y ~ forest_ed, family = binomial(), data = dat)
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
  expect_named(final$scl_params$mean, "forest_ed")
})


test_that("kernel_scale.raster supports standalone fixed-radius landscape scale_vars without sigma", {
  forest <- landscape_test_raster(910, 0:1, name = "forest")
  pts <- landscape_test_points(forest)
  radius <- 80
  vars <- msr_vars(
    forest_ai = landscape_var("forest", metric = "ai", radius = radius),
    forest_ed = landscape_var("forest", metric = "ed", radius = radius)
  )

  cached_ai <- landscape_cached_metric(
    raster = forest,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_adjacency_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        n_cols = terra::ncol(forest),
        metric = "ai"
      )[[1]]
    }
  )

  cached_ed <- landscape_cached_metric(
    raster = forest,
    pts = pts,
    radius = radius,
    metric_fun = function(d, values, cells) {
      .landscape_edge_by_buffer(
        d = d,
        r_stack.df = values,
        cells = cells,
        radius = radius,
        resolution = terra::res(forest)[[1]],
        n_cols = terra::ncol(forest),
        metric = "ed"
      )[[1]]
    }
  )

  projected <- kernel_scale.raster(
    raster_stack = forest,
    scale_vars = vars,
    verbose = FALSE
  )

  expect_named(projected, c("forest_ai", "forest_ed"))
  expect_lt(
    max(abs(terra::extract(projected[["forest_ai"]], pts)[, 2] - cached_ai),
        na.rm = TRUE),
    2.5
  )
  expect_equal(terra::extract(projected[["forest_ed"]], pts)[, 2],
               cached_ed,
               tolerance = 25)
})


test_that("composition helpers validate inputs and handle empty buffers", {
  expect_equal(.landscape_shdi(c(1, 1, 2, 2)), log(2))
  expect_equal(.landscape_composition_metric(c(1, 1, 2, 2), "sidi"), 0.5)
  expect_true(is.na(.landscape_shdi(numeric(0))))

  expect_error(
    .landscape_composition_by_buffer(
      d = c(1, 2),
      r_stack.df = matrix(1:3, ncol = 1),
      radius = 1,
      metric = "shdi"
    ),
    "one row for each distance"
  )

  out <- .landscape_composition_by_buffer(
    d = c(10, 20),
    r_stack.df = data.frame(x = c(1, 2)),
    radius = 1,
    metric = "shdi"
  )
  expect_named(out, "x")
  expect_true(is.na(out[["x"]]))

  expect_error(
    .landscape_composition_by_buffer(
      d = c(1, 2),
      r_stack.df = data.frame(x = c(1, 1.5)),
      radius = 2,
      metric = "shdi"
    ),
    "integer-like class values"
  )
})


test_that("edge and adjacency helpers validate inputs and handle empty buffers", {
  local_matrix <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
  expect_equal(.landscape_edge_density(local_matrix, resolution = 10), 500)
  expect_equal(.landscape_edge_metric(local_matrix, 10, "te"), 20)
  expect_equal(.landscape_adjacency_metric(local_matrix, "ai"), 75)
  expect_equal(.landscape_adjacency_metric(local_matrix, "pladj"), 50)

  expect_error(
    .landscape_ed_by_buffer(
      d = c(1, 2),
      r_stack.df = data.frame(x = c(1, 0)),
      cells = 1,
      radius = 1,
      resolution = 10,
      n_cols = 2
    ),
    "same length as `d`"
  )

  out <- .landscape_ed_by_buffer(
    d = c(10, 20),
    r_stack.df = data.frame(x = c(1, 0)),
    cells = c(1, 2),
    radius = 1,
    resolution = 10,
    n_cols = 2
  )
  expect_named(out, "x")
  expect_true(is.na(out[["x"]]))
})


test_that("adjacency projections enforce class-count ceilings", {
  landcover <- terra::rast(
    nrows = 8,
    ncols = 8,
    xmin = 0,
    xmax = 80,
    ymin = 0,
    ymax = 80,
    crs = "EPSG:3857"
  )
  terra::values(landcover) <- seq_len(terra::ncell(landcover))
  names(landcover) <- "landcover"

  expect_error(
    .landscape_adjacency_raster_fft(landcover, radius = 20, metric = "contag"),
    "exceeds the current supported ceiling"
  )
})


test_that("kernel preparation rejects unscalable landscape covariates", {
  landcover <- terra::rast(
    nrows = 20,
    ncols = 20,
    xmin = 0,
    xmax = 200,
    ymin = 0,
    ymax = 200,
    crs = "EPSG:3857"
  )
  terra::values(landcover) <- 1
  names(landcover) <- "landcover"
  pts <- terra::vect(cbind(c(55, 85, 115), c(55, 85, 115)),
                     crs = terra::crs(landcover))
  vars <- msr_vars(
    landcover_shdi = landscape_var("landcover", metric = "shdi", radius = 30)
  )

  expect_error(
    kernel_prep(
      pts = pts,
      raster_stack = landcover,
      max_D = 30,
      scale_vars = vars,
      verbose = FALSE
    ),
    "zero variance"
  )
})


test_that("fragile landscape metrics give metric-specific scale-variation guidance", {
  landcover_two <- terra::rast(
    nrows = 20,
    ncols = 20,
    xmin = 0,
    xmax = 200,
    ymin = 0,
    ymax = 200,
    crs = "EPSG:3857"
  )
  terra::values(landcover_two) <- rep(c(1L, 2L), length.out = terra::ncell(landcover_two))
  names(landcover_two) <- "landcover"
  pts <- terra::vect(cbind(c(55, 85, 115), c(55, 85, 115)),
                     crs = terra::crs(landcover_two))

  expect_error(
    kernel_prep(
      pts = pts,
      raster_stack = landcover_two,
      max_D = 30,
      scale_vars = msr_vars(
        landcover_iji = landscape_var("landcover", metric = "iji", radius = 30)
      ),
      verbose = FALSE
    ),
    "Metric `iji` is fragile"
  )

  landcover_single <- terra::rast(landcover_two)
  terra::values(landcover_single) <- 1L
  names(landcover_single) <- "landcover"

  expect_error(
    kernel_prep(
      pts = pts,
      raster_stack = landcover_single,
      max_D = 30,
      scale_vars = msr_vars(
        landcover_siei = landscape_var("landcover", metric = "siei", radius = 30)
      ),
      verbose = FALSE
    ),
    "Metric `siei` can lose variation"
  )

  expect_error(
    kernel_prep(
      pts = pts,
      raster_stack = landcover_single,
      max_D = 30,
      scale_vars = msr_vars(
        landcover_msiei = landscape_var("landcover", metric = "msiei", radius = 30)
      ),
      verbose = FALSE
    ),
    "Metric `msiei` can lose variation"
  )
})


# Whole-landscape helpers: a buffer large enough to cover the entire raster makes
# the per-window co-occurrence equal to the whole-landscape value, so the direct
# C++ path can be compared to landscapemetrics at machine precision.
landscape_whole_adjacency <- function(raster, metric, base = 2, focal = NULL) {
  v <- terra::values(raster)[, 1]
  n <- length(v)
  .landscape_adjacency_by_buffer(
    d = rep(0, n),
    r_stack.df = data.frame(x = v),
    cells = seq_len(n),
    radius = 1e6,
    n_cols = terra::ncol(raster),
    metric = metric,
    base = base,
    focal_class = focal,
    validate = FALSE
  )[[1]]
}

landscape_whole_composition <- function(raster, metric, focal = NULL) {
  v <- terra::values(raster)[, 1]
  n <- length(v)
  .landscape_composition_by_buffer(
    d = rep(0, n),
    r_stack.df = data.frame(x = v),
    radius = 1e6,
    metric = metric,
    resolution = terra::res(raster)[[1]],
    focal_class = focal,
    validate = FALSE
  )[[1]]
}


test_that("IJI and information-theory metrics match landscapemetrics", {
  skip_if_not_installed("landscapemetrics")

  lsm_funs <- list(
    iji = landscapemetrics::lsm_l_iji,
    ent = landscapemetrics::lsm_l_ent,
    condent = landscapemetrics::lsm_l_condent,
    joinent = landscapemetrics::lsm_l_joinent,
    mutinf = landscapemetrics::lsm_l_mutinf,
    relmutinf = landscapemetrics::lsm_l_relmutinf
  )

  for (seed in c(11, 27)) {
    landcover <- landscape_test_raster(seed, 1:4)
    for (metric in names(lsm_funs)) {
      # base = 2 to match the log2 convention landscapemetrics uses
      got <- landscape_whole_adjacency(landcover, metric, base = 2)
      reference <- lsm_funs[[metric]](landcover)$value
      expect_equal(got, reference, tolerance = 1e-6,
                   info = paste0(metric, " (seed ", seed, ")"))
    }
  }
})


test_that("class-level CLUMPY, PLAND, and CA match landscapemetrics", {
  skip_if_not_installed("landscapemetrics")

  landcover <- landscape_test_raster(33, 1:4)
  lm_clumpy <- landscapemetrics::lsm_c_clumpy(landcover)
  lm_pland <- landscapemetrics::lsm_c_pland(landcover)
  lm_ca <- landscapemetrics::lsm_c_ca(landcover)

  for (k in 1:4) {
    expect_equal(landscape_whole_adjacency(landcover, "clumpy", focal = k),
                 lm_clumpy$value[lm_clumpy$class == k], tolerance = 1e-6,
                 info = paste("clumpy class", k))
    expect_equal(landscape_whole_composition(landcover, "pland", focal = k),
                 lm_pland$value[lm_pland$class == k], tolerance = 1e-6,
                 info = paste("pland class", k))
    expect_equal(landscape_whole_composition(landcover, "ca", focal = k),
                 lm_ca$value[lm_ca$class == k], tolerance = 1e-6,
                 info = paste("ca class", k))
  }
})


test_that("FFT projections of the new metrics agree with windowed landscapemetrics", {
  skip_if_not_installed("landscapemetrics")

  landcover <- landscape_test_raster(606, 1:4)
  pts <- landscape_test_points(landcover)
  radius <- 80
  n_cols <- terra::ncol(landcover)

  windowed_reference <- function(metric_call) {
    landscape_cached_metric(
      raster = landcover, pts = pts, radius = radius,
      metric_fun = function(d, values, cells) {
        keep <- d <= radius
        local_matrix <- .landscape_cells_to_matrix(values[[1]][keep], cells[keep], n_cols)
        metric_call(local_matrix)
      }
    )
  }

  # landscape-level (base = 2 for the information-theory metrics)
  level_funs <- list(
    iji = landscapemetrics::lsm_l_iji,
    ent = landscapemetrics::lsm_l_ent,
    mutinf = landscapemetrics::lsm_l_mutinf,
    relmutinf = landscapemetrics::lsm_l_relmutinf
  )
  for (metric in names(level_funs)) {
    projected <- terra::extract(
      .landscape_adjacency_raster_fft(landcover, radius = radius, base = 2, metric = metric),
      pts
    )[, 2]
    reference <- windowed_reference(function(m) {
      suppressWarnings(level_funs[[metric]](m)$value)
    })
    expect_equal(projected, reference, tolerance = 1,
                 info = paste("FFT", metric))
  }

  # class-level: clumpy (adjacency) and pland (composition) for one class
  projected_clumpy <- terra::extract(
    .landscape_adjacency_raster_fft(landcover, radius = radius, metric = "clumpy", focal_class = 2),
    pts
  )[, 2]
  reference_clumpy <- windowed_reference(function(m) {
    res <- suppressWarnings(landscapemetrics::lsm_c_clumpy(terra::rast(m)))
    out <- res$value[res$class == 2]
    if (length(out) == 0) NA_real_ else out
  })
  expect_equal(projected_clumpy, reference_clumpy, tolerance = 0.05,
               info = "FFT clumpy")

  projected_pland <- terra::extract(
    .landscape_composition_raster_fft(landcover, radius = radius, metric = "pland", focal_class = 2),
    pts
  )[, 2]
  reference_pland <- windowed_reference(function(m) {
    res <- suppressWarnings(landscapemetrics::lsm_c_pland(terra::rast(m)))
    out <- res$value[res$class == 2]
    if (length(out) == 0) 0 else out
  })
  expect_equal(projected_pland, reference_pland, tolerance = 1,
               info = "FFT pland")
})


test_that("new metrics handle invariant and degenerate landscapes", {
  single <- landscape_test_raster(5, 1L)          # one class everywhere
  two <- terra::rast(single)
  terra::values(two) <- rep(c(1L, 2L), terra::ncell(single) / 2)

  # IJI needs at least three classes
  expect_true(is.na(landscape_whole_adjacency(single, "iji")))
  expect_true(is.na(landscape_whole_adjacency(two, "iji")))

  # A single class carries no information: entropy 0, relative MI defined as 1
  expect_equal(landscape_whole_adjacency(single, "ent", base = 2), 0)
  expect_equal(landscape_whole_adjacency(single, "mutinf", base = 2), 0)
  expect_equal(landscape_whole_adjacency(single, "relmutinf"), 1)

  # CLUMPY is NA where the focal class fills the landscape or is absent
  expect_true(is.na(landscape_whole_adjacency(single, "clumpy", focal = 1)))
  expect_true(is.na(landscape_whole_adjacency(single, "clumpy", focal = 9)))

  # PLAND and CA: full cover for the present class, zero for an absent class
  expect_equal(landscape_whole_composition(single, "pland", focal = 1), 100)
  expect_equal(landscape_whole_composition(single, "pland", focal = 9), 0)
  expect_equal(landscape_whole_composition(single, "ca", focal = 9), 0)
})


test_that("landscape_var validates the class argument", {
  expect_error(landscape_var("x", "clumpy"), "`class` is required")
  expect_error(landscape_var("x", "pland"), "`class` is required")
  expect_error(landscape_var("x", "ca"), "`class` is required")
  expect_error(landscape_var("x", "iji", class = 2), "applies only to the class-level")
  expect_error(landscape_var("x", "shdi", class = 2), "applies only to the class-level")

  spec <- landscape_var("x", "clumpy", class = 3)
  expect_equal(spec$metric, "clumpy")
  expect_equal(spec$class, 3)

  vars <- msr_vars(
    cover3 = landscape_var("x", "pland", radius = 100, class = 3),
    juxta = landscape_var("x", "iji", radius = 100)
  )
  expect_true("class" %in% names(vars))
  expect_equal(vars$class, c(3, NA_real_))
})

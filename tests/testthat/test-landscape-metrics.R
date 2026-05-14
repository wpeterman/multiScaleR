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
    c("pladj", "contag"),
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
  names(cached) <- c("pladj", "contag")

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
    what = "lsm_l_pladj",
    verbose = FALSE,
    progress = FALSE
  )

  expect_equal(cached$pladj, landscape_reference_values(reference, "pladj"),
               tolerance = 5)
  expect_equal(cached$contag, contag_reference, tolerance = 0.01)
})


test_that("FFT adjacency projections agree with cached point metrics", {
  landcover <- landscape_test_raster(606, 1:4)
  pts <- landscape_test_points(landcover)
  radius <- 80

  for (metric in c("pladj", "contag")) {
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
    expect_equal(projected, cached, tolerance = if (metric == "pladj") 1.5 else 0.25)
  }
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
    forest_ed = landscape_var("forest", metric = "ed", radius = radius)
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

  expect_named(projected, "forest_ed")
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
})


test_that("edge and adjacency helpers validate inputs and handle empty buffers", {
  local_matrix <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
  expect_equal(.landscape_edge_density(local_matrix, resolution = 10), 500)
  expect_equal(.landscape_edge_metric(local_matrix, 10, "te"), 20)
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

# Exploratory spike for fixed-buffer landscape metrics in multiScaleR.
#
# This script is intentionally not package API. It checks whether a
# multiScaleR-style cached cell extraction can reproduce simple fixed-buffer
# landscape metrics and whether selected metrics can be projected with the
# existing FFT convolution machinery.

suppressPackageStartupMessages({
  library(exactextractr)
  library(fields)
  library(landscapemetrics)
  library(pkgload)
  library(sf)
  library(terra)
})

pkgload::load_all(quiet = TRUE)

if (!exists("fft_convolution", mode = "function")) {
  source("R/fft_smooth.R")
}
if (!exists("validate_scalar_numeric", mode = "function")) {
  source("R/validation_helpers.R")
}
if (!exists(".landscape_composition_by_buffer", mode = "function")) {
  source("R/landscape_metrics.R")
}

set.seed(101)

landcover <- terra::rast(
  nrows = 80,
  ncols = 80,
  xmin = 0,
  xmax = 800,
  ymin = 0,
  ymax = 800,
  crs = "EPSG:3857"
)
terra::values(landcover) <- sample(1:5, terra::ncell(landcover), replace = TRUE)
names(landcover) <- "landcover"

points <- terra::vect(
  cbind(c(165, 265, 425, 615, 705), c(165, 585, 425, 245, 675)),
  crs = terra::crs(landcover)
)
points_sf <- sf::st_as_sf(points)
point_xy <- sf::st_coordinates(points_sf)
radius <- 80

value_column <- function(x) {
  x[, setdiff(names(x), c("x", "y", "cell", "coverage_fraction"))[[1]],
    drop = FALSE]
}

cached_metric <- function(raster, radius, metric_fun) {
  extracted <- exactextractr::exact_extract(
    raster,
    sf::st_buffer(points_sf, dist = radius),
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

      metric_fun(d, value_column(extracted[[i]]), extracted[[i]]$cell)
    },
    numeric(1)
  )
}

cat("\n--- Composition/diversity metrics ---\n")
composition_metrics <- c("shdi", "shei", "sidi", "siei", "msidi", "msiei",
                         "pr", "prd", "rpr", "ta")
composition_check <- do.call(
  rbind,
  lapply(composition_metrics, function(metric) {
    cached <- cached_metric(
      landcover,
      radius,
      function(d, values, cells) {
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
    projected <- terra::extract(
      .landscape_composition_raster_fft(
        landcover,
        radius = radius,
        metric = metric,
        classes_max = 5
      ),
      points
    )[, 2]

    data.frame(
      metric = metric,
      max_abs_fft_diff = max(abs(projected - cached), na.rm = TRUE)
    )
  })
)
print(composition_check)

cat("\n--- Edge metrics ---\n")
binary_landcover <- landcover == 1
names(binary_landcover) <- "binary_landcover"
edge_metrics <- c("ed", "te", "lsi")
edge_check <- do.call(
  rbind,
  lapply(edge_metrics, function(metric) {
    cached <- cached_metric(
      binary_landcover,
      radius,
      function(d, values, cells) {
        .landscape_edge_by_buffer(
          d = d,
          r_stack.df = values,
          cells = cells,
          radius = radius,
          resolution = terra::res(binary_landcover)[[1]],
          n_cols = terra::ncol(binary_landcover),
          metric = metric
        )[[1]]
      }
    )
    projected <- if (metric %in% c("ed", "te")) {
      terra::extract(
        .landscape_edge_raster_fft(binary_landcover, radius, metric = metric),
        points
      )[, 2]
    } else if (metric == "lsi") {
      terra::extract(
        .landscape_edge_raster_fft(binary_landcover, radius, metric = metric),
        points
      )[, 2]
    } else {
      rep(NA_real_, length(cached))
    }

    data.frame(
      metric = metric,
      cached_range = paste(round(range(cached, na.rm = TRUE), 3),
                           collapse = " - "),
      max_abs_fft_diff = if (all(is.na(projected))) {
        NA_real_
      } else {
        max(abs(projected - cached), na.rm = TRUE)
      }
    )
  })
)
print(edge_check)

cat("\n--- Adjacency metrics ---\n")
adjacency_metrics <- c("pladj", "contag")
adjacency_check <- do.call(
  rbind,
  lapply(adjacency_metrics, function(metric) {
    cached <- cached_metric(
      landcover,
      radius,
      function(d, values, cells) {
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
    projected <- terra::extract(
      .landscape_adjacency_raster_fft(landcover, radius, metric = metric),
      points
    )[, 2]

    data.frame(
      metric = metric,
      cached_range = paste(round(range(cached, na.rm = TRUE), 3),
                           collapse = " - "),
      max_abs_fft_diff = max(abs(projected - cached), na.rm = TRUE)
    )
  })
)
print(adjacency_check)

cat("\n--- Takeaways ---\n")
cat("* Composition/diversity metrics reuse class-count kernels and project cleanly by FFT.\n")
cat("* ED, TE, and LSI share the same local edge-count primitives for cached and FFT workflows.\n")
cat("* PLADJ and CONTAG are adjacency-table metrics; FFT projections are feasible but cost grows with class count.\n")

spatial_source_fixture <- function() {
  raster <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 80, ymin = 0, ymax = 80,
    crs = "EPSG:26917")
  terra::values(raster) <- seq_len(64L)
  names(raster) <- "habitat"
  points <- sf::st_as_sf(data.frame(id = c("a", "b"), east = c(5, 40), north = c(5, 40)),
    coords = c("east", "north"), crs = 26917)
  rownames(points) <- points$id
  inputs <- kernel_prep(points, raster, max_D = 20, sigma = 8, nbins = 16,
    store_cell_data = FALSE, verbose = FALSE)
  data <- data.frame(key = c("b", "a", "b"), east = c(40, 5, 40), north = c(40, 5, 40))
  list(points = points, raster = raster, inputs = inputs, data = data)
}

test_that("source reports align explicit projected sources without changing model inputs", {
  f <- spatial_source_fixture()
  report <- multiScaleR:::.msr_spatial_source_report
  before <- serialize(f$inputs, NULL)
  set.seed(713404); rng <- .Random.seed
  x <- report(f$points, f$raster, f$inputs, f$data, c("east", "north"), "key", "id")
  expect_s3_class(x, "msr_spatial_source_report")
  expect_identical(.Random.seed, rng)
  expect_identical(serialize(f$inputs, NULL), before)
  expect_identical(x$alignment$observations, 3L)
  expect_identical(x$alignment$distinct_observed_points, 2L)
  expect_equal(x$alignment$maximum_coordinate_difference, 0)
  expect_equal(x$alignment$coordinate_tolerance, 0)
  expect_identical(x$source$linear_unit, sf::st_crs(f$points)$units_gdal)
  expect_identical(x$source$raster_layers, "habitat")
  expect_identical(x$extent$points$point_id, c("a", "b"))
  expect_equal(x$extent$points$minimum_clearance, c(5, 40))
  expect_identical(x$extent$points$point_inside_extent, c(TRUE, TRUE))
  expect_identical(x$extent$points$buffer_within_extent, c(FALSE, TRUE))
  expect_equal(x$distance_scaling$value, f$inputs$unit_conv)
  expect_match(x$distance_scaling$meaning, "not a conversion")
  expect_match(x$source$historical_raster_lineage, "unverified")
  expect_match(x$source$historical_distance_unit_agreement, "unverified")
  expect_match(x$extent$unit_assumption, "unverified")
  expect_identical(x$extent$valid_cell_coverage, "not_evaluated")
  expect_identical(x$extent$kernel_weight_coverage, "not_evaluated")
  expect_false(x$fit_eligibility)
  expect_identical(x$calibration, "not_established")
  expect_identical(report(f$points[2:1, ], f$raster, f$inputs, f$data,
    c("east", "north"), "key", "id"), x)
  expect_identical(report(terra::vect(f$points), f$raster, f$inputs, f$data,
    c("east", "north"), "key", "id"), x)
  d <- f$data[2:1, ]; rownames(d) <- d$key
  expect_equal(report(f$points, f$raster, f$inputs, d, c("east", "north"))$alignment$observations, 2)
  # A raster hole can leave this metadata-only report identical. Extent interior
  # therefore never certifies actual valid-cell or kernel-weight coverage.
  hole <- f$raster; hole[29] <- NA_real_
  expect_identical(report(f$points, hole, f$inputs, f$data, c("east", "north"), "key", "id"), x)
  near <- f$data; near$east[1] <- near$east[1] + .01
  expect_error(report(f$points, f$raster, f$inputs, near, c("east", "north"), "key", "id"), "disagree")
  tolerant <- report(f$points, f$raster, f$inputs, near, c("east", "north"), "key", "id", .02)
  expect_equal(tolerant$alignment$coordinate_tolerance, .02)
  expect_equal(tolerant$alignment$maximum_coordinate_difference, .01)
  spatial_data <- sf::st_as_sf(f$data, coords = c("east", "north"), crs = 26917, remove = FALSE)
  expect_identical(report(f$points, f$raster, f$inputs, spatial_data,
    c("east", "north"), "key", "id"), x)
  outside <- f$points; sf::st_geometry(outside)[[1L]] <- sf::st_point(c(-5, 5))
  outside_data <- f$data; outside_data$east[outside_data$key == "a"] <- -5
  edge <- report(outside, f$raster, f$inputs, outside_data, c("east", "north"), "key", "id")
  expect_equal(edge$extent$points$minimum_clearance[1], -5)
  expect_false(edge$extent$points$point_inside_extent[1])
  expect_false(edge$fit_eligibility)
})

test_that("source reports reject unverified coordinate and row interpretations", {
  f <- spatial_source_fixture()
  call_report <- function(points = f$points, raster = f$raster, inputs = f$inputs,
      data = f$data, coordinates = c("east", "north"), row_id = "key", point_id = "id",
      coordinate_tolerance = 0) multiScaleR:::.msr_spatial_source_report(points, raster,
        inputs, data, coordinates, row_id, point_id, coordinate_tolerance)
  expect_error(call_report(points = matrix(1)), "points")
  expect_error(call_report(points = sf::st_set_crs(f$points, NA)), "defined CRS")
  geographic <- sf::st_transform(f$points, 4326)
  expect_error(call_report(points = geographic), "projected")
  expect_error(call_report(raster = matrix(1)), "SpatRaster")
  wrong_crs <- f$raster; terra::crs(wrong_crs) <- "EPSG:26918"
  expect_error(call_report(raster = wrong_crs), "same defined projected CRS")
  no_crs <- f$raster; terra::crs(no_crs) <- ""
  expect_error(call_report(raster = no_crs), "same defined projected CRS")
  polygon <- sf::st_buffer(f$points, 1)
  expect_error(call_report(points = polygon), "POINT")
  empty <- f$points[FALSE, ]
  expect_error(call_report(points = empty), "nonempty")
  empty <- f$points; sf::st_geometry(empty)[[1L]] <- sf::st_point()
  expect_error(call_report(points = empty), "nonempty")
  three_d <- sf::st_as_sf(data.frame(id = c("a", "b"), east = c(5, 40), north = c(5, 40), z = 1),
    coords = c("east", "north", "z"), crs = 26917)
  expect_error(call_report(points = three_d), "two-dimensional")
  duplicate <- f$points; duplicate$id[] <- "a"
  expect_error(call_report(points = duplicate), "uniquely match")
  missing <- f$points; missing$id[1] <- NA_character_
  expect_error(call_report(points = missing), "nonmissing")
  expect_error(call_report(point_id = "absent"), "point_id")
  expect_error(call_report(row_id = "absent"), "row_id")
  bad <- f$data; bad$key[1] <- "unknown"
  expect_error(call_report(data = bad), "Observation IDs")
  bad <- f$data; bad$key[1] <- NA_character_
  expect_error(call_report(data = bad), "nonmissing")
  bad <- f$data; bad$east[1] <- Inf
  expect_error(call_report(data = bad), "finite numeric")
  expect_error(call_report(coordinates = c("east", "east")), "coordinates")
  bad <- f$data; bad$east <- bad$east / 1000
  expect_error(call_report(data = bad), "disagree")
  for (tol in list(-1, NA_real_, Inf, c(0, 1), "0", matrix(0), 1i))
    expect_error(call_report(coordinate_tolerance = tol), "coordinate_tolerance")
  bad <- f$inputs; bad$binned <- NULL
  expect_error(call_report(inputs = bad), "binned")
  bad <- f$inputs; bad$binned$point_ids[] <- "a"
  expect_error(call_report(inputs = bad), "Binned point IDs")
  renamed <- f$raster; names(renamed) <- "other"
  expect_error(call_report(raster = renamed), "source layers")
  bad <- f$inputs; bad$binned$vsum$habitat <- bad$binned$vsum$habitat[1, , drop = FALSE]
  expect_error(call_report(inputs = bad), "Binned point IDs")
  bad <- sf::st_as_sf(f$data, coords = c("east", "north"), crs = 4326, remove = FALSE)
  expect_error(call_report(data = bad), "Spatial observation data")
})

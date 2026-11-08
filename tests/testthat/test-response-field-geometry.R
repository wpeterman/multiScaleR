test_that("field geometry binds explicit projected rows and rejects changed contexts", {
  skip_if_not_installed("glmmTMB")
  f <- native_response_field_fixture("poisson")
  xy <- terra::crds(f$points)
  crs <- terra::crs(f$points)
  geometry <- response_field_geometry(f$object, xy, crs)
  expect_s3_class(geometry, "multiScaleR_field_geometry")
  expect_identical(geometry$position, "field_pos")
  expect_identical(geometry$group, "field_group")
  expect_identical(geometry$group_value, "shared")
  expect_identical(geometry$training_locations, 36L)
  expect_identical(geometry$training_rows, 108L)
  expect_identical(geometry$provider_version, as.character(getNamespaceVersion("glmmTMB")))
  expect_true(nzchar(geometry$units))
  expect_false(any(vapply(geometry, is.matrix, logical(1))))
  expect_equal(nchar(geometry$reference_sha256), 64L)
  expect_error(response_field_geometry(f$object, xy[-1L, ], crs), "every retained")
  expect_error(response_field_geometry(f$object, xy[, 2:1], crs), "row order and units")
  expect_error(response_field_geometry(f$object, xy[c(4:108, 1:3), ], crs), "row order and units")
  expect_error(response_field_geometry(f$object, xy / 1000, crs), "row order and units")
  bad <- xy; bad[1L, 1L] <- Inf
  expect_error(response_field_geometry(f$object, bad, crs), "finite coordinate pairs")
  expect_error(response_field_geometry(f$object, xy, 4326), "known projected CRS")
  expect_error(response_field_geometry(f$object, xy, NA), "known projected CRS")
  for (definition in c(4978, 5703, 7405))
    expect_error(response_field_geometry(f$object, xy, definition), "known projected CRS")
  bound <- sf::st_crs("+proj=utm +zone=17 +ellps=clrk66 +towgs84=-8,160,176 +units=m +no_defs")
  expect_s3_class(.msr_response_field_crs(bound), "crs")
  expect_error(response_field_geometry(f$object, xy, crs, batch_size = 0), "batch_size")
  expect_error(response_field_geometry(f$object, xy, crs, batch_size = 1.5), "batch_size")
  expect_error(response_field_geometry(f$object, xy, crs, max_memory_mb = Inf), "finite")
  expect_error(response_field_geometry(refit_spatial_fixture()$object, xy, crs), "native glmmTMB")
  path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
  saveRDS(geometry, path)
  expect_identical(readRDS(path), geometry)
})

test_that("response-field cells retain native values and SEs across serial batches", {
  skip_if_not_installed("glmmTMB")
  for (family in c("poisson", "gaussian")) {
    f <- native_response_field_fixture(family)
    model <- f$model
    xy <- terra::crds(f$points)
    geometry <- response_field_geometry(f$object, xy, terra::crs(f$points), batch_size = 3L)
    whole <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
      scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
      site_values = list(z = .2, field_pos = as.character(f$data$field_pos[1L]),
        field_group = "shared"), offset_values = c(exposure = 2))
    # Crop already computed cells; their values, transformation and descriptor
    # remain exact. This is test preparation, not reconstruction from a file.
    raster <- terra::crop(whole, terra::ext(40, 80, 40, 80))
    attr(raster, "multiScaleR_projection") <- attr(whole, "multiScaleR_projection")
    cells <- seq_len(terra::ncell(raster))
    coords <- terra::xyFromCell(raster, cells)
    data <- terra::as.data.frame(raster)
    data$field_pos <- glmmTMB::numFactor(coords)
    data$field_group <- factor("shared", levels = levels(f$data$field_group))
    for (zero in c(FALSE, TRUE)) {
      controls <- list(allow.new.levels = TRUE)
      if (zero) controls$re.form <- NA
      expected <- do.call(predict, c(list(object = model, newdata = data,
        type = "response", se.fit = TRUE), controls))
      output <- predict(f$object, raster, type = "response", se.fit = TRUE,
        uncertainty = "conditional", native_args = controls,
        raster_args = list(field_geometry = geometry))
      expect_equal(as.numeric(terra::values(output)[, 1L]), as.numeric(expected$fit), tolerance = 1e-7)
      expect_equal(as.numeric(terra::values(output)[, 2L]), as.numeric(expected$se.fit), tolerance = 1e-7)
      contract <- attr(output, "prediction_contract")$response_spatial_field
      expect_identical(contract$geometry$coordinate_source, "projected_raster_cell_centers")
      expect_identical(contract$geometry$units, geometry$units)
      expect_false(contract$kernel_uncertainty_included)
      expect_identical(contract$conditioning, if (zero)
        "all_response_random_effects_set_to_zero" else "native_conditional_spatial_prediction")
    }
    # Repeated/permuted/new cell rows are checked independently of terra's order.
    thin <- list(object = list(opt_mod = model), type = "link", se.fit = TRUE,
      native_args = list(allow.new.levels = TRUE), logical_names = character())
    rows <- c(1, 7, 1, 16)
    data <- data[rows, , drop = FALSE]
    data$.msr_field_x <- coords[rows, 1L]; data$.msr_field_y <- coords[rows, 2L]
    for (batch in c(1L, 4L)) {
      thin$field_geometry <- response_field_geometry(f$object, xy,
        terra::crs(f$points), batch_size = batch)
      actual <- .msr_native_raster_eval(thin, data)
      expected <- predict(model, data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
      expect_equal(actual$fit, as.numeric(expected$fit), tolerance = 1e-7)
      expect_equal(actual$se.fit, as.numeric(expected$se.fit), tolerance = 1e-7)
    }
    known <- f$data[c(1, 4, 1, 40), , drop = FALSE]
    known$.msr_field_x <- xy[c(1, 4, 1, 40), 1L]
    known$.msr_field_y <- xy[c(1, 4, 1, 40), 2L]
    known_thin <- thin; known_thin$native_args <- list()
    known_expected <- predict(model, known, type = "link", se.fit = TRUE)
    known_actual <- .msr_native_raster_eval(known_thin, known)
    expect_equal(known_actual$fit, as.numeric(known_expected$fit), tolerance = 1e-7)
    expect_equal(known_actual$se.fit, as.numeric(known_expected$se.fit), tolerance = 1e-7)
    if (family == "poisson") {
      changed <- geometry; changed$batch_size <- 1L
      expect_error(predict(f$object, raster, raster_args = list(field_geometry = changed)), "Recreate")
      changed_fit <- f$object; changed_fit$opt_mod$fit$par[1L] <- changed_fit$opt_mod$fit$par[1L] + .01
      expect_error(predict(changed_fit, raster, raster_args = list(field_geometry = geometry)), "fitted context")
      wrong <- raster; terra::crs(wrong) <- "EPSG:32617"
      expect_error(predict(f$object, wrong, raster_args = list(field_geometry = geometry)), "CRS and linear units")
      # Revalidate a saved descriptor admitted by the initial constructor.
      legacy <- geometry; legacy$crs <- sf::st_crs(4978)$wkt
      legacy$definition_sha256 <- .msr_gdr_hash(unclass(legacy)[
        setdiff(names(legacy), "definition_sha256")])
      wrong <- raster; terra::crs(wrong) <- legacy$crs
      expect_error(.msr_response_field_preflight(legacy, f$object, wrong, 1L, NULL),
        "known projected CRS")
      expect_error(predict(f$object, raster, raster_args = list(field_geometry = geometry, cores = 2)), "cores = 1")
      tiny <- response_field_geometry(f$object, xy, terra::crs(f$points), max_memory_mb = .001)
      path <- tempfile(fileext = ".tif"); on.exit(unlink(path), add = TRUE)
      expect_error(predict(f$object, raster, raster_args = list(field_geometry = tiny, filename = path)), "exceeds max_memory_mb")
      expect_false(file.exists(path))
      expect_error(predict(f$object, raster, raster_args = list(field_geometry = geometry,
        const = data.frame(.msr_field_x = 0))), "Reserved field coordinate names")
      expect_error(predict(f$object, raster, type = "zprob", raster_args = list(field_geometry = geometry)), "response, link or conditional")
      bad <- data; bad$field_group <- "different"
      expect_error(.msr_native_raster_eval(thin, bad), "verified shared field group")
      bad <- data; bad$.msr_field_x[1L] <- Inf
      expect_error(.msr_native_raster_eval(thin, bad), "finite projected field coordinates")
      expect_equal(dim(.msr_native_raster_eval(thin, data[FALSE, , drop = FALSE])), c(0L, 2L))
      expect_error(predict(f$object, raster, raster_args = list(field_geometry = geometry)),
        "New response-field positions require")
      # Verify warning forwarding separately from the provider's unsupported
      # default new-position contrast path, using a native-boundary warning.
      (function() {
        warn_thin <- thin
        local_mocked_bindings(.msr_native_prediction_rows = function(object, newdata, ...) {
          warning("native field warning", call. = FALSE)
          data.frame(fit = rep(1, nrow(newdata)), se.fit = rep(.1, nrow(newdata)))
        }, .package = "multiScaleR", .env = environment())
        expect_warning(.msr_native_raster_eval(warn_thin, data), "native field warning")
      })()
      file <- tempfile(fileext = ".tif"); on.exit(unlink(file), add = TRUE)
      masked <- raster; masked[1L] <- NA
      attr(masked, "multiScaleR_projection") <- attr(raster, "multiScaleR_projection")
      disk <- predict(f$object, masked, native_args = list(allow.new.levels = TRUE),
        raster_args = list(field_geometry = geometry, filename = file, wopt = list(steps = 2)))
      expect_true(file.exists(file))
      expect_true(is.na(terra::values(disk)[1L, 1L]))
      point_rows <- terra::as.data.frame(raster)
      point_rows$field_pos <- glmmTMB::numFactor(coords)
      point_rows$field_group <- factor("shared", levels = levels(f$data$field_group))
      reference <- predict(model, point_rows, allow.new.levels = TRUE, type = "response")
      expect_equal(as.numeric(terra::values(disk)[-1L, 1L]),
        as.numeric(reference[-1L]), tolerance = 1e-7)
      native_se <- predict(model, point_rows, allow.new.levels = TRUE,
        type = "response", se.fit = TRUE)
      # An empty probe row must not collapse the requested fit/SE schema.
      holes <- raster; holes[5:8] <- NA
      attr(holes, "multiScaleR_projection") <- attr(raster, "multiScaleR_projection")
      hole_map <- predict(f$object, holes, se.fit = TRUE, uncertainty = "conditional",
        native_args = list(allow.new.levels = TRUE), raster_args = list(field_geometry = geometry))
      expect_equal(as.numeric(terra::values(hole_map)[-(5:8), 1L]),
        as.numeric(native_se$fit[-(5:8)]), tolerance = 1e-7)
      expect_equal(as.numeric(terra::values(hole_map)[-(5:8), 2L]),
        as.numeric(native_se$se.fit[-(5:8)]), tolerance = 1e-7)
      blank <- raster; blank[] <- NA
      attr(blank, "multiScaleR_projection") <- attr(raster, "multiScaleR_projection")
      empty_map <- predict(f$object, blank, se.fit = TRUE, uncertainty = "conditional",
        native_args = list(allow.new.levels = TRUE), raster_args = list(field_geometry = geometry))
      expect_equal(terra::nlyr(empty_map), 2L)
      expect_identical(names(empty_map), c("fit", "se.fit"))
      expect_true(all(is.na(terra::values(empty_map))))
      scenarios <- predict(f$object, raster, se.fit = TRUE, uncertainty = "conditional",
        native_args = list(allow.new.levels = TRUE),
        raster_args = list(field_geometry = geometry, const = data.frame(scenario_id = 1:2)))
      expect_equal(unname(terra::values(scenarios)),
        cbind(as.numeric(native_se$fit), as.numeric(native_se$se.fit),
          as.numeric(native_se$fit), as.numeric(native_se$se.fit)), tolerance = 1e-7)
      expect_identical(names(scenarios), c("fit", "se.fit", "fit.1", "se.fit.1"))
    }
  }
  cat("RESPONSE_FIELD_GEOMETRY families=2 no_year=1 projected=1 serial_batches=1 joint_guard=retained\n")
})

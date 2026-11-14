# Response-field maps and new-position predictions from a real glmmTMB field
# fit made by multiScale_optim(). The default refit backend stores no
# `backend` element, so fixtures that set one by hand hid a rejection of
# real fits by response_field_geometry().

# One small field fit, shared by the tests in this file.
field_fit_env <- new.env()
small_field_fit <- function() {
  if (!is.null(field_fit_env$fit)) return(field_fit_env$fit)
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  terra::crs(land) <- "EPSG:32617"
  land <- land[["land2"]]
  pts <- terra::vect(surv_pts)
  terra::crs(pts) <- "EPSG:32617"
  set.seed(3)
  pts <- pts[sort(sample(nrow(pts), 50))]
  xy <- terra::crds(pts)
  truth <- kernel_scale.raster(land, sigma = 400, kernel = "gaussian", verbose = FALSE)
  x2 <- as.numeric(scale(terra::extract(truth, pts, ID = FALSE)[, 1]))
  field <- drop(t(chol(0.5^2 * exp(-as.matrix(stats::dist(xy)) / 800))) %*%
                  stats::rnorm(nrow(xy)))
  inputs <- kernel_prep(pts = pts, raster_stack = land, max_D = 1200,
                        kernel = "gaussian", verbose = FALSE)
  df <- data.frame(counts = stats::rpois(nrow(xy), exp(1 + 0.6 * x2 + field)),
                   inputs$kernel_dat)
  df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])
  df$grp <- factor("all")
  # Explicit starts on meter coordinates: log field SD and log range.
  m0 <- suppressWarnings(glmmTMB::glmmTMB(
    counts ~ land2 + exp(pos + 0 | grp), family = poisson(), data = df,
    start = list(theta = c(log(0.5), log(500)))))
  fit <- suppressWarnings(multiScale_optim(m0, inputs, n_cores = 1,
                                           verbose = FALSE, par = 400 / 1200))
  field_fit_env$fit <- list(fit = fit, land = land, xy = xy, df = df)
  field_fit_env$fit
}

test_that("response_field_geometry() accepts a real default-backend field fit", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  f <- small_field_fit()
  expect_null(f$fit$backend)
  expect_true(f$fit$opt_mod$sdr$pdHess)

  geometry <- response_field_geometry(f$fit, f$xy, 32617, batch_size = 64L)
  expect_s3_class(geometry, "multiScaleR_field_geometry")
  expect_identical(geometry$training_locations, 50L)

  # Map a small window with the field. Cells are new field positions.
  site <- list(pos = as.character(f$df$pos[1]), grp = "all")
  whole <- kernel_scale.raster(f$land, multiScaleR = f$fit, scale_center = TRUE,
                               match_buffer = TRUE, clamp = TRUE, verbose = FALSE,
                               site_values = site)
  layers <- terra::crop(whole, terra::ext(1500, 1700, 1500, 1700))
  attr(layers, "multiScaleR_projection") <- attr(whole, "multiScaleR_projection")
  map <- predict(f$fit, layers, native_args = list(allow.new.levels = TRUE),
                 raster_args = list(field_geometry = geometry))
  values <- terra::values(map)[, 1]
  expect_true(all(is.finite(values)))

  # Each cell equals glmmTMB's own prediction at that cell's coordinates.
  cell_xy <- terra::xyFromCell(layers, seq_len(terra::ncell(layers)))
  native <- predict(f$fit$opt_mod, type = "response", allow.new.levels = TRUE,
                    newdata = data.frame(land2 = terra::values(layers[["land2"]])[, 1],
                                         pos = glmmTMB::numFactor(cell_xy[, 1], cell_xy[, 2]),
                                         grp = factor("all", levels = "all")))
  expect_equal(values, as.numeric(native), tolerance = 1e-8)
})

test_that("response_field_geometry() still rejects analytic-backend results", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  f <- small_field_fit()
  other <- f$fit
  other$backend <- "profile_glm"
  expect_error(response_field_geometry(other, f$xy, 32617),
               "native glmmTMB multiScaleR refit result")
})

test_that("predicting a field at new positions names allow.new.levels and re.form", {
  skip_if_not_installed("glmmTMB")
  f <- small_field_fit()
  new_xy <- cbind(c(1234, 2345), c(1111, 2222))
  new_inputs <- kernel_prep(pts = terra::vect(new_xy, crs = "EPSG:32617"),
                            raster_stack = f$land, max_D = 1200,
                            kernel = "gaussian", verbose = FALSE)
  nd <- data.frame(pos = glmmTMB::numFactor(new_xy[, 1], new_xy[, 2]),
                   grp = factor("all"),
                   row.names = rownames(new_inputs$kernel_dat))
  expect_error(predict(f$fit, nd, new_inputs),
               "allow.new.levels = TRUE.*re.form = NA")

  # Both named routes work.
  with_field <- predict(f$fit, nd, new_inputs, allow.new.levels = TRUE)
  zero_field <- predict(f$fit, nd, new_inputs, re.form = NA)
  expect_true(all(is.finite(as.numeric(with_field))))
  expect_true(all(is.finite(as.numeric(zero_field))))

  # Fitted rows need neither control.
  expect_true(all(is.finite(as.numeric(predict(f$fit)))))
})

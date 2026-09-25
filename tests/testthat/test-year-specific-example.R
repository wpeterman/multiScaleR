test_that("group-specific bins and final model use the assigned annual map", {
  example_file <- system.file("examples", "year_specific_rasters.R",
                              package = "multiScaleR")
  if (!nzchar(example_file)) {
    example_file <- test_path("..", "..", "inst", "examples",
                              "year_specific_rasters.R")
  }
  env <- new.env(parent = globalenv())
  invisible(capture.output(sys.source(example_file, envir = env)))

  # Prepare each map at all points only as an independent reference.
  ref_2024 <- kernel_prep(env$pts, env$habitat_2024, max_D = 400,
                          sigma = 100, bin = FALSE, verbose = FALSE)
  ref_2025 <- kernel_prep(env$pts, env$habitat_2025, max_D = 400,
                          sigma = 100, bin = FALSE, verbose = FALSE)
  manual_covariate <- function(sigma_m) {
    vapply(seq_len(nrow(env$obs)), function(i) {
      annual <- if (env$obs$year[i] == "2025") ref_2025 else ref_2024
      d <- annual$d_list[[i]] * annual$unit_conv
      val <- as.matrix(annual$raw_cov[[i]])[, "habitat"]
      weighted.mean(val, exp(-0.5 * (d / sigma_m)^2))
    }, numeric(1))
  }
  for (sigma_m in c(60, 100, 180)) {
    actual <- .msr_kernel_cov_w(
      par = sigma_m / env$pooled$unit_conv,
      d_list = NULL, cov_df = NULL,
      kernel = "gaussian", opt_context = env$fit$opt_context
    )
    expect_lt(max(abs(actual[, "habitat"] - manual_covariate(sigma_m))),
              0.001)
  }
  expect_equal(
    as.numeric(env$pooled$kernel_dat$habitat),
    as.numeric(scale(manual_covariate(100))), tolerance = 1e-10
  )
  direct_data <- env$model_data
  binned_at_fit <- .msr_kernel_cov_w(
    par = env$fit$scale_est["habitat", "Mean"] / env$pooled$unit_conv,
    d_list = NULL, cov_df = NULL,
    kernel = "gaussian", opt_context = env$fit$opt_context
  )
  direct_data$habitat <- as.numeric(scale(binned_at_fit[, "habitat"]))
  direct_fit <- glm(survived ~ habitat + year, data = direct_data,
                    family = binomial())
  expect_equal(nobs(env$fit$opt_mod), nrow(env$obs))
  expect_equal(coef(env$fit$opt_mod), coef(direct_fit), tolerance = 1e-8)
  expect_null(env$pooled$raw_cov)
  expect_identical(unname(env$pooled$raster_group),
                   as.character(env$obs$year))

  project_year <- function(map, year) {
    surface <- suppressWarnings(kernel_scale.raster(
      map, multiScaleR = env$fit, scale_center = TRUE,
      verbose = FALSE
    ))
    terra::predict(surface, env$fit$opt_mod, type = "response",
                   fun = function(model, data, ...) {
                     data$year <- factor(year, levels = levels(env$obs$year))
                     predict(model, newdata = data, ...)
                   })
  }
  pred_2024 <- project_year(env$habitat_2024, "2024")
  pred_2025 <- project_year(env$habitat_2025, "2025")
  expect_true(any(is.finite(terra::values(pred_2024))))
  expect_true(any(is.finite(terra::values(pred_2025))))
  expect_false(isTRUE(all.equal(terra::values(pred_2024),
                                terra::values(pred_2025))))
})

test_that("group preparation pools constant annual maps and tracks changing NA", {
  r1 <- terra::rast(nrows = 18, ncols = 21, xmin = 0, xmax = 210,
                    ymin = 0, ymax = 180, crs = "EPSG:26915")
  terra::values(r1) <- 0
  names(r1) <- "habitat"
  r2 <- r1
  terra::values(r2) <- 1
  r2[terra::cellFromXY(r2, matrix(c(50, 50), ncol = 2))] <- NA_real_
  xy <- cbind(c(50, 80, 100, 50, 80, 100),
              c(50, 80, 100, 50, 80, 100))
  pts <- terra::vect(xy, type = "points", crs = terra::crs(r1))
  years <- rep(c("a", "b"), 3)
  out <- kernel_prep_by_group(pts, list(a = r1, b = r2), years,
                              max_D = 25, store_cell_data = FALSE,
                              verbose = FALSE)
  expect_identical(unname(out$raster_group), years)
  expect_null(out$raw_cov)
  expect_equal(unname(as.numeric(out$scl_params$mean)), 0.5)
  expect_equal(unique(out$kernel_dat$habitat[c(1, 3, 5)]),
               -unique(out$kernel_dat$habitat[c(2, 4, 6)]))
  expect_true(any(out$binned$csum$habitat[2, ] > 0))
  expect_lt(sum(out$binned$csum$habitat[4, ]),
            sum(out$binned$csum$habitat[1, ]))
  specs <- msr_vars(mean_habitat = kernel_var("habitat"),
                    cover = landscape_var("habitat", metric = "pland",
                                          class = 1, radius = 25))
  mixed <- kernel_prep_by_group(pts, list(a = r1, b = r2), years,
                                max_D = 25, scale_vars = specs,
                                store_cell_data = FALSE, verbose = FALSE)
  expect_true(mixed$cell_data_stored)
  expect_identical(mixed$binned$covariates, "mean_habitat")
  expect_true(all(is.finite(as.matrix(mixed$kernel_dat))))
  expect_error(
    kernel_prep_by_group(pts, list(a = r1, b = r2),
                         c(years[-1], "missing"), max_D = 25),
    "assign every point"
  )
  misaligned <- terra::shift(r2, dx = 10)
  expect_error(
    kernel_prep_by_group(pts, list(a = r1, b = misaligned), years,
                         max_D = 25),
    "identical grid"
  )
})

test_that("one unchanged map reproduces ordinary kernel preparation", {
  r <- terra::rast(nrows = 20, ncols = 23, xmin = 0, xmax = 230,
                   ymin = 0, ymax = 200, crs = "EPSG:26915")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- sin(xy[, 1] / 32) + cos(xy[, 2] / 46)
  names(r) <- "habitat"
  pts <- terra::vect(cbind(c(40, 65, 90, 115, 140, 165),
                           c(45, 75, 105, 135, 85, 55)),
                     type = "points", crs = terra::crs(r))
  ordinary <- kernel_prep(pts, r, max_D = 30, sigma = 12,
                          verbose = FALSE)
  grouped <- kernel_prep_by_group(pts, list(a = r, b = r),
                                  rep(c("a", "b"), 3), max_D = 30,
                                  sigma = 12, verbose = FALSE)
  expect_equal(grouped$kernel_dat, ordinary$kernel_dat)
  expect_equal(grouped$binned, ordinary$binned)
  exact <- kernel_prep_by_group(pts, list(a = r, b = r),
                                rep(c("a", "b"), 3), max_D = 30,
                                sigma = 12, bin = FALSE, verbose = FALSE)
  expect_null(exact$binned)
  expect_equal(exact$kernel_dat, ordinary$kernel_dat)
})

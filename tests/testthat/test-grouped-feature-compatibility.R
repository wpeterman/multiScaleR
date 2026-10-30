test_that("grouped mean, weighted SQ, and fixed metrics share pooled prediction normalization", {
  f <- refit_spatial_fixture()
  early <- f$rasters
  late <- early
  late[["hab"]] <- early[["hab"]] + 1
  late[["elev"]] <- early[["elev"]] * 1.5
  p <- sf::st_as_sf(f$points)
  points <- rbind(p, p)
  rownames(points) <- paste0("site-year-", seq_len(nrow(points)))
  year <- rep(c("early", "late"), each = nrow(p))
  x <- kernel_prep_by_group(points, list(early = early, late = late), year,
    max_D = 40, scale_vars = f$inputs$scale_vars, bin = FALSE, verbose = FALSE)
  covs <- c("mean", "rough")
  par <- c(19, 29) / x$unit_conv
  context <- list(covs = covs, n_covs = 2L, scale_vars = x$scale_vars,
    unit_conv = x$unit_conv, resolution = x$resolution, n_cols = x$n_cols,
    binned = x$binned)
  raw <- .msr_kernel_cov_w(par, x$d_list, x$raw_cov, x$kernel, context)
  scaled <- scale(raw)
  dat <- x$kernel_dat
  dat[covs] <- as.data.frame(scaled)
  dat$year <- factor(year)
  dat$y <- rep(f$data$y, 2L)
  model <- glm(y ~ mean + rough + fixed + year, dat, family = poisson())
  object <- structure(list(opt_mod = model, kernel_inputs = x,
    scale_est = data.frame(Mean = c(19, 29), SE = NA_real_, row.names = covs),
    optim_results = list(par = par), min_D = x$min_D, max_D = x$max_D,
    scl_params = .msr_merge_scl_params(list(mean = attr(scaled, "scaled:center"),
      sd = attr(scaled, "scaled:scale")), x$scl_params)), class = "multiScaleR")
  expect_identical(unname(x$raster_group), year)
  expect_identical(names(x$raster_group), rownames(dat))
  prepared <- .msr_refit_prediction_data(object, dat, x, include_fixed = TRUE)
  expect_equal(prepared[covs], dat[covs], tolerance = 1e-12)
  expect_equal(prepared$fixed, dat$fixed, tolerance = 1e-12)
  expect_equal(unname(predict(object, dat, x, include_fixed = TRUE)),
    unname(predict(model, dat, type = "response")), ignore_attr = TRUE)
  # Annual maps use the pooled fitting normalization, rather than a fresh
  # standardization within the projected year.
  rasters <- kernel_scale.raster(late, multiScaleR = object, scale_center = TRUE,
    include_fixed = TRUE, match_buffer = TRUE, site_values = list(year = "late"),
    verbose = FALSE)
  selected <- dat[year == "late", , drop = FALSE]
  expect_equal(terra::extract(rasters[[covs]], f$points, ID = FALSE),
    unname(selected[covs]), ignore_attr = TRUE, tolerance = 1e-8)
  output <- predict(object, rasters, se.fit = TRUE, uncertainty = "conditional")
  expected <- predict(model, selected, type = "response", se.fit = TRUE)
  actual <- terra::extract(output, f$points, ID = FALSE)
  expect_equal(actual$fit, unname(expected$fit), tolerance = 1e-8)
  expect_equal(actual$se.fit, unname(expected$se.fit), tolerance = 1e-8)
  snapshot <- prepare_prediction(object, dat, x, include_fixed = TRUE)
  expect_equal(predict(snapshot), predict(object, dat, x, include_fixed = TRUE))
})

test_that("analytic raster design cannot borrow a missing predictor from a formula environment", {
  z <- 2
  model <- list(analytic = list(contract = list(terms = terms(y ~ z + rough))))
  expect_error(.msr_native_raster_eval(model, data.frame(rough = 1)),
    "missing variables: z")
})

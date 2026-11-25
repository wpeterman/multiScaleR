# Regression tests for simulation CRS handling, the kernel_prep() start-scale
# hint, and model-list validation in aic_tab()/bic_tab() (audit findings
# U-7, B-7, and M-9).

test_that("sim_dat_unmarked() points carry the raster CRS without changing values", {
  rs <- terra::subset(sim_rast(dim = 50, resolution = 10, user_seed = 99),
                      c("bin1", "cont1"))
  args <- list(alpha = 0.5, beta = c(0.5, -0.25), kernel = "gaussian",
    sigma = c(40, 60), n_points = 15, n_surv = 4, det = 0.6, type = "count",
    max_D = 250, user_seed = 99)
  # Without a CRS the result is the reference simulation.
  plain <- do.call(sim_dat_unmarked, c(args, list(raster_stack = rs)))
  labeled_rs <- rs
  terra::crs(labeled_rs) <- "EPSG:32617"
  warnings <- character()
  labeled <- withCallingHandlers(
    do.call(sim_dat_unmarked, c(args, list(raster_stack = labeled_rs))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_false(any(grepl("No CRS specified", warnings)))
  expect_true(sf::st_crs(labeled$pts) == sf::st_crs(terra::crs(labeled_rs)))
  expect_true(is.na(sf::st_crs(plain$pts)))
  # Same seed, same locations and simulated values; only the CRS label differs.
  expect_identical(labeled$y, plain$y)
  expect_identical(labeled$df, plain$df)
  expect_equal(sf::st_coordinates(labeled$pts), sf::st_coordinates(plain$pts))
  expect_identical(labeled$pts$obs, plain$pts$obs)
})

test_that("sim_dat() points carry the raster CRS without changing values", {
  rs <- terra::subset(sim_rast(dim = 50, resolution = 10, user_seed = 42), "cont1")
  args <- list(alpha = 0.5, beta = 0.75, kernel = "gaussian", sigma = 40,
    type = "count", n_points = 15, max_D = 250, user_seed = 42)
  plain <- do.call(sim_dat, c(args, list(raster_stack = rs)))
  labeled_rs <- rs
  terra::crs(labeled_rs) <- "EPSG:32617"
  expect_no_warning(labeled <- do.call(sim_dat, c(args, list(raster_stack = labeled_rs))))
  expect_true(sf::st_crs(labeled$pts) == sf::st_crs(terra::crs(labeled_rs)))
  expect_identical(labeled$obs, plain$obs)
  expect_identical(labeled$df, plain$df)
  expect_equal(sf::st_coordinates(labeled$pts), sf::st_coordinates(plain$pts))
})

test_that("a metric constant at the kernel_prep() starting scale suggests sigma", {
  # Three classes in 40 m blocks on a 10 m grid. Every 100 m buffer (the
  # default starting radius, max_D / 2) holds all three classes, so patch
  # richness is constant there, while 15 m buffers hold one or two.
  rc <- expand.grid(col = 0:39, row = 0:39)
  cover <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 400, ymin = 0,
    ymax = 400, crs = "EPSG:32617", names = "cover",
    vals = ((rc$row %/% 4) + (rc$col %/% 4)) %% 3 + 1)
  set.seed(4)
  xy <- cbind(stats::runif(25, 120, 280), stats::runif(25, 120, 280))
  pts <- terra::vect(xy, crs = "EPSG:32617")
  specs <- msr_vars(rich = landscape_var("cover", metric = "pr"))
  err <- tryCatch(kernel_prep(pts = pts, raster_stack = cover, max_D = 200,
    scale_vars = specs, verbose = FALSE), error = conditionMessage)
  expect_match(err, "covariate `rich` has zero variance", fixed = TRUE)
  expect_match(err, "constant across the sample points at the starting scale", fixed = TRUE)
  expect_match(err, "`sigma =`", fixed = TRUE)
  # The suggested route gets past the refusal.
  prepared <- kernel_prep(pts = pts, raster_stack = cover, max_D = 200,
    sigma = 15, scale_vars = specs, verbose = FALSE)
  expect_s3_class(prepared, "multiScaleR_data")
  # A fixed radius has no starting scale, so the message stays generic.
  err_fixed <- tryCatch(kernel_prep(pts = pts, raster_stack = cover, max_D = 200,
    scale_vars = msr_vars(rich = landscape_var("cover", metric = "pr", radius = 100)),
    verbose = FALSE), error = conditionMessage)
  expect_match(err_fixed, "has zero variance", fixed = TRUE)
  expect_no_match(err_fixed, "starting scale", fixed = TRUE)
  # Optimization-time checks keep the generic message as well.
  specs <- msr_vars(rich = landscape_var("cover", metric = "pr"))
  expect_identical(.msr_start_scale_hint("rich", specs, "optimized covariates"), "")
})

test_that("aic_tab() and bic_tab() name malformed model lists", {
  fix <- make_core_fixture()
  for (fn in c("aic_tab", "bic_tab")) {
    tab_fun <- get(fn)
    expect_error(tab_fun(list()), "`mod_list` is empty")
    expect_error(tab_fun(fix$opt), paste0("wrap a single fit in list\\(\\)"))
    err <- tryCatch(tab_fun(list(a = fix$opt, b = NULL)), error = conditionMessage)
    expect_match(err, "Not a model: `b` (NULL)", fixed = TRUE)
    expect_no_match(err, "FUN(X", fixed = TRUE)
    err <- tryCatch(tab_fun(list(fix$opt, "fit2")), error = conditionMessage)
    expect_match(err, "element 2 (a `character`)", fixed = TRUE)
    # A Bayesian element keeps its existing refusal.
    bayes <- structure(list(), class = "multiScaleR_bayes")
    expect_error(tab_fun(list(fix$opt, bayes)), "maximized likelihoods")
  }
  # Valid lists, including a plain native model, are unaffected.
  expect_s3_class(aic_tab(list(fix$opt, fix$plain_mod), AICc = FALSE), "data.frame")
})

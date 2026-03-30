test_that("kernel_scale.raster supports both fft and focal smoothing", {
  r <- terra::rast(matrix(runif(20 * 20), 20, 20))
  names(r) <- "x"

  fft_r <- kernel_scale.raster(
    raster_stack = r,
    sigma = 3,
    kernel = "gaussian",
    fft = TRUE,
    verbose = FALSE
  )
  focal_r <- kernel_scale.raster(
    raster_stack = r,
    sigma = 3,
    kernel = "gaussian",
    fft = FALSE,
    verbose = FALSE
  )

  fft_mat <- terra::as.matrix(fft_r, wide = TRUE)
  focal_mat <- terra::as.matrix(focal_r, wide = TRUE)

  expect_true(inherits(fft_r, "SpatRaster"))
  expect_equal(names(fft_r), "x")
  expect_equal(fft_mat[10, 10], focal_mat[10, 10], tolerance = 1e-8)
})

test_that("kernel_scale.raster adds numeric site covariates as dummy rasters", {
  fix <- make_core_fixture()

  scaled <- kernel_scale.raster(
    raster_stack = fix$rs,
    multiScaleR = fix$opt,
    scale_center = TRUE,
    clamp = TRUE,
    verbose = FALSE
  )

  expect_true(all(c("cont1", "site") %in% names(scaled)))

  site_vals <- terra::values(scaled[["site"]], mat = FALSE)
  expect_true(all(site_vals[is.finite(site_vals)] == 0))

  cont1_vals <- terra::values(scaled[["cont1"]], mat = FALSE)
  train_range <- range(fix$opt$opt_mod$model$cont1, na.rm = TRUE)

  expect_gte(min(cont1_vals, na.rm = TRUE), train_range[1] - 1e-8)
  expect_lte(max(cont1_vals, na.rm = TRUE), train_range[2] + 1e-8)
})

test_that("kernel_scale.raster warns and skips categorical site covariates", {
  fix <- make_core_fixture()

  expect_warning(
    scaled <- kernel_scale.raster(
      raster_stack = fix$rs,
      multiScaleR = fix$opt_factor,
      scale_center = FALSE,
      verbose = FALSE
    ),
    "categorical"
  )

  expect_false("group" %in% names(scaled))
})

test_that("plot helpers return ggplot objects without error", {
  fix <- make_core_fixture()

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  kernel_plot <- plot_kernel(sigma = 30, kernel = "gaussian")
  opt_plots <- plot(fix$opt, prob = 0.95)
  marginal_plots <- plot_marginal_effects(fix$opt)

  expect_s3_class(kernel_plot, "ggplot")
  expect_length(opt_plots, 1)
  expect_length(marginal_plots, 2)
  expect_true(all(vapply(opt_plots, inherits, logical(1), what = "ggplot")))
  expect_true(all(vapply(marginal_plots, inherits, logical(1), what = "ggplot")))
})

test_that("plot_marginal_effects covers unmarked and zeroinfl models", {
  unmark <- make_unmarked_fixture()
  zi <- make_zeroinfl_fixture()

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  unmark_plots <- plot_marginal_effects(unmark$obj, type = "state", length.out = 5)
  zi_plots <- plot_marginal_effects(zi$obj, length.out = 5)

  expect_length(unmark_plots, 2)
  expect_length(zi_plots, 1)
  expect_error(plot_marginal_effects(unmark$obj, type = NULL, length.out = 5), "must be specified")
})

test_that("plot_marginal_effects mocked branches cover HLfit, missing CI, and numeric predictions", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  numeric_plots <- with_mocked_bindings(
    {
      fake <- structure(
        list(opt_mod = structure(list(), class = "fake_numeric"),
             scl_params = list(mean = c(x = 0), sd = c(x = 1))),
        class = "multiScaleR"
      )
      plot_marginal_effects(fake, length.out = 5, link = FALSE)
    },
    namespace = function(x) invisible("fakepkg"),
    find_predictors = function(x) list(c("x")),
    get_data = function(x, ...) data.frame(y = 1:5, x = seq(1, 5)),
    link_inverse = function(x) NULL,
    predict = function(object, newdata, se.fit = TRUE) seq_len(nrow(newdata)),
    .package = "multiScaleR"
  )

  hlfit_plots <- with_mocked_bindings(
    {
      fake <- structure(
        list(opt_mod = structure(list(), class = "HLfit"),
             scl_params = list(mean = c(x = 0), sd = c(x = 1))),
        class = "multiScaleR"
      )
      plot_marginal_effects(fake, length.out = 5, link = FALSE)
    },
    namespace = function(x) invisible("fakepkg"),
    find_predictors = function(x) list(c("x")),
    get_data = function(x, ...) data.frame(y = 1:5, x = seq(1, 5)),
    link_inverse = function(x) NULL,
    predict = function(object, newdata, variances = NULL, re.form = NULL) {
      out <- seq_len(nrow(newdata))
      attr(out, "fixefVar") <- rep(0.25, nrow(newdata))
      out
    },
    .package = "multiScaleR"
  )

  zeroinfl_plots <- with_mocked_bindings(
    {
      fake <- structure(
        list(opt_mod = structure(list(), class = "zeroinfl"),
             scl_params = list(mean = c(x = 0), sd = c(x = 1))),
        class = "multiScaleR"
      )
      plot_marginal_effects(fake, length.out = 5, link = FALSE)
    },
    namespace = function(x) invisible("fakepkg"),
    find_predictors = function(x) list(c("x")),
    get_data = function(x, ...) data.frame(y = 1:5, x = seq(1, 5)),
    get_predicted = function(x, data) data.frame(Predicted = seq_len(nrow(data))),
    .package = "multiScaleR"
  )

  expect_length(numeric_plots, 1)
  expect_length(hlfit_plots, 1)
  expect_length(zeroinfl_plots, 1)
})

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

test_that("kernel_scale.raster ignores unused default scale specs in fitted models", {
  r <- terra::rast(list(
    a = terra::rast(matrix(runif(20 * 20), 20, 20)),
    b = terra::rast(matrix(runif(20 * 20), 20, 20)),
    c = terra::rast(matrix(runif(20 * 20), 20, 20))
  ))
  names(r) <- c("a", "b", "c")

  dat <- data.frame(y = rnorm(10), a = rnorm(10), b = rnorm(10))
  mod <- glm(y ~ a + b, data = dat)

  obj <- structure(
    list(
      scale_est = data.frame(Mean = c(a = 25, b = 40),
                             SE = c(a = 5, b = 6)),
      shape_est = NULL,
      kernel_inputs = list(
        kernel = "gaussian",
        scale_vars = multiScaleR:::.msr_default_scale_vars(r)
      ),
      opt_mod = mod
    ),
    class = "multiScaleR"
  )

  scaled <- kernel_scale.raster(
    raster_stack = r,
    multiScaleR = obj,
    scale_center = FALSE,
    verbose = FALSE
  )

  expect_true(inherits(scaled, "SpatRaster"))
  expect_equal(names(scaled), c("a", "b"))
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

test_that("plot.multiScaleR includes the updated corner annotation label", {
  fix <- make_core_fixture()

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  opt_plots <- plot(fix$opt, prob = 0.95, add_label = TRUE)
  text_layer <- opt_plots[[1]]$layers[[3]]

  expect_true(inherits(text_layer$geom, "GeomText"))
  expect_match(text_layer$aes_params$label, "95% density")
  expect_match(text_layer$aes_params$label, "95% CI")
})

test_that("plot_marginal_effects covers unmarked and zeroinfl models", {
  unmark <- make_unmarked_fixture()
  zi <- make_zeroinfl_fixture()

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  unmark_plots <- plot_marginal_effects(unmark$obj, type = "state", length.out = 5)
  zi_plots <- expect_no_warning(plot_marginal_effects(zi$obj, length.out = 5))

  expect_length(unmark_plots, 2)
  expect_length(zi_plots, 1)
  expect_error(plot_marginal_effects(unmark$obj, type = NULL, length.out = 5), "must be specified")
})

test_that("plot_marginal_effects uses nested analysis models for wrapped clogit fits", {
  testthat::skip_if_not_installed("survival")

  fix <- make_core_fixture()
  n <- nrow(fix$kernel_inputs$kernel_dat)
  dat <- data.frame(
    case_ = rep(c(TRUE, FALSE, FALSE), length.out = n),
    cont1 = fix$kernel_inputs$kernel_dat$cont1,
    site = seq_len(n) / 10,
    stratum = factor(rep(seq_len(ceiling(n / 3)), each = 3, length.out = n))
  )
  coxph <- survival::coxph
  strata <- survival::strata

  inner <- suppressWarnings(
    survival::clogit(
      case_ ~ cont1 + site + strata(stratum),
      data = dat,
      model = TRUE
    )
  )
  wrapped <- structure(
    list(model = inner, sl_ = NULL, ta_ = NULL, more = NULL),
    class = c("fit_clogit", "list")
  )
  obj <- structure(
    list(
      opt_mod = wrapped,
      scl_params = fix$kernel_inputs$scl_params
    ),
    class = "multiScaleR"
  )

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  plots <- expect_no_warning(plot_marginal_effects(obj, length.out = 5))

  expect_named(plots, c("cont1", "site"))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
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

test_that("plot_marginal_effects reports prediction failures clearly", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_error(
    with_mocked_bindings(
      {
        fake <- structure(
          list(
            opt_mod = structure(list(), class = "fake_predict_error"),
            scl_params = list(mean = c(x = 0), sd = c(x = 1))
          ),
          class = "multiScaleR"
        )
        plot_marginal_effects(fake, length.out = 5, link = FALSE)
      },
      namespace = function(x) invisible("fakepkg"),
      find_predictors = function(x) list(c("x")),
      get_data = function(x, ...) data.frame(y = 1:5, x = seq(1, 5)),
      link_inverse = function(x) NULL,
      predict = function(object, newdata, se.fit = TRUE) stop("mock prediction failure"),
      .package = "multiScaleR"
    ),
    "Failed to compute marginal effects for covariate 'x'"
  )
})

test_that("plot_marginal_effects annotates polynomial and interaction terms", {
  set.seed(1)
  dat <- data.frame(
    y = rpois(40, lambda = 4),
    x = seq(-1, 1, length.out = 40),
    z = rep(seq(-0.5, 0.5, length.out = 10), each = 4)
  )
  dat$x2 <- dat$x^2

  mod <- glm(y ~ x + I(x^2) + z + x:z, family = poisson(), data = dat)
  obj <- structure(
    list(
      opt_mod = mod,
      scl_params = list(mean = c(x = 0, z = 0), sd = c(x = 1, z = 1))
    ),
    class = "multiScaleR"
  )

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_message(
    plots <- plot_marginal_effects(obj, length.out = 6),
    "Interaction term\\(s\\) detected"
  )

  expect_named(plots, c("x", "z"))
  expect_match(plots$x$labels$subtitle, "Includes: I\\(x\\^2\\)")
  expect_match(plots$x$labels$subtitle, "At mean of: z")
  expect_match(plots$z$labels$subtitle, "At mean of: x")
})

test_that("profile_sigma returns complete profiles and validates inputs", {
  fix <- make_core_fixture()

  prof <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt, n_pts = 5, metric = "AICc", verbose = FALSE)
    )
  )
  prof_ll <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt, n_pts = 4, metric = "LL", verbose = FALSE)
    )
  )

  expect_s3_class(prof, "sigma_profile")
  expect_equal(prof$metric, "AICc")
  expect_equal(sort(unique(prof$profiles$variable)), rownames(fix$opt$scale_est))
  expect_equal(nrow(prof$profiles), 5 * nrow(fix$opt$scale_est))
  expect_true(all(c("variable", "sigma", "LL", "AICc") %in% names(prof$profiles)))
  expect_true(all(is.finite(prof$profiles$sigma)))
  expect_equal(prof_ll$metric, "LL")

  expect_error(profile_sigma(structure(list(), class = "not_multiScaleR")), "multiScaleR")
  expect_error(profile_sigma(fix$opt, n_pts = 2, verbose = FALSE), "n_pts")
  expect_error(profile_sigma(fix$opt, n_pts = 5, verbose = NA), "verbose")
  expect_error(profile_sigma(fix$opt, n_cores = 0, verbose = FALSE), "n_cores")
  expect_error(profile_sigma(fix$opt, n_cores = 1.5, verbose = FALSE), "n_cores")
})

test_that("profile_sigma supports linear and custom sigma grids", {
  fix <- make_core_fixture()

  prof_linear <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt,
                    n_pts = 3,
                    spacing = "linear",
                    sigma_range = c(10, 30),
                    verbose = FALSE)
    )
  )
  prof_custom <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt,
                    sigma_values = c(30, 10, 20, 20),
                    verbose = FALSE)
    )
  )

  expect_equal(prof_linear$spacing, "linear")
  expect_equal(prof_linear$sigma_grid, c(10, 20, 30))
  expect_equal(unique(prof_linear$profiles$sigma), c(10, 20, 30))

  expect_equal(prof_custom$spacing, "custom")
  expect_equal(prof_custom$sigma_grid, c(10, 20, 30))
  expect_equal(unique(prof_custom$profiles$sigma), c(10, 20, 30))

  expect_error(profile_sigma(fix$opt, spacing = "equal", verbose = FALSE),
               "should be one of")
  expect_error(profile_sigma(fix$opt, sigma_range = c(10, 10), verbose = FALSE),
               "distinct")
  expect_error(profile_sigma(fix$opt, sigma_values = c(10, 10, 20),
                             verbose = FALSE),
               "at least 3 unique")
})

test_that("profile_sigma parallel output matches serial output", {
  cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(cores) || cores < 2) {
    skip("Parallel parity test requires at least 2 physical cores.")
  }

  fix <- make_core_fixture()

  prof_serial <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt, n_pts = 4, verbose = FALSE)
    )
  )
  prof_parallel <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt, n_pts = 4, n_cores = 2, verbose = FALSE)
    )
  )

  expect_identical(prof_parallel$metric, prof_serial$metric)
  expect_identical(prof_parallel$spacing, prof_serial$spacing)
  expect_equal(prof_parallel$sigma_grid, prof_serial$sigma_grid)
  expect_equal(prof_parallel$opt_sigma, prof_serial$opt_sigma)
  expect_equal(prof_parallel$profiles, prof_serial$profiles, tolerance = 1e-10)
})

test_that("profile_sigma uses the same parameter count logic as model selection", {
  fix <- make_core_fixture()
  opt_multi <- fix$opt
  opt_multi$scale_est <- data.frame(
    Mean = c(40, 80),
    SE = c(1, 1),
    row.names = c("cont1", "site")
  )

  prof <- with_mocked_bindings(
    profile_sigma(opt_multi, n_pts = 3, verbose = FALSE),
    kernel_scale_fn = function(...) 10,
    .package = "multiScaleR"
  )

  expected_k <- nrow(insight::get_parameters(opt_multi$opt_mod)) +
    nrow(opt_multi$scale_est)
  expected_aic <- -2 * (-10) + 2 * expected_k
  expected_aicc <- expected_aic +
    (2 * expected_k * (expected_k + 1)) /
    (.msr_model_nobs(opt_multi$opt_mod) - expected_k - 1)

  expect_equal(unique(prof$profiles$AICc), expected_aicc)
})

test_that("plot.sigma_profile returns named ggplots and rejects bad input", {
  fix <- make_core_fixture()
  prof <- suppressWarnings(
    suppressMessages(
      profile_sigma(fix$opt, n_pts = 4, metric = "LL", verbose = FALSE)
    )
  )

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  plots <- plot(prof)

  expect_length(plots, nrow(fix$opt$scale_est))
  expect_named(plots, rownames(fix$opt$scale_est))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_error(plot.sigma_profile(structure(list(), class = "not_sigma_profile")), "sigma_profile")
})

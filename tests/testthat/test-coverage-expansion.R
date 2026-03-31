test_that("simulation helpers cover alternate response types, spatial inputs, and validations", {
  fix <- make_simulation_fixture()
  rs1 <- terra::subset(fix$rs, "cont1")

  pts_sf <- sf::st_as_sf(
    data.frame(x = c(100, 200, 300), y = c(100, 200, 300)),
    coords = c("x", "y"),
    crs = terra::crs(rs1)
  )
  pts_sv <- terra::vect(pts_sf)

  sim_auto <- sim_dat(
    alpha = 0.25,
    beta = 0.5,
    kernel = "gaussian",
    sigma = 25,
    type = "count",
    n_points = pts_sv,
    raster_stack = rs1,
    max_D = NULL,
    user_seed = 12
  )
  sim_nb <- sim_dat(
    alpha = 0.25,
    beta = 0.5,
    kernel = "gaussian",
    sigma = 25,
    type = "count_nb",
    n_points = pts_sv,
    raster_stack = rs1,
    max_D = 120,
    user_seed = 13
  )
  sim_occ <- sim_dat(
    alpha = -0.25,
    beta = 0.5,
    kernel = "gaussian",
    sigma = 25,
    type = "occ",
    n_points = pts_sv,
    raster_stack = rs1,
    max_D = 120,
    user_seed = 14
  )
  sim_gaussian <- sim_dat(
    alpha = 0.25,
    beta = 0.5,
    kernel = "gaussian",
    sigma = 25,
    type = "gaussian",
    StDev = 0.2,
    n_points = pts_sv,
    raster_stack = rs1,
    max_D = 120,
    user_seed = 15
  )

  expect_equal(nrow(sim_auto$df), 3)
  expect_true(all(sim_nb$obs >= 0))
  expect_true(all(sim_occ$obs %in% 0:1))
  expect_true(is.numeric(sim_gaussian$obs))
  expect_true(multiScaleR:::is_spatial(pts_sf))
  expect_true(multiScaleR:::is_spatial(sf::st_geometry(pts_sf)))
  expect_true(multiScaleR:::is_spatial(pts_sv))
  expect_false(multiScaleR:::is_spatial(data.frame(x = 1)))

  expect_error(
    sim_dat(
      beta = c(0.5, 0.25),
      sigma = 25,
      raster_stack = rs1,
      max_D = 120,
      n_points = 3
    ),
    "beta"
  )
  expect_error(
    sim_dat(
      beta = 0.5,
      sigma = 25,
      kernel = "expow",
      raster_stack = rs1,
      max_D = 120,
      n_points = 3
    ),
    "Shape parameter"
  )
  expect_error(
    sim_dat(
      beta = 0.5,
      sigma = 25,
      raster_stack = rs1,
      max_D = 120,
      n_points = "bad"
    ),
    "spatVector or sf"
  )
})

test_that("sim_dat_unmarked covers alternate response types and validations", {
  fix <- make_simulation_fixture()

  sim_nb <- sim_dat_unmarked(
    alpha = 0.25,
    beta = c(0.4, -0.2),
    kernel = "gaussian",
    sigma = c(20, 30),
    type = "count_nb",
    n_points = 6,
    n_surv = 3,
    det = 0.5,
    raster_stack = fix$rs_two,
    user_seed = 16
  )
  sim_occ <- sim_dat_unmarked(
    alpha = -0.5,
    beta = c(0.4, -0.2),
    kernel = "gaussian",
    sigma = c(20, 30),
    type = "occ",
    n_points = 6,
    n_surv = 3,
    det = 0.5,
    raster_stack = fix$rs_two,
    user_seed = 17
  )

  expect_equal(dim(sim_nb$y), c(6, 3))
  expect_true(all(sim_nb$y >= 0))
  expect_true(all(sim_occ$y %in% 0:1))

  expect_error(
    sim_dat_unmarked(
      beta = c(0.5, 0.25, 0.1),
      sigma = c(20, 30),
      raster_stack = fix$rs_two,
      n_points = 5
    ),
    "beta"
  )
  expect_error(
    sim_dat_unmarked(
      beta = c(0.5, 0.25),
      sigma = c(20, 30),
      kernel = "expow",
      raster_stack = fix$rs_two,
      n_points = 5
    ),
    "Shape parameter"
  )
})

test_that("helper math functions cover validation and alternate branches", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted <- sim_rast(dim = 15, resolution = 5, plot = TRUE, user_seed = NULL)
  expect_true(inherits(plotted, "SpatRaster"))

  expect_error(
    fft_convolution(matrix(1:4, 2, 2), array(1, dim = c(1, 1, 1))),
    "kernel"
  )
  expect_error(
    fft_convolution(matrix("a", 2, 2), matrix(1, 1, 1)),
    "numeric"
  )
  expect_error(
    fft_convolution(matrix(1:4, 2, 2), matrix(1, 1, 1), fun = "median"),
    "either 'mean' or 'sum'"
  )
  expect_error(
    fft_convolution(matrix(1:4, 2, 2), matrix(0, 1, 1), fun = "mean"),
    "Sum of kernel is zero"
  )

  x <- matrix(c(5, 1), ncol = 2, dimnames = list("a", c("Estimate", "SE")))
  ci_no_min <- multiScaleR:::ci_func_r(x, df = 20, min_D = NULL, names = "a")
  expect_true(all(is.finite(ci_no_min[1, c("2.5%", "97.5%")])))

  d <- c(0, 1, 2)
  rs_df <- data.frame(x = c(1, 2, 3))
  expect_equal(length(scale_type_r(d, kernel = "exp", sigma = 2, r_stack.df = rs_df)), 1)
  expect_equal(length(scale_type_r(d, kernel = "fixed", sigma = 2, r_stack.df = rs_df)), 1)
  expect_equal(length(scale_type_r(d, kernel = "expow", sigma = 2, shape = 2, r_stack.df = rs_df)), 1)
  expect_error(scale_type(d, kernel = "gaussian", sigma = 2), "r_stack_df")
  expect_error(multiScaleR:::k_dist(1, prob = 1, kernel = "gaussian"), "between 0 and 1")
  expect_error(multiScaleR:::k_dist(1, prob = 0.9, kernel = "expow"), "beta must be specified")
  expect_error(multiScaleR:::k_dist(1, prob = 0.9, kernel = "expow", beta = -1), "beta must be positive")
})

test_that("scale_center_raster covers multiScaleR_data, gls, unmarked, and fallback branches", {
  skip_if_not_installed("nlme")

  fix <- make_core_fixture()
  ki_scaled <- multiScaleR:::scale_center_raster(
    r = fix$rs,
    multiScaleR = fix$kernel_inputs,
    clamp = FALSE
  )

  raw_vals <- terra::values(fix$rs, mat = FALSE)
  scaled_vals <- terra::values(ki_scaled, mat = FALSE)
  expect_true(inherits(ki_scaled, "SpatRaster"))
  expect_equal(
    scaled_vals[1],
    unname(
      (raw_vals[1] - fix$kernel_inputs$scl_params$mean["cont1"]) /
        fix$kernel_inputs$scl_params$sd["cont1"]
    )
  )

  gls_mod <- nlme::gls(y ~ cont1 + site, data = fix$df)
  gls_obj <- structure(
    list(opt_mod = gls_mod, scl_params = fix$opt$scl_params),
    class = "multiScaleR"
  )
  gls_scaled <- multiScaleR:::scale_center_raster(
    r = fix$rs,
    multiScaleR = gls_obj,
    clamp = FALSE
  )
  expect_true(inherits(gls_scaled, "SpatRaster"))

  unmark <- make_unmarked_fixture()
  unmark_scaled <- multiScaleR:::scale_center_raster(
    r = unmark$rs_two,
    multiScaleR = unmark$obj,
    clamp = FALSE
  )
  expect_equal(names(unmark_scaled), c("bin1", "cont1"))

  fake_r <- terra::rast(matrix(1:4, 2, 2))
  names(fake_r) <- "x"
  fake_obj <- structure(
    list(
      opt_mod = structure(list(), class = "fakeScaleCenter"),
      scl_params = list(mean = c(x = 2), sd = c(x = 4))
    ),
    class = "multiScaleR"
  )

  fallback_scaled <- with_mocked_bindings(
    multiScaleR:::scale_center_raster(
      r = fake_r,
      multiScaleR = fake_obj,
      clamp = FALSE
    ),
    get_data = function(mod, effects = "all") NULL,
    extract_model_data = function(mod) data.frame(y = 1:4, x = 1:4),
    .package = "multiScaleR"
  )
  expect_equal(terra::values(fallback_scaled, mat = FALSE), c(-0.25, 0.25, 0, 0.5))
})

test_that("kernel_scale.raster and site-covariate helpers cover validation and fallback branches", {
  fix <- make_core_fixture()

  expect_warning(
    scaled <- kernel_scale.raster(
      raster_stack = fix$rs,
      scale_opt = fix$kernel_inputs,
      scale_center = TRUE,
      verbose = FALSE
    ),
    "deprecated"
  )
  expect_equal(names(scaled), "cont1")

  expect_output(
    kernel_scale.raster(
      raster_stack = fix$rs,
      sigma = 3,
      kernel = "gaussian",
      verbose = TRUE
    ),
    "Smoothing spatRaster"
  )

  bad_r <- fix$rs
  names(bad_r) <- "other"

  expect_error(
    kernel_scale.raster(bad_r, multiScaleR = fix$opt, verbose = FALSE),
    "optimized covariate"
  )
  expect_error(
    kernel_scale.raster(bad_r, multiScaleR = fix$kernel_inputs, verbose = FALSE),
    "optimized covariate"
  )
  expect_error(kernel_scale.raster(fix$rs, sigma = NULL, verbose = FALSE), "sigma values")
  expect_error(
    kernel_scale.raster(data.frame(x = 1:3), sigma = 1, verbose = FALSE),
    "SpatRaster"
  )
  expect_error(
    kernel_scale.raster(fix$rs, sigma = 1, na.rm = 1, verbose = FALSE),
    "na.rm"
  )
  expect_error(
    kernel_scale.raster(fix$rs, sigma = 1, fft = 1, verbose = FALSE),
    "fft"
  )
  expect_warning(
    kernel_scale.raster(
      raster_stack = terra::subset(fix$rs_all, c("bin1", "cont1")),
      sigma = 1,
      kernel = "gaussian",
      verbose = FALSE
    ),
    "Number of sigma values"
  )

  base_r <- terra::rast(matrix(1:4, 2, 2))
  names(base_r) <- "x"

  plain_out <- multiScaleR:::.add_site_covariate_rasters(
    smooth_stack = base_r,
    multiScaleR = structure(list(), class = "not_multiScaleR"),
    raster_covs = "x"
  )
  expect_equal(names(plain_out), "x")

  null_mod_out <- multiScaleR:::.add_site_covariate_rasters(
    smooth_stack = base_r,
    multiScaleR = structure(list(opt_mod = NULL), class = "multiScaleR"),
    raster_covs = "x"
  )
  expect_equal(names(null_mod_out), "x")

  assign("model.frame.fakeSiteError", function(formula, ...) stop("boom"), envir = .GlobalEnv)
  assign(
    "model.frame.fakeListSite",
    function(formula, ...) data.frame(y = 1:3, weird = I(list(1, 2, 3))),
    envir = .GlobalEnv
  )
  on.exit(
    rm(list = c("model.frame.fakeSiteError", "model.frame.fakeListSite"), envir = .GlobalEnv),
    add = TRUE
  )

  expect_warning(
    err_out <- multiScaleR:::.add_site_covariate_rasters(
      smooth_stack = base_r,
      multiScaleR = structure(
        list(opt_mod = structure(list(), class = "fakeSiteError")),
        class = "multiScaleR"
      ),
      raster_covs = "x"
    ),
    "Could not inspect fitted model frame"
  )
  expect_equal(names(err_out), "x")

  no_site_mod <- glm(am ~ wt, family = binomial(), data = mtcars)
  no_site_r <- terra::rast(matrix(runif(4), 2, 2))
  names(no_site_r) <- "wt"
  no_site_out <- multiScaleR:::.add_site_covariate_rasters(
    smooth_stack = no_site_r,
    multiScaleR = structure(list(opt_mod = no_site_mod), class = "multiScaleR"),
    raster_covs = "wt"
  )
  expect_equal(names(no_site_out), "wt")

  expect_warning(
    weird_out <- multiScaleR:::.add_site_covariate_rasters(
      smooth_stack = base_r,
      multiScaleR = structure(
        list(opt_mod = structure(list(), class = "fakeListSite")),
        class = "multiScaleR"
      ),
      raster_covs = "x"
    ),
    "not numeric/integer/logical"
  )
  expect_equal(names(weird_out), "x")
})

test_that("kernel_dist, summary, print, and plotting helpers cover alternate branches", {
  skip_if_not_installed("nlme")

  fix <- make_core_fixture()
  unmark <- make_unmarked_fixture()

  expect_error(kernel_dist(structure(list(), class = "not_multiScaleR")), "Provide a fitted")
  expect_error(kernel_dist(prob = 2, sigma = 1, kernel = "gaussian"), "within \\(0, 1\\)")
  expect_error(kernel_dist(), "Parameters not correctly specified")
  expect_error(kernel_dist(sigma = NULL, kernel = "gaussian"), "sigma")
  expect_error(kernel_dist(sigma = 1, kernel = NULL), "kernel")
  expect_error(kernel_dist(sigma = 1, kernel = "expow"), "shape")
  expect_error(kernel_dist(sigma = 1, kernel = "expow", beta = -1), "beta")

  gls_mod <- nlme::gls(y ~ cont1 + site, data = fix$df)
  gls_obj <- structure(
    list(
      opt_mod = gls_mod,
      scale_est = matrix(c(10, 1), ncol = 2, dimnames = list("cont1", c("Estimate", "SE"))),
      shape_est = NULL,
      min_D = 1,
      kernel_inputs = list(kernel = "gaussian"),
      warn_message = c(1, 2),
      call = quote(gls_call())
    ),
    class = "multiScaleR"
  )
  expect_true(is.data.frame(kernel_dist(gls_obj)))

  um_obj <- structure(
    list(
      opt_mod = unmark$mod,
      scale_est = matrix(
        c(10, 1, 12, NaN),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("bin1", "cont1"), c("Estimate", "SE"))
      ),
      shape_est = matrix(
        c(2, 0.25, 3, 0.5),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("bin1", "cont1"), c("Estimate", "SE"))
      ),
      min_D = 1,
      kernel_inputs = list(kernel = "expow"),
      warn_message = c(1, 2, 3),
      call = quote(unmarked_call())
    ),
    class = "multiScaleR"
  )

  expect_warning(dist_um <- kernel_dist(um_obj, sigma = 99), "Ignoring specified")
  nan_dist <- with_mocked_bindings(
    kernel_dist(um_obj),
    ci_func = function(x, df, min_D, names) {
      out <- matrix(
        c(10, NaN, NaN, NaN,
          12, 1, 10, 14),
        ncol = 4,
        byrow = TRUE
      )
      colnames(out) <- c("Estimate", "SE", "2.5%", "97.5%")
      rownames(out) <- c("bin1", "cont1")
      out
    },
    .package = "multiScaleR"
  )
  expect_true(is.nan(nan_dist["bin1", "Estimate"]))

  gls_summary <- summary(gls_obj, prob = 0.8)
  um_summary <- summary(um_obj, prob = 0.7)

  expect_s3_class(gls_summary, "summary_multiScaleR")
  expect_s3_class(um_summary, "summary_multiScaleR")
  expect_output(print(gls_summary), "Fitted Model Summary")
  expect_output(print(um_summary), "Optimized Kernel Shape")
  expect_output(print(um_obj), "Optimized Kernel Shape Parameter")
  expect_output(print(um_obj), "WARNING")

  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  low_prob_plots <- plot(fix$opt, prob = 0.7)
  no_scale_plots <- plot(fix$opt, scale_dist = FALSE)
  low_prob_kernel <- plot_kernel(prob = 0.7, sigma = 30, kernel = "gaussian")
  no_scale_kernel <- plot_kernel(sigma = 30, kernel = "gaussian", scale_dist = FALSE)

  expect_length(low_prob_plots, 1)
  expect_length(no_scale_plots, 1)
  expect_s3_class(low_prob_kernel, "ggplot")
  expect_s3_class(no_scale_kernel, "ggplot")
  expect_error(plot(fix$opt, prob = 2), "decimal between 0 and 1")
  expect_error(plot_kernel(sigma = NULL, kernel = "gaussian"), "sigma")
  expect_error(plot_kernel(sigma = 1, kernel = NULL), "kernel")
})

test_that("new front-end validation catches malformed user inputs early", {
  fix <- make_core_fixture()
  sim_fix <- make_simulation_fixture()

  expect_error(
    sim_dat(
      beta = 0.5,
      sigma = 25,
      raster_stack = NULL,
      max_D = 120,
      n_points = 3
    ),
    "raster_stack"
  )
  expect_error(
    sim_dat(
      beta = 0.5,
      sigma = 25,
      raster_stack = terra::subset(sim_fix$rs, "cont1"),
      max_D = 120,
      n_points = 0
    ),
    "n_points"
  )
  expect_error(
    sim_dat_unmarked(
      beta = c(0.5, 0.25),
      sigma = c(20, 30),
      raster_stack = sim_fix$rs_two,
      n_points = 5,
      n_surv = 0
    ),
    "n_surv"
  )
  expect_error(
    sim_dat_unmarked(
      beta = c(0.5, 0.25),
      sigma = c(20, 30),
      raster_stack = sim_fix$rs_two,
      n_points = 5,
      det = 2
    ),
    "det"
  )

  expect_error(plot_marginal_effects(structure(list(), class = "not_multiScaleR")), "multiScaleR")
  expect_error(plot_marginal_effects(fix$opt, length.out = 0), "length.out")
  expect_error(plot_marginal_effects(fix$opt, link = NA), "link")

  expect_error(
    multiScaleR:::scale_center_raster(
      r = data.frame(x = 1),
      multiScaleR = fix$kernel_inputs
    ),
    "SpatRaster"
  )
  expect_error(
    multiScaleR:::scale_center_raster(
      r = fix$rs,
      multiScaleR = structure(list(), class = "not_multiScaleR")
    ),
    "multiScaleR"
  )

  expect_warning(
    kernel_scale.raster(
      raster_stack = fix$rs,
      sigma = 1,
      scale_center = TRUE,
      verbose = FALSE
    ),
    "scale_center = TRUE"
  )
  expect_warning(
    kernel_scale.raster(
      raster_stack = fix$rs,
      sigma = 1,
      clamp = TRUE,
      verbose = FALSE
    ),
    "clamp = TRUE"
  )
  expect_error(
    kernel_scale.raster(
      raster_stack = fix$rs,
      sigma = 1,
      pct_wt = 1,
      verbose = FALSE
    ),
    "pct_wt"
  )
  expect_error(
    kernel_scale.raster(
      raster_stack = fix$rs,
      sigma = 1,
      pct_mx = 2,
      verbose = FALSE
    ),
    "pct_mx"
  )
})

test_that("marginal effects and model selection cover remaining alternate branches", {
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)

  unmark <- make_unmarked_fixture()
  mod_nuisance <- unmarked::pcount(~1 ~ bin1 + nuisance, data = unmark$umf, K = 50)
  nuisance_obj <- structure(
    list(opt_mod = mod_nuisance, scl_params = unmark$kernel_inputs$scl_params),
    class = "multiScaleR"
  )
  nuisance_plots <- plot_marginal_effects(nuisance_obj, type = "state", length.out = 5)
  expect_equal(names(nuisance_plots), c("bin1", "nuisance"))

  link_plots <- with_mocked_bindings(
    {
      fake <- structure(
        list(
          opt_mod = structure(list(), class = "fake_link"),
          scl_params = list(mean = c(x = 0), sd = c(x = 1))
        ),
        class = "multiScaleR"
      )
      plot_marginal_effects(fake, length.out = 5, link = TRUE)
    },
    namespace = function(x) invisible("fakepkg"),
    find_predictors = function(x) list(c("x")),
    get_data = function(x, ...) data.frame(y = 1:5, x = seq(1, 5)),
    link_inverse = function(x) exp,
    predict = function(object, newdata, se.fit = TRUE) seq_len(nrow(newdata)),
    .package = "multiScaleR"
  )
  expect_length(link_plots, 1)

  s4_obj <- methods::new("dgCMatrix")
  s4_pred <- with_mocked_bindings(
    multiScaleR:::safe_predict(s4_obj, data.frame(x = 1)),
    extract_namespace = function(x) NULL,
    getS3method = function(...) NULL,
    predict = function(object, newdata, se.fit = TRUE) list(fit = 1, se.fit = 0),
    .package = "multiScaleR"
  )
  s4_ns <- with_mocked_bindings(
    multiScaleR:::namespace(s4_obj),
    extract_namespace = function(x) NULL,
    getS3method = function(...) NULL,
    .package = "multiScaleR"
  )
  expect_true(is.list(s4_pred))
  expect_equal(s4_ns, "Matrix")

  cl_s4 <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl_s4), add = TRUE)
  expect_equal(multiScaleR:::cluster_prep(s4_obj, cl_s4), "Matrix")

  ms_obj <- structure(
    list(opt_mod = unmark$mod, kernel_inputs = list(kernel = "gaussian")),
    class = "multiScaleR"
  )

  aic_verbose <- with_mocked_bindings(
    capture.output(
      aic_tab(
        list(ms_obj, unmark$mod),
        AICc = FALSE,
        mod_names = c("opt", "plain"),
        verbose = TRUE
      )
    ),
    n_obs = function(x) 15,
    find_formula = function(x) list(conditional = c("~", "y", "x")),
    get_parameters = function(x) data.frame(Parameter = 1:2),
    logLik = function(x) structure(-10, class = "logLik"),
    aictabCustom = function(logL, K, modnames, second.ord, nobs, sort) {
      data.frame(Modnames = modnames, K = K, AIC = seq_along(K))
    },
    .package = "multiScaleR"
  )
  bic_verbose <- with_mocked_bindings(
    capture.output(
      bic_tab(
        list(ms_obj, unmark$mod),
        mod_names = c("opt", "plain"),
        verbose = TRUE
      )
    ),
    n_obs = function(x) 15,
    find_formula = function(x) list(conditional = c("~", "y", "x")),
    get_parameters = function(x) data.frame(Parameter = 1:2),
    logLik = function(x) structure(-10, class = "logLik"),
    bictabCustom = function(logL, K, modnames, nobs, sort) {
      data.frame(Modnames = modnames, K = K, BIC = seq_along(K))
    },
    .package = "multiScaleR"
  )

  aicc_verbose <- with_mocked_bindings(
    capture.output(
      aic_tab(
        list(ms_obj, unmark$mod),
        AICc = TRUE,
        mod_names = c("opt", "plain"),
        verbose = TRUE
      )
    ),
    n_obs = function(x) 15,
    find_formula = function(x) list(conditional = c("~", "y", "x")),
    get_parameters = function(x) data.frame(Parameter = 1:2),
    logLik = function(x) structure(-10, class = "logLik"),
    aictabCustom = function(logL, K, modnames, second.ord, nobs, sort) {
      data.frame(Modnames = modnames, K = K, AICc = seq_along(K))
    },
    .package = "multiScaleR"
  )

  expect_gt(length(aic_verbose), 0)
  expect_gt(length(bic_verbose), 0)
  expect_gt(length(aicc_verbose), 0)

  smaller_mod <- glm(y ~ site, family = poisson(), data = make_core_fixture()$df[-1, ])
  expect_error(bic_tab(list(make_core_fixture()$opt, smaller_mod)), "different number of sample locations")
})

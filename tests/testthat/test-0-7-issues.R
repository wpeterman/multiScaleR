test_that("adjacency projections preserve non-square raster geometry", {
  set.seed(604)
  r <- terra::rast(matrix(sample(1:3, 29 * 37, replace = TRUE),
                          nrow = 29, ncol = 37))
  names(r) <- "landcover"
  for (metric in c("ai", "pladj", "contag", "iji", "ent", "condent",
                   "joinent", "mutinf", "relmutinf", "clumpy")) {
    spec <- if (metric == "clumpy") {
      landscape_var("landcover", metric = metric, radius = 3, class = 1)
    } else {
      landscape_var("landcover", metric = metric, radius = 3)
    }
    projection <- kernel_scale.raster(
      r, scale_vars = msr_vars(metric_value = spec), verbose = FALSE
    )
    expect_equal(terra::nrow(projection), terra::nrow(r), info = metric)
    expect_equal(terra::ncol(projection), terra::ncol(r), info = metric)
    expect_true(any(is.finite(terra::values(projection))), info = metric)
  }
})

test_that("undefined candidate covariates cannot improve likelihood by dropping rows", {
  dat <- data.frame(y = c(1, 2, 4, 3), x = c(1, 2, 3, 4))
  mod <- lm(y ~ x, data = dat)
  context <- build_opt_context(mod, cov_df = lapply(seq_len(4),
                                                    function(i) data.frame(x = i)))
  expect_equal(
    with_mocked_bindings(
      kernel_scale_fn(0.5, d_list = as.list(rep(0.1, 4)),
                      cov_df = lapply(seq_len(4), function(i) data.frame(x = i)),
                      kernel = "gaussian", fitted_mod = mod,
                      opt_context = context),
      .msr_kernel_cov_w = function(...) {
        matrix(c(1, 2, NA, 4), ncol = 1,
               dimnames = list(as.character(1:4), "x"))
      },
      .package = "multiScaleR"
    ),
    1e6^10
  )
})

test_that("optimization reports a sample already reduced by the starting model", {
  fix <- make_core_fixture()
  dat <- fix$df
  dat$y[1] <- NA_real_
  initial <- glm(y ~ cont1 + site, family = poisson(), data = dat)
  expect_warning(
    fit <- multiScale_optim(initial, fix$kernel_inputs,
                            par = 40 / fix$kernel_inputs$unit_conv,
                            verbose = FALSE),
    "uses 14 of 15 prepared observations"
  )
  expect_true(diagnostics(fit)$sample_size$triggered)
  expect_equal(diagnostics(fit)$sample_size$fitted_n, 14)
})

test_that("hard-radius summaries omit Wald uncertainty and identify radius", {
  fix <- make_core_fixture()
  opt <- fix$opt
  opt$kernel_inputs$scale_vars$type[] <- "landscape"
  summary_out <- summary(opt)
  expect_identical(summary_out$hard_radius, "cont1")
  expect_true(is.na(summary_out$opt_scale$SE))
  expect_true(all(is.na(summary_out$opt_scale[, c("2.5%", "97.5%")])))
  expect_equal(summary_out$opt_dist[1, 1],
               round(opt$scale_est$Mean, 2))
  printed <- capture.output(print(summary_out))
  expect_true(any(grepl("Radius", printed)))
  expect_true(any(grepl("profile = TRUE", printed, fixed = TRUE)))
  print_fit <- capture.output(print(opt))
  expect_true(any(grepl("Radius", print_fit)))
  opt$kernel_inputs$scale_vars$weighted <- NULL
  expect_identical(summary(opt)$hard_radius, "cont1")
})

test_that("fitted hard-radius metrics store no Hessian standard error", {
  set.seed(625)
  r <- terra::rast(matrix(sample(1:3, 30 * 30, replace = TRUE),
                          nrow = 30, ncol = 30))
  terra::crs(r) <- "EPSG:26915"
  names(r) <- "landcover"
  pts <- terra::vect(cbind(runif(22, 6, 24), runif(22, 6, 24)),
                     type = "points", crs = terra::crs(r))
  vars <- msr_vars(edge = landscape_var("landcover", metric = "ed"))
  prepared <- kernel_prep(pts, r, max_D = 5, sigma = 3,
                          scale_vars = vars, bin = FALSE, verbose = FALSE)
  dat <- data.frame(y = rnorm(22) + prepared$kernel_dat$edge,
                    prepared$kernel_dat)
  initial <- lm(y ~ edge, data = dat)
  fit <- suppressWarnings(multiScale_optim(
    initial, prepared, par = 3 / prepared$unit_conv,
    n_cores = 1, verbose = FALSE
  ))
  expect_true(is.na(fit$scale_est["edge", "SE"]))
  expect_true(is.finite(kernel_dist(fit)["edge", "Mean"]))
  expect_true(is.null(diagnostics(fit)$sigma_precision))
})

test_that("projection uses input offset variables and scenario values", {
  fix <- make_core_fixture()
  opt <- fix$opt
  dat <- fix$df
  dat$n_weeks <- rep(2, nrow(dat))
  opt$opt_mod <- glm(y ~ cont1 + offset(log(n_weeks)),
                     data = dat, family = poisson())
  one <- kernel_scale.raster(fix$rs, multiScaleR = opt,
                             scale_center = TRUE, verbose = FALSE)
  four <- kernel_scale.raster(fix$rs, multiScaleR = opt,
                              scale_center = TRUE,
                              offset_values = c(n_weeks = 4), verbose = FALSE)
  expect_true("n_weeks" %in% names(one))
  expect_false("offset(log(n_weeks))" %in% names(one))
  expect_equal(unique(stats::na.omit(as.vector(terra::values(one$n_weeks)))), 1)
  expect_equal(unique(stats::na.omit(as.vector(terra::values(four$n_weeks)))), 4)
  pred_one <- terra::predict(one, opt$opt_mod, type = "response")
  pred_four <- terra::predict(four, opt$opt_mod, type = "response")
  expect_equal(terra::values(pred_four), 4 * terra::values(pred_one),
               tolerance = 1e-7)
  expect_error(kernel_scale.raster(fix$rs, multiScaleR = opt,
                                   offset_values = 4, verbose = FALSE),
               "named finite")
})

test_that("all-penalty profiles report failure", {
  fix <- make_core_fixture()
  expect_error(
    with_mocked_bindings(
      profile_sigma(fix$opt, n_pts = 3, verbose = FALSE),
      kernel_scale_fn = function(...) 1e6^10,
      .package = "multiScaleR"
    ),
    "Could not profile"
  )
  broken <- fix$opt
  broken$opt_context$refit_fn <- function(...) {
    stop("formula object missing")
  }
  expect_error(profile_sigma(broken, n_pts = 3, verbose = FALSE),
               "formula object missing")
  cache_key <- profile_scale_cache_key(broken, broken$min_D,
                                       rownames(broken$scale_est))
  if (exists(cache_key, envir = .profile_scale_cache, inherits = FALSE)) {
    rm(list = cache_key, envir = .profile_scale_cache)
  }
  expect_error(summary(broken, profile = TRUE), "formula object missing")
})

test_that("negative-binomial dispersion contributes to model selection K", {
  skip_if_not_installed("MASS")
  set.seed(11)
  dat <- data.frame(x = seq_len(90) / 90)
  dat$y <- MASS::rnegbin(90, mu = exp(1 + dat$x), theta = 2)
  nb <- MASS::glm.nb(y ~ x, data = dat)
  expect_equal(.msr_parameter_count(nb),
               as.integer(attr(logLik(nb), "df")))
  expect_equal(.msr_parameter_count(nb), length(coef(nb)) + 1L)
})

test_that("msr_vars errors name all supported constructors", {
  expect_error(msr_vars(bad = list()), "surface_var")
})

# Fixed metrics must reproduce an ordinary, already standardized numeric column.
# These bounded fixtures assess implementation contracts, not scale recovery.
fixed_inference_fixture <- function(kernel = "gaussian", bin = FALSE) {
  f <- refit_spatial_fixture(kernel = kernel, bin = bin)
  f$data$key <- rownames(f$data)
  f$reference <- f$inputs
  f$reference$scale_vars <- f$inputs$scale_vars[
    f$inputs$scale_vars$covariate != "fixed", , drop = FALSE]
  f
}

test_that("legacy binned inputs without scale_vars retain their original contract", {
  f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
  empty <- list(covariates = character(), center = numeric(), spread = numeric())
  expect_identical(.msr_fixed_reference_contract(f$inputs, f$data,
    seq_len(nrow(f$data)), "habitat"), empty)
  native <- glm(count ~ temp, data = f$data, family = poisson())
  wrapper <- structure(list(opt_mod = native, kernel_inputs = f$inputs), class = "multiScaleR")
  expect_identical(.msr_fixed_training_contract(wrapper), empty)
  expect_identical(.msr_attach_fixed_reference(wrapper), wrapper)
  expect_identical(.msr_fixed_prediction_data(wrapper, f$data, f$inputs,
    seq_len(nrow(f$data)), list(), numeric(), NULL), f$data)
  legacy <- prepare_bayesian_fixture(f, family = "gaussian")
  declared <- f
  declared$inputs$scale_vars <- msr_vars(habitat = kernel_var("habitat"))
  reference <- prepare_bayesian_fixture(declared, family = "gaussian")
  expect_identical(legacy$specification, reference$specification)
  expect_identical(legacy$parameters, reference$parameters)
  skip_if_not_installed("RTMB")
  a <- RTMB::MakeADFun(.msr_bayes_objective(legacy$specification, TRUE), legacy$parameters, silent = TRUE)
  b <- RTMB::MakeADFun(.msr_bayes_objective(reference$specification, TRUE), reference$parameters, silent = TRUE)
  expect_equal(a$fn(a$par), b$fn(b$par), tolerance = 0)
  expect_equal(a$gr(a$par), b$gr(b$par), tolerance = 0)
  expect_equal(a$he(a$par), b$he(b$par), tolerance = 0)
})

test_that("fixed features leave analytic objectives and scores unchanged", {
  skip_if_not_installed("MASS")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- fixed_inference_fixture(kernel, bin = TRUE)
    dat <- f$data
    dat$y <- if (family == "binomial") rep(0:1, 18) else rep(c(0, 1, 0, 9), 9)
    model <- if (family == "nbinom2") MASS::glm.nb(y ~ mean + rough + fixed + z, dat) else
      glm(y ~ mean + rough + fixed + z, dat, family = switch(family,
        gaussian = gaussian(), poisson = poisson(), binomial = binomial()))
    actual <- .msr_profile_contract(model, f$inputs)
    reference <- .msr_profile_contract(model, f$reference)
    expect_identical(actual$covariates, c("mean", "rough"))
    expect_identical(actual$fixed_columns, reference$fixed_columns)
    expect_identical(actual$objective$parameter_names, reference$objective$parameter_names)
    par <- rep(log(.7), length(actual$objective$parameter_names))
    if (family == "nbinom2") par[length(par)] <- log(3)
    expect_equal(actual$objective$fn(par), reference$objective$fn(par), tolerance = 0)
    expect_equal(actual$objective$gr(par), reference$objective$gr(par), tolerance = 0)
  }
})

test_that("fixed Bayesian columns preserve priors, target, derivatives, and repeated-row predictions", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- fixed_inference_fixture(kernel)
    dat <- f$data
    if (family == "binomial") dat$y <- rep(0:1, 18)
    form <- y ~ rough + fixed + z
    scales <- list(rough = ~1 + fixed)
    actual <- multiScale_bayes(form, f$inputs, dat, scales, family,
      row_id = "key", priors = list(scale_median = 25), method = "prepare")
    reference <- multiScale_bayes(form, f$reference, dat, scales, family,
      row_id = "key", priors = list(scale_median = 25), method = "prepare")
    a <- actual$specification; r <- reference$specification
    expect_identical(a$covariates, "rough")
    expect_identical(a$priors, r$priors)
    expect_identical(a$parameters, r$parameters)
    expect_equal(a$X_fixed, r$X_fixed, tolerance = 0)
    expect_equal(a$center, r$center, tolerance = 0)
    expect_equal(a$spread, r$spread, tolerance = 0)
    expect_equal(a$fixed_reference$center, f$inputs$scl_params$mean["fixed"], tolerance = 0)
    p <- a$parameters
    p$beta[] <- .2
    ad <- RTMB::MakeADFun(.msr_bayes_objective(a, TRUE), p, silent = TRUE)
    rd <- RTMB::MakeADFun(.msr_bayes_objective(r, TRUE), p, silent = TRUE)
    expect_equal(ad$fn(ad$par), rd$fn(rd$par), tolerance = 0)
    expect_equal(ad$gr(ad$par), rd$gr(rd$par), tolerance = 0)
    expect_equal(ad$he(ad$par), rd$he(rd$par), tolerance = 0)
    actual <- fake_bayesian_draws(actual)
    actual$draws$beta[,] <- .2
    actual$draws$scale_beta[, 2L] <- .15
    old <- predict(actual, summary = FALSE)
    ids <- c(8, 2, 2, 5)
    nd <- dat[ids, , drop = FALSE]
    nd$fixed <- NULL
    got <- predict(actual, nd, kernel_inputs = f$inputs, include_fixed = TRUE, summary = FALSE)
    expect_equal(as.numeric(got), as.numeric(old[, ids]), tolerance = 1e-12)
    expect_error(predict(actual, nd, kernel_inputs = f$inputs), "fixed")
    explicit <- dat[ids, , drop = FALSE]
    expect_equal(simulate(actual, 2, seed = 9, newdata = explicit, kernel_inputs = f$inputs),
      simulate(actual, 2, seed = 9, newdata = nd, kernel_inputs = f$inputs, include_fixed = TRUE))
  }
})

test_that("fixed preparation proof rejects altered data and definitions", {
  f <- fixed_inference_fixture()
  dat <- f$data
  model <- glm(y ~ rough + fixed + z, dat, family = poisson())
  bad <- model
  bad$model$fixed <- bad$model$fixed * 2
  expect_error(.msr_profile_contract(bad, f$inputs), "does not match")
  dat$fixed <- dat$fixed * 2
  expect_error(multiScale_bayes(y ~ rough + fixed + z, f$inputs, dat,
    family = "poisson", method = "prepare"), "does not match")
  reference <- .msr_fixed_reference_contract(f$inputs, f$data, seq_len(nrow(f$data)), "fixed")
  changed <- f$inputs
  changed$scale_vars$radius[changed$scale_vars$covariate == "fixed"] <- 26
  expect_error(.msr_fixed_reference_values(reference, changed, 1:3), "definitions")
  changed <- f$inputs
  changed$scale_vars <- NULL
  expect_error(.msr_fixed_reference_values(reference, changed, 1:3), "definitions")
  changed <- f$inputs
  changed$unit_conv <- changed$unit_conv * 2
  changed$d_list <- lapply(changed$d_list, function(x) x / 2)
  expect_equal(.msr_fixed_reference_values(reference, changed, c(3, 1, 1)),
    as.matrix(f$data[c(3, 1, 1), "fixed", drop = FALSE]), ignore_attr = TRUE, tolerance = 0)
  changed <- f$inputs; changed$max_D <- 15
  expect_error(.msr_fixed_reference_values(reference, changed, 1:3), "buffer extent")
})

test_that("fixed automatic prediction binds finite buffer support and preserves older metadata", {
  f <- fixed_inference_fixture()
  short <- kernel_prep(f$points, f$rasters, max_D = 15,
    scale_vars = f$inputs$scale_vars, bin = FALSE, verbose = FALSE)
  fit <- fake_bayesian_draws(multiScale_bayes(y ~ rough + fixed + z, f$inputs,
    f$data, family = "poisson", method = "prepare", priors = list(scale_median = 25)))
  fit$draws$beta[,] <- .2
  nd <- f$data; nd$fixed <- NULL
  expected <- predict(fit, nd, kernel_inputs = f$inputs, include_fixed = TRUE, summary = FALSE)
  older <- fit
  older$specification$fixed_reference$max_D <- NULL
  expect_equal(predict(older, nd, kernel_inputs = f$inputs, include_fixed = TRUE,
    summary = FALSE), expected, tolerance = 0)
  expect_error(predict(fit, nd, kernel_inputs = short, include_fixed = TRUE), "buffer extent")
  set.seed(446); rng <- .Random.seed
  expect_error(simulate(fit, nsim = 1, seed = 17, newdata = nd, kernel_inputs = short,
    include_fixed = TRUE), "buffer extent")
  expect_identical(.Random.seed, rng)
})

test_that("all finite varying fixed surface and landscape metrics enter numeric designs", {
  f <- fixed_inference_fixture()
  cover <- f$rasters[[1]]
  xy <- terra::xyFromCell(cover, seq_len(terra::ncell(cover)))
  terra::values(cover) <- 1L + (floor(xy[, 1] / 20) + floor(xy[, 2] / 20)) %% 2L
  patch <- (xy[, 1] - 145)^2 + (xy[, 2] - 145)^2 < 18^2
  cover[which(patch)] <- 3L
  names(cover) <- "cover"
  rasters <- c(f$rasters, cover)
  surface <- .msr_surface_metrics()
  landscape <- .msr_landscape_metrics()
  metrics <- c(paste0("surface_", surface), paste0("landscape_", landscape))
  status <- setNames(rep("pending", length(metrics)), metrics)
  for (id in metrics) {
    metric <- sub("^[^_]+_", "", id)
    fixed <- if (startsWith(id, "surface_")) surface_var("elev", metric, radius = 65) else
      landscape_var("cover", metric, radius = 65,
        classes_max = 3,
        class = if (metric %in% .msr_class_metrics()) 1 else NULL)
    prep <- tryCatch(kernel_prep(f$points, rasters, max_D = 90, sigma = 29,
      scale_vars = msr_vars(rough = surface_var("elev", "sq", weighted = TRUE),
        fixed = fixed), bin = FALSE, verbose = FALSE), error = identity)
    if (inherits(prep, "error")) {
      # Undefined or constant columns remain preparation failures, not fitted features.
      expect_match(conditionMessage(prep), "constant|finite|variance|variation|standard deviation|all NA")
      status[id] <- "invalid_preparation"
      next
    }
    dat <- prep$kernel_dat
    dat$y <- seq_len(nrow(dat)) / nrow(dat) + sin(seq_len(nrow(dat)))
    complete <- complete.cases(dat)
    if (any(!complete)) {
      expect_error(multiScale_bayes(y ~ rough + fixed, prep, dat, family = "gaussian",
        method = "prepare", priors = list(scale_median = 25)), "missing values")
      message(id, ": explicitly retained ", sum(complete), "/", nrow(dat), " finite rows")
      dat <- dat[complete, , drop = FALSE]
    }
    model <- glm(y ~ rough + fixed, data = dat, na.action = na.fail)
    reference <- prep
    reference$scale_vars <- reference$scale_vars[1L, , drop = FALSE]
    a <- .msr_profile_contract(model, prep)
    b <- .msr_profile_contract(model, reference)
    expect_equal(a$objective$fn(log(.7)), b$objective$fn(log(.7)), tolerance = 0, info = id)
    expect_equal(a$objective$gr(log(.7)), b$objective$gr(log(.7)), tolerance = 0, info = id)
    bayes <- multiScale_bayes(y ~ rough + fixed, prep, dat, family = "gaussian",
      method = "prepare", priors = list(scale_median = 25))
    expect_equal(bayes$specification$X_fixed[, "fixed"], dat$fixed,
      tolerance = 0, ignore_attr = TRUE, info = id)
    status[id] <- "verified_numeric_design"
  }
  expect_false(any(status == "pending"))
  expect_true(all(status[paste0("surface_", surface)] == "verified_numeric_design"))
  expect_gte(sum(status == "verified_numeric_design"), 20L)
  message(paste(names(status), status, sep = ": ", collapse = "\n"))
})

test_that("public analytic fixed features share point, snapshot, raster, and K contracts", {
  f <- fixed_inference_fixture()
  dat <- f$data
  set.seed(442)
  dat$y <- rpois(nrow(dat), exp(3 + 1.5 * dat$rough + .1 * dat$z))
  model <- glm(y ~ rough + fixed + z, dat, family = poisson())
  warnings <- character()
  fit <- withCallingHandlers(suppressMessages(multiScale_optim(model, f$inputs,
    backend = "profile_glm", verbose = FALSE)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  if (length(warnings)) {
    message(paste(warnings, collapse = "\n"))
    expect_true(all(grepl("Joint Wald uncertainty unavailable", warnings, fixed = TRUE)))
  }
  expect_identical(rownames(fit$scale_est), "rough")
  expect_equal(fit$scl_params$mean["fixed"], f$inputs$scl_params$mean["fixed"], tolerance = 0)
  expect_equal(.msr_parameter_count(fit$opt_mod), attr(logLik(fit$opt_mod), "df"))
  expect_equal(.msr_selection_parameter_count(fit),
    attr(logLik(fit$opt_mod), "df") + 1L)
  nd <- dat; nd$z <- .2; nd$fixed <- NULL
  point <- predict(fit, nd, f$inputs, include_fixed = TRUE, se.fit = TRUE,
    uncertainty = "conditional")
  cached <- prepare_prediction(fit, nd, f$inputs, include_fixed = TRUE)
  expect_equal(predict(cached, se.fit = TRUE, uncertainty = "conditional"), point)
  changed <- f$inputs; changed$max_D <- 15
  expect_error(predict(fit, nd, changed, include_fixed = TRUE), "buffer extent")
  expect_error(prepare_prediction(fit, nd, changed, include_fixed = TRUE), "buffer extent")
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = fit, scale_center = TRUE,
    include_fixed = TRUE, match_buffer = TRUE, site_values = list(z = .2), verbose = FALSE)
  map <- predict(fit, inputs, se.fit = TRUE, uncertainty = "conditional")
  got <- terra::extract(map, f$points, ID = FALSE)
  expect_equal(got$fit, as.numeric(point$fit), tolerance = 1e-8)
  expect_equal(got$se.fit, as.numeric(point$se.fit), tolerance = 1e-8)
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(fit, length.out = 5)
  expect_equal(plots$fixed$data$x,
    seq(min(dat$fixed), max(dat$fixed), length.out = 5) *
      f$inputs$scl_params$sd[["fixed"]] + f$inputs$scl_params$mean[["fixed"]],
    tolerance = 1e-12)
})

test_that("fixed mixed predictors retain all starts, cold targets, and native prediction", {
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  f <- fixed_inference_fixture()
  dat <- f$data
  dat$site <- factor(rep(1:6, each = 6))
  set.seed(443)
  effects <- c(-.8, .4, -.2, .9, -.5, .2)
  dat$y <- rpois(nrow(dat), exp(3 + 1.5 * dat$rough + .2 * dat$fixed +
    .1 * dat$z + effects[dat$site]))
  native <- glmmTMB::glmmTMB(y ~ rough + fixed + z + (1 | site), dat, family = poisson())
  fit <- suppressMessages(multiScale_optim(native, f$inputs, backend = "joint_mixed",
    verbose = FALSE, n_cores = 1,
    profile_control = list(starts = matrix(c(15, 25, 35) / f$inputs$unit_conv, ncol = 1))))
  expect_length(fit$mixed_fit$start_results, 3L)
  selected <- fit$mixed_fit$start_results[[fit$mixed_fit$selected]]
  expect_true(selected$eligible)
  expect_true(selected$design_unchanged)
  expect_lt(selected$native_delta, 1e-5)
  expect_identical(rownames(fit$scale_est), "rough")
  expect_equal(fit$fixed_reference$center, f$inputs$scl_params$mean["fixed"], tolerance = 0)
  nd <- dat; nd$fixed <- NULL
  got <- predict(fit, nd, f$inputs, include_fixed = TRUE, se.fit = TRUE,
    uncertainty = "conditional")
  expect_equal(got, predict(fit, dat, f$inputs, se.fit = TRUE, uncertainty = "conditional"))
  cached <- prepare_prediction(fit, nd, f$inputs, include_fixed = TRUE)
  expect_equal(predict(cached, se.fit = TRUE, uncertainty = "conditional"), got)
  expect_false(fit$uncertainty$available)
})

test_that("fixed conditional-logistic predictors share verified point and snapshot preparation", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- fixed_inference_fixture()
  dat <- f$data
  dat$choice <- rep(1:12, each = 3)
  set.seed(444)
  dat$case_ <- 0L
  dat$case_[(0:11) * 3L + sample.int(3L, 12, replace = TRUE)] <- 1L
  form <- as.formula("case_ ~ rough + fixed + z + strata(choice)",
    env = asNamespace("survival"))
  native <- survival::clogit(form, dat, model = TRUE, x = TRUE, y = TRUE)
  fit <- suppressMessages(multiScale_optim(native, f$inputs, backend = "conditional_logit",
    verbose = FALSE, n_cores = 1,
    profile_control = list(starts = matrix(c(15, 25, 35) / f$inputs$unit_conv, ncol = 1))))
  nd <- dat; nd$fixed <- NULL
  expected <- predict(fit, dat, f$inputs)
  expect_equal(predict(fit, nd, f$inputs, include_fixed = TRUE), expected)
  expect_equal(predict(prepare_prediction(fit, nd, f$inputs, include_fixed = TRUE)), expected)
})

test_that("Bayesian fixed-feature curves retain original preparation units", {
  f <- fixed_inference_fixture()
  fit <- fake_bayesian_draws(multiScale_bayes(y ~ rough + fixed + z, f$inputs,
    f$data, family = "gaussian", method = "prepare", priors = list(scale_median = 25)))
  fit$draws$beta[,] <- .2
  open_null_device()
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(fit, length.out = 5)
  expect_equal(plots$fixed$data$x,
    seq(min(f$data$fixed), max(f$data$fixed), length.out = 5) *
      f$inputs$scl_params$sd[["fixed"]] + f$inputs$scl_params$mean[["fixed"]],
    tolerance = 1e-12)
  expect_match(plots$fixed$labels$caption, "Hypothetical|reference")
})

test_that("fixed NB2 predictors share analytic point and snapshot information", {
  skip_if_not_installed("MASS")
  f <- fixed_inference_fixture()
  dat <- f$data
  set.seed(445)
  dat$y <- rnbinom(nrow(dat), mu = exp(3 + 1.5 * dat$rough + .1 * dat$z), size = 4)
  native <- MASS::glm.nb(y ~ rough + fixed + z, data = dat)
  warnings <- character()
  fit <- withCallingHandlers(suppressMessages(multiScale_optim(native, f$inputs,
    backend = "joint_nb2", verbose = FALSE, n_cores = 1)), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
  if (length(warnings)) {
    message(paste(warnings, collapse = "\n"))
    expect_true(all(grepl("Joint Wald uncertainty unavailable", warnings, fixed = TRUE)))
  }
  expect_identical(rownames(fit$scale_est), "rough")
  nd <- dat; nd$fixed <- NULL
  expected <- predict(fit, dat, f$inputs, se.fit = TRUE, uncertainty = "conditional")
  expect_equal(predict(fit, nd, f$inputs, include_fixed = TRUE,
    se.fit = TRUE, uncertainty = "conditional"), expected)
  expect_equal(predict(prepare_prediction(fit, nd, f$inputs, include_fixed = TRUE),
    se.fit = TRUE, uncertainty = "conditional"), expected)
})

# Small native conformance fixtures, not a scale-recovery or timing panel.
test_that("raw native refits preserve metric units, factors, SEs and spatial prediction", {
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    f <- refit_spatial_fixture(kernel)
    dat <- f$data
    dat[c("mean", "rough")] <- as.data.frame(f$raw)
    initial <- stats::glm(y ~ mean*group + rough + fixed + z + offset(log(exposure)),
      data = dat, family = stats::gaussian(), contrasts = list(group = stats::contr.sum(3)))
    context <- build_opt_context(initial, f$inputs$raw_cov, scale_vars = f$inputs$scale_vars,
      unit_conv = f$inputs$unit_conv, resolution = f$inputs$resolution,
      n_cols = f$inputs$n_cols, binned = f$inputs$binned, d_list = f$inputs$d_list)
    par <- f$object$optim_results$par
    legacy <- kernel_scale_fn(par, f$inputs$d_list, f$inputs$raw_cov, kernel, initial,
      opt_context = context, mod_return = TRUE)
    expect_equal(as.numeric(unlist(legacy$mod$model[c("mean", "rough")])),
      as.numeric(scale(f$raw)), tolerance = 1e-12)
    context$covariate_transform <- "raw"
    final <- kernel_scale_fn(par, f$inputs$d_list, f$inputs$raw_cov, kernel, initial,
      opt_context = context, mod_return = TRUE)
    direct <- stats::update(initial, data = dat)
    expect_equal(as.numeric(unlist(final$mod$model[c("mean", "rough")])),
      as.numeric(f$raw), tolerance = 0)
    expect_identical(rownames(final$mod$model), rownames(dat))
    expect_identical(final$scl_params$mean, c(mean = 0, rough = 0))
    expect_identical(final$scl_params$sd, c(mean = 1, rough = 1))
    expect_equal(unname(stats::coef(final$mod)), unname(stats::coef(direct)), tolerance = 1e-10)
    expect_equal(unname(stats::vcov(final$mod)), unname(stats::vcov(direct)), tolerance = 1e-10)
    expect_identical(final$mod$contrasts, direct$contrasts)
    expect_equal(.msr_parameter_count(final$mod), .msr_parameter_count(direct))
    expect_equal(kernel_scale_fn(par, f$inputs$d_list, f$inputs$raw_cov, kernel, initial,
      opt_context = context), -as.numeric(stats::logLik(direct)), tolerance = 1e-10)
    object <- f$object
    object$opt_mod <- final$mod
    object$fitted_mod_original <- initial
    object$opt_context <- context
    object$scl_params <- .msr_merge_scl_params(final$scl_params, f$inputs$scl_params)
    expect_equal(object$scl_params$mean['fixed'], f$inputs$scl_params$mean['fixed'])
    expect_equal(object$scl_params$sd['fixed'], f$inputs$scl_params$sd['fixed'])
    rows <- c(19, 3, 28, 9)
    nd <- dat[rows, , drop = FALSE]
    expected <- stats::predict(direct, nd, type = "response", se.fit = TRUE)
    nd[c("mean", "rough", "fixed")] <- -999
    actual <- stats::predict(object, nd, f$inputs, include_fixed = TRUE,
      se.fit = TRUE, uncertainty = "conditional")
    expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-10)
    expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-10)
    expect_identical(object$scl_params, .msr_merge_scl_params(final$scl_params, f$inputs$scl_params))
    if (kernel == "gaussian") {
      maps <- kernel_scale.raster(f$rasters, multiScaleR = object, scale_center = TRUE,
        include_fixed = TRUE, match_buffer = TRUE, verbose = FALSE,
        site_values = list(group = "b", z = .2), offset_values = c(exposure = 2))
      scenario <- dat
      scenario$group <- factor("b", levels = levels(dat$group))
      scenario$z <- .2; scenario$exposure <- 2
      reference <- stats::predict(direct, scenario, type = "response", se.fit = TRUE)
      mapped <- stats::predict(object, maps, se.fit = TRUE, uncertainty = "conditional")
      extracted <- terra::extract(mapped, f$points, ID = FALSE)
      expect_equal(extracted$fit, as.numeric(reference$fit), tolerance = 1e-8)
      expect_equal(extracted$se.fit, as.numeric(reference$se.fit), tolerance = 1e-8)
      expect_identical(attr(mapped, "prediction_contract")$uncertainty$kernel_parameters,
        "fixed_at_estimate")
      legacy_object <- object; legacy_object$opt_context$covariate_transform <- NULL
      expect_false(identical(profile_scale_cache_key(object, object$min_D, rownames(object$scale_est)),
        profile_scale_cache_key(legacy_object, object$min_D, rownames(object$scale_est))))
      open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
      curves <- plot_marginal_effects(object, length.out = 7, reference_data = dat)
      expect_equal(range(curves$mean$data$x), range(dat$mean), tolerance = 1e-12)
      expect_equal(range(curves$rough$data$x), range(dat$rough), tolerance = 1e-12)
    }
  }
})

test_that("fixed raw spline knots and hard radii keep native contracts", {
  f <- refit_spatial_fixture()
  dat <- f$data; dat[c("mean", "rough")] <- as.data.frame(f$raw)
  model <- stats::glm(y ~ splines::ns(mean, knots = c(0, 1), Boundary.knots = c(-3, 5))*group +
    rough + fixed + z + offset(log(exposure)), data = dat, family = stats::gaussian())
  context <- build_opt_context(model, f$inputs$raw_cov, scale_vars = f$inputs$scale_vars,
    unit_conv = f$inputs$unit_conv, resolution = f$inputs$resolution,
    n_cols = f$inputs$n_cols, binned = f$inputs$binned, d_list = f$inputs$d_list)
  context$covariate_transform <- "raw"
  final <- kernel_scale_fn(f$object$optim_results$par, f$inputs$d_list, f$inputs$raw_cov,
    f$inputs$kernel, model, opt_context = context, mod_return = TRUE)
  native <- stats::update(model, data = dat)
  expect_equal(stats::model.matrix(final$mod), stats::model.matrix(native), tolerance = 0)
  expect_equal(unname(stats::vcov(final$mod)), unname(stats::vcov(native)), tolerance = 1e-10)
  for (backend in c("profile_glm", "joint_nb2", "joint_mixed", "conditional_logit"))
    expect_error(multiScale_optim(NULL, NULL, backend = backend, covariate_transform = "raw"),
      "require backend = 'refit'")
  f <- refit_spatial_fixture("fixed")
  dat <- f$data; dat[c("mean", "rough")] <- as.data.frame(f$raw)
  model <- stats::glm(y ~ mean + rough + fixed + group, data = dat, family = stats::gaussian())
  fitted <- multiScale_optim(model, f$inputs, par = c(.475, .725), n_cores = 1L,
    verbose = FALSE, covariate_transform = "raw")
  expect_s3_class(fitted, "multiScaleR")
  expect_identical(fitted$opt_context$covariate_transform, "raw")
  expect_true(all(is.na(fitted$scale_est$SE)))
  expect_true(all(is.finite(diag(stats::vcov(fitted$opt_mod)))))
  expect_identical(fitted$scl_params$mean[c("mean", "rough")], c(mean = 0, rough = 0))
  expect_equal(.msr_selection_parameter_count(fitted, 0L), 2L)
})

test_that("raw native C GDR callbacks preserve observations, controls and final SEs", {
  skip_if_not_installed("unmarked")
  skip_if(as.character(packageVersion("unmarked")) != "1.5.2", "Native GDR provider needs qualification")
  f <- gdr_refit_fixture("P")
  frame <- f$model@data
  n <- nrow(frame@siteCovs)
  distances <- rep(list(c(.1, .4, .9)), n)
  cells <- lapply(frame@siteCovs$hab, function(x) matrix(x + c(-.2, .1, .4), ncol = 1,
    dimnames = list(NULL, "hab")))
  names(cells) <- names(distances) <- rownames(frame@siteCovs)
  callback <- make_gdr_refitter(f$model, f$starts, gdr_test_audit_path("raw-native-C-"))
  context <- build_opt_context(f$model, cells, refit_fn = callback)
  context$covariate_transform <- "raw"
  w <- exp(-distances[[1]]^2/(2*.3^2))
  raw <- vapply(cells, function(x) sum(w*x[, 1])/sum(w), 0)
  result <- kernel_scale_fn(.3, distances, cells, "gaussian", f$model,
    opt_context = context, mod_return = TRUE)
  expect_equal(unname(result$mod@data@siteCovs$hab), unname(raw), tolerance = 1e-12)
  expect_identical(result$mod@data@yDistance, frame@yDistance)
  expect_identical(result$mod@data@yRemoval, frame@yRemoval)
  expect_identical(result$mod@K, f$model@K)
  expect_equal(.msr_parameter_count(result$mod), .msr_parameter_count(f$model))
  expect_identical(.msr_gdr_variance_status(result$mod, TRUE)$status, "available")
  expect_identical(result$scl_params, list(mean = c(hab = 0), sd = c(hab = 1)))
})

test_that("raw spatial CV fits and scores all supported GLM families in training units", {
  # Seven bounded fold fits: six raw routes and one existing candidate route.
  candidates <- list()
  cases <- list(gaussian = list(family = "gaussian"), poisson = list(family = "poisson"),
    binomial = list(family = "binomial"), nb = list(family = "nbinom2", theta = 2),
    gamma = list(family = "Gamma", dispersion = .4), sq = list(family = "gaussian"))
  for (label in names(cases)) {
    case <- cases[[label]]
    f <- cv_formula_spatial_inputs(label == "sq")
    f$data$y <- f$data[[switch(case$family, gaussian = "continuous", poisson = "count",
      binomial = "binary", nbinom2 = "count", Gamma = "positive")]]
    held <- f$folds[[1]]$assessment; analysis <- f$folds[[1]]$analysis
    f$data[held, "z"] <- f$data[held, "z"] + 3
    form <- y ~ hab*group + I(base::scale(z)^2) + offset(base::scale(z)) + offset(log(exposure))
    cv <- do.call(cross_validate, c(list(formula = form, data = f$data, kernel_inputs = f$inputs,
      folds = f$folds[1], backend = "refit", keep_fits = TRUE, keep_selection = TRUE,
      optim_args = list(par = .3, covariate_transform = "raw")), case))
    expect_true(cv$complete, info = paste(label, paste(unlist(cv$errors), collapse = "; ")))
    if (!cv$complete) next
    fit <- cv$fits[[1]]
    expect_identical(fit$scl_params$mean['hab'], c(hab = 0))
    expect_identical(fit$scl_params$sd['hab'], c(hab = 1))
    expect_identical(cv$training[[1]]$initial_normalization$mean['hab'], c(hab = 0))
    expect_identical(cv$training[[1]]$initial_normalization$sd['hab'], c(hab = 1))
    sigma <- fit$optim_results$par[1]
    reference <- cv$training[[1]]$bin_reference
    raw <- vapply(rownames(f$data), function(id) {
      d <- f$inputs$d_list[[id]]; v <- f$inputs$raw_cov[[id]][, 1]
      if (!is.null(reference)) d <- reference$bin_rep[.bincode(d, reference$breaks, include.lowest = TRUE)]
      w <- exp(-d^2/(2*sigma^2)); mean <- sum(w*v)/sum(w)
      if (label == "sq") sqrt(sum(w*(v - mean)^2)/sum(w)) else mean
    }, 0)
    dat <- f$data; dat$hab <- raw
    dat$scaled_z <- (dat$z - mean(dat[analysis, "z"]))/sd(dat[analysis, "z"])
    dat$squared_z <- dat$scaled_z^2
    # Keep the native family's SE convention. A fixed NB size does not change
    # stats::glm's separate Pearson-dispersion covariance into unit dispersion.
    family <- fit$opt_mod$family
    direct <- stats::glm(y ~ hab*group + squared_z + offset(scaled_z) + offset(log(exposure)),
      data = dat[analysis, ], family = family)
    expect_equal(cv$predictions$prediction,
      as.numeric(stats::predict(direct, dat[held, ], type = "response")), tolerance = 1e-7, info = label)
    expect_equal(unname(stats::vcov(fit$opt_mod)), unname(stats::vcov(direct)), tolerance = 1e-7, info = label)
    expect_equal(cv$selection_records[[1]]$data$ecological_K, 1L)
    expect_equal(cv$selection_records[[1]]$data$K, cv$selection_records[[1]]$data$native_K + 1L)
    expect_identical(cv$settings$parameter_uncertainty, "not_integrated")
    compact <- cv; compact$fits <- NULL
    selected <- select_cv(list(raw = compact, copy = compact), criterion = "AIC", ties = "first")
    expect_true(selected$complete)
    expect_true(all(selected$selection$chosen == "raw"))
    if (label == "gaussian") {
      candidates$raw <- cv
      candidates$candidate <- cross_validate(formula = form, data = f$data, folds = f$folds[1],
        family = "gaussian", kernel_inputs = f$inputs,
        backend = "refit", optim_args = list(par = .3), keep_fits = TRUE, keep_selection = TRUE)
    }
  }
  expect_true(candidates$candidate$complete)
  expect_equal(nrow(cv_tab(candidates, "rmse")), 2L)
  expect_true(select_cv(candidates, criterion = "AIC", ties = "first")$complete)
  expect_error(.msr_cv_covariate_transform(list(covariate_transform = "raw"), "profile_glm"),
    "require spatial CV")
  cat("\nRAW_REFIT_CONTRACT kernels=4 native_C_GDR=1 GLM_CV_families=5 guards=retained\n")
})

test_that("matched-set native CV preserves raw units and conditional-choice scoring", {
  skip_if_not_installed("survival")
  strata <- survival::strata; coxph <- survival::coxph
  f <- cv_formula_choice_inputs(); d <- f$data
  form <- case_ ~ hab + group + I(base::scale(z)^2) + offset(base::scale(z)) + strata(choice)
  cv <- cross_validate_clogit(form, d, f$folds, f$inputs, backend = "refit",
    optim_args = list(par = .25, covariate_transform = "raw"), keep_fits = TRUE, keep_selection = TRUE)
  expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = "; "))
  if (!cv$complete) return(invisible(NULL))
  fit <- cv$fits[[1]]; analysis <- cv$folds[[1]]$analysis; held <- cv$folds[[1]]$assessment
  reference <- cv$training[[1]]$bin_reference; sigma <- fit$optim_results$par[1]
  d$hab <- vapply(rownames(d), function(id) {
    distance <- f$inputs$d_list[[id]]; v <- f$inputs$raw_cov[[id]][, 1]
    if (!is.null(reference)) distance <- reference$bin_rep[.bincode(distance, reference$breaks, include.lowest = TRUE)]
    w <- exp(-distance^2/(2*sigma^2)); sum(w*v)/sum(w)
  }, 0)
  d$scaled_z <- (d$z - mean(d[analysis, "z"]))/sd(d[analysis, "z"])
  d$squared_z <- d$scaled_z^2
  native <- survival::clogit(case_ ~ hab + group + squared_z + offset(scaled_z) + strata(choice),
    data = d[analysis, ], method = "exact", model = TRUE)
  assessment <- d[held, ]; assessment$choice <- d[analysis[1], "choice"]
  eta <- as.numeric(stats::predict(native, assessment, type = "lp", reference = "zero"))
  probability <- numeric(length(eta))
  for (rows in split(seq_along(eta), cv$predictions$choice))
    probability[rows] <- exp(eta[rows] - max(eta[rows]))/sum(exp(eta[rows] - max(eta[rows])))
  expect_equal(cv$predictions$prediction, probability, tolerance = 1e-7)
  expect_equal(unname(stats::vcov(fit$opt_mod)), unname(stats::vcov(native)), tolerance = 1e-7)
  expect_identical(fit$scl_params$mean['hab'], c(hab = 0))
  expect_identical(fit$scl_params$sd['hab'], c(hab = 1))
  expect_equal(cv$selection_records[[1]]$data$ecological_K, 1L)
  expect_identical(cv$settings$prediction_target, "within_set_choice_probability")
})

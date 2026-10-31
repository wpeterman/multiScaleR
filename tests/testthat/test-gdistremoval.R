test_that("GDR defaults and snapshots retain native component predictions", {
  f <- gdr_fixture()
  expect_equal(as.numeric(predict(f$object)),
    unmarked::predict(f$object$opt_mod, "lambda")$Predicted)
  for (type in c("lambda", "phi", "dist", "rem")) {
    native <- unmarked::predict(f$object$opt_mod, type)
    prepared <- prepare_prediction(f$object, component = type)
    expect_equal(predict(prepared), native)
    expect_equal(predict(unserialize(serialize(prepared, NULL))), native)
  }
  expect_error(predict(f$object, type = "state"), "active native component")
  expect_error(prepare_prediction(f$object, component = "det"), "active native component")
})

test_that("inactive GDR formulas cannot introduce spurious scale parameters", {
  model <- gdr_fixture()$model
  model@data@numPrimary <- 1
  model@keyfun <- "uniform"
  model@formlist$phi <- ~unused_phi
  model@formlist$dist <- ~unused_distance
  expect_setequal(.unmarked_model_predictors(model), c("hab", "tod"))
  expect_error(.msr_gdr_prediction_type(model, "phi"), "active native component")
  expect_error(.msr_gdr_prediction_type(model, "dist"), "active native component")
})

test_that("GDR selection compares both histories and their observation definitions", {
  model <- gdr_fixture()$model
  # A distance-bin swap preserves every survey total but changes its likelihood.
  changed <- model
  row <- which(changed@data@yDistance[, 1] != changed@data@yDistance[, 2])[1]
  changed@data@yDistance[row, 1:2] <- rev(changed@data@yDistance[row, 1:2])
  expect_error(.msr_assert_selection_targets(list(model, changed)), "response")
  for (slot in c("dist.breaks", "unitsIn", "period.lengths", "numPrimary")) {
    changed <- model
    methods::slot(changed@data, slot) <- switch(slot,
      dist.breaks = c(0, 30, 60, 100), unitsIn = "km",
      period.lengths = c(2, 2, 2), numPrimary = 1)
    expect_error(.msr_assert_selection_targets(list(model, changed)), "observation definitions")
  }
  expect_no_error(.msr_assert_selection_targets(list(model, model)))
  # Different hypotheses may change the mixture or key, not the observations.
  changed <- model
  changed@keyfun <- "exp"
  expect_no_error(.msr_assert_selection_targets(list(model, changed)))
  changed <- model
  changed@data <- model@data[nrow(model@data@yDistance):1, ]
  # S4 subsetting in some unmarked versions preserves old site row names.
  expect_no_error(.msr_assert_selection_targets(list(model, changed)))
})

test_that("GDR site predictors cannot be shadowed by visit-level covariates", {
  f <- gdr_fixture()
  model <- f$model
  model@formlist$dist <- ~hab
  model@data@yearlySiteCovs$hab <- rep(model@data@siteCovs$hab, each = 2)
  expect_error(build_opt_context(model, f$inputs$raw_cov), "shadowed")
})

test_that("one raster source can have separate scales in all GDR components", {
  f <- gdr_fixture()
  vars <- msr_vars(ha = kernel_var("hab"), hp = kernel_var("hab"),
    hd = kernel_var("hab"), hr = kernel_var("hab"))
  inputs <- kernel_prep(f$points, f$raster, max_D = 80, sigma = rep(30, 4),
    scale_vars = vars, bin = FALSE, verbose = FALSE)
  frame <- f$frame
  for (nm in vars$covariate) frame@siteCovs[[nm]] <- frame@siteCovs$hab
  model <- unmarked::update(f$model, data = frame, lambdaformula = ~ha,
    phiformula = ~hp, distanceformula = ~hd + visit, removalformula = ~hr + tod)
  context <- build_opt_context(model, inputs$raw_cov, scale_vars = vars,
    unit_conv = inputs$unit_conv)
  expect_identical(context$covs, c("ha", "hp", "hd", "hr"))
  evaluated <- kernel_scale_fn(c(22, 30, 45, 60) / inputs$unit_conv,
    inputs$d_list, inputs$raw_cov, inputs$kernel, model,
    opt_context = context, mod_return = TRUE)
  values <- evaluated$mod@data@siteCovs[vars$covariate]
  expect_equal(nrow(unique(t(as.matrix(values)))), 4L)
  expect_identical(evaluated$mod@data@yDistance, frame@yDistance)
  expect_identical(evaluated$mod@data@yRemoval, frame@yRemoval)
  object <- f$object
  object$opt_mod <- evaluated$mod
  object$kernel_inputs <- inputs
  object$scl_params <- evaluated$scl_params
  object$scale_est <- data.frame(Mean = c(22, 30, 45, 60), SE = NA_real_,
    row.names = c("ha", "hp", "hd", "hr"))
  scaled <- kernel_scale.raster(f$raster, multiScaleR = object,
    scale_vars = vars[1, , drop = FALSE], scale_center = TRUE,
    match_buffer = TRUE, verbose = FALSE)
  expect_identical(names(scaled), "ha")
  projected <- predict(object, scaled, type = "lambda")
  native <- unmarked::predict(evaluated$mod, "lambda",
    newdata = as.data.frame(scaled, na.rm = FALSE))$Predicted
  expect_equal(as.numeric(terra::values(projected)), native)
  expect_error(predict(object, scaled, type = "rem"), "projected layer")
})

test_that("combined GDR means agree with independent analytical probabilities", {
  model <- gdr_oracle()
  out <- predict_gdistremoval(model)
  sigma <- rep(c(50, 100), nrow(model@data@yDistance))
  pd <- 2 * sigma^2 / 100^2 * (1 - exp(-100^2 / (2 * sigma^2)))
  pa <- 1 - (1 - 0.3)^5
  expect_equal(out$phi, rep(0.7, nrow(out)))
  expect_equal(out$availability, rep(pa, nrow(out)))
  expect_equal(out$perceptibility, pd, tolerance = 1e-12)
  expect_equal(out$abundance, rep(8, nrow(out)))
  expect_equal(out$density, rep(8 / pi, nrow(out)))
  expect_equal(out$detection, 0.7 * pa * pd)
  expect_equal(out$expected_count, 8 * 0.7 * pa * pd)
  expect_identical(attr(out, "density_units"), "ha")
  expect_equal(attr(out, "survey_area"), pi)
  expect_identical(attr(out, "prediction_contract")$uncertainty$intervals, "not_provided")
  expect_gt(out$perceptibility[2], out$perceptibility[1])
})

test_that("GDR density units and ZIP mean use their actual model definitions", {
  model <- gdr_oracle("ZIP", "density")
  out <- predict_gdistremoval(model)
  expect_equal(out$lambda, rep(8, nrow(out)))
  expect_equal(out$zero_probability, rep(0.4, nrow(out)))
  expect_equal(out$density, rep(8 * 0.6, nrow(out)))
  expect_equal(out$abundance, rep(8 * 0.6 * pi, nrow(out)))
  model@unitsOut <- "kmsq"
  out <- predict_gdistremoval(model)
  expect_equal(out$abundance, rep(8 * 0.6 * pi / 100, nrow(out)))
  km <- model
  km@data@dist.breaks <- km@data@dist.breaks / 1000
  km@data@unitsIn <- "km"
  km@estimates@estimates$dist@estimates[1] <- log(0.05)
  converted <- predict_gdistremoval(km)
  expect_equal(converted$perceptibility, out$perceptibility)
  expect_equal(converted$abundance, out$abundance)
})

test_that("distance keys and zero removal probabilities remain valid", {
  for (key in c("uniform", "exp", "hazard")) {
    model <- gdr_oracle(keyfun = key)
    out <- predict_gdistremoval(model)
    expected <- if (key == "uniform") 1 else vapply(rep(c(50, 100), 24), function(scale) {
      stats::integrate(function(r) 2 * r / 100^2 * if (key == "exp")
        exp(-r / scale) else (1 - exp(-(r / scale)^-2)), 0, 100)$value
    }, numeric(1))
    expect_equal(out$perceptibility, rep(expected, length.out = nrow(out)), tolerance = 1e-7)
  }
  p <- rbind(c(0, 0.5, 0), c(1, 0, 1), c(0, 0, 0))
  expect_equal(.msr_gdr_removal(p, c(1, 2, 2)), rbind(c(0, 0.75, 0), c(1, 0, 0), c(0, 0, 0)))
  expect_equal(.msr_gdr_distance("halfnorm", c(0, Inf, NA), c(0, 100), 3),
    matrix(c(0, 1, NA), ncol = 1))
})

test_that("new GDR protocols preserve native rows, fitted normalization and inputs", {
  f <- gdr_fixture()
  frame <- f$object$opt_mod@data
  original <- frame
  frame@period.lengths <- c(2, 2, 3)
  keys <- data.frame(site = seq_len(nrow(frame@yDistance)))
  manual <- frame
  manual@siteCovs <- prepare_unmarked_data(f$object, frame@siteCovs, f$inputs, join_by = keys)
  got <- predict_gdistremoval(f$object, frame, f$inputs, join_by = keys)
  expected <- predict_gdistremoval(f$object$opt_mod, manual)
  expect_equal(as.matrix(got), as.matrix(expected))
  expect_identical(attr(got, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
  expect_equal(frame@siteCovs, original@siteCovs)
  expect_equal(frame@period.lengths, c(2, 2, 3))
  expect_gt(got$availability[1], predict_gdistremoval(f$object)$availability[1])
  expect_error(predict_gdistremoval(f$object, frame, f$inputs), "join_by")
  expect_error(predict_gdistremoval(f$model, newdata = data.frame(x = 1)), "unmarkedFrameGDR")
  changed <- frame
  changed@dist.breaks <- c(0, 20, 60, 100)
  expect_error(predict_gdistremoval(f$model, changed), "survey footprint")
  changed <- frame
  changed@siteCovs$hab[2] <- NA_real_
  predicted <- predict_gdistremoval(f$model, changed)
  expect_true(all(is.na(predicted$abundance[3:4])))
  expect_equal(predicted$perceptibility[-(3:4)], got$perceptibility[-(3:4)])
})

test_that("GDR scale optimization and checked raster prediction complete", {
  f <- gdr_fixture()
  fit <- suppressMessages(suppressWarnings(multiScale_optim(f$model, f$inputs,
    join_by = data.frame(site = seq_len(nrow(f$frame@yDistance))),
    par = 30 / f$inputs$unit_conv, verbose = FALSE)))
  expect_s4_class(fit$opt_mod, "unmarkedFitGDR")
  expect_identical(fit$opt_mod@data@yDistance, f$frame@yDistance)
  expect_identical(fit$opt_mod@data@yRemoval, f$frame@yRemoval)
  expect_identical(fit$opt_mod@data@obsCovs, f$frame@obsCovs)
  expect_identical(fit$opt_mod@data@yearlySiteCovs, f$frame@yearlySiteCovs)
  expect_identical(fit$opt_mod@formlist, f$model@formlist)
  expect_true(is.finite(fit$scale_est$Mean))
  scaled <- kernel_scale.raster(f$raster, multiScaleR = fit,
    scale_center = TRUE, match_buffer = TRUE, verbose = FALSE)
  mapped <- predict(fit, scaled, type = "lambda")
  native <- unmarked::predict(fit$opt_mod, "lambda",
    newdata = terra::as.data.frame(scaled, na.rm = FALSE))$Predicted
  expect_equal(as.numeric(terra::values(mapped)), native)
  expect_error(predict(fit, scaled, type = "state"), "active native component")
  expect_identical(attr(mapped, "prediction_contract")$native_component$units,
    "individuals_per_survey_circle")
})

test_that("single-visit and uniform GDR fits need no inactive component", {
  f <- gdr_fixture()
  frame <- f$frame[, 1]
  model <- unmarked::gdistremoval(~hab, ~1, ~tod, ~1, data = frame,
    keyfun = "uniform", K = 35, se = FALSE)
  out <- predict_gdistremoval(model)
  expect_equal(out$phi, rep(1, nrow(out)))
  expect_equal(out$perceptibility, rep(1, nrow(out)))
  expect_identical(.msr_gdr_components(model), c("lambda", "rem"))
})

test_that("GDR native random effects preserve explicit prediction policies", {
  skip_if_not_installed("TMB")
  f <- gdr_fixture()
  frame <- f$frame
  frame@siteCovs$route <- factor(rep(seq_len(6), each = 4))
  model <- suppressWarnings(unmarked::update(f$model, data = frame,
    lambdaformula = ~hab + (1 | route), engine = "TMB", se = FALSE))
  expect_false(is.null(model@TMB))
  refitted <- refit_gdistremoval(model, frame, list(return_model = FALSE))
  expect_false(is.null(refitted@TMB))
  expect_equal(unmarked::coef(refitted), unmarked::coef(model), tolerance = 1e-8)
  expect_equal(as.numeric(unmarked::logLik(refitted)),
    as.numeric(unmarked::logLik(model)), tolerance = 1e-8)
  for (policy in list(NULL, NA)) {
    expected <- unmarked::predict(model, "lambda", newdata = frame@siteCovs,
      re.form = policy, level = NULL)$Predicted
    got <- predict_gdistremoval(model, re.form = policy)
    expect_equal(got$lambda, rep(expected, each = frame@numPrimary))
  }
  # A zero group deviation is an explicit prediction scenario for new groups.
  unseen <- frame
  unseen@siteCovs$route <- factor(rep("unsampled", nrow(frame@yDistance)))
  expected <- unmarked::predict(model, "lambda", newdata = unseen@siteCovs,
    re.form = NA, level = NULL)$Predicted
  got <- predict_gdistremoval(model, unseen, re.form = NA)
  expect_equal(got$lambda, rep(expected, each = frame@numPrimary))
  expect_identical(attr(got, "prediction_contract")$effects$response,
    "zero_deviation_reference")
})

test_that("GDR simulation reproduces paired totals and conditional means", {
  model <- gdr_oracle("ZIP")
  # Four sites suffice for checking Monte Carlo moments without a slow refit.
  model@data <- model@data[1:4, ]
  set.seed(81)
  simulated <- simulate_gdistremoval(model, nsim = 600)
  set.seed(81)
  expect_identical(simulate_gdistremoval(model, nsim = 2),
    structure(simulated[1:2], prediction_contract = attr(simulated, "prediction_contract")))
  expect_true(all(vapply(simulated, function(frame) {
    totals <- function(x) cbind(rowSums(x[, 1:3, drop = FALSE]), rowSums(x[, 4:6, drop = FALSE]))
    identical(totals(frame@yDistance), totals(frame@yRemoval)) && identical(frame@y, frame@yRemoval)
  }, logical(1))))
  empirical <- Reduce(`+`, lapply(simulated, function(frame) cbind(
    rowSums(frame@yRemoval[, 1:3, drop = FALSE]), rowSums(frame@yRemoval[, 4:6, drop = FALSE])))) / length(simulated)
  expected <- matrix(predict_gdistremoval(model)$expected_count, ncol = 2, byrow = TRUE)
  expect_lt(max(abs(empirical - expected)), 0.45)
  expect_true(attr(simulated, "prediction_contract")$uncertainty$observation_noise)
  expect_error(simulate_gdistremoval(model, nsim = 0), "nsim")
  missing <- model@data
  missing@yDistance[2, 4] <- NA
  simulated <- simulate_gdistremoval(model, newdata = missing)[[1]]
  expect_true(all(is.na(simulated@yDistance[2, 4:6])))
  expect_true(all(is.na(simulated@yRemoval[2, 4:6])))
})

test_that("GDR NB simulations retain overdispersion and shared site abundance", {
  model <- gdr_oracle("NB")
  model@data <- model@data[1, ]
  expected <- predict_gdistremoval(model)
  set.seed(482)
  simulated <- simulate_gdistremoval(model, nsim = 1500)
  counts <- t(vapply(simulated, function(frame)
    c(sum(frame@yRemoval[1, 1:3]), sum(frame@yRemoval[1, 4:6])), numeric(2)))
  mean <- expected$expected_count
  theoretical_variance <- mean + mean^2 / 2
  theoretical_covariance <- expected$detection[1] * expected$detection[2] * (8 + 8^2 / 2)
  expect_lt(max(abs(colMeans(counts) - mean)), 0.35)
  expect_lt(max(abs(apply(counts, 2, var) / theoretical_variance - 1)), 0.25)
  expect_lt(abs(cov(counts)[1, 2] / theoretical_covariance - 1), 0.3)
})

test_that("GDR effect curves use the correct component covariate layout", {
  f <- gdr_fixture()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (type in c("lambda", "dist", "rem")) {
    curves <- plot_marginal_effects(f$object, type = type, length.out = 5)
    variable <- switch(type, lambda = "hab", dist = "visit", rem = "tod")
    expect_named(curves, variable)
    nd <- stats::setNames(data.frame(curves[[variable]]$data$x), variable)
    if (variable == "hab") nd$hab <- (nd$hab - f$object$scl_params$mean["hab"]) /
      f$object$scl_params$sd["hab"]
    expected <- unmarked::predict(f$object$opt_mod, type, newdata = nd, re.form = NA)
    expect_equal(curves[[variable]]$data$fit, expected$Predicted)
  }
  expect_named(plot_marginal_effects(f$object, length.out = 5), "hab")
})

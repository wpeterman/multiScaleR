test_that("selected lme4 effects share native point, snapshot and raster rows", {
  skip_if_not_installed("lme4")
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(791)
  dat$other <- factor(rep(1:4, 9))
  dat$only_other <- 1.5 + cos(seq_len(nrow(dat)))
  dat$y <- dat$mean + dat$rough + dat$z + 2 * as.integer(dat$site) +
    as.integer(dat$other) * dat$only_other + rnorm(nrow(dat), sd = .2)
  model <- lme4::lmer(y ~ mean + rough + fixed + z + (1 | site) +
    (0 + only_other | other), data = dat)
  object <- f$object
  object$opt_mod <- model
  nd <- dat[c(3, 21, 14), ]
  nd[c("other", "only_other")] <- NULL
  control <- stats::as.formula("~(1|site)", env = new.env(parent = baseenv()))
  expected <- predict(model, nd, re.form = control)
  actual <- predict(object, nd, f$inputs, re.form = control, include_fixed = TRUE)
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)
  expect_identical(attr(actual, "prediction_contract")$effects$response,
    "selected_conditional_effects")
  expect_identical(attr(actual, "prediction_contract")$integration$random_effects,
    "native_method")
  snap <- prepare_prediction(object, nd, f$inputs, include_fixed = TRUE,
    native_args = list(re.form = control))
  expect_equal(predict(snap), actual)
  expect_match(snap$prediction_contract$response$effects$selection, "site")
  bad <- nd
  bad$site <- NULL
  expect_error(predict(object, bad, f$inputs, re.form = control, include_fixed = TRUE),
    "formula variables")
  bad <- nd
  bad$z <- NULL
  expect_error(predict(object, bad, f$inputs, re.form = control, include_fixed = TRUE),
    "formula variables")
  expect_error(predict(object, nd, f$inputs, include_fixed = TRUE), "formula variables")
  expect_error(predict(object, nd, f$inputs, re.form = ~(1 | other),
    include_fixed = TRUE), "formula variables")
  expect_equal(.msr_selection_parameter_count(object),
    attr(logLik(model), "df") + nrow(object$scale_est))

  rasters <- suppressWarnings(kernel_scale.raster(f$rasters, multiScaleR = object,
    scale_center = TRUE, include_fixed = TRUE, match_buffer = TRUE,
    site_values = list(z = .2, site = levels(dat$site)[1]), verbose = FALSE))
  expect_false("other" %in% names(rasters))
  map <- predict(object, rasters, re.form = control)
  reference <- dat
  reference$z <- .2
  reference$site <- factor(levels(dat$site)[1], levels = levels(dat$site))
  expect_equal(terra::extract(map, f$points, ID = FALSE)[[1]],
    as.numeric(predict(model, reference, re.form = control)), tolerance = 1e-8)
})

test_that("GAM term selections are explicit and keep checked input requirements", {
  skip_if_not_installed("mgcv")
  f <- refit_spatial_fixture()
  dat <- f$data
  set.seed(793)
  dat$y <- dat$mean + dat$rough + sin(dat$z * 2) + rnorm(nrow(dat), sd = .2)
  object <- f$object
  model <- object$opt_mod <- mgcv::gam(y ~ mean + rough + fixed + s(z, k = 4),
    data = dat, method = "ML")
  nd <- dat[c(3, 21, 14), ]
  controls <- list(terms = "s(z)", iterms.type = 2)
  point <- predict(object, nd, f$inputs, type = "iterms", native_args = controls)
  snap <- prepare_prediction(object, nd, f$inputs, component = "iterms",
    native_args = controls)
  expect_equal(predict(snap), point)
  expect_identical(snap$prediction_contract$iterms$term_selection,
    list(terms = "s(z)", exclude = NULL, iterms.type = 2))
  expected <- predict(model, nd, type = "iterms", terms = "s(z)", iterms.type = 2)
  expect_equal(as.numeric(point), as.numeric(expected), tolerance = 1e-10)
  omitted <- nd
  omitted$z <- NULL
  expect_error(predict(object, omitted, f$inputs,
    native_args = list(exclude = "s(z)")), "formula variables")
  expect_error(prepare_prediction(object, omitted, f$inputs,
    native_args = list(exclude = "s(z)")), "formula variables")
  info <- .msr_native_prediction_contract(NULL, "gam", "terms",
    list(exclude = "s(z)", unconditional = TRUE))
  expect_identical(info$term_selection$exclude, "s(z)")
  expect_identical(info$uncertainty$smoothing_parameters, "native_unconditional_requested")
  expect_identical(info$uncertainty$kernel_parameters, "fixed_at_estimate")
})

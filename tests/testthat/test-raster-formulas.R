test_that("analytic raster products retain native factors, learned bases and covariance", {
  f <- refit_spatial_fixture()
  dat <- f$data
  # Identical response generation to the established analytic SQ raster test.
  set.seed(442)
  dat$y <- rpois(nrow(dat), exp(3 + 1.5 * dat$rough + .1 * dat$z))
  model <- glm(y ~ rough * group + poly(z, 2) + offset(log(exposure)),
    dat, family = poisson())
  warnings <- character()
  object <- withCallingHandlers(multiScale_optim(model, f$inputs,
    backend = "profile_glm", n_cores = 1, verbose = FALSE), warning = function(w) {
      text <- conditionMessage(w)
      if (grepl("Joint Wald uncertainty unavailable", text, fixed = TRUE)) {
        warnings <<- c(warnings, text)
        invokeRestart("muffleWarning")
      }
    })
  if (length(warnings)) message(paste(warnings, collapse = "\n"))
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = object,
    scale_center = TRUE, match_buffer = TRUE,
    site_values = list(group = "b", z = .2), offset_values = c(exposure = 2),
    verbose = FALSE)
  nd <- dat; nd$group <- factor("b", levels = levels(dat$group))
  nd$z <- .2; nd$exposure <- 2
  point <- predict(object, nd, f$inputs, se.fit = TRUE, uncertainty = "conditional")
  cached <- prepare_prediction(object, nd, f$inputs)
  expect_equal(predict(cached, se.fit = TRUE, uncertainty = "conditional"), point)
  map <- predict(object, inputs, se.fit = TRUE, uncertainty = "conditional")
  actual <- terra::extract(map, f$points, ID = FALSE)
  expect_equal(actual$fit, as.numeric(point$fit), tolerance = 1e-8)
  expect_equal(actual$se.fit, as.numeric(point$se.fit), tolerance = 1e-8)
  link <- predict(object, inputs, type = "link")
  expect_equal(terra::extract(link, f$points, ID = FALSE)[[1]],
    as.numeric(predict(object, nd, f$inputs, type = "link")), tolerance = 1e-8)
  expect_error(predict(object, inputs, se.fit = TRUE, uncertainty = "joint"), "conditional SEs")
  expect_identical(attr(map, "prediction_contract")$uncertainty$standard_errors, "provided")
  expect_equal(.msr_selection_parameter_count(object),
    as.numeric(attr(logLik(object$opt_mod), "df")) + 1L)
})

test_that("private spatial responses match cell-weight oracles across kernels, families and SQ", {
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
    for (metric in c("mean", "sq")) {
      x <- spatial_response_fixture(kernel, family, metric); b <- x$object; p <- x$parameters
      before <- list(b, p, .Random.seed)
      data <- x$fixture$data; expected <- spatial_response_oracle(x, data)
      actual <- .msr_spatial_response_predict(b, parameters = p, type = "link")
      expect_equal(as.numeric(actual), expected, tolerance = 1e-10)
      rows <- c(7, 2, 2, 5); nd <- data[rows, ]; nd$count <- NULL
      nd$east <- nd$east + .13; nd$north <- nd$north - .07
      expected <- spatial_response_oracle(x, nd)
      link <- .msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = p, type = "link")
      response <- .msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = p)
      expect_equal(as.numeric(link), expected, tolerance = 1e-10)
      mu <- switch(family, gaussian = expected, binomial = plogis(expected), exp(expected))
      expect_equal(as.numeric(response), mu, tolerance = 1e-10)
      expect_identical(response[2], response[3])
      expect_identical(list(b, p, .Random.seed), before)
    }
  # These metadata flags share one branch across the numerical case matrix.
  contract <- attr(response, "prediction_contract")
  expect_identical(contract$uncertainty$population_parameters, "fixed_parameters")
  expect_false(contract$uncertainty$observation_noise)
  expect_false(contract$metadata_establishes_eligibility)
  expect_false(contract$effects$spatial_scale$new_location_innovation)
  expect_identical(contract$effects$spatial_scale$basis, "frozen_training_projection_and_normalization")
})

test_that("spatial joint response predictions retain vectors, diagnostics and zero-field distinctions", {
  x <- spatial_response_fixture(metric = "sq"); b <- x$object
  b$parameters <- x$parameters; b <- fake_bayesian_draws(b, 4)
  b$draws$scale_z <- matrix(seq(-.4, .5, length.out = 12), 4)
  b$draws$beta[, 3] <- seq(.2, .8, length.out = 4)
  b$inference <- "joint Gaussian approximation"
  b$inference_metadata <- list(integration = list(method = "joint_laplace"))
  b$diagnostics <- list(passed = FALSE, original_failure = "retained")
  nd <- x$fixture$data[c(3, 3, 1), ]; nd$count <- NULL
  answer <- .msr_spatial_response_predict(b, nd, x$fixture$inputs)
  expected <- do.call(rbind, lapply(1:4, function(i)
    as.numeric(.msr_spatial_response_predict(b, nd, x$fixture$inputs,
      parameters = .msr_bayes_draw(b, i)))))
  expect_identical(dim(answer), c(4L, 3L))
  expect_equal(as.numeric(answer), as.numeric(expected), tolerance = 1e-12)
  expect_identical(answer[, 1], answer[, 2])
  contract <- attr(answer, "prediction_contract")
  expect_identical(contract$diagnostics_status, "failed")
  expect_identical(contract$method, b$inference)
  expect_identical(contract$integration$method, "joint_laplace")
  expect_true(contract$joint_parameter_map_complete)
  expect_false(contract$metadata_establishes_eligibility)
  expect_identical(contract$uncertainty$intervals, "not_provided")
  subset <- .msr_spatial_response_predict(b, nd, x$fixture$inputs, draws = 2)
  expect_equal(as.numeric(subset), as.numeric(expected[c(1, 4), ]))
  p <- x$parameters; p$log_scale_sd[] <- -Inf
  zero <- .msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = p)
  reference <- .msr_spatial_response_predict(b, nd, x$fixture$inputs,
    parameters = x$parameters, scale_re.form = NA)
  expect_equal(as.numeric(zero), as.numeric(reference), tolerance = 1e-12)
  expect_false(attr(reference, "prediction_contract")$effects$spatial_scale$included)
  no_coordinates <- nd; no_coordinates$east <- no_coordinates$north <- NULL
  expect_equal(as.numeric(.msr_spatial_response_predict(b, no_coordinates,
    x$fixture$inputs, parameters = x$parameters, scale_re.form = NA)), as.numeric(reference))
  bad <- b; bad$draws$scale_z <- bad$draws$scale_z[, -1, drop = FALSE]
  expect_error(.msr_spatial_response_predict(bad), "complete finite joint")
  bad <- b; bad$draws$scale_z[1, 1] <- Inf
  expect_error(.msr_spatial_response_predict(bad), "complete finite joint")
  expect_error(.msr_spatial_response_predict(b, parameters = p, draws = 2), "not both")
})

test_that("spatial response inputs keep fixed transforms, units and existing guards", {
  for (metric in c("mean", "sq")) {
    x <- spatial_response_fixture(metric = metric); b <- x$object; p <- x$parameters
    nd <- x$fixture$data[c(4, 1, 4), ]; nd$count <- NULL
    inputs <- x$fixture$inputs
    reference <- .msr_spatial_response_predict(b, nd, inputs, parameters = p)
    units <- inputs; units$unit_conv <- units$unit_conv / 10
    if (metric == "mean") units$binned$bin_rep <- units$binned$bin_rep * 10 else
      units$d_list <- lapply(units$d_list, function(d) d * 10)
    expect_equal(as.numeric(.msr_spatial_response_predict(b, nd, units, parameters = p)),
      as.numeric(reference), tolerance = 1e-12)
    future <- nd; future$year <- c("any", "other", "label")
    expect_identical(.msr_spatial_response_predict(b, future, inputs, parameters = p), reference)
    missing <- nd; missing$east <- NULL
    expect_error(.msr_spatial_response_predict(b, missing, inputs, parameters = p), "coordinates")
    invalid <- nd; invalid$north[1] <- NA_real_
    expect_error(.msr_spatial_response_predict(b, invalid, inputs, parameters = p), "finite numeric")
    bad <- inputs; bad$kernel <- "exp"
    expect_error(.msr_spatial_response_predict(b, nd, bad, parameters = p), "fitted kernel")
    invalid <- nd; invalid$key[1] <- "not_a_prepared_ID"
    expect_error(.msr_spatial_response_predict(b, invalid, inputs, parameters = p), "IDs do not match")
    expect_error(.msr_bayes_prediction_prepare(b, nd, inputs, "key", "error", "error", TRUE, TRUE),
      "Gaussian projection")
    extra <- b; extra$specification$response$random$unsupported <- b$specification$scales$habitat$random$.spatial_field
    expect_error(.msr_spatial_response_predict(extra, nd, inputs, parameters = p), "Gaussian projection")
    bad <- b; bad$target <- "unknown"
    expect_error(.msr_spatial_response_predict(bad, parameters = p), "Unknown spatial prototype")
    if (metric == "sq") {
      bad <- inputs; bad$scale_vars$source <- "different_source"
      expect_error(.msr_spatial_response_predict(b, nd, bad, parameters = p), "specifications must match")
    }
  }
})

test_that("new ordinary group coefficients append after the spatial field and remain coherent", {
  x <- spatial_response_fixture(groups = TRUE); b <- x$object; p <- x$parameters
  nd <- x$fixture$data[rep(1, 3), ]; nd$count <- NULL
  nd$site <- "new_site"; nd$year <- "new_year"
  expect_error(.msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = p), "New scale group")
  context <- .msr_spatial_response_context(b, nd, x$fixture$inputs, "key",
    "sample", "sample", TRUE, TRUE, FALSE)
  field <- context$spec$scales$habitat$random$.spatial_field
  expect_identical(field$z_index, b$specification$scales$habitat$random$.spatial_field$z_index)
  expect_equal(tail(context$spec$scales$habitat$random$year$z_index, 1), length(p$scale_z) + 1)
  before <- list(b, p, .Random.seed)
  result <- .msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = p,
    new_scale_levels = "sample", new_response_levels = "sample", seed = 319)
  expect_equal(as.numeric(result), rep(as.numeric(result)[1], 3))
  expect_identical(result, .msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = p,
    new_scale_levels = "sample", new_response_levels = "sample", seed = 319))
  expect_identical(list(b, p, .Random.seed), before)
  expect_identical(attr(result, "new_groups")$action, rep("sample", 2))
  bad <- p; bad$beta[1] <- NA_real_
  expect_error(.msr_spatial_response_predict(b, nd, x$fixture$inputs, parameters = bad,
    new_scale_levels = "sample", new_response_levels = "sample", seed = 319), "numeric parameter blocks")
  expect_identical(.Random.seed, before[[3]])
})

test_that("spatial responses reevaluate fixed metrics at their frozen reference scales", {
  skip_if_not_installed("terra")
  for (kernel in c("gaussian", "exp")) for (field in c("mean", "rough")) {
    f <- refit_spatial_fixture(kernel = kernel, bin = TRUE)
    f$data$key <- rownames(f$data)
    xy <- terra::crds(f$points)
    f$data$east <- xy[, 1]; f$data$north <- xy[, 2]
    b <- .msr_spatial_scale_prepare(y ~ mean + rough + fixed + z + offset(log(exposure)),
      f$inputs, f$data, scales = list(mean = ~1, rough = ~1), family = "poisson",
      row_id = "key", priors = list(scale_median = 20), covariate = field,
      coordinates = c("east", "north"), knots = rbind(c(60, 60), c(150, 60), c(100, 150)),
      range = 50, field_sd_prior = .4)
    p <- b$parameters; p$beta[] <- c(.2, .1, -.2, .4, .5)
    p$scale_z[] <- c(.3, -.2, .4)
    rows <- c(8, 2, 8); nd <- f$data[rows, ]; nd$y <- nd$fixed <- NULL
    before <- list(b, f$inputs, p)
    actual <- .msr_spatial_response_predict(b, nd, f$inputs, parameters = p,
      type = "link", include_fixed = TRUE)
    expected <- .msr_bayes_components(b$specification, p)$eta[rows]
    expect_equal(as.numeric(actual), unname(expected), tolerance = 1e-10)
    context <- .msr_spatial_response_context(b, nd, f$inputs, "key",
      "error", "error", TRUE, TRUE, TRUE)
    expect_identical(context$spec$center, b$specification$center)
    expect_identical(context$spec$spread, b$specification$spread)
    expect_identical(context$spec$fixed_reference, b$specification$fixed_reference)
    expect_identical(list(b, f$inputs, p), before)
    expect_error(.msr_spatial_response_predict(b, nd, f$inputs, parameters = p), "fixed")
  }
})

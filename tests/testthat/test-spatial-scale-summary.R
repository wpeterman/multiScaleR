spatial_summary_fixture <- function(kernel = "gaussian", metric = "mean") {
  x <- spatial_response_fixture(kernel = kernel, metric = metric)
  b <- fake_bayesian_draws(x$object, 4)
  b$draws$scale_z <- matrix(seq(-.4, .5, length.out = 12), 4)
  b$draws$scale_beta[, 1] <- log(c(20, 200, 1200, 400) / b$specification$unit_conv)
  b$draws$log_scale_sd[, 1] <- log(c(.15, .35, .5, .2))
  b
}

# Independent radial/projection algebra, not the prediction or summary helper.
spatial_summary_log_oracle <- function(b, data) {
  nm <- b$prototype$covariate; sm <- b$specification$scales[[nm]]
  field <- sm$random$.spatial_field$gaussian_block; basis <- field$basis
  xy <- as.matrix(data[field$coordinates])
  squared <- outer(xy[, 1], basis$knots[, 1], "-")^2 +
    outer(xy[, 2], basis$knots[, 2], "-")^2
  design <- (exp(-.5 * squared / basis$range^2) -
    matrix(1, nrow(xy), 1) %*% basis$fixed_projection) / basis$normalizer
  result <- b$draws$scale_z %*% t(design %*% field$chol)
  result * as.numeric(exp(b$draws$log_scale_sd[, 1])) +
    as.numeric(b$draws$scale_beta[, 1]) + log(b$specification$unit_conv)
}

test_that("location scale summaries preserve radial algebra, row dependence and reference limits", {
  for (kernel in c("gaussian", "exp")) for (metric in c("mean", "sq")) {
    b <- spatial_summary_fixture(kernel, metric)
    nd <- data.frame(east = c(.2, -.1, .2), north = c(.3, .4, .3))
    before <- list(b, .Random.seed)
    expected_log <- spatial_summary_log_oracle(b, nd)
    raw_log <- .msr_spatial_scale_summary(b, nd, scale = "log", summary = FALSE)
    raw <- .msr_spatial_scale_summary(b, nd, summary = FALSE)
    expect_identical(dim(raw), c(4L, 3L))
    expect_equal(as.numeric(raw_log), as.numeric(expected_log), tolerance = 1e-12)
    expect_equal(as.numeric(raw), as.numeric(exp(expected_log)), tolerance = 1e-12)
    expect_identical(raw[, 1], raw[, 3])
    actual <- .msr_spatial_scale_summary(b, nd, level = .8)
    values <- exp(expected_log)
    manual <- cbind(mean = colMeans(values), median = apply(values, 2, median),
      sd = apply(values, 2, sd), lower = apply(values, 2, quantile, probs = .1, names = FALSE),
      upper = apply(values, 2, quantile, probs = .9, names = FALSE))
    expect_equal(as.numeric(as.matrix(actual[names(as.data.frame(manual))])),
      as.numeric(manual), tolerance = 1e-12)
    expect_equal(actual$p_below_min_D, colMeans(values < b$specification$support[["min_D"]]))
    expect_equal(actual$p_above_max_D, colMeans(values > b$specification$support[["max_D"]]))
    expect_identical(actual$row, 1:3)
    expect_equal(as.numeric(as.matrix(actual[c("x", "y")])), as.numeric(as.matrix(nd)))
    expect_identical(attr(actual, "level"), .8)
    subset <- .msr_spatial_scale_summary(b, nd, draws = 2, summary = FALSE)
    expect_equal(as.numeric(subset), as.numeric(values[c(1, 4), ]), tolerance = 1e-12)
    training <- .msr_spatial_scale_summary(b, summary = FALSE)
    expect_equal(as.numeric(training), as.numeric(exp(spatial_summary_log_oracle(b, b$specification$data))),
      tolerance = 1e-12)
    permuted <- .msr_spatial_scale_summary(b, nd[c(3, 1, 2), ], summary = FALSE)
    expect_equal(as.numeric(permuted), as.numeric(values[, c(3, 1, 2)]), tolerance = 1e-12)
    nd$year <- c("irrelevant", "different", "labels")
    expect_identical(.msr_spatial_scale_summary(b, nd, summary = FALSE), raw)
    expect_identical(list(b, .Random.seed), before)
  }
})

test_that("private field summaries retain approximation provenance and reject unsupported layouts", {
  b <- spatial_summary_fixture(metric = "sq")
  b$inference <- "joint Gaussian approximation"
  b$inference_metadata <- list(integration = list(method = "joint_laplace"))
  b$importance <- list(resample_ids = 1:4); b$diagnostics <- list(passed = FALSE)
  nd <- data.frame(east = .2, north = .3)
  answer <- .msr_spatial_scale_summary(b, nd)
  contract <- attr(answer, "prediction_contract")
  expect_identical(contract$quantity, "kernel_scale")
  expect_identical(contract$diagnostics_status, "failed")
  expect_identical(contract$method, b$inference)
  expect_identical(contract$integration$method, "joint_laplace")
  expect_identical(contract$integration$prediction_vectors, "stored_joint_resamples")
  expect_true(contract$joint_parameter_map_complete)
  expect_false(contract$uncertainty$observation_noise)
  expect_false(contract$metadata_establishes_eligibility)
  expect_false(contract$effects$spatial_scale$new_location_innovation)
  expect_identical(attr(answer, "coordinate_columns"), c("east", "north"))
  raw <- .msr_spatial_scale_summary(b, nd, summary = FALSE)
  expect_identical(attr(raw, "prediction_contract")$uncertainty$intervals, "not_provided")
  expect_error(.msr_spatial_scale_summary(fake_bayesian_draws(spatial_response_fixture(groups = TRUE)$object)),
    "no ordinary scale-group")
  bad <- b; bad$draws$scale_z <- bad$draws$scale_z[, -1, drop = FALSE]
  expect_error(.msr_spatial_scale_summary(bad, nd), "complete finite joint")
  bad <- b; bad$draws$beta[2, 1] <- NA_real_
  expect_error(.msr_spatial_scale_summary(bad, nd, draws = 2), "complete finite joint")
  expect_error(.msr_spatial_scale_summary(b, nd, level = 1), "level")
  expect_error(.msr_spatial_scale_summary(b, nd[FALSE, ]), "nonempty data frame")
  bad_data <- nd; bad_data$east <- Inf
  expect_error(.msr_spatial_scale_summary(b, bad_data), "finite numeric")
  expect_error(predict_scales(b, nd), "Gaussian projection")
  expect_error(posterior_scales(b), "multiScaleR_bayes fit")
  before <- .Random.seed
  expect_identical(.msr_spatial_scale_summary(b, nd, seed = 319), answer)
  expect_error(.msr_spatial_scale_summary(b, bad_data, seed = 319), "finite numeric")
  expect_identical(.Random.seed, before)
})

test_that("private scale maps retain rows and original diagnostic state", {
  b <- spatial_summary_fixture()
  b$diagnostics <- list(passed = FALSE)
  nd <- data.frame(east = c(.2, -.1, .2), north = c(.3, .4, .3))
  render <- function(statistic, scale) {
    path <- tempfile(fileext = ".pdf")
    grDevices::pdf(path, width = 7, height = 6)
    on.exit(grDevices::dev.off())
    margin <- graphics::par("mar")
    table <- .msr_spatial_scale_plot(b, nd, statistic = statistic, scale = scale)
    list(table = table, path = path, margin_restored = identical(graphics::par("mar"), margin))
  }
  before <- list(b, .Random.seed)
  for (statistic in c("median", "sd", "lower", "upper")) {
    output <- render(statistic, "log")
    expect_identical(output$table, .msr_spatial_scale_summary(b, nd, scale = "log"))
    expect_gt(file.info(output$path)$size, 1000)
  }
  expect_true(output$margin_restored)
  b$draws <- lapply(b$draws, function(x) x[rep(1, 4), , drop = FALSE])
  b$draws$scale_z[,] <- 0
  constant <- render("median", "distance")
  expect_true(all(is.finite(constant$table$median)))
  expect_identical(length(unique(constant$table$median)), 1L)
  expect_error(.msr_spatial_scale_plot(b, nd, statistic = "sd", draws = 1), "Nonfinite spatial summary colors")
  expect_identical(.Random.seed, before[[2]])
})

test_that("response effects do not enter the ecological-scale summary", {
  x <- spatial_response_fixture()
  b <- .msr_spatial_scale_prepare(count ~ habitat + temp + (1 | site),
    x$fixture$inputs, x$fixture$data, scales = list(habitat = ~1), family = "gaussian",
    row_id = "key", priors = list(scale_median = 200), covariate = "habitat",
    coordinates = c("east", "north"), knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2),
    range = .7, field_sd_prior = .4)
  b <- fake_bayesian_draws(b, 4)
  nd <- data.frame(east = c(.2, -.1), north = c(.3, .4))
  original <- .msr_spatial_scale_summary(b, nd, summary = FALSE)
  before <- .Random.seed
  b$draws$response_z[,] <- 100; b$draws$beta[,] <- -100
  changed <- .msr_spatial_scale_summary(b, nd, summary = FALSE)
  expect_identical(changed, original)
  expect_identical(attr(changed, "prediction_contract")$effects$response$known, "not_applicable")
  expect_identical(.Random.seed, before)
})

spatial_contrast_fixture <- function(kernel = "gaussian", metric = "mean") {
  b <- fake_bayesian_draws(spatial_response_fixture(kernel = kernel, metric = metric)$object, 4)
  b$draws$scale_z <- matrix(seq(-.4, .5, length.out = 12), 4)
  b$draws$scale_beta[, 1] <- log(c(20, 200, 1200, 400) / b$specification$unit_conv)
  b$draws$log_scale_sd[, 1] <- log(c(.15, .35, .5, .2))
  b
}

# Independent radial projection and joint-vector algebra for the selected field.
spatial_contrast_design <- function(b, data) {
  field <- b$specification$scales$habitat$random$.spatial_field$gaussian_block
  basis <- field$basis; xy <- as.matrix(data[field$coordinates])
  squared <- outer(xy[, 1], basis$knots[, 1], "-")^2 +
    outer(xy[, 2], basis$knots[, 2], "-")^2
  projected <- (exp(-.5 * squared / basis$range^2) -
    matrix(1, nrow(xy), 1) %*% basis$fixed_projection) / basis$normalizer
  projected %*% field$chol
}

test_that("paired scales propagate the same joint field at both locations", {
  nd <- data.frame(east = c(.2, -.1, .2, 99), north = c(.3, .4, .3, 88))
  pairs <- rbind(c(1, 2), c(1, 3), c(1, 1), c(2, 1), c(1, 2))
  rownames(pairs) <- c("forward", "equal coordinates", "same row", "reverse", "forward")
  for (kernel in c("gaussian", "exp")) for (metric in c("mean", "sq")) {
    b <- spatial_contrast_fixture(kernel, metric); before <- list(b, .Random.seed)
    design <- spatial_contrast_design(b, nd)
    expected <- b$draws$scale_z %*% t(design[pairs[, 2], , drop = FALSE] -
      design[pairs[, 1], , drop = FALSE]) * as.numeric(exp(b$draws$log_scale_sd[, 1]))
    actual <- .msr_spatial_scale_contrasts(b, pairs, nd, summary = FALSE)
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
    expect_identical(actual[, 2:3], matrix(0, 4, 2,
      dimnames = list(NULL, rownames(pairs)[2:3])))
    expect_identical(actual[, 1], actual[, 5])
    expect_equal(actual[, 4], -actual[, 1], tolerance = 1e-12)
    expect_identical(attr(actual, "location_pairs"), pairs)
    expect_identical(attr(actual, "prediction_locations")$row, 1:3)
    expect_identical(colnames(actual), rownames(pairs))
    ratio <- .msr_spatial_scale_contrasts(b, pairs, nd, scale = "ratio", summary = FALSE)
    expect_equal(as.numeric(ratio), as.numeric(exp(expected)), tolerance = 1e-12)
    for (scale in c("log_ratio", "ratio")) {
      values <- if (scale == "ratio") exp(expected) else expected
      table <- .msr_spatial_scale_contrasts(b, pairs, nd, scale, level = .8)
      manual <- data.frame(mean = colMeans(values), median = apply(values, 2, median),
        sd = apply(values, 2, sd), lower = apply(values, 2, quantile, probs = .1, names = FALSE),
        upper = apply(values, 2, quantile, probs = .9, names = FALSE), row.names = NULL)
      expect_equal(table[names(manual)], manual, tolerance = 1e-12)
      expect_identical(table$from_row, as.integer(pairs[, 1]))
      expect_identical(table$to_row, as.integer(pairs[, 2]))
      expect_identical(attr(table, "level"), .8)
    }
    subset <- .msr_spatial_scale_contrasts(b, pairs, nd, draws = 2, summary = FALSE)
    expect_equal(as.numeric(subset), as.numeric(expected[c(1, 4), ]), tolerance = 1e-12)
    expect_identical(list(b, .Random.seed), before)
    b$specification$unit_conv <- 2 * b$specification$unit_conv
    b$draws$scale_beta <- b$draws$scale_beta - log(2)
    expect_equal(as.numeric(.msr_spatial_scale_contrasts(b, pairs, nd, summary = FALSE)),
      as.numeric(actual), tolerance = 1e-12)
  }
})

test_that("contrasts retain failed approximate provenance and scale-only conditioning", {
  b <- spatial_contrast_fixture(metric = "sq")
  b$inference <- "joint Gaussian approximation"
  b$inference_metadata <- list(integration = list(method = "joint_laplace"))
  b$importance <- list(resample_ids = 1:4); b$diagnostics <- list(passed = FALSE)
  nd <- data.frame(east = c(.2, -.1), north = c(.3, .4)); pairs <- matrix(1:2, 1)
  before <- list(b, .Random.seed)
  d <- .msr_spatial_scale_contrasts(b, pairs, nd)
  contract <- attr(d, "prediction_contract")
  expect_identical(contract$quantity, "kernel_scale_contrast")
  expect_identical(contract$output_scale, "log_ratio")
  expect_identical(contract$target_population, "specified_location_pairs")
  expect_identical(contract$method, b$inference)
  expect_identical(contract$integration$prediction_vectors, "stored_joint_resamples")
  expect_identical(contract$integration$method, "joint_laplace")
  expect_identical(contract$diagnostics_status, "failed")
  expect_false(contract$metadata_establishes_eligibility)
  expect_false(contract$uncertainty$observation_noise)
  expect_identical(contract$uncertainty$calibration, "not_established")
  expect_identical(contract$uncertainty$intervals, "conditional_posterior_contrast_quantiles")
  expect_identical(contract$component$direction, "to_over_from")
  expect_false(contract$effects$spatial_scale$new_location_innovation)
  expect_identical(contract$effects$response$known, "not_applicable")
  raw <- .msr_spatial_scale_contrasts(b, pairs, nd, summary = FALSE)
  expect_identical(attr(raw, "prediction_contract")$uncertainty$intervals, "not_provided")
  expect_identical(.msr_spatial_scale_contrasts(b, pairs, nd, seed = 129), d)
  expect_identical(list(b, .Random.seed), before)
  b$draws$scale_z[,] <- 0
  expect_identical(as.numeric(.msr_spatial_scale_contrasts(b, pairs, nd, summary = FALSE)), rep(0, 4))
  expect_identical(as.numeric(.msr_spatial_scale_contrasts(b, pairs, nd, scale = "ratio", summary = FALSE)), rep(1, 4))
  training <- .msr_spatial_scale_contrasts(b, matrix(c(1, 2), 1), summary = FALSE)
  expect_identical(dim(training), c(4L, 1L))
})

test_that("response effects do not enter a scale contrast", {
  x <- spatial_response_fixture(); f <- x$fixture
  b <- .msr_spatial_scale_prepare(count ~ habitat + temp + offset(log(effort)) + (1 | site),
    f$inputs, f$data, scales = list(habitat = ~1), family = "gaussian", row_id = "key",
    priors = list(scale_median = 200), covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7, field_sd_prior = .4)
  b <- fake_bayesian_draws(b, 4); b$draws$scale_z <- spatial_contrast_fixture()$draws$scale_z
  nd <- data.frame(east = c(.2, -.1), north = c(.3, .4)); pairs <- matrix(1:2, 1)
  d <- .msr_spatial_scale_contrasts(b, pairs, nd)
  b$draws$response_z <- b$draws$response_z + 1000; b$draws$beta <- b$draws$beta + 1000
  expect_identical(.msr_spatial_scale_contrasts(b, pairs, nd), d)
  expect_identical(attr(d, "prediction_contract")$effects$response$known, "not_applicable")
})

test_that("ratio summaries and numeric limits preserve their requested scale", {
  b <- spatial_contrast_fixture(); nd <- data.frame(east = c(.2, -.1), north = c(.3, .4))
  design <- spatial_contrast_design(b, nd); loading <- design[2, ] - design[1, ]
  b$draws$log_scale_sd[,] <- 0
  b$draws$scale_z <- outer(c(-1, 0, 1, 3), loading / sum(loading^2))
  pairs <- matrix(1:2, 1)
  log_result <- .msr_spatial_scale_contrasts(b, pairs, nd)
  ratio <- .msr_spatial_scale_contrasts(b, pairs, nd, scale = "ratio")
  expect_equal(log_result$median, .5, tolerance = 1e-12)
  expect_equal(ratio$median, (1 + exp(1)) / 2, tolerance = 1e-12)
  expect_gt(abs(ratio$median - exp(log_result$median)), .1)
  b$draws$scale_z <- outer(rep(1000, 4), loading / sum(loading^2))
  expect_true(all(is.finite(.msr_spatial_scale_contrasts(b, pairs, nd, summary = FALSE))))
  expect_error(.msr_spatial_scale_contrasts(b, pairs, nd, scale = "ratio"), "use scale = 'log_ratio'")
  expect_error(.msr_spatial_scale_contrasts(b, pairs[, 2:1, drop = FALSE], nd, scale = "ratio"),
    "use scale = 'log_ratio'")
})

test_that("invalid pairs and incomplete source vectors retain explicit errors", {
  b <- spatial_contrast_fixture(); nd <- data.frame(east = c(.2, -.1, 1), north = c(.3, .4, 1))
  unsupported <- spatial_response_fixture(groups = TRUE)$object
  pairs <- matrix(1:2, 1); before <- .Random.seed
  for (bad in list(NULL, 1:2, data.frame(from = 1, to = 2), matrix(numeric(), 0, 2),
      matrix(1, 1, 3), matrix(c(0, 1), 1), matrix(c(-1, 1), 1), matrix(c(1, 4), 1),
      matrix(c(1, NA_real_), 1), matrix(c(1, Inf), 1), matrix(c(1, 1.5), 1), matrix(1i, 1, 2)))
    expect_error(.msr_spatial_scale_contrasts(b, bad, nd), "pairs must have two columns")
  expect_error(.msr_spatial_scale_contrasts(list(), pairs, nd), "private spatial")
  expect_error(.msr_spatial_scale_contrasts(b, pairs, nd[FALSE, ]), "nonempty data frame")
  bad <- nd; bad$east[3] <- Inf
  expect_error(.msr_spatial_scale_contrasts(b, pairs, bad), "finite numeric")
  expect_error(.msr_spatial_scale_contrasts(b, pairs, nd, level = 1), "level")
  expect_error(.msr_spatial_scale_contrasts(b, pairs, nd, draws = 0), "draws")
  expect_error(.msr_spatial_scale_contrasts(b, pairs, nd, seed = 0), "seed")
  bad <- b; bad$draws$beta[2, 1] <- NA_real_
  expect_error(.msr_spatial_scale_contrasts(bad, pairs, nd, draws = 2), "complete finite joint")
  expect_error(.msr_spatial_scale_contrasts(unsupported, pairs, nd),
    "no ordinary scale-group")
  expect_identical(.Random.seed, before)
})

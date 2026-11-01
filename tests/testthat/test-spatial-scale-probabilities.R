# Synthetic arrays exercise software contracts only, without fitting or sampling.
spatial_probability_fixture <- function(kernel = "gaussian", family = "gaussian",
    metric = "mean", iterations = 80L, chains = 3L) {
  x <- spatial_response_fixture(kernel, family, metric)
  b <- x$object; b$sampler <- "deterministic software array"
  b$diagnostics <- list(passed = FALSE, marker = "retain original failure")
  columns <- unlist(lapply(names(b$parameters), function(nm) {
    n <- length(b$parameters[[nm]])
    if (!n) character() else if (n == 1L) nm else paste0(nm, "[", seq_len(n), "]")
  }))
  n <- iterations * chains
  set.seed(478621)
  values <- matrix(rep(unlist(b$parameters), each = n), n) +
    matrix(rnorm(n * length(columns), sd = .3), n)
  list(object = b, array = array(values, c(iterations, chains, length(columns)),
    dimnames = list(NULL, NULL, columns)))
}

spatial_probability_log_oracle <- function(b, a, data) {
  sm <- b$specification$scales$habitat
  field <- sm$random$.spatial_field$gaussian_block; basis <- field$basis
  xy <- as.matrix(data[field$coordinates])
  squared <- outer(xy[, 1], basis$knots[, 1], "-")^2 +
    outer(xy[, 2], basis$knots[, 2], "-")^2
  design <- (exp(-.5 * squared / basis$range^2) -
    matrix(1, nrow(xy), 1) %*% basis$fixed_projection) / basis$normalizer
  result <- array(NA_real_, c(dim(a)[1:2], nrow(data)))
  for (chain in seq_len(dim(a)[2])) for (row in seq_len(nrow(data))) {
    result[, chain, row] <- a[, chain, "scale_beta"] + log(b$specification$unit_conv) +
      exp(a[, chain, "log_scale_sd"]) * drop(a[, chain, paste0("scale_z[", 1:3, "]")] %*%
        as.numeric(design[row, , drop = FALSE] %*% field$chol))
  }
  result
}

test_that("location events retain original chain order and independent radial algebra", {
  skip_if_not_installed("posterior")
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
    for (metric in c("mean", "sq")) {
      f <- spatial_probability_fixture(kernel, family, metric); b <- f$object; a <- f$array
      # Deliberately contradictory permuted vectors must never supply chains.
      b <- fake_bayesian_draws(b, 3); b$draws$scale_beta[,] <- 1000
      nd <- data.frame(east = c(.2, -.1, .2), north = c(.3, .4, .3))
      thresholds <- c(200, 100, 400)
      expected <- spatial_probability_log_oracle(b, a, nd)
      local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
      before <- serialize(b, NULL); rng <- .Random.seed
      d <- .msr_spatial_scale_probabilities(b, thresholds, nd, mcse_target = .1)
      manual <- do.call(rbind, lapply(1:3, function(row) do.call(rbind,
        lapply(thresholds, function(threshold) {
          indicator <- 1 * (expected[, , row] <= log(threshold))
          data.frame(estimate = mean(indicator), mcse = posterior::mcse_mean(indicator),
            ess_mean = posterior::ess_mean(indicator), rhat = posterior::rhat(indicator))
        }))))
      expect_equal(d$probabilities[names(manual)], manual, tolerance = 1e-12)
      expected_chains <- unlist(lapply(1:3, function(row) lapply(thresholds,
        function(threshold) colMeans(expected[, , row] <= log(threshold)))))
      expect_equal(d$chain_probabilities$estimate, expected_chains)
      expect_identical(d$probabilities$row, rep(1:3, each = 3))
      expect_identical(list(serialize(b, NULL), .Random.seed), list(before, rng))
      log_report <- .msr_spatial_scale_probabilities(b, log(thresholds), nd, scale = "log", mcse_target = .1)
      expect_identical(log_report$probabilities[names(manual)], d$probabilities[names(manual)])
    }
})

test_that("spatial probability reports preserve provenance, repeated rows and zero-field scalar events", {
  skip_if_not_installed("posterior")
  f <- spatial_probability_fixture(); b <- f$object; a <- f$array
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  nd <- data.frame(east = c(.2, -.1, .2), north = c(.3, .4, .3))
  d <- .msr_spatial_scale_probabilities(b, c(100, 200), nd)
  expect_true(d$available)
  expect_identical(d$raw_diagnostics, b$diagnostics)
  expect_false(d$metadata_establishes_eligibility)
  expect_false("passed" %in% names(d))
  expect_false(any(c("draws", "sampler", "parameters") %in% names(d)))
  expect_identical(d$spatial$basis, "frozen_training_projection_and_normalization")
  expect_false(d$spatial$new_location_innovation)
  expect_identical(d$chain_source, "original_unpermuted_post_warmup_NUTS")
  expect_identical(d$method, b$inference)
  expect_identical(d$coordinate_columns, c("east", "north"))
  expect_equal(d$probabilities$x, rep(nd$east, each = 2))
  columns <- c("estimate", "mcse", "ess_mean", "rhat", "reason")
  expect_identical(d$probabilities[1:2, columns],
    `rownames<-`(d$probabilities[5:6, columns], NULL))
  nd$year <- c("ignored", "never fitted", "ignored")
  expect_identical(.msr_spatial_scale_probabilities(b, c(100, 200), nd), d)
  original <- a
  dimnames(a)[[3]] <- ifelse(grepl("[", dimnames(a)[[3]], fixed = TRUE),
    dimnames(a)[[3]], paste0(dimnames(a)[[3]], "[1]"))
  expect_identical(.msr_spatial_scale_probabilities(b, c(100, 200), nd), d)
  a <- original[, , rev(seq_len(dim(original)[3])), drop = FALSE]
  expect_identical(.msr_spatial_scale_probabilities(b, c(100, 200), nd), d)
  a <- original
  b$specification$unit_conv <- 2 * b$specification$unit_conv
  a[, , "scale_beta"] <- a[, , "scale_beta"] - log(2)
  expect_identical(.msr_spatial_scale_probabilities(b, c(100, 200), nd), d)
  b$specification$unit_conv <- b$specification$unit_conv / 2
  a <- original
  training <- .msr_spatial_scale_probabilities(b, 200)
  expect_equal(nrow(training$probabilities), nrow(b$specification$data))
  a[, , paste0("scale_z[", 1:3, "]")] <- 0
  private <- .msr_spatial_scale_probabilities(b, c(100, 200), nd)
  scalar <- b; scalar$specification$scales$habitat$random$.spatial_field <- NULL
  class(scalar) <- "multiScaleR_bayes"
  reference <- diagnose_bayes_probabilities(scalar, list(habitat = c(100, 200)))
  expect_identical(private$probabilities[1:2, columns], reference$probabilities[columns])
  expect_error(diagnose_bayes_probabilities(b, list(habitat = 200)), "multiScaleR_bayes")
  class(b) <- c(class(b), "multiScaleR_bayes")
  expect_match(diagnose_bayes_probabilities(b, list(habitat = 200))$reason, "dedicated diagnostic adapter")
})

test_that("unavailable and invalid spatial chains never acquire precision", {
  f <- spatial_probability_fixture(); b <- f$object; a <- f$array
  nd <- data.frame(east = .2, north = .3)
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  missing <- b; missing$sampler <- NULL
  d <- .msr_spatial_scale_probabilities(missing, 200, nd)
  expect_false(d$available)
  expect_match(d$reason, "not retained")
  expect_equal(nrow(d$probabilities), 0)
  expect_identical(d$probabilities$precision_target_met, logical())
  for (mode in c("laplace", "marginal", "importance")) {
    approx <- b; approx[[mode]] <- list(passed = TRUE)
    expect_match(.msr_spatial_scale_probabilities(approx, 200, nd)$reason, "approximate draws")
  }
  expect_error(.msr_spatial_scale_probabilities(list(), 200, nd), "private spatial")
  expect_error(.msr_spatial_scale_probabilities(spatial_response_fixture(groups = TRUE)$object, 200, nd),
    "no ordinary scale-group")
  for (threshold in list(NULL, numeric(), c(1, 1), 0, -1, NA_real_, Inf, "200", matrix(200), 1i))
    expect_error(.msr_spatial_scale_probabilities(b, threshold, nd), "thresholds")
  for (target in list(0, 1, NA_real_, Inf, c(.1, .2), ".1", matrix(.1), .1i))
    expect_error(.msr_spatial_scale_probabilities(b, 200, nd, mcse_target = target), "mcse_target")
  expect_error(.msr_spatial_scale_probabilities(b, 200, nd[FALSE, ]), "nonempty")
  bad <- nd; bad$east <- Inf
  expect_error(.msr_spatial_scale_probabilities(b, 200, bad), "finite numeric")
  original <- a; a[1, 1, "beta[1]"] <- NA_real_
  expect_error(.msr_spatial_scale_probabilities(b, 200, nd), "complete finite original")
  a <- original[, , -1, drop = FALSE]
  expect_error(.msr_spatial_scale_probabilities(b, 200, nd), "original parameter layout")
  a <- original; dimnames(a)[[3]][1] <- dimnames(a)[[3]][2]
  expect_error(.msr_spatial_scale_probabilities(b, 200, nd), "named iterations")
  a <- array(c(original, original[, , "scale_beta"]), c(dim(original)[1:2], dim(original)[3] + 1L),
    dimnames = list(NULL, NULL, c(dimnames(original)[[3]], "scale_beta[1]")))
  expect_error(.msr_spatial_scale_probabilities(b, 200, nd), "ambiguous original")
  a <- original; a[, , "scale_beta"] <- 0; a[, , paste0("scale_z[", 1:3, "]")] <- 0
  d <- .msr_spatial_scale_probabilities(b, c(-1, 0, 1) + log(b$specification$unit_conv), nd,
    scale = "log", mcse_target = .1)
  expect_identical(d$probabilities$estimate, c(0, 1, 1))
  expect_true(all(d$probabilities$constant))
  expect_false(any(d$probabilities$precision_available))
  expect_true(all(is.na(d$probabilities$mcse)))
  expect_true(all(is.na(d$probabilities$precision_target_met)))
  a[, , "log_scale_sd"] <- 1000
  d <- .msr_spatial_scale_probabilities(b, 200, nd)
  expect_match(d$probabilities$reason, "no iterations were dropped")
  expect_true(is.na(d$probabilities$estimate))
  a <- original[, 1, , drop = FALSE]
  expect_match(.msr_spatial_scale_probabilities(b, 200, nd)$probabilities$reason, "two chains")
  a <- original[1:3, , , drop = FALSE]
  expect_match(.msr_spatial_scale_probabilities(b, 200, nd)$probabilities$reason, "four iterations")
  a <- original
  a[, , "scale_beta"] <- a[, , "scale_beta"] - b$parameters$scale_beta + 1000
  d <- .msr_spatial_scale_probabilities(b, 1000 + log(b$specification$unit_conv), nd, scale = "log")
  expect_true(d$probabilities$precision_available)
  expect_true(is.finite(d$probabilities$estimate))
})

test_that("private probability mean precision remains distinct from auxiliary Rhat", {
  skip_if_not_installed("posterior")
  f <- spatial_probability_fixture(); b <- f$object; a <- f$array
  a[, , paste0("scale_z[", 1:3, "]")] <- 0
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  local_mocked_bindings(rhat = function(x, ...) Inf, .package = "posterior")
  d <- .msr_spatial_scale_probabilities(b, 200, data.frame(east = .2, north = .3), mcse_target = .5)
  expect_true(d$probabilities$precision_available)
  expect_true(d$probabilities$precision_target_met)
  expect_identical(d$probabilities$rhat, Inf)
  expect_false(d$probabilities$rhat_available)
  expect_match(d$probabilities$rhat_reason, "nonfinite")
  expect_identical(d$raw_diagnostics, b$diagnostics)
})

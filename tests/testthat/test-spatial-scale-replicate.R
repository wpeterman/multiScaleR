# Small deterministic family-RNG oracles, not scientific calibration.
spatial_replicate_oracle <- function(object, data, inputs, nsim, seed, trials = NULL, include_fixed = FALSE) {
  saved <- .Random.seed
  on.exit(assign(".Random.seed", saved, envir = .GlobalEnv))
  set.seed(seed)
  selected <- sample(seq_len(nrow(object$draws$beta)), nsim, replace = TRUE)
  ans <- matrix(NA_real_, nrow(data), nsim)
  for (j in seq_len(nsim)) {
    p <- .msr_bayes_draw(object, selected[j])
    eta <- as.numeric(.msr_spatial_response_predict(object, data, inputs, parameters = p, type = "link", include_fixed = include_fixed))
    ans[, j] <- switch(object$specification$family,
      gaussian = rnorm(length(eta), eta, exp(p$log_dispersion)),
      poisson = rpois(length(eta), exp(eta)),
      binomial = rbinom(length(eta), if (is.null(trials)) 1 else trials, plogis(eta)),
      nbinom2 = rnbinom(length(eta), mu = exp(eta), size = exp(p$log_dispersion)))
  }
  ans
}

test_that("private spatial replicates match independent family RNG and preserve joint metadata", {
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
    for (metric in c("mean", "sq")) {
      x <- spatial_response_fixture(kernel, family, metric); b <- x$object
      b$parameters <- x$parameters; b <- fake_bayesian_draws(b, 4)
      b$draws$beta[, 3] <- seq(.2, .8, length.out = 4)
      b$draws$scale_z <- matrix(seq(-.4, .5, length.out = 12), 4)
      b$diagnostics <- list(passed = FALSE)
      nd <- x$fixture$data[c(4, 1, 4), ]; nd$count <- NULL
      before <- list(b, .Random.seed)
      expected <- spatial_replicate_oracle(b, nd, x$fixture$inputs, 3, 718)
      actual <- .msr_spatial_response_simulate(b, 3, seed = 718, newdata = nd,
        kernel_inputs = x$fixture$inputs)
      expect_identical(dim(actual), c(3L, 3L))
      expect_equal(unname(as.matrix(actual)), expected, tolerance = 0)
      expect_identical(list(b, .Random.seed), before)
      expect_identical(attr(actual, "point_ids"), as.character(nd$key))
    }
  # The metadata path is shared by all family-specific noise calculations.
  contract <- attr(actual, "prediction_contract")
  expect_identical(contract$quantity, "replicate_response")
  expect_true(contract$uncertainty$observation_noise)
  expect_identical(contract$diagnostics_status, "failed")
  expect_false(contract$metadata_establishes_eligibility)
  expect_false(contract$effects$spatial_scale$new_location_innovation)
})

test_that("private replicate zero-field references match existing scalar simulation", {
  for (metric in c("mean", "sq")) {
    x <- spatial_response_fixture(metric = metric); b <- x$object
    b$parameters <- x$parameters; b <- fake_bayesian_draws(b, 4)
    scalar <- multiScale_bayes(count ~ habitat + temp + offset(log(effort)),
      x$fixture$inputs, x$fixture$data, scales = list(habitat = ~1),
      family = "gaussian", row_id = "key", priors = list(scale_median = 200), method = "prepare")
    scalar <- fake_bayesian_draws(scalar, 4); scalar$draws$beta <- b$draws$beta
    actual <- .msr_spatial_response_simulate(b, 3, seed = 719, scale_re.form = NA)
    expected <- simulate(scalar, 3, seed = 719)
    expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), tolerance = 0)
    expect_false(attr(actual, "prediction_contract")$effects$spatial_scale$included)
    expect_error(simulate.multiScaleR_bayes(b, 3, seed = 719), "Gaussian projection")
  }
})

test_that("private replicate trial rules keep grouped-binomial targets explicit", {
  x <- spatial_response_fixture(family = "binomial")
  data <- x$fixture$data; data$failures <- 5L - data$count
  b <- .msr_spatial_scale_prepare(cbind(count, failures) ~ habitat + temp + offset(log(effort)),
    x$fixture$inputs, data, scales = list(habitat = ~1), family = "binomial", row_id = "key",
    priors = list(scale_median = 200), covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7, field_sd_prior = .4)
  b$parameters <- x$parameters; b <- fake_bayesian_draws(b, 4)
  nd <- data[c(3, 1, 3), ]; nd$count <- nd$failures <- NULL
  before <- .Random.seed
  expect_error(.msr_spatial_response_simulate(b, 3, seed = 718, newdata = nd,
    kernel_inputs = x$fixture$inputs), "trials explicitly")
  actual <- .msr_spatial_response_simulate(b, 3, seed = 718, newdata = nd,
    kernel_inputs = x$fixture$inputs, trials = c(2, 5, 7))
  expected <- spatial_replicate_oracle(b, nd, x$fixture$inputs, 3, 718, c(2, 5, 7))
  expect_equal(unname(as.matrix(actual)), expected, tolerance = 0)
  expect_true(all(as.matrix(actual) >= 0 & as.matrix(actual) <= c(2, 5, 7)))
  expect_error(.msr_spatial_response_simulate(b, trials = 0, seed = 718), "positive integer totals")
  expect_identical(.Random.seed, before)
})

test_that("spatial replicates share new group effects but draw separate observation noise", {
  x <- spatial_response_fixture(groups = TRUE); b <- x$object
  b$parameters <- x$parameters; b <- fake_bayesian_draws(b, 4)
  nd <- x$fixture$data[rep(1, 3), ]; nd$count <- NULL
  nd$site <- "new_site"; nd$year <- "new_year"
  context <- .msr_spatial_response_context(b, nd, x$fixture$inputs, "key",
    "sample", "sample", TRUE, TRUE, FALSE)
  manual <- function() {
    saved <- .Random.seed
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv))
    set.seed(718); selected <- sample(1:4, 3, replace = TRUE)
    result <- matrix(NA_real_, 3, 3)
    for (j in 1:3) {
      p <- .msr_bayes_draw(b, selected[j])
      p$scale_z <- c(p$scale_z, rnorm(context$n_scale))
      p$response_z <- c(p$response_z, rnorm(context$n_response))
      eta <- .msr_bayes_components(context$spec, p)$eta
      result[, j] <- rnorm(3, eta, exp(p$log_dispersion))
    }
    result
  }
  before <- list(b, .Random.seed); expected <- manual()
  actual <- .msr_spatial_response_simulate(b, 3, seed = 718, newdata = nd,
    kernel_inputs = x$fixture$inputs, new_scale_levels = "sample",
    new_response_levels = "sample")
  expect_equal(unname(as.matrix(actual)), expected, tolerance = 0)
  expect_identical(list(b, .Random.seed), before)
  expect_identical(attr(actual, "new_groups")$action, rep("sample", 2))
})

test_that("private spatial replicates reject incomplete vectors and restore RNG on failure", {
  x <- spatial_response_fixture(family = "poisson"); b <- x$object
  b$parameters <- x$parameters; b <- fake_bayesian_draws(b, 4)
  before <- .Random.seed
  bad <- b; bad$draws$scale_z <- bad$draws$scale_z[, -1, drop = FALSE]
  expect_error(.msr_spatial_response_simulate(bad, seed = 718), "complete finite joint")
  bad <- b; bad$draws$scale_z[1, 1] <- Inf
  expect_error(.msr_spatial_response_simulate(bad, seed = 718), "complete finite joint")
  bad <- b; bad$target <- "marginal_likelihood"
  expect_error(.msr_spatial_response_simulate(bad, seed = 718), "posterior target")
  expect_error(.msr_spatial_response_simulate(b, trials = 2, seed = 718), "only to binomial")
  bad <- b; bad$draws$beta[, 1] <- 1000
  expect_warning(expect_error(.msr_spatial_response_simulate(bad, seed = 718),
    "Nonfinite spatial response replicate"), "NAs produced")
  expect_identical(.Random.seed, before)
})

test_that("spatial replicates retain fixed references and approximation provenance", {
  skip_if_not_installed("terra")
  f <- refit_spatial_fixture(bin = TRUE)
  f$data$key <- rownames(f$data); xy <- terra::crds(f$points)
  f$data$east <- xy[, 1]; f$data$north <- xy[, 2]
  b <- .msr_spatial_scale_prepare(y ~ mean + rough + fixed + z + offset(log(exposure)),
    f$inputs, f$data, scales = list(mean = ~1, rough = ~1), family = "poisson",
    row_id = "key", priors = list(scale_median = 20), covariate = "rough",
    coordinates = c("east", "north"), knots = rbind(c(60, 60), c(150, 60), c(100, 150)),
    range = 50, field_sd_prior = .4)
  b$parameters$beta[] <- .2
  b$parameters$scale_z[] <- c(.3, -.2, .4)
  b <- fake_bayesian_draws(b, 4)
  b$inference <- "joint Gaussian approximation"
  b$inference_metadata <- list(integration = list(method = "joint_laplace"))
  b$importance <- list(resample_ids = 1:4)
  nd <- f$data[c(8, 2, 8), ]; nd$y <- nd$fixed <- NULL
  before <- list(b, f$inputs, .Random.seed)
  expected <- spatial_replicate_oracle(b, nd, f$inputs, 3, 718, include_fixed = TRUE)
  actual <- .msr_spatial_response_simulate(b, 3, seed = 718, newdata = nd,
    kernel_inputs = f$inputs, include_fixed = TRUE)
  expect_equal(unname(as.matrix(actual)), expected, tolerance = 0)
  expect_identical(list(b, f$inputs, .Random.seed), before)
  contract <- attr(actual, "prediction_contract")
  expect_identical(contract$method, b$inference)
  expect_identical(contract$integration$method, "joint_laplace")
  expect_identical(contract$integration$prediction_vectors, "stored_joint_resamples")
  expect_identical(contract$uncertainty$calibration, "not_established")
  expect_error(.msr_spatial_response_simulate(b, newdata = nd, kernel_inputs = f$inputs), "fixed")
})

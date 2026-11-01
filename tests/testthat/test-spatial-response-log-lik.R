spatial_density_oracle <- function(x, p = x$parameters) {
  eta <- spatial_response_oracle(x, x$object$specification$data, p)
  s <- x$object$specification
  switch(s$family,
    gaussian = dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE),
    poisson = dpois(s$y, exp(eta), log = TRUE),
    binomial = dbinom(s$y, s$trials, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(s$y, mu = exp(eta), size = exp(p$log_dispersion), log = TRUE))
}

test_that("conditional spatial response densities retain constants, priors and original rows", {
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
    for (metric in c("mean", "sq")) {
      x <- spatial_response_fixture(kernel, family, metric)
      b <- x$object; p <- x$parameters; s <- b$specification
      before <- list(b, p, .Random.seed)
      expected <- spatial_density_oracle(x)
      d <- .msr_spatial_response_log_lik(b, p)
      expect_identical(dim(d), c(1L, length(s$y)))
      expect_equal(as.numeric(d), expected, tolerance = 1e-12)
      expect_identical(attr(d, "observation_ids"), s$observation_ids)
      expect_identical(attr(d, "training_rows"), s$rows)
      # Independently account for every normalized prior and coordinate Jacobian.
      pr <- s$priors; sm <- s$scales$habitat
      log_prior <- sum(dnorm(p$beta, 0, pr$beta_sd, log = TRUE)) +
        sum(dnorm(p$scale_beta, sm$prior_mean, sm$prior_sd, log = TRUE)) +
        sum(dnorm(p$scale_z, 0, 1, log = TRUE)) +
        sum(dnorm(exp(p$log_scale_sd), 0, pr$scale_group_sd, log = TRUE) + log(2) + p$log_scale_sd) +
        sum(dnorm(p$log_dispersion, pr$dispersion_log_mean, pr$dispersion_log_sd, log = TRUE))
      expect_equal(.msr_spatial_scale_objective(b)(p), -sum(d) - log_prior, tolerance = 1e-10)
      zero <- p; zero$scale_z[] <- 0
      expect_equal(as.numeric(.msr_spatial_response_log_lik(b, zero)),
        spatial_density_oracle(x, zero), tolerance = 1e-12)
      zero$beta[3] <- 0
      expect_equal(as.numeric(.msr_spatial_response_log_lik(b, zero)),
        spatial_density_oracle(x, zero), tolerance = 1e-12)
      expect_identical(list(b, p, .Random.seed), before)
    }
})

test_that("conditional densities propagate complete vectors without claiming marginal selection", {
  x <- spatial_response_fixture(metric = "sq"); b <- fake_bayesian_draws(x$object, 4)
  vectors <- lapply(1:4, function(i) {
    p <- x$parameters; p$beta <- p$beta + i / 10; p$scale_z <- p$scale_z * i
    p$log_scale_sd <- p$log_scale_sd + i / 20; p
  })
  b$draws <- lapply(names(b$parameters), function(nm)
    matrix(unlist(lapply(vectors, `[[`, nm)), nrow = 4L,
      ncol = length(vectors[[1]][[nm]]), byrow = TRUE))
  names(b$draws) <- names(b$parameters)
  b$inference <- "joint Gaussian approximation"
  b$inference_metadata <- list(integration = list(method = "joint_laplace"))
  b$importance <- list(resample_ids = c(8L, 2L, 8L, 5L))
  b$diagnostics <- list(passed = FALSE, marker = "original numerical failure")
  expected <- t(vapply(vectors, function(p) spatial_density_oracle(x, p), numeric(nrow(x$fixture$data))))
  before <- list(b, .Random.seed)
  d <- .msr_spatial_response_log_lik(b)
  expect_equal(as.numeric(d), as.numeric(expected), tolerance = 1e-12)
  subset <- .msr_spatial_response_log_lik(b, draws = 2)
  expect_equal(as.numeric(subset), as.numeric(expected[c(1, 4), ]), tolerance = 1e-12)
  expect_identical(attr(subset, "draw_ids"), c(1L, 4L))
  contract <- attr(d, "likelihood_contract")
  expect_identical(contract$target, "conditional_training_response_log_density")
  expect_identical(contract$parameter_source, "stored_joint_resamples")
  expect_identical(contract$method, b$inference)
  expect_identical(contract$vector_integration$method, "joint_laplace")
  expect_identical(contract$raw_diagnostics, b$diagnostics)
  expect_identical(contract$population_dimension, 6L)
  expect_identical(contract$latent_dimension, 3L)
  expect_identical(contract$parameter_blocks, lengths(b$specification$parameters))
  expect_identical(contract$marginalization, "none")
  expect_false(contract$metadata_establishes_eligibility)
  expect_false(contract$IC_eligible)
  expect_identical(contract$uncertainty$intervals, "not_provided")
  fixed <- attr(.msr_spatial_response_log_lik(b, vectors[[1]]), "likelihood_contract")
  expect_identical(fixed$parameter_source, "supplied_fixed_vector")
  expect_identical(fixed$uncertainty$parameters, "fixed_at_supplied_vector")
  expect_identical(list(b, .Random.seed), before)
  bad <- b; bad$draws$beta[2, 1] <- NA_real_
  expect_error(.msr_spatial_response_log_lik(bad, draws = 2), "complete finite joint")
  bad <- b; bad$draws$scale_z <- bad$draws$scale_z[, -1, drop = FALSE]
  expect_error(.msr_spatial_response_log_lik(bad), "complete finite joint")
})

test_that("grouped trials and shared response groups retain their conditional density", {
  x <- spatial_response_fixture(family = "binomial", metric = "sq")
  f <- x$fixture; f$data$success <- seq_len(nrow(f$data)) %% 3
  total <- rep(c(3, 4), length.out = nrow(f$data)); f$data$failure <- total - f$data$success
  b <- .msr_spatial_scale_prepare(cbind(success, failure) ~ habitat + temp + offset(log(effort)),
    f$inputs, f$data, scales = list(habitat = ~1), family = "binomial", row_id = "key",
    priors = list(scale_median = 200), covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7, field_sd_prior = .4)
  x$object <- b
  expect_equal(as.numeric(.msr_spatial_response_log_lik(b, x$parameters)), spatial_density_oracle(x), tolerance = 1e-12)
  expect_identical(b$specification$trials, as.numeric(total))
  x <- spatial_response_fixture(family = "gaussian")
  f <- x$fixture
  b <- .msr_spatial_scale_prepare(count ~ habitat + temp + offset(log(effort)) + (1 | site),
    f$inputs, f$data, scales = list(habitat = ~1), family = "gaussian", row_id = "key",
    priors = list(scale_median = 200), covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7, field_sd_prior = .4)
  p <- b$parameters; p$beta <- x$parameters$beta
  p$scale_z <- x$parameters$scale_z; p$log_scale_sd <- x$parameters$log_scale_sd
  p$response_z[] <- seq(-.4, .3, length.out = length(p$response_z)); p$log_response_sd[] <- log(.2)
  x$object <- b; eta <- spatial_response_oracle(x, f$data, p)
  term <- b$specification$response$random[[1]]
  eta <- eta + .2 * p$response_z[term$z_index][term$index]
  d <- .msr_spatial_response_log_lik(b, p)
  expect_equal(as.numeric(d), dnorm(b$specification$y, eta, exp(p$log_dispersion), log = TRUE), tolerance = 1e-12)
  expect_identical(attr(d, "likelihood_contract")$population_dimension, 7L)
  expect_identical(attr(d, "likelihood_contract")$latent_dimension, 11L)
  expect_false("year" %in% names(b$specification$data))
})

test_that("density and selection boundaries remain explicit", {
  x <- spatial_response_fixture(); b <- x$object; p <- x$parameters
  expect_error(.msr_spatial_response_log_lik(list(), p), "private spatial")
  expect_error(.msr_spatial_response_log_lik(b), "no posterior draws")
  expect_error(.msr_spatial_response_log_lik(b, p, draws = 2), "either parameters or draws")
  for (bad in list(p[-1], c(p, list(extra = 0)), c(p, p[1])))
    expect_error(.msr_spatial_response_log_lik(b, bad), "complete finite original")
  bad <- p; bad$beta[1] <- NA_real_
  expect_error(.msr_spatial_response_log_lik(b, bad), "complete finite original")
  bad <- p; bad$scale_z <- matrix(bad$scale_z)
  expect_error(.msr_spatial_response_log_lik(b, bad), "complete finite original")
  expect_error(.msr_spatial_response_log_lik(spatial_response_fixture(groups = TRUE)$object, p),
    "no ordinary scale-group")
  for (fun in list(aic_tab, bic_tab)) {
    expect_error(fun(list(b)), "Private spatial prototypes")
    likelihood <- b; likelihood$target <- "marginal_likelihood"; likelihood$specification$priors <- NULL
    expect_error(fun(list(likelihood)), "Private spatial prototypes")
  }
  x <- spatial_response_fixture(family = "poisson"); p <- x$parameters; p$beta <- c(1000, 0, 0)
  d <- .msr_spatial_response_log_lik(x$object, p)
  expect_true(all(d == -Inf))
  expect_identical(attr(d, "likelihood_contract")$zero_density_count, length(x$object$specification$y))
  b$specification$y <- .msr_bayes_components(b$specification, x$parameters)$eta
  p <- x$parameters; p$log_dispersion <- -1000
  expect_error(.msr_spatial_response_log_lik(b, p), "Invalid conditional response log density")
})

test_that("training densities preserve mixed fixed metrics and frozen references", {
  skip_if_not_installed("terra")
  for (kernel in c("gaussian", "exp")) for (field in c("mean", "rough")) {
    f <- refit_spatial_fixture(kernel = kernel, bin = TRUE)
    f$data$key <- rownames(f$data); xy <- terra::crds(f$points)
    f$data$east <- xy[, 1]; f$data$north <- xy[, 2]
    b <- .msr_spatial_scale_prepare(y ~ mean + rough + fixed + z + offset(log(exposure)),
      f$inputs, f$data, scales = list(mean = ~1, rough = ~1), family = "poisson",
      row_id = "key", priors = list(scale_median = 20), covariate = field,
      coordinates = c("east", "north"), knots = rbind(c(60, 60), c(150, 60), c(100, 150)),
      range = 50, field_sd_prior = .4)
    p <- b$parameters; p$beta[] <- c(.2, .1, -.2, .4, .5); p$scale_z[] <- c(.3, -.2, .4)
    before <- list(b, f$inputs, p, .Random.seed)
    d <- .msr_spatial_response_log_lik(b, p)
    expected <- dpois(b$specification$y,
      exp(.msr_bayes_components(b$specification, p)$eta), log = TRUE)
    expect_equal(as.numeric(d), expected, tolerance = 1e-12)
    expect_identical(attr(d, "likelihood_contract")$population_dimension, 8L)
    expect_identical(list(b, f$inputs, p, .Random.seed), before)
  }
})

test_that("repeated training observations and unpenalized targets keep separate meanings", {
  x <- spatial_response_fixture(); f <- x$fixture; p <- x$parameters
  rows <- c(1:8, 4L, 2L, 4L); f$data <- f$data[rows, ]
  b <- .msr_spatial_scale_prepare(count ~ habitat + temp + offset(log(effort)),
    f$inputs, f$data, scales = list(habitat = ~1), family = "gaussian", row_id = "key",
    priors = list(scale_median = 200), covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7, field_sd_prior = .4)
  x$object <- b; x$fixture <- f
  d <- .msr_spatial_response_log_lik(b, p)
  expect_identical(ncol(d), 11L)
  expect_equal(as.numeric(d), spatial_density_oracle(x, p), tolerance = 1e-12)
  expect_identical(attr(d, "observation_ids"), as.character(f$data$key))
  expect_identical(as.numeric(d[1, c(9, 11)]), rep(as.numeric(d[1, 4]), 2))
  f <- spatial_response_fixture()$fixture
  m <- .msr_spatial_scale_prepare(count ~ habitat + temp + offset(log(effort)),
    f$inputs, f$data, scales = list(habitat = ~1), family = "gaussian", row_id = "key",
    covariate = "habitat", coordinates = c("east", "north"),
    knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2), range = .7,
    target = "marginal_likelihood", reference_scale = 200)
  d <- .msr_spatial_response_log_lik(m, p)
  expect_equal(.msr_spatial_scale_objective(m)(p), -sum(d) - sum(dnorm(p$scale_z, log = TRUE)), tolerance = 1e-10)
  expect_null(m$specification$priors)
  expect_identical(attr(d, "likelihood_contract")$source_target, "marginal_likelihood")
  expect_false(attr(d, "likelihood_contract")$IC_eligible)
  expect_error(.msr_spatial_response_log_lik(fake_bayesian_draws(m, 2)), "posterior target")
})

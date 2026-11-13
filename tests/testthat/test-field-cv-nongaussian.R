test_that("Poisson Laplace block scores approximate independent quadrature", {
  skip_if_not_installed("glmmTMB")
  yt <- c(0, 2); ya <- c(0, 1)
  eta_t <- c(.3, -.2); eta_a <- c(-.1, .4)
  xy_t <- matrix(0, nrow = 2, ncol = 2)
  xy_a <- cbind(c(1, 1), 0)
  field_sd <- .25; range <- 3
  rho <- exp(-1/range)
  weight <- function(z) stats::dnorm(z) *
    stats::dpois(yt[1L], exp(eta_t[1L] + field_sd * z)) *
    stats::dpois(yt[2L], exp(eta_t[2L] + field_sd * z))
  denominator <- stats::integrate(weight, -9, 9, rel.tol = 1e-10)$value
  integration <- function(z1) {
    stats::integrate(function(z2) {
      latent <- field_sd * (rho * z1 + sqrt(1 - rho^2) * z2)
      stats::dnorm(z2) * stats::dpois(ya[1L], exp(eta_a[1L] + latent)) *
        stats::dpois(ya[2L], exp(eta_a[2L] + latent))
    }, -9, 9, rel.tol = 1e-10)$value
  }
  numerator <- stats::integrate(function(z) weight(z) * vapply(z, integration, 0),
    -9, 9, rel.tol = 1e-10)$value
  exact_moment <- function(order) exp(order * eta_a + order^2 * field_sd^2 * (1 - rho^2)/2) *
    stats::integrate(function(z) weight(z) * exp(order * field_sd * rho * z),
      -9, 9, rel.tol = 1e-10)$value / denominator
  mean <- exact_moment(1)
  variance <- mean + exact_moment(2) - mean^2
  score <- function(rows, tilted = FALSE) .msr_field_laplace_block(yt, ya[rows], eta_t, eta_a[rows],
    xy_t, xy_a[rows, , drop = FALSE], "exponential", field_sd, range, "poisson", 1,
    tilted_log_moments = tilted)
  joint <- score(1:2)
  tilted <- score(1:2, TRUE)
  expect_identical(tilted$joint, joint$joint)
  expect_equal(tilted$prediction, mean, tolerance = 1e-4)
  expect_equal(tilted$observation_predictive_SE^2, variance, tolerance = 1e-4)
  # A small-field reference checks approximation error to a thousandth of a nat.
  expect_lt(abs(joint$joint$neg_log_score + log(numerator / denominator)), 1e-3)
  expect_lt(max(abs(joint$prediction - mean)), .005)
  marginal <- sum(vapply(1:2, function(i) score(i)$joint$neg_log_score, 0))
  expect_gt(abs(marginal - joint$joint$neg_log_score), 1e-4)
})

test_that("non-Gaussian field blocks recover independent normalized likelihoods at zero SD", {
  skip_if_not_installed("glmmTMB")
  train <- cbind(c(0, 1, 1), 0)
  test <- cbind(c(3, 3), 0)
  eta_t <- c(.1, -.2, .3); eta_a <- c(-.4, .2)
  for (family in c("poisson", "binomial", "nbinom2", "Gamma")) {
    y <- switch(family, poisson = c(0, 1, 5, 2, 1), binomial = c(0, 1, 0, 1, 0),
      nbinom2 = c(0, 1, 5, 2, 1), Gamma = c(.1, 1, 5, .2, 1))
    dispersion <- 2.8
    score <- function(sd) .msr_field_laplace_block(y[1:3], y[4:5], eta_t, eta_a,
      train, test, "gaussian", sd, 1.8, family, dispersion, tilted_log_moments = TRUE)
    result <- score(0)
    mean <- if (family == "binomial") stats::plogis(eta_a) else exp(eta_a)
    variance <- switch(family, poisson = mean, binomial = mean * (1 - mean),
      nbinom2 = mean + mean^2/dispersion, Gamma = mean^2/dispersion)
    loglik <- switch(family, poisson = stats::dpois(y[4:5], mean, log = TRUE),
      binomial = stats::dbinom(y[4:5], 1, mean, log = TRUE),
      nbinom2 = stats::dnbinom(y[4:5], mu = mean, size = dispersion, log = TRUE),
      Gamma = stats::dgamma(y[4:5], shape = dispersion, scale = mean/dispersion, log = TRUE))
    expect_equal(result$joint$neg_log_score, -sum(loglik), tolerance = 1e-12)
    expect_equal(result$prediction, mean, tolerance = 1e-12)
    expect_equal(result$observation_predictive_SE^2, variance, tolerance = 1e-12)
    expect_equal(result$fitted_training_field_variance + result$new_location_field_variance, c(0, 0))
    expect_equal(score(1e-8)$joint$neg_log_score, result$joint$neg_log_score, tolerance = 1e-10)
  }
})

test_that("field moments retain representable variance in response tails", {
  binary <- .msr_field_response_moments("binomial", c(-40, 40), c(0, 0), 1)
  tail <- stats::plogis(-40)
  expect_identical(binary$mean, c(tail, 1))
  expect_true(all(binary$variance > 0))
  expect_equal(binary$variance / tail, c(1, 1), tolerance = 1e-14)
  expect_equal(binary$noise / binary$variance, c(1, 1), tolerance = 1e-14)

  poisson <- .msr_field_response_moments("poisson", 400, 0, 1)
  expect_identical(poisson$variance, exp(400))
  expect_identical(poisson$variance, poisson$noise)
  small_mixing <- .msr_field_log_response_moments("poisson", 355, 710, 1,
    log_ratio = -1e-6)
  expect_equal(log(small_mixing$variance), 710 + log(-expm1(-1e-6)),
    tolerance = 1e-13)
  for (family in c("nbinom2", "Gamma")) {
    answer <- .msr_field_response_moments(family, 400, 0, exp(600))
    expected <- exp(200) + if (family == "nbinom2") exp(400) else 0
    expect_equal(log(answer$variance), log(expected), tolerance = 1e-13)
    expect_identical(answer$variance, answer$noise)
  }
})

test_that("a discrete field block cannot report mass above one", {
  skip_if_not_installed("glmmTMB")
  local_mocked_bindings(.msr_field_laplace = function(y, ...) {
    list(log_density = length(y), mode = numeric(), iterations = 0L, max_score = 0)
  })
  expect_error(.msr_field_laplace_block(0, 0, 0, 0, matrix(c(0, 0), 1L),
    matrix(c(1, 0), 1L), "exponential", 0, 1, "poisson", 1),
    "impossible discrete block probability")
})

test_that("native non-Gaussian CV uses frozen training parameters for joint assessment", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  set.seed(412)
  locations <- as.matrix(expand.grid(x = seq(0, 70, by = 10), y = seq(0, 60, by = 10)))
  n <- nrow(locations)
  lower <- t(chol(1.3^2 * exp(-as.matrix(stats::dist(locations)) / 25)))
  field <- drop(lower %*% stats::rnorm(n))
  xy <- locations[rep(seq_len(n), each = 3), , drop = FALSE]
  d <- data.frame(x = stats::rnorm(nrow(xy)), exposure = exp(stats::runif(nrow(xy), -.3, .3)))
  ids <- sprintf("r%03d", seq_len(nrow(d)))
  rownames(d) <- rownames(xy) <- ids
  eta <- .4 + .5 * d$x + log(d$exposure) + rep(field, each = 3)
  assessment <- ids[c(1:9, 97:105)]
  buffered <- utils::tail(ids, 3L)
  analysis <- setdiff(ids, c(assessment, buffered))
  folds <- list(first = list(analysis = analysis, assessment = assessment))
  form <- y ~ x + offset(log(exposure))
  run <- function(family, data = d, splits = folds, covariance = "exponential")
    cross_validate_field(form, data, xy, splits,
      audit_dir = gdr_test_audit_path(paste0("field-cv-", family, "-")),
      crs = 32617, family = family, covariance = covariance, keep_fits = TRUE)
  for (family in c("poisson", "binomial", "nbinom2", "Gamma")) {
    d$y <- switch(family, poisson = stats::rpois(nrow(d), exp(eta)),
      binomial = stats::rbinom(nrow(d), 1, stats::plogis(eta)),
      nbinom2 = stats::rnbinom(nrow(d), mu = exp(eta), size = 4),
      Gamma = stats::rgamma(nrow(d), shape = 4, scale = exp(eta)/4))
    cv <- run(family)
    expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = " "))
    model <- cv$fits$first
    joint <- cv$joint_scores$first
    expect_equal(cv$selection_records$first$K, if (family %in% c("nbinom2", "Gamma")) 5L else 4L)
    expect_equal(joint$analysis_log_density, as.numeric(stats::logLik(model)), tolerance = 1e-7)

    # The native objective on both blocks gets exactly the training parameters.
    native_data <- d[c(analysis, assessment), ]
    native_data$position <- glmmTMB::numFactor(xy[rownames(native_data), ])
    native_data$group <- factor("shared")
    native_family <- stats::family(model)
    design <- glmmTMB::glmmTMB(y ~ x + offset(log(exposure)) +
      exp(position + 0 | group), family = native_family, data = native_data, doFit = FALSE)
    joint_objective <- glmmTMB::fitTMB(design, doOptim = FALSE)
    expect_identical(names(joint_objective$par), names(model$fit$par))
    expect_equal(joint$combined_log_density, -as.numeric(joint_objective$fn(model$fit$par)), tolerance = 1e-7)
    expect_equal(cv_tab(list(field = cv), "neg_log_score")$Score,
      joint$neg_log_score/length(assessment))
    contract <- cv$prediction_contracts$first
    expect_identical(contract$integration$assessment, "joint_Laplace_marginal_likelihood_ratio")
    expect_identical(contract$integration$field, if (family == "binomial")
      "training_posterior_Laplace_Gaussian_moments" else
        "training_posterior_exponentially_tilted_Laplace_moments")
    expect_identical(contract$response_spatial_field$field_variance_scale, "linear_predictor")
    expect_true(all(cv$predictions$observation_predictive_SE > 0))
    if (family == "poisson") {
      expect_equal(cv$predictions$observation_noise_variance, cv$predictions$prediction)
      altered <- d
      altered[assessment, "y"] <- altered[assessment, "y"] + 5L
      altered[buffered, "y"] <- altered[buffered, "y"] + 200L
      changed <- run(family, altered)
      expect_true(changed$complete)
      expect_identical(changed$training$first$parameters, cv$training$first$parameters)
      expect_identical(changed$predictions$prediction, cv$predictions$prediction)
      expect_gt(changed$joint_scores$first$neg_log_score, joint$neg_log_score)
    }
    if (family == "binomial")
      expect_equal(cv$predictions$observation_predictive_SE^2,
        cv$predictions$prediction * (1 - cv$predictions$prediction))
    if (family == "nbinom2") {
      count <- list(poisson = run("poisson"), nbinom2 = cv)
      expect_equal(nrow(cv_tab(count, "neg_log_score")), 2L)
      selected <- select_cv(count, criterion = "AIC")
      expect_true(selected$complete)
      expect_identical(selected$selection$chosen[["first"]], "nbinom2")
      expect_identical(selected$settings$family, "count")
      expect_identical(selected$settings$fold_settings$first$family, "nbinom2")
      expect_equal(cv_tab(list(selected = selected), "neg_log_score")$Score,
        joint$neg_log_score/length(assessment))
    }
  }
})


test_that("Poisson field CV accepts Gaussian correlation and refits ecological bandwidths", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("TMB")
  f <- native_response_field_fixture("poisson")
  d <- f$data
  xy <- terra::crds(f$points)
  rownames(xy) <- rownames(d)
  crs <- terra::crs(f$points)
  folds <- spatial_cv_folds(xy, block_size = 40, k = 2L, crs = crs, ids = rownames(d))[1L]
  form <- y ~ mean + z + offset(log(exposure))
  run <- function(covariance, inputs = NULL) cross_validate_field(form, d, xy, folds,
    crs = crs, family = "poisson", covariance = covariance, kernel_inputs = inputs,
    keep_fits = TRUE, audit_dir = gdr_test_audit_path("field-cv-poisson-kernel-"))
  gaussian <- run("gaussian")
  expect_true(gaussian$complete, info = paste(unlist(gaussian$errors), collapse = " "))
  expect_identical(gaussian$prediction_contracts$fold1$response_spatial_field$covariance, "gaussian")
  scaled <- run("exponential", f$inputs)
  expect_true(scaled$complete, info = paste(unlist(scaled$errors), collapse = " "))
  expect_s3_class(scaled$fits$fold1, "multiScaleR")
  expect_identical(rownames(scaled$fits$fold1$kernel_inputs$kernel_dat), folds$fold1$analysis)
  expect_equal(scaled$selection_records$fold1$ecological_K, 1L)
  expect_true(is.finite(scaled$training$fold1$scale_est["mean", "Mean"]))
  expect_gt(scaled$training$fold1$refits, 0L)
})

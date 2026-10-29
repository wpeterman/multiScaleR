test_that("response effects and hierarchical scales have separate formula designs", {
  f <- bayesian_scale_fixture()
  fit <- prepare_bayesian_fixture(f)
  s <- fit$specification
  expect_s3_class(fit, "multiScaleR_bayes")
  expect_identical(fit$inference, "prepared")
  expect_identical(s$beta_names, c("(Intercept)", "temp", "habitat"))
  expect_identical(colnames(s$scales$habitat$X), "(Intercept)")
  expect_length(s$scales$habitat$random$year$levels, 4)
  expect_length(s$response$random$site$levels, 8)
  expect_equal(s$offset, log(f$data$effort))
  expect_equal(s$parameters$scale_beta, log(0.2))
  fixed <- prepare_bayesian_fixture(f, scales = list(habitat = ~year))
  expect_equal(ncol(fixed$specification$scales$habitat$X), 4)
  expect_length(fixed$parameters$scale_z, 0)
  shared <- multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", method = "prepare")
  expect_length(shared$parameters$scale_beta, 1)
  no_intercept <- multiScale_bayes(count ~ 0 + habitat + (1 | site), f$inputs,
    f$data, row_id = "key", method = "prepare")
  expect_identical(no_intercept$specification$beta_names, "habitat")
})

test_that("stacked observations preserve annual and static landscape identities", {
  f <- bayesian_scale_fixture()
  original <- prepare_bayesian_fixture(f)
  perm <- rev(seq_len(nrow(f$data)))
  f$data <- f$data[perm, ]
  shuffled <- prepare_bayesian_fixture(f)
  expect_equal(shuffled$specification$V$habitat,
                original$specification$V$habitat[perm, ])
  expect_equal(shuffled$specification$offset, original$specification$offset[perm])
  static <- bayesian_scale_fixture(annual = FALSE)
  reused <- prepare_bayesian_fixture(static)
  expect_equal(length(unique(reused$specification$rows)), 8)
  expect_equal(length(reused$specification$rows), 32)
  expect_equal(reused$specification$V$habitat[1, ],
                reused$specification$V$habitat[9, ])
  bad <- f
  bad$data$key[1] <- "unknown"
  expect_error(prepare_bayesian_fixture(bad), "Observation IDs")
  bad <- f
  bad$inputs$binned$point_ids[2] <- bad$inputs$binned$point_ids[1]
  expect_error(prepare_bayesian_fixture(bad), "uniquely")
  bad <- f
  bad$data$temp[1] <- NA
  expect_error(prepare_bayesian_fixture(bad), "missing")
})

test_that("kernel means match bin sums and retain fixed pooled standardization", {
  for (kernel in c("gaussian", "exp")) {
    f <- bayesian_scale_fixture(kernel = kernel)
    fit <- prepare_bayesian_fixture(f)
    s <- fit$specification
    p <- fit$parameters
    p$scale_z <- seq(-0.7, 0.7, length.out = length(p$scale_z))
    ans <- multiScaleR:::.msr_bayes_components(s, p)
    raw <- vapply(seq_len(nrow(f$data)), function(i) {
      ds <- f$inputs$binned$bin_rep / (ans$scales$habitat[i] / 1000)
      w <- if (kernel == "gaussian") exp(-ds^2 / 2) else exp(-ds)
      sum(w * s$V$habitat[i, ]) / sum(w * s$C$habitat[i, ])
    }, 0)
    expect_equal(ans$kernel$habitat, (raw - s$center[1]) / s$spread[1],
                   ignore_attr = TRUE, tolerance = 1e-10)
    expect_equal(length(unique(ans$scales$habitat)), 4)
    small <- p
    small$scale_beta[1] <- log(1e-6)
    expect_true(all(is.finite(multiScaleR:::.msr_bayes_components(s, small)$kernel$habitat)))
    expect_identical(s$center, fit$specification$center)
  }
})

test_that("full posterior includes likelihood constants, priors and SD Jacobians", {
  f <- bayesian_scale_fixture()
  fit <- prepare_bayesian_fixture(f)
  s <- fit$specification
  p <- fit$parameters
  p$beta <- c(0.3, -0.1, 0.7)
  p$scale_z <- seq(-0.2, 0.3, length.out = length(p$scale_z))
  p$response_z <- seq(-0.4, 0.4, length.out = length(p$response_z))
  mu <- exp(multiScaleR:::.msr_bayes_components(s, p)$eta)
  expected <- -sum(dnbinom(s$y, mu = mu, size = exp(p$log_dispersion), log = TRUE)) -
    sum(dnorm(p$beta, 0, 2.5, log = TRUE)) -
    dnorm(p$scale_beta, log(.2), 1, log = TRUE) -
    sum(dnorm(p$scale_z, log = TRUE)) - sum(dnorm(p$response_z, log = TRUE)) -
    sum(dnorm(exp(p$log_scale_sd), 0, .5, log = TRUE) + log(2) + p$log_scale_sd) -
    sum(dnorm(exp(p$log_response_sd), 0, 1, log = TRUE) + log(2) + p$log_response_sd) -
    sum(dnorm(p$log_dispersion, 0, 1, log = TRUE))
  expect_equal(multiScaleR:::.msr_bayes_objective(s)(p), expected, tolerance = 1e-11)
  centered <- prepare_bayesian_fixture(f, priors = list(scale_median = c(habitat = 350)))
  expect_equal(centered$parameters$scale_beta, log(.35))
  expect_error(prepare_bayesian_fixture(f, priors = list(typo = 1)), "documented")
})

test_that("unsupported designs fail explicitly instead of changing the model", {
  f <- bayesian_scale_fixture()
  expect_error(multiScale_bayes(count ~ habitat * year, f$inputs, f$data,
    row_id = "key", method = "prepare"), "bare additive")
  expect_error(multiScale_bayes(count ~ habitat + (temp | site), f$inputs, f$data,
    row_id = "key", method = "prepare"), "random intercepts")
  expect_error(prepare_bayesian_fixture(f, scales = list(habitat = ~0 + temp)), "factor levels")
  expect_error(prepare_bayesian_fixture(f, scales = list(other = ~1)), "every kernel")
  expect_error(prepare_bayesian_fixture(f, scales = list(habitat = ~habitat)), "cannot use")
  expect_error(prepare_bayesian_fixture(f, scales = list(habitat = ~1 + (1 | temp))),
    NA)
  f$inputs$kernel <- "fixed"
  expect_error(prepare_bayesian_fixture(f), "Gaussian or exponential")
})

test_that("posterior summaries distinguish scales, response effects and fitted means", {
  f <- bayesian_scale_fixture()
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f))
  fit$draws$beta[, 1] <- 0.3
  fit$draws$beta[, 3] <- 0.6
  fit$draws$response_z[,] <- 0.5
  fit$draws$scale_z[,] <- rep(seq(-.5, .5, length.out = 4), each = 8)
  tab <- posterior_scales(fit)
  expect_equal(nrow(tab), 4)
  expect_equal(tab$group, paste0("year=y", 1:4))
  expect_true(all(tab$lower == tab$upper))
  pred <- predict(fit, summary = FALSE)
  zero <- predict(fit, re.form = NA, summary = FALSE)
  expect_equal(log(pred / zero), matrix(exp(fit$parameters$log_response_sd) * .5, 8, 32))
  expect_equal(predict(fit)$mean, colMeans(pred))
  expect_equal(nrow(predict(fit, draws = 2)), 32)
  expect_error(predict(fit, newdata = f$data), "matching binned kernel_inputs")
  expect_error(predict(fit, unknown = 1), "Unused")
  expect_error(predict(fit, re.form = ~1), "NULL or NA")
  expect_error(posterior_scales(prepare_bayesian_fixture(f)), "no posterior")
  expect_equal(nrow(summary(fit)$coefficients), 3)
  expect_equal(summary(fit)$group_sd$log_scale_group_sd$term, "habitat:year")
  restored <- unserialize(serialize(fit, NULL))
  expect_equal(predict(restored), predict(fit))
})

test_that("AD posterior values and gradients agree with an independent numeric calculation", {
  skip_if_not_installed("RTMB")
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- bayesian_scale_fixture(n_site = 4, n_year = 3)
    if (family == "binomial") f$data$count <- rep(0:1, length.out = nrow(f$data))
    fit <- prepare_bayesian_fixture(f, family)
    s <- fit$specification
    obj <- RTMB::MakeADFun(multiScaleR:::.msr_bayes_objective(s, ad = TRUE),
                            fit$parameters, silent = TRUE)
    p <- obj$par + seq(-.03, .03, length.out = length(obj$par))
    unpack <- function(x) {
      ends <- cumsum(lengths(fit$parameters))
      starts <- ends - lengths(fit$parameters)
      lapply(seq_along(ends), function(j) x[starts[j] + seq_len(ends[j] - starts[j])])
    }
    fn <- function(x) {
      par <- unpack(x)
      names(par) <- names(fit$parameters)
      multiScaleR:::.msr_bayes_objective(s)(par)
    }
    expect_equal(obj$fn(p), unname(fn(p)), tolerance = 1e-9, info = family)
    h <- 1e-5
    fd <- vapply(seq_along(p), function(i) {
      hi <- lo <- p
      hi[i] <- hi[i] + h
      lo[i] <- lo[i] - h
      (fn(hi) - fn(lo)) / (2 * h)
    }, 0)
    expect_equal(as.numeric(obj$gr(p)), unname(fd), tolerance = 1e-5, info = family)
  }
})

test_that("posterior predictive simulation includes observation variation and preserves seeded RNG", {
  f <- bayesian_scale_fixture()
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f))
  set.seed(812)
  previous <- .Random.seed
  a <- simulate(fit, nsim = 10, seed = 456)
  expect_identical(.Random.seed, previous)
  expect_identical(a, simulate(fit, nsim = 10, seed = 456))
  expect_equal(dim(a), c(32, 10))
  expect_true(all(as.matrix(a) >= 0))
  expect_true(any(apply(a, 1, stats::var) > 0))
  expect_true(is.function(getS3method("predict", "multiScaleR_bayes", optional = TRUE)))
  expect_true(is.function(getS3method("simulate", "multiScaleR_bayes", optional = TRUE)))
})

test_that("multiple predictors retain independent scale hierarchies and parameter labels", {
  f <- bayesian_scale_fixture()
  f$inputs$binned$vsum$forest <- f$inputs$binned$vsum$habitat[, 15:1]
  f$inputs$binned$csum$forest <- f$inputs$binned$csum$habitat[, 15:1]
  fit <- multiScale_bayes(count ~ habitat + forest + temp, f$inputs, f$data,
    scales = list(habitat = ~1, forest = ~1 + (1 | year)),
    row_id = "key", method = "prepare")
  fit <- fake_bayesian_draws(fit)
  expect_equal(nrow(posterior_scales(fit)), 5)
  ans <- summary(fit)
  expect_equal(ans$group_sd$log_scale_group_sd$term, "forest:year")
  expect_equal(ans$scale_coefficients$habitat$median, log(200))
  expect_equal(names(ans$scale_coefficients), c("habitat", "forest"))
})

test_that("grouped binomial encodings preserve trials and posterior prediction units", {
  f <- bayesian_scale_fixture()
  f$data$success <- rep(1:3, length.out = nrow(f$data))
  f$data$failure <- rep(2:4, length.out = nrow(f$data))
  fit <- multiScale_bayes(cbind(success, failure) ~ habitat, f$inputs, f$data,
    family = "binomial", row_id = "key", method = "prepare")
  expect_equal(fit$specification$trials, f$data$success + f$data$failure)
  fit <- fake_bayesian_draws(fit)
  expect_equal(predict(fit)$mean, rep(.5, nrow(f$data)))
  reps <- simulate(fit, 4, seed = 4)
  expect_true(all(as.matrix(reps) <= fit$specification$trials))
})

test_that("saved model specifications do not retain unrelated calling-frame data", {
  f <- bayesian_scale_fixture()
  prepared <- local({
    unrelated_payload <- raw(2e6)
    multiScale_bayes(count ~ habitat, f$inputs, f$data,
      row_id = "key", method = "prepare")
  })
  expect_lt(length(serialize(prepared, NULL)), 250000)
  restored <- unserialize(serialize(prepared, NULL))
  expect_equal(multiScaleR:::.msr_bayes_objective(restored$specification)(restored$parameters),
    multiScaleR:::.msr_bayes_objective(prepared$specification)(prepared$parameters))
})

test_that("independent annual scale priors are invariant to the reference year", {
  f <- bayesian_scale_fixture()
  fit <- prepare_bayesian_fixture(f, scales = list(habitat = ~0 + year))
  sm <- fit$specification$scales$habitat
  expect_equal(sm$prior_mean, rep(log(.2), 4))
  expect_equal(sm$prior_sd, rep(1, 4))
  expect_length(fit$parameters$scale_z, 0)
  fit <- fake_bayesian_draws(fit)
  fit$draws$scale_beta[,] <- rep(log(c(.1, .2, .3, .4)), each = 8)
  expect_equal(posterior_scales(fit)$median, c(100, 200, 300, 400))
  expect_equal(summary(fit)$scale_coefficients$habitat$median,
                 log(c(100, 200, 300, 400)))
  g <- f
  g$data$year <- stats::relevel(g$data$year, "y3")
  other <- prepare_bayesian_fixture(g, scales = list(habitat = ~0 + year))
  other$parameters$scale_beta <- log(c(.3, .1, .2, .4))
  a <- multiScaleR:::.msr_bayes_objective(other$specification)(other$parameters)
  original <- prepare_bayesian_fixture(f, scales = list(habitat = ~0 + year))
  original$parameters$scale_beta <- log(c(.1, .2, .3, .4))
  b <- multiScaleR:::.msr_bayes_objective(original$specification)(original$parameters)
  expect_equal(a, b)
})

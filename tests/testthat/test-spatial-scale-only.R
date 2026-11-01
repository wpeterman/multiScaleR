# Spatial-only software fixtures: no calendar column or annual/response effect.
spatial_only_fixture <- function(kernel, family, metric, target = "posterior") {
  f <- bayesian_scale_fixture(n_site = 8, n_year = 1, annual = FALSE, kernel = kernel)
  f$data$year <- NULL
  angle <- 2 * pi * (as.integer(f$data$site) - 1) / 8
  f$data$east <- cos(angle); f$data$north <- sin(angle)
  if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
  if (metric == "sq") {
    d <- f$inputs$binned$bin_rep; ids <- f$inputs$binned$point_ids
    f$inputs$d_list <- rep(list(d), length(ids))
    f$inputs$raw_cov <- lapply(seq_along(ids), function(i)
      data.frame(elev = sin(d * (i + 1)) + i * d^2))
    f$inputs$raw_cov[[1L]]$elev[] <- 2
    for (i in 2:4) {
      f$inputs$d_list[[i]][2L] <- d[1L]
      f$inputs$raw_cov[[i]]$elev[1:2] <- c(0, i)
    }
    f$inputs$kernel_dat <- data.frame(row.names = ids)
    f$inputs$scale_vars <- msr_vars(habitat = surface_var("elev", "sq", weighted = TRUE))
    f$inputs$binned <- NULL
  }
  args <- list(formula = count ~ habitat + temp, inputs = f$inputs, data = f$data,
    scales = list(habitat = ~1), family = family, row_id = "key", covariate = "habitat",
    coordinates = c("east", "north"), knots = matrix(c(-.5, 0, .5, 0, .5, -.5), 3, 2),
    range = .7, target = target)
  if (target == "posterior") {
    args$priors <- list(scale_median = 200); args$field_sd_prior <- .4
  } else args$reference_scale <- 200
  list(fixture = f, args = args, object = do.call(.msr_spatial_scale_prepare, args))
}

test_that("spatial-only fields preserve target, derivatives, priors and Gaussian-block guards", {
  skip_if_not_installed("RTMB")
  cases <- 0L
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
    for (metric in c("mean", "sq")) {
      x <- spatial_only_fixture(kernel, family, metric)
      b <- x$object; spec <- b$specification; sm <- spec$scales$habitat
      m <- if (metric == "mean") spatial_only_fixture(kernel, family, metric, "marginal_likelihood")$object else NULL
      cases <- cases + 1L; label <- paste(kernel, family, metric)
      expect_false("year" %in% names(spec$data), info = label)
      expect_identical(names(sm$random), ".spatial_field")
      expect_identical(sm$random[[1L]]$sd_index, 1L)
      expect_identical(sm$random[[1L]]$z_index, 1:3)
      expect_identical(lengths(spec$parameters)[c("log_scale_sd", "scale_z", "log_response_sd", "response_z")],
        c(log_scale_sd = 1L, scale_z = 3L, log_response_sd = 0L, response_z = 0L))
      expect_identical(spec$priors$scale_group_sd, .4)
      p <- b$parameters; p$beta <- c(.2, -.1, .5)
      p$scale_z <- c(.4, -.2, .3); p$log_scale_sd <- log(.35)
      before <- list(b, m, p, .Random.seed)
      pr <- spec$priors
      penalty <- -sum(dnorm(p$beta, 0, pr$beta_sd, log = TRUE)) -
        sum(dnorm(p$scale_beta, sm$prior_mean, sm$prior_sd, log = TRUE)) -
        sum(dnorm(exp(p$log_scale_sd), 0, pr$scale_group_sd, log = TRUE) + log(2) + p$log_scale_sd) -
        sum(dnorm(p$log_dispersion, pr$dispersion_log_mean, pr$dispersion_log_sd, log = TRUE))
      neutral_value <- if (metric == "mean") .msr_spatial_scale_objective(m)(p) else
        sum(.msr_marginal_response_nll(spec, p)) - sum(dnorm(p$scale_z, log = TRUE))
      if (metric == "sq") expect_error(spatial_only_fixture(kernel, family, metric, "marginal_likelihood"), "Weighted SQ")
      expect_equal(.msr_spatial_scale_objective(b)(p) - neutral_value,
        penalty, tolerance = 1e-10, info = label)
      joint <- if (metric == "mean") m else b; joint$parameters <- p
      ad <- .msr_spatial_scale_adfun(joint)
      expect_equal(as.numeric(ad$fn(ad$par)), .msr_spatial_scale_objective(joint)(p), tolerance = 1e-10)
      step <- 1e-5
      gradient <- vapply(seq_along(ad$par), function(j) {
        hi <- lo <- ad$par; hi[j] <- hi[j] + step; lo[j] <- lo[j] - step
        (ad$fn(hi) - ad$fn(lo)) / (2 * step)
      }, 0)
      expect_equal(as.numeric(ad$gr(ad$par)), gradient, tolerance = 2e-5, info = label)
      expect_equal(unname(ad$he(ad$par)), unname(optimHess(ad$par, ad$fn, ad$gr)), tolerance = 2e-5, info = label)
      # One field SD remains a coupled Gaussian block, not an IID year term.
      expect_error(.msr_scale_check_capability(spec, "marginal_likelihood", "factorized"),
        "Gaussian projection|Weighted SQ")
      expect_error(.msr_scale_check_capability(spec, "posterior", "marginal_laplace"), "Gaussian projection")
      zero <- p; zero$log_scale_sd <- -Inf
      base <- .msr_bayes_prepare(x$args$formula, x$args$inputs, x$args$data, x$args$scales,
        family, "key", x$args$priors)
      q <- zero; q$log_scale_sd <- q$scale_z <- numeric()
      expect_equal(.msr_marginal_response_nll(spec, zero), .msr_marginal_response_nll(base, q), tolerance = 1e-10)
      expect_identical(list(b, m, p, .Random.seed), before)
    }
  expect_identical(cases, 16L)
  cat("\nSPATIAL_ONLY_TARGET cases=16 annual_sd=0 field_sd=1 latent_dimension=3\n")
})

test_that("spatial-only prediction needs no year and propagates a coherent joint field", {
  for (metric in c("mean", "sq")) {
    x <- spatial_only_fixture("gaussian", "gaussian", metric); b <- x$object
    p <- b$parameters; p$scale_z <- c(.4, -.2, .3); p$log_scale_sd <- log(.35)
    spec <- b$specification; sm <- spec$scales$habitat; block <- sm$random[[1L]]$gaussian_block
    expected <- drop(sm$X %*% p$scale_beta) + exp(p$log_scale_sd) * drop(block$design %*% block$chol %*% p$scale_z) +
      log(spec$unit_conv)
    expect_equal(.msr_spatial_scale_predict(b, x$fixture$data, parameters = p, scale = "log"), expected, tolerance = 1e-12)
    future <- data.frame(east = c(.15, .15), north = c(.2, .2))
    before <- list(b, p, .Random.seed)
    predicted <- .msr_spatial_scale_predict(b, future, parameters = p, new_levels = "sample")
    expect_identical(unname(predicted[1L]), unname(predicted[2L]))
    with_year <- future; with_year$year <- c("arbitrary", "different")
    expect_identical(.msr_spatial_scale_predict(b, with_year, parameters = p, new_levels = "sample"), predicted)
    expect_identical(list(b, p, .Random.seed), before)
    b$draws <- lapply(p, function(v) matrix(rep(v, each = 4L), nrow = 4L))
    b$draws$scale_z <- matrix(seq(-.4, .5, length.out = 12L), 4L)
    joint <- .msr_spatial_scale_predict(b, future, scale = "log")
    expected_draws <- do.call(rbind, lapply(1:4, function(i)
      .msr_spatial_scale_predict(b, future, parameters = .msr_bayes_draw(b, i), scale = "log")))
    expect_identical(dim(joint), c(4L, 2L))
    expect_equal(joint, expected_draws, tolerance = 1e-12)
    expect_equal(joint[, 1L], joint[, 2L], tolerance = 1e-12)
  }
})

test_that("spatial-only zero-slope latent integration is normalized and SD remains unidentifiable", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (metric in c("mean", "sq")) {
    target <- if (metric == "mean") "marginal_likelihood" else "posterior"
    x <- spatial_only_fixture(kernel, "gaussian", metric, target)
    b <- x$object; b$parameters$beta <- c(.2, -.1, 0)
    b$parameters$log_scale_sd <- log(.5)
    ad <- .msr_spatial_scale_adfun(b, integrate_latent = TRUE)
    reference <- if (metric == "mean") sum(.msr_marginal_response_nll(b$specification, b$parameters)) else
      .msr_spatial_scale_objective(b)(b$parameters) + sum(dnorm(b$parameters$scale_z, log = TRUE))
    expect_equal(as.numeric(ad$fn(ad$par)), reference, tolerance = 1e-8)
    at <- names(ad$par) == "log_scale_sd"
    prior_gradient <- if (metric == "mean") 0 else exp(2 * b$parameters$log_scale_sd) / .4^2 - 1
    expect_equal(as.numeric(ad$gr(ad$par))[at], prior_gradient, tolerance = 1e-8)
    hi <- lo <- ad$par; hi[at] <- hi[at] + .5; lo[at] <- lo[at] - .5
    # Uninformative response: ML is flat; posterior variation is exactly its
    # declared half-normal field-SD prior with the original log-SD Jacobian.
    difference <- if (metric == "mean") 0 else {
      value <- function(log_sd) .5 * exp(2 * log_sd) / .4^2 - log_sd
      value(hi[at]) - value(lo[at])
    }
    expect_equal(as.numeric(ad$fn(hi) - ad$fn(lo)), as.numeric(difference), tolerance = 1e-8)
    if (metric == "mean") expect_error(.msr_spatial_scale_laplace(b), "posterior target")
  }
})

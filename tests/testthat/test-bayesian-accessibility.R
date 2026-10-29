test_that("setup checks describe method-specific requirements without fitting", {
  p <- check_bayes_setup("prepare")
  expect_s3_class(p, "multiScaleR_bayes_setup")
  expect_true(p$ready)
  expect_equal(nrow(p$packages), 0L)
  expect_output(print(p), "no optional Bayesian")
  expect_error(check_bayes_setup("sample", "importance"), "only")
  expect_error(check_bayes_setup("other"), "arg")
  expect_error(print(p, verbose = TRUE), "Unused")
  local_mocked_bindings(.msr_bayes_package_status = function(package, minimum)
    data.frame(package = package, minimum = unname(minimum), version = NA_character_,
      loadable = FALSE, meets_version = FALSE, passed = FALSE, message = "Missing"),
    .package = "multiScaleR")
  s <- check_bayes_setup("sample")
  expect_false(s$ready)
  expect_identical(s$packages$package, c("RTMB", "tmbstan", "rstan", "posterior"))
  expect_identical(check_bayes_setup("laplace")$packages$package, "RTMB")
  expect_identical(check_bayes_setup("marginal", "importance")$packages$package,
    c("RTMB", "posterior"))
})

test_that("package status distinguishes a missing dependency", {
  s <- multiScaleR:::.msr_bayes_package_status("multiScaleRnonexistentpackage", "1.0")
  expect_false(s$passed)
  expect_false(s$loadable)
  expect_match(s$message, "missing|loaded")
})

test_that("scale plots preserve posterior summaries, labels and RNG state", {
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(f), n = 5)
  fit$draws$scale_beta[, 1] <- log(c(.1, .2, .3, .4, .5))
  fit$draws$scale_z[, 1] <- c(-.3, -.2, 0, .2, .3)
  fit$diagnostics <- list(passed = TRUE)
  original <- serialize(fit, NULL); rng <- .Random.seed
  p <- plot(fit, level = .8)
  expect_s3_class(p, "ggplot")
  expected <- posterior_scales(fit, level = .8)
  expect_equal(p$data$median, expected$median)
  expect_equal(p$data$lower, expected$lower)
  expect_identical(as.character(p$data$group), expected$group)
  expect_equal(plot(fit, scale = "log")$data$median, log(plot(fit)$data$median))
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L)
  expect_identical(serialize(fit, NULL), original)
  expect_identical(.Random.seed, rng)
  size <- length(serialize(plot(fit), NULL))
  fit$unrelated_payload <- rep(1, 100000)
  expect_lte(length(serialize(plot(fit), NULL)), size + 1024)
  expect_error(plot(fit, covariate = "unknown"), "covariate")
  expect_error(plot(fit, level = 1), "level")
  expect_error(plot(fit, nonsense = 1), "Unused")
  expect_error(plot(prepare_bayesian_fixture(f)), "no posterior draws")
  fit$diagnostics$passed <- FALSE
  expect_match(plot(fit)$labels$subtitle, "Diagnostics need attention")
})

test_that("group-SD plots distinguish population variation from annual uncertainty", {
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(bayesian_scale_fixture()), n = 5)
  fit$draws$log_scale_sd[, 1] <- log(c(.1, .2, .3, .4, .5))
  p <- plot(fit, type = "group_sd", level = .8)
  expect_equal(p$data$median, c(.3, .5))
  expect_equal(p$data$lower[1], unname(quantile(c(.1, .2, .3, .4, .5), .1)))
  expect_identical(as.character(p$data$term), c("habitat:year", "site"))
  expect_identical(p$data$component, c("Log-scale group SD", "Response group SD"))
  expect_equal(nrow(ggplot2::ggplot_build(p)$data[[1]]), 2L)
  expect_error(plot(fit, type = "group_sd", covariate = "habitat"), "only")
  f <- bayesian_scale_fixture()
  shared <- fake_bayesian_draws(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", method = "prepare"))
  expect_error(plot(shared, type = "group_sd"), "no between-group")
})

test_that("printed evidence distinguishes raw checks from derived scale precision", {
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(bayesian_scale_fixture()), n = 5)
  fit$diagnostics <- list(passed = TRUE)
  original <- serialize(fit, NULL)
  expect_output(print(fit), "Raw-parameter sampling checks: passed")
  expect_output(print(fit), "diagnose_bayes_scales")
  expect_identical(serialize(fit, NULL), original)
  for (flag in list(NA, NULL)) {
    fit$diagnostics <- list(passed = flag)
    expect_output(print(fit), "checks: unavailable")
    expect_match(plot(fit)$labels$subtitle, "Stored numerical checks unavailable")
  }
  fit$diagnostics <- NULL
  expect_output(print(fit), "checks: unavailable")
  expect_match(plot(fit)$labels$subtitle, "Stored numerical checks unavailable")
  fit$diagnostics <- list(passed = FALSE)
  expect_output(print(fit), "checks: need attention")
  expect_match(plot(fit)$labels$subtitle, "Diagnostics need attention")
  fit$marginal <- list()
  expect_output(print(fit), "Marginal Laplace numerical checks")
  fit$importance <- list()
  expect_output(print(fit), "Marginal importance diagnostics")
})

test_that("multiple predictors omit absent scale SDs and retain selectable panels", {
  f <- bayesian_scale_fixture(n_site = 5, n_year = 3)
  b <- f$inputs$binned
  f$inputs$binned$vsum$forest <- cos(b$vsum$habitat / pmax(1, b$csum$habitat)) * b$csum$habitat
  f$inputs$binned$csum$forest <- b$csum$habitat
  p <- multiScale_bayes(count ~ habitat + forest, f$inputs, f$data,
    scales = list(habitat = ~1, forest = ~1 + (1 | year)), row_id = "key", method = "prepare")
  fit <- fake_bayesian_draws(p)
  expect_identical(as.character(plot(fit, type = "group_sd")$data$term), "forest:year")
  expect_identical(unique(plot(fit, covariate = "forest")$data$covariate), "forest")
  expect_equal(nrow(plot(fit)$data), 4L)
})

test_that("importance plots use direct weighted candidates by default", {
  fit <- fake_bayesian_draws(prepare_bayesian_fixture(bayesian_scale_fixture()), n = 4)
  fit$draws$scale_beta[, 1] <- log(c(.1, .2, .3, .4))
  fit$draws$log_scale_sd[, 1] <- log(c(.1, .2, .3, .4))
  fit$importance <- list(candidates = do.call(cbind, fit$draws), weights = c(.1, .1, .1, .7))
  expected <- posterior_scales(fit, estimator = "weighted")
  expect_equal(plot(fit)$data$median, expected$median)
  expect_equal(plot(fit, type = "group_sd")$data$median[1], .4)
  expect_equal(plot(fit, type = "group_sd", estimator = "draws")$data$median[1], .25)
  fit$importance$weights[1] <- -1
  expect_error(plot(fit, type = "group_sd"), "Invalid retained")
  fit$importance <- NULL
  expect_error(plot(fit, estimator = "weighted"), "importance-corrected")
  expect_error(plot(fit, type = "group_sd", estimator = "weighted"), "retained")
})

test_that("the bundled annual example is portable and uses eligible actual fits", {
  path <- system.file("extdata", "bayesian-annual-example.rds", package = "multiScaleR")
  expect_true(nzchar(path))
  e <- readRDS(path)
  expect_identical(e$format_version, 1L)
  expect_true(all(e$fit_info$passed))
  expect_true(all(e$fit_info$max_rhat <= 1.01, e$fit_info$min_tail_ess >= 400,
    e$fit_info$min_bulk_ess >= 400, e$fit_info$divergences == 0,
    e$fit_info$depth_hits == 0, e$fit_info$min_ebfmi >= .3))
  contains_pointer <- function(x) {
    if (typeof(x) == "externalptr" || inherits(x, "SpatRaster") || inherits(x, "stanfit")) return(TRUE)
    if (is.list(x)) return(any(vapply(x, contains_pointer, TRUE)))
    FALSE
  }
  expect_false(contains_pointer(e))
  expect_equal(as.integer(table(e$data$year)), e$truth$observations)
  expect_gt(length(unique(e$truth$scale)), 1L)
  f <- e$fits$pooled
  expect_null(f$sampler)
  expect_equal(nrow(f$draws$beta), 2000L)
  expect_gt(f$example_provenance$original_draws, 2000L)
  prepared <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
    e$inputs, e$data, scales = list(habitat = ~1 + (1 | year)),
    row_id = "site_year", priors = list(scale_median = 250), method = "prepare")
  for (nm in c("V", "C", "delta", "center", "spread", "X_fixed", "y"))
    expect_identical(prepared$specification[[nm]], f$specification[[nm]])
  expect_s3_class(plot(f), "ggplot")
  expect_s3_class(plot(f, type = "group_sd"), "ggplot")
  expect_equal(nrow(e$comparison), 15L)
  expect_identical(e$sensitivity$prior_sd, c(.25, .5, .75))
  before <- .Random.seed
  a <- simulate(f, nsim = 3, seed = 94)
  expect_identical(.Random.seed, before)
  expect_equal(dim(a), c(nrow(e$data), 3L))
  expect_equal(a, simulate(f, nsim = 3, seed = 94))
})

test_that("the computation example retains failures and matches its reference", {
  path <- system.file("extdata", "bayesian-computation-example.rds", package = "multiScaleR")
  expect_true(nzchar(path))
  e <- readRDS(path)
  reference <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
    package = "multiScaleR"))
  expect_identical(e$format_version, 1L)
  expect_setequal(e$availability$method, c("laplace", "marginal", "nested"))
  expect_identical(e$availability$status[e$availability$method == "marginal"], "rejected")
  expect_false(e$diagnostics$marginal$passed)
  expect_false(e$diagnostics$marginal$grid_converged)
  expect_true(e$diagnostics$nested$passed)
  expect_identical(e$provenance$priors, reference$fits$pooled$specification$priors)
  expect_false(any(grepl("^[A-Za-z]:|^/", names(e$provenance$engine_md5))))
  expect_lt(file.info(path)$size, 100000L)
  # The teaching artifact is a table, not another retained numerical population.
  expect_null(e$draws)
  expect_null(e$sampler)
  expect_null(e$candidates)
  expect_true(all(is.finite(as.matrix(e$quantiles[c("lower", "median", "upper")]))))
  expect_true(all(e$quantiles$lower <= e$quantiles$median,
    e$quantiles$median <= e$quantiles$upper))
  sd <- subset(e$quantiles, method == "nuts" & type == "scale_sd")
  expect_equal(unname(unlist(sd[c("lower", "median", "upper")])),
    unname(quantile(exp(reference$fits$pooled$draws$log_scale_sd[, 1]),
      c(.025, .5, .975))))
})

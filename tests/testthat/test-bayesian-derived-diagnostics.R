diagnostic_chain_fixture <- function(fit, iterations = 1500L, chains = 4L) {
  set.seed(61591)
  columns <- unlist(lapply(names(fit$parameters), function(nm) {
    n <- length(fit$parameters[[nm]])
    if (!n) character() else if (n == 1L) nm else paste0(nm, "[", seq_len(n), "]")
  }))
  n <- iterations * chains
  values <- matrix(rep(unlist(fit$parameters), each = n), nrow = n)
  values <- values + matrix(rnorm(length(values), sd = .15), nrow = n)
  array(values, c(iterations, chains, ncol(values)), dimnames = list(NULL, NULL, columns))
}

test_that("derived diagnostics require original chains and preserve unavailable evidence", {
  fit <- prepare_bayesian_fixture(bayesian_scale_fixture(n_site = 3, n_year = 2))
  fit$diagnostics <- list(passed = TRUE, marker = "retained")
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) stop("must not extract"),
    .package = "multiScaleR")
  d <- diagnose_bayes_scales(fit)
  expect_false(d$available)
  expect_match(d$reason, "not retained")
  expect_identical(d$raw_diagnostics, fit$diagnostics)
  expect_true(is.na(d$passed))
  expect_true(is.na(d$derived_passed))
  expect_equal(nrow(d$parameters), 0L)
  expect_equal(nrow(d$quantiles), 0L)
  fit <- fake_bayesian_draws(fit)
  expect_false(diagnose_bayes_scales(fit)$available)
  for (nm in c("laplace", "marginal", "importance")) {
    f <- fit; f[[nm]] <- list(passed = TRUE)
    expect_match(diagnose_bayes_scales(f)$reason, "approximate draws")
  }
  expect_error(diagnose_bayes_scales(list()), "multiScaleR_bayes")
  for (p in list(numeric(), 0, 1, NA_real_, Inf, c(.5, .5), ".5", matrix(.5)))
    expect_error(diagnose_bayes_scales(fit, p), "probs")
  e <- readRDS(system.file("extdata", "bayesian-annual-example.rds", package = "multiScaleR"))
  expect_false(diagnose_bayes_scales(e$fits$pooled)$available)
})

test_that("derived reports use unpermuted chains, log map units and quantile MCSE", {
  skip_if_not_installed("posterior")
  fit <- prepare_bayesian_fixture(bayesian_scale_fixture(n_site = 3, n_year = 2))
  a <- diagnostic_chain_fixture(fit)
  fit$sampler <- "mock extraction only"; fit$diagnostics <- list(passed = TRUE)
  # Deliberately unrelated permuted blocks must never affect chain diagnostics.
  fit <- fake_bayesian_draws(fit, n = 3)
  fit$draws$scale_beta[,] <- 100
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  before <- serialize(fit, NULL); rng <- .Random.seed
  d <- diagnose_bayes_scales(fit)
  expect_true(d$available)
  expect_true(d$derived_passed)
  expect_true(d$passed)
  expect_identical(d$chains, 4L)
  expect_identical(d$iterations, 1500L)
  expect_identical(serialize(fit, NULL), before)
  expect_identical(.Random.seed, rng)
  expect_identical(d$parameters$quantity,
    c("log_distance", "log_distance", "scale_group_sd", "response_group_sd"))
  sm <- fit$specification$scales$habitat
  expected <- a[, , "scale_beta"] + exp(a[, , "log_scale_sd"]) * a[, , "scale_z[1]"] +
    log(fit$specification$unit_conv)
  expect_identical(d$parameters$group[1:2], sm$labels)
  expect_identical(d$parameters$row[1:2], sm$representatives)
  expect_equal(d$parameters$rhat[1], posterior::rhat(expected))
  expect_equal(d$parameters$ess_bulk[1], posterior::ess_bulk(expected))
  expect_equal(d$parameters$ess_tail[1], posterior::ess_tail(expected))
  expect_equal(d$quantiles$estimate[1:3], unname(quantile(expected, c(.025, .5, .975))))
  expect_equal(d$quantiles$mcse[1:3], unname(posterior::mcse_quantile(expected, c(.025, .5, .975))))
  expect_identical(d$parameters$units,
    c("log(map units)", "log(map units)", "log distance", "response linear predictor"))
  q <- diagnose_bayes_scales(fit, probs = c(.9, .1))$quantiles
  expect_equal(q$probability[1:2], c(.9, .1))
  expect_equal(q$estimate[1:2], unname(quantile(expected, c(.9, .1))))
  # Scalar blocks can also use Stan's bracketed naming convention.
  dimnames(a)[[3]] <- ifelse(grepl("[", dimnames(a)[[3]], fixed = TRUE),
    dimnames(a)[[3]], paste0(dimnames(a)[[3]], "[1]"))
  expect_equal(diagnose_bayes_scales(fit)$parameters, d$parameters)
})

test_that("stuck derived chains and raw failures remain failures with unchanged gates", {
  skip_if_not_installed("posterior")
  fit <- prepare_bayesian_fixture(bayesian_scale_fixture(n_site = 3, n_year = 2))
  a <- diagnostic_chain_fixture(fit)
  fit$sampler <- "mock"; fit$diagnostics <- list(passed = FALSE, marker = "raw failure")
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  d <- diagnose_bayes_scales(fit)
  expect_true(d$derived_passed)
  expect_false(d$passed)
  expect_identical(d$raw_diagnostics, fit$diagnostics)
  fit$diagnostics$passed <- TRUE
  a[, 4, "scale_beta"] <- a[, 4, "scale_beta"] + 3
  d <- diagnose_bayes_scales(fit)
  expect_false(d$derived_passed)
  expect_false(d$passed)
  expect_gt(d$parameters$rhat[1], 1.01)
  expect_true(fit$diagnostics$passed)
  a[, , "scale_beta"] <- NA_real_
  d <- diagnose_bayes_scales(fit)
  expect_false(d$passed)
  expect_true(is.na(d$parameters$rhat[1]))
  expect_true(all(is.na(d$quantiles$mcse[d$quantiles$quantity == "log_distance"])))
  a <- diagnostic_chain_fixture(fit, chains = 1L)
  expect_false(diagnose_bayes_scales(fit)$derived_passed)
  a <- diagnostic_chain_fixture(fit)
  a[, , "log_scale_sd"] <- log(.3)
  d <- diagnose_bayes_scales(fit)
  expect_false(d$derived_passed)
  expect_true(is.na(d$parameters$rhat[d$parameters$quantity == "scale_group_sd"]))
  a <- a[, , -1, drop = FALSE]
  expect_error(diagnose_bayes_scales(fit), "parameter layout")
  a <- matrix(1, 2, 2)
  expect_error(diagnose_bayes_scales(fit), "dimensions")
})

test_that("derived diagnostics retain multiple predictors and crossed SD mappings", {
  skip_if_not_installed("posterior")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  counts <- f$inputs$binned$csum$habitat
  f$inputs$binned$vsum$forest <- cos(f$inputs$binned$vsum$habitat / pmax(1, counts)) * counts
  f$inputs$binned$csum$forest <- f$inputs$binned$csum$habitat
  fit <- multiScale_bayes(count ~ habitat + forest + (1 | site) + (1 | year),
    f$inputs, f$data, row_id = "key", method = "prepare", scales = list(
      habitat = ~1 + temp + (1 | year) + (1 | site), forest = ~0 + year))
  a <- diagnostic_chain_fixture(fit)
  fit$sampler <- "mock"; fit$diagnostics <- list(passed = TRUE)
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  d <- diagnose_bayes_scales(fit)
  expect_true(d$passed)
  p <- d$parameters
  expect_equal(sum(p$quantity == "log_distance" & p$covariate == "habitat"), 6L)
  expect_equal(sum(p$quantity == "log_distance" & p$covariate == "forest"), 2L)
  expect_identical(p$group[p$quantity == "scale_group_sd"], c("year", "site"))
  expect_true(all(p$covariate[p$quantity == "scale_group_sd"] == "habitat"))
  expect_identical(p$group[p$quantity == "response_group_sd"], c("site", "year"))
  for (group in c("year", "site")) {
    sm <- fit$specification$scales$habitat$random[[group]]
    expected <- exp(a[, , paste0("log_scale_sd[", sm$sd_index, "]")])
    q <- d$quantiles[d$quantiles$quantity == "scale_group_sd" & d$quantiles$group == group, ]
    expect_equal(q$estimate, unname(quantile(expected, c(.025, .5, .975))))
  }
})

test_that("a shared scale without random effects still has one chain target", {
  skip_if_not_installed("posterior")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key", method = "prepare")
  a <- diagnostic_chain_fixture(fit)
  fit$sampler <- "mock"; fit$diagnostics <- list(passed = TRUE)
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  d <- diagnose_bayes_scales(fit)
  expect_true(d$passed)
  expect_equal(nrow(d$parameters), 1L)
  expect_identical(d$parameters$group, "shared")
  expect_equal(nrow(d$quantiles), 3L)
})

probability_chain_fixture <- function(fit, iterations = 600L, chains = 4L) {
  set.seed(82917)
  columns <- unlist(lapply(names(fit$parameters), function(nm) {
    n <- length(fit$parameters[[nm]])
    if (!n) character() else if (n == 1L) nm else paste0(nm, "[", seq_len(n), "]")
  }))
  n <- iterations * chains
  values <- matrix(rep(unlist(fit$parameters), each = n), nrow = n)
  values <- values + matrix(rnorm(length(values), sd = .25), nrow = n)
  array(values, c(iterations, chains, ncol(values)), dimnames = list(NULL, NULL, columns))
}

test_that("probability diagnostics validate requests without extracting unavailable chains", {
  fit <- prepare_bayesian_fixture(bayesian_scale_fixture(n_site = 3, n_year = 2))
  fit$diagnostics <- list(passed = FALSE, marker = "original failure")
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) stop("must not extract"),
    .package = "multiScaleR")
  d <- diagnose_bayes_probabilities(fit, list(habitat = 100))
  expect_false(d$available)
  expect_match(d$reason, "not retained")
  expect_identical(d$raw_diagnostics, fit$diagnostics)
  expect_equal(nrow(d$probabilities), 0L)
  expect_equal(nrow(d$chain_probabilities), 0L)
  expect_identical(d$probabilities$precision_target_met, logical())
  expect_false(d$metadata_establishes_eligibility)
  expect_false("passed" %in% names(d))
  fit <- fake_bayesian_draws(fit)
  expect_false(diagnose_bayes_probabilities(fit, list(habitat = 100))$available)
  for (nm in c("laplace", "marginal", "importance")) {
    f <- fit; f[[nm]] <- list(passed = TRUE); f$sampler <- "mock"
    expect_match(diagnose_bayes_probabilities(f, list(habitat = 100))$reason, "approximate draws")
  }
  expect_error(diagnose_bayes_probabilities(list(), list(habitat = 100)), "multiScaleR_bayes")
  for (x in list(NULL, list(), 100, list(100), list(unknown = 100),
      list(habitat = numeric()), list(habitat = NA_real_), list(habitat = Inf),
      list(habitat = c(100, 100)), list(habitat = "100"), list(habitat = matrix(100)),
      list(habitat = 0), list(habitat = -1), list(habitat = 1 + 1i), setNames(list(100), NA_character_),
      setNames(list(100), ""), setNames(list(100, 200), c("habitat", "habitat"))))
    expect_error(diagnose_bayes_probabilities(fit, x), "thresholds")
  for (target in list(0, 1, NA_real_, Inf, c(.01, .02), ".01", matrix(.01), .01 + .01i))
    expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100), mcse_target = target), "mcse_target")
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100), scale = "radius"), "arg")
  expect_false(diagnose_bayes_probabilities(fit, list(habitat = c(-2, 0)), scale = "log")$available)
  fit$sampler <- "mock"
  local_mocked_bindings(.msr_has_gaussian_blocks = function(spec) TRUE, .package = "multiScaleR")
  expect_match(diagnose_bayes_probabilities(fit, list(habitat = 100))$reason, "dedicated diagnostic adapter")
})

test_that("probabilities use original ordered indicators and preserve fit and RNG", {
  skip_if_not_installed("posterior")
  fit <- prepare_bayesian_fixture(bayesian_scale_fixture(n_site = 3, n_year = 2))
  a <- probability_chain_fixture(fit)
  fit$sampler <- "mock"; fit$diagnostics <- list(passed = FALSE, marker = "retain raw failure")
  fit <- fake_bayesian_draws(fit, n = 3)
  fit$draws$scale_beta[,] <- 1000
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  thresholds <- c(250, 100, 400)
  before <- serialize(fit, NULL); rng <- .Random.seed
  d <- diagnose_bayes_probabilities(fit, list(habitat = thresholds), mcse_target = .03)
  expect_true(d$available)
  expect_identical(d$chains, 4L)
  expect_identical(d$iterations, 600L)
  expect_identical(d$raw_diagnostics, fit$diagnostics)
  expect_identical(serialize(fit, NULL), before)
  expect_identical(.Random.seed, rng)
  expect_false("passed" %in% names(d))
  expect_false(d$metadata_establishes_eligibility)
  expect_false(any(c("draws", "sampler", "parameters") %in% names(d)))
  sm <- fit$specification$scales$habitat
  expect_identical(d$probabilities$group, rep(sm$labels, each = 3L))
  expect_identical(d$probabilities$row, rep(sm$representatives, each = 3L))
  expect_identical(d$probabilities$threshold, rep(thresholds, 2L))
  # Reconstruct each chain independently, without the function's flattening code.
  for (j in seq_along(sm$representatives)) {
    expected <- a[, , "scale_beta"] + exp(a[, , "log_scale_sd"]) *
      a[, , paste0("scale_z[", j, "]")] + log(fit$specification$unit_conv)
    for (k in seq_along(thresholds)) {
      indicator <- 1 * (expected <= log(thresholds[k]))
      row <- (j - 1L) * 3L + k
      expect_equal(d$probabilities$estimate[row], mean(indicator))
      expect_equal(d$probabilities$mcse[row], posterior::mcse_mean(indicator))
      expect_equal(d$probabilities$ess_mean[row], posterior::ess_mean(indicator))
      expect_equal(d$probabilities$rhat[row], posterior::rhat(indicator))
      expect_true(d$probabilities$precision_available[row])
      expect_identical(d$probabilities$precision_target_met[row],
        posterior::mcse_mean(indicator) <= .03)
      rows <- (row - 1L) * 4L + seq_len(4L)
      expect_equal(d$chain_probabilities$estimate[rows], colMeans(indicator))
      expect_identical(d$chain_probabilities$chain[rows], seq_len(4L))
    }
  }
  log_report <- diagnose_bayes_probabilities(fit, list(habitat = log(thresholds)), scale = "log")
  columns <- c("estimate", "mcse", "ess_mean", "rhat", "precision_available")
  expect_equal(log_report$probabilities[columns], d$probabilities[columns])
  expect_true(all(is.na(log_report$probabilities$precision_target_met)))
  dimnames(a)[[3L]] <- ifelse(grepl("[", dimnames(a)[[3L]], fixed = TRUE),
    dimnames(a)[[3L]], paste0(dimnames(a)[[3L]], "[1]"))
  expect_equal(diagnose_bayes_probabilities(fit, list(habitat = thresholds),
    mcse_target = .03), d)
  # Equivalent fitted map units must give the same events.
  fit$specification$unit_conv <- fit$specification$unit_conv * 2
  a[, , "scale_beta[1]"] <- a[, , "scale_beta[1]"] - log(2)
  expect_equal(diagnose_bayes_probabilities(fit, list(habitat = thresholds),
    mcse_target = .03), d)
})

test_that("constant, nonfinite and undersized chains retain explicit precision limits", {
  skip_if_not_installed("posterior")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key", method = "prepare")
  fit$sampler <- "mock"; fit$specification$unit_conv <- 1
  a <- probability_chain_fixture(fit)
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  a[, , "scale_beta"] <- 0
  d <- diagnose_bayes_probabilities(fit, list(habitat = c(-1, 0, 1)), scale = "log", mcse_target = .001)
  expect_identical(d$probabilities$estimate, c(0, 1, 1))
  expect_true(all(d$probabilities$constant))
  expect_true(all(d$chain_probabilities$constant))
  expect_false(any(d$probabilities$precision_available))
  expect_true(all(is.na(d$probabilities$precision_target_met)))
  expect_true(all(is.na(as.matrix(d$probabilities[c("mcse", "ess_mean", "rhat")]))))
  expect_true(all(grepl("Constant indicator", d$probabilities$reason)))
  a[1, 1, "scale_beta"] <- NA_real_
  d <- diagnose_bayes_probabilities(fit, list(habitat = 1))
  expect_true(d$available)
  expect_match(d$probabilities$reason, "no iterations were dropped")
  expect_true(is.na(d$probabilities$estimate))
  expect_true(all(is.na(d$chain_probabilities$estimate)))
  a <- probability_chain_fixture(fit, chains = 1L)
  a[, , "scale_beta"] <- a[, , "scale_beta"] - fit$parameters$scale_beta
  d <- diagnose_bayes_probabilities(fit, list(habitat = 0), scale = "log")
  expect_match(d$probabilities$reason, "two chains")
  expect_true(is.na(d$probabilities$mcse))
  a <- probability_chain_fixture(fit, iterations = 3L)
  a[, , "scale_beta"] <- a[, , "scale_beta"] - fit$parameters$scale_beta
  d <- diagnose_bayes_probabilities(fit, list(habitat = 0), scale = "log")
  expect_match(d$probabilities$reason, "four iterations")
  # Finite log scales beyond exp()'s range remain valid event comparisons.
  a <- probability_chain_fixture(fit)
  a[, , "scale_beta"] <- a[, , "scale_beta"] - fit$parameters$scale_beta + 1000
  d <- diagnose_bayes_probabilities(fit, list(habitat = 1000), scale = "log")
  expect_true(is.finite(d$probabilities$estimate))
  expect_true(d$probabilities$precision_available)
  # A computable MCSE does not turn poor chain mixing into a model pass.
  a[, 4, "scale_beta"] <- a[, 4, "scale_beta"] + 1
  d <- diagnose_bayes_probabilities(fit, list(habitat = 1000), scale = "log", mcse_target = .99)
  expect_gt(d$probabilities$rhat, 1.01)
  expect_false("passed" %in% names(d))
})

test_that("multiple scale designs and crossed group maps preserve event identities", {
  skip_if_not_installed("posterior")
  f <- bayesian_scale_fixture(n_site = 3, n_year = 2)
  counts <- f$inputs$binned$csum$habitat
  f$inputs$binned$vsum$forest <- cos(f$inputs$binned$vsum$habitat / pmax(1, counts)) * counts
  f$inputs$binned$csum$forest <- counts
  fit <- multiScale_bayes(count ~ habitat + forest + (1 | site) + (1 | year),
    f$inputs, f$data, row_id = "key", method = "prepare", scales = list(
      habitat = ~1 + temp + (1 | year) + (1 | site), forest = ~0 + year))
  fit$sampler <- "mock"
  a <- probability_chain_fixture(fit)
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  d <- diagnose_bayes_probabilities(fit, list(forest = 250, habitat = c(100, 250)))
  expect_identical(d$probabilities$covariate, c(rep("forest", 2), rep("habitat", 12)))
  for (i in seq_len(nrow(d$probabilities))) {
    target <- d$probabilities[i, ]; sm <- fit$specification$scales[[target$covariate]]
    expected <- matrix(0, dim(a)[1L], dim(a)[2L])
    for (k in seq_along(sm$beta_index))
      expected <- expected + a[, , paste0("scale_beta[", sm$beta_index[k], "]")] * sm$X[target$row, k]
    for (term in sm$random)
      expected <- expected + exp(a[, , paste0("log_scale_sd[", term$sd_index, "]")]) *
        a[, , paste0("scale_z[", term$z_index[term$index[target$row]], "]")]
    expected <- expected + log(fit$specification$unit_conv)
    expect_equal(target$estimate, mean(expected <= log(target$threshold)))
  }
})

test_that("malformed original parameter layouts fail before event calculation", {
  skip_if_not_installed("posterior")
  fit <- prepare_bayesian_fixture(bayesian_scale_fixture(n_site = 3, n_year = 2))
  fit$sampler <- "mock"
  a <- probability_chain_fixture(fit)
  original <- a
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  a <- matrix(1, 2, 2)
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "dimensions")
  a <- original; dimnames(a) <- NULL
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "dimensions")
  a <- original; dimnames(a)[[3L]][1L] <- NA_character_
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "dimensions")
  a <- original; dimnames(a)[[3L]][1L] <- ""
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "dimensions")
  a <- original; dimnames(a)[[3L]][1L] <- dimnames(a)[[3L]][2L]
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "dimensions")
  a <- original[, , -1L, drop = FALSE]
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "original parameter layout")
  a <- array(c(original, original[, , "scale_beta"]), c(dim(original)[1:2], dim(original)[3] + 1L),
    dimnames = list(NULL, NULL, c(dimnames(original)[[3L]], "scale_beta[1]")))
  expect_error(diagnose_bayes_probabilities(fit, list(habitat = 100)), "ambiguous original parameter layout")
})

test_that("balanced binary events separate mean precision from indicator Rhat", {
  skip_if_not_installed("posterior")
  f <- bayesian_scale_fixture(n_site = 3L, n_year = 2L)
  fit <- multiScale_bayes(count ~ habitat, f$inputs, f$data, row_id = "key", method = "prepare")
  fit$specification$unit_conv <- 1
  fit$sampler <- "deterministic chain fixture"
  fit$diagnostics <- list(passed = FALSE, marker = "raw failure retained")
  a <- probability_chain_fixture(fit, iterations = 400L)
  set.seed(917613)
  a[, , "scale_beta"] <- sample(rep(c(-1, 1), each = 800L))
  indicator <- 1 * (a[, , "scale_beta"] <= 0)
  expected_mcse <- posterior::mcse_mean(indicator)
  expected_ess <- posterior::ess_mean(indicator)
  expected_rhat <- posterior::rhat(indicator)
  local_mocked_bindings(.msr_bayes_sampler_array = function(sampler) a, .package = "multiScaleR")
  before <- serialize(fit, NULL); rng <- .Random.seed
  d <- diagnose_bayes_probabilities(fit, list(habitat = 0), scale = "log", mcse_target = .1)
  p <- d$probabilities
  expect_identical(p$estimate, .5)
  expect_equal(p$mcse, expected_mcse)
  expect_equal(p$ess_mean, expected_ess)
  expect_equal(p$rhat, expected_rhat)
  expect_true(p$precision_available)
  expect_identical(p$precision_target_met, expected_mcse <= .1)
  expect_identical(p$rhat_available, is.finite(expected_rhat))
  if (is.finite(expected_rhat)) expect_identical(p$rhat_reason, "") else
    expect_match(p$rhat_reason, "undefined")
  expect_identical(p$reason, "")
  expect_identical(d$raw_diagnostics, fit$diagnostics)
  expect_false(d$metadata_establishes_eligibility)
  expect_false("passed" %in% names(d))
  expect_identical(serialize(fit, NULL), before)
  expect_identical(.Random.seed, rng)
  unrequested <- diagnose_bayes_probabilities(fit, list(habitat = 0), scale = "log")
  expect_true(unrequested$probabilities$precision_available)
  expect_true(is.na(unrequested$probabilities$precision_target_met))
  constants <- diagnose_bayes_probabilities(fit, list(habitat = c(-2, 2)), scale = "log", mcse_target = .1)
  expect_true(all(constants$probabilities$constant))
  expect_false(any(constants$probabilities$precision_available))
  expect_false(any(constants$probabilities$rhat_available))
  expect_true(all(is.na(constants$probabilities$mcse)))
  expect_true(all(is.na(constants$probabilities$rhat)))
  expect_true(all(is.na(constants$probabilities$precision_target_met)))
  # Preserve a nonfinite auxiliary Rhat as reported; never manufacture a finite
  # replacement or erase the independently computed mean precision.
  local_mocked_bindings(rhat = function(x, ...) Inf, .package = "posterior")
  nonfinite <- diagnose_bayes_probabilities(fit, list(habitat = 0), scale = "log", mcse_target = .1)
  expect_identical(nonfinite$probabilities$rhat, Inf)
  expect_false(nonfinite$probabilities$rhat_available)
  expect_match(nonfinite$probabilities$rhat_reason, "nonfinite")
  expect_true(nonfinite$probabilities$precision_available)
  expect_equal(nonfinite$probabilities$mcse, expected_mcse)
  expect_true(nonfinite$probabilities$precision_target_met)
})

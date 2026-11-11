test_that("Bayesian setup reports the optional parallel dependency", {
  expect_false("processx" %in% check_bayes_setup("slice")$packages$package)
  expect_identical(check_bayes_setup("slice", cores = 2)$packages$package, "processx")
  expect_true(check_bayes_setup("prepare", cores = 2)$ready)
  expect_error(check_bayes_setup("slice", cores = 0), "cores")
  f <- bayesian_scale_fixture(n_site = 4L, n_year = 1L)
  expect_error(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", method = "prepare", cores = 1.5), "cores")
  single <- suppressWarnings(multiScale_bayes(count ~ habitat, f$inputs, f$data,
    row_id = "key", family = "gaussian", method = "slice", chains = 1L, iter = 10L,
    cores = 2L, seed = 129L))
  expect_identical(single$sampling$cores, 1L)
})

test_that("parallel slice sampling keeps the same ordered draws and caller RNG", {
  skip_on_cran()
  skip_if_not_installed("processx")
  f <- bayesian_scale_fixture(n_site = 5L, n_year = 2L)
  run <- function(cores) suppressWarnings(multiScale_bayes(count ~ habitat,
    f$inputs, f$data, row_id = "key", family = "gaussian", method = "slice",
    cores = cores, chains = 2L, iter = 24L, seed = 312L))
  before <- .Random.seed
  serial <- run(1L)
  parallel <- run(2L)
  expect_identical(parallel$sampling$cores, 2L)
  expect_identical(parallel$chains, serial$chains)
  expect_identical(parallel$draws, serial$draws)
  expect_false(identical(parallel$chains[, 1L, ], parallel$chains[, 2L, ]))
  expect_identical(.Random.seed, before)
})

test_that("Bayesian workers preserve task order and can be stopped after an error", {
  skip_on_cran()
  skip_if_not_installed("processx")
  pool <- .msr_bayes_pool(2L, identity, 10L)
  on.exit(.msr_bayes_pool_stop(pool))
  task <- function(job, context) job + context
  environment(task) <- baseenv()
  expect_identical(.msr_bayes_pool_map(pool, as.list(1:5), task), as.list(11:15))
  failure <- function(job, context) stop("intentional worker error")
  environment(failure) <- baseenv()
  expect_error(.msr_bayes_pool_map(pool, list(1L, 2L), failure), "intentional worker error")
  .msr_bayes_pool_stop(pool)
  expect_false(dir.exists(pool$path))
  expect_false(any(vapply(pool$workers, function(w) w$process$is_alive(), logical(1))))
})

test_that("parallel NUTS preserves the ordered Stan chains", {
  skip_on_cran()
  skip_if_not_installed("processx")
  skip_if_not_installed("RTMB", "2.0")
  skip_if_not_installed("tmbstan", "1.2.1")
  skip_if_not_installed("posterior", "1.6.0")
  f <- bayesian_scale_fixture(n_site = 5L, n_year = 2L)
  run <- function(cores) suppressWarnings(multiScale_bayes(count ~ habitat,
    f$inputs, f$data, row_id = "key", family = "gaussian", method = "sample",
    cores = cores, chains = 2L, iter = 40L, seed = 182L))
  serial <- run(1L)
  fit <- run(2L)
  expect_identical(fit$sampling$cores, 2L)
  ordered <- rstan::extract(fit$sampler, permuted = FALSE)
  expect_equal(ordered, rstan::extract(serial$sampler, permuted = FALSE))
  expect_equal(length(rstan::get_sampler_params(fit$sampler)), 2L)
  expect_true(diagnose_bayes_scales(fit)$available)
})

test_that("parallel Laplace starts retain the fitted mode and marginal weights", {
  skip_on_cran()
  skip_if_not_installed("processx")
  skip_if_not_installed("RTMB", "2.0")
  f <- bayesian_scale_fixture(n_site = 6L, n_year = 3L)
  for (method in c("laplace", "marginal")) {
    control <- if (method == "laplace") list(starts = 2L) else
      list(starts = 2L, max_grid = 33L)
    run <- function(cores) suppressWarnings(multiScale_bayes(count ~ habitat,
      f$inputs, f$data, row_id = "key", family = "gaussian", method = method,
      scales = list(habitat = ~1 + (1 | year)), control = control,
      laplace_draws = 10L, cores = cores, seed = 91L))
    serial <- run(1L)
    fit <- run(2L)
    expect_identical(fit[[method]]$cores, 2L)
    if (method == "laplace") expect_equal(fit$laplace$mode, serial$laplace$mode, tolerance = 1e-6)
    else expect_equal(fit$marginal$weights, serial$marginal$weights, tolerance = 1e-6)
    expect_equal(fit$draws, serial$draws, tolerance = 1e-5)
  }
})

test_that("worker pool controls are opt-in and restricted to joint NB2", {
  expect_null(.msr_profile_control(list(), TRUE)$worker_pool)
  expect_error(.msr_profile_control(list(worker_pool = 2), TRUE), "PSOCK")
  expect_error(.msr_profile_control(list(worker_pool = list()), FALSE), "joint_nb2")
  expect_error(testthat::with_mocked_bindings(.msr_joint_pool_ram(list(x = 1), 2),
    .msr_total_ram_bytes = function() 1024^3), "RAM budget")
})

test_that("reusable workers preserve all starts, predictions and uncertainty", {
  skip_on_cran()
  joint_nb2_fit <- function(f, ...) suppressWarnings(multiScale_optim(f$model, f$inputs, backend = "joint_nb2", n_cores = 1, verbose = FALSE, ...))
  f <- profile_glm_fixture(n = 100)
  pool <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(pool), add = TRUE)
  # Test installed code, as workers must use the same build as the controller.
  bootstrap <- function(paths) {
    .libPaths(paths)
    f <- get0(".msr_joint_pool_signature", asNamespace("multiScaleR"), inherits = FALSE)
    if (is.function(f)) f() else NULL
  }
  environment(bootstrap) <- baseenv()
  compatible <- parallel::clusterCall(pool, bootstrap, .libPaths())
  expect_identical(compatible[[1]], .msr_joint_pool_signature())
  ctl <- list(n_starts = 3, coefficient_start = "warm")
  serial <- joint_nb2_fit(f, profile_control = ctl)
  ctl$worker_pool <- pool
  for (i in 1:2) {
    fit <- joint_nb2_fit(f, profile_control = ctl)
    a <- serial$profile_fit$start_results; b <- fit$profile_fit$start_results
    a$seconds <- b$seconds <- NULL
    rownames(a) <- rownames(b) <- NULL
    expect_equal(b, a, tolerance = 1e-10)
    expect_identical(fit$profile_fit$starts, serial$profile_fit$starts)
    expect_identical(fit$profile_fit$winner, serial$profile_fit$winner)
    expect_equal(predict(fit, se.fit = TRUE), predict(serial, se.fit = TRUE), tolerance = 1e-10)
    expect_equal(fit$uncertainty, serial$uncertainty, tolerance = 1e-10)
    expect_null(fit$profile_fit$control$worker_pool)
    expect_length(fit$profile_fit$worker_coefficient_audits, 2)
    saved <- unserialize(serialize(fit, NULL))
    expect_equal(predict(saved), predict(serial), tolerance = 1e-10)
  }
  expect_length(parallel::clusterCall(pool, Sys.getpid), 2)
  expect_error(joint_nb2_fit(f, profile_control = list(worker_pool = pool,
    starts = matrix(100, 1, 1))), "bounds")
  expect_length(parallel::clusterCall(pool, Sys.getpid), 2)
  # A stale worker must fail explicitly, without silently fitting fewer starts.
  parallel::clusterCall(pool[1], function() {
    ns <- asNamespace("multiScaleR")
    unlockBinding(".msr_joint_pool_signature", ns)
    assign(".msr_joint_pool_signature", function() "stale", ns)
    lockBinding(".msr_joint_pool_signature", ns)
    NULL
  })
  expect_error(joint_nb2_fit(f, profile_control = ctl), "differs")
  expect_length(parallel::clusterCall(pool, Sys.getpid), 2)
})


test_that("estimate_multiscale_ram returns a structured recommendation", {
  fix <- make_core_fixture()

  ram <- with_mocked_bindings(
    estimate_multiscale_ram(
      kernel_inputs = fix$kernel_inputs,
      fitted_mod = fix$fitted_mod,
      n_cores = 3,
      PSOCK = TRUE,
      safety_factor = 1.25,
      ram_fraction = 0.8,
      leave_free = 2
    ),
    .msr_detect_cores = function(logical = FALSE) {
      if (isTRUE(logical)) 16L else 8L
    },
    .msr_total_ram_bytes = function() 8 * 1024^3,
    .package = "multiScaleR"
  )

  expect_s3_class(ram, "multiScaleR_ram")
  expect_equal(ram$requested_n_cores, 3L)
  expect_equal(ram$max_cores_by_cpu, 6L)
  expect_true(is.numeric(ram$max_cores_by_ram))
  expect_true(ram$recommended_n_cores <= ram$max_cores_by_cpu)
  expect_equal(ram$system$physical_cores, 8L)
  expect_equal(ram$system$logical_cores, 16L)
  expect_equal(ram$backend, "PSOCK")
  expect_true(all(
    c("kernel_inputs_total", "worker_bundle", "peak_parallel_estimate") %in%
      ram$component_bytes$component
  ))
  expect_equal(nrow(ram$point_summary), 1)
  expect_true(ram$point_summary$n_points > 0)
  expect_false(any(grepl("No fitted model was supplied", ram$notes, fixed = TRUE)))
})

test_that("estimate_multiscale_ram falls back cleanly when RAM is unavailable", {
  fix <- make_core_fixture()

  ram <- with_mocked_bindings(
    estimate_multiscale_ram(
      kernel_inputs = fix$kernel_inputs,
      fitted_mod = NULL,
      n_cores = NULL,
      PSOCK = TRUE
    ),
    .msr_detect_cores = function(logical = FALSE) {
      if (isTRUE(logical)) 8L else 4L
    },
    .msr_total_ram_bytes = function() NA_real_,
    .package = "multiScaleR"
  )

  expect_equal(ram$max_cores_by_cpu, 2L)
  expect_true(is.na(ram$max_cores_by_ram))
  expect_equal(ram$recommended_n_cores, 2L)
  expect_equal(ram$requested_n_cores, 2L)
  expect_true(any(grepl("No fitted model was supplied", ram$notes, fixed = TRUE)))
  expect_true(any(grepl("Total system RAM could not be detected", ram$notes, fixed = TRUE)))
})

test_that("RAM helper utilities compute conservative bounds", {
  expect_equal(multiScaleR:::.msr_max_cores_by_cpu(physical_cores = 8, leave_free = 2), 6L)
  expect_equal(multiScaleR:::.msr_max_cores_by_cpu(physical_cores = NA, leave_free = 2), 1L)

  expect_equal(
    multiScaleR:::.msr_max_cores_by_ram(
      main_bundle_bytes = 100,
      worker_bundle_bytes = 50,
      usable_ram_bytes = 1000,
      safety_factor = 2
    ),
    9L
  )

  expect_equal(
    multiScaleR:::.msr_peak_parallel_bytes(
      n_cores = 3,
      main_bundle_bytes = 100,
      worker_bundle_bytes = 50,
      safety_factor = 2,
      PSOCK = TRUE
    ),
    400
  )

  expect_equal(
    multiScaleR:::.msr_recommended_cores(
      max_cores_by_cpu = 6L,
      max_cores_by_ram = 4L
    ),
    4L
  )
})

test_that("print.multiScaleR_ram reports headline recommendations", {
  fix <- make_core_fixture()

  ram <- with_mocked_bindings(
    estimate_multiscale_ram(
      kernel_inputs = fix$kernel_inputs,
      fitted_mod = NULL,
      n_cores = 2
    ),
    .msr_detect_cores = function(logical = FALSE) {
      if (isTRUE(logical)) 8L else 4L
    },
    .msr_total_ram_bytes = function() 4 * 1024^3,
    .package = "multiScaleR"
  )

  expect_output(print(ram), "Recommended max workers")
  expect_output(print(ram), "Peak estimate at 2 worker")
})

test_that("parallel RAM guard reports the budget and actionable alternatives", {
  inputs <- structure(list(), class = "multiScaleR_data")
  for (backend in c("PSOCK", "fork")) {
    err <- tryCatch(with_mocked_bindings(
      .msr_stop_unsafe_parallel_ram(inputs, fitted_mod = NULL, join_by = NULL,
                                    refit_fn = NULL, n_cores = 3,
                                    PSOCK = backend == "PSOCK"),
      estimate_multiscale_ram = function(..., n_cores, PSOCK) {
        expect_identical(n_cores, 3)
        expect_identical(PSOCK, backend == "PSOCK")
        list(max_cores_by_ram = 1L, backend = backend,
             component_bytes = data.frame(component = "peak_parallel_estimate",
                                           pretty = "233.9 Gb"))
      }, .package = "multiScaleR"), error = identity)
    expect_s3_class(err, "error")
    msg <- conditionMessage(err)
    for (phrase in c("Parallel optimization was not started", "n_cores = 3",
                     "233.9 Gb", backend, "worker limit: 1", "n_cores = 1",
                     "bin = TRUE, store_cell_data = FALSE", "landscape or surface metrics",
                     "estimate_multiscale_ram")) {
      expect_match(msg, phrase, fixed = TRUE)
    }
    expect_null(conditionCall(err))
  }
})

test_that("serial requests bypass the RAM estimator", {
  inputs <- structure(list(), class = "multiScaleR_data")
  for (n_cores in list(NULL, 1L)) {
    result <- with_mocked_bindings(
      .msr_stop_unsafe_parallel_ram(inputs, fitted_mod = NULL, join_by = NULL,
                                    refit_fn = NULL, n_cores = n_cores, PSOCK = TRUE),
      estimate_multiscale_ram = function(...) stop("Estimator must not be called"),
      .package = "multiScaleR")
    expect_null(result)
  }
})

test_that("RAM guard allows requests within budget and ignores unavailable estimates", {
  inputs <- structure(list(), class = "multiScaleR_data")
  for (limit in c(3L, 4L, NA_integer_)) {
    result <- with_mocked_bindings(
      .msr_stop_unsafe_parallel_ram(inputs, fitted_mod = NULL, join_by = NULL,
                                    refit_fn = NULL, n_cores = 3L, PSOCK = TRUE),
      estimate_multiscale_ram = function(...) {
        list(max_cores_by_ram = limit, recommended_n_cores = 1L)
      }, .package = "multiScaleR")
    expect_null(result)
  }
  result <- with_mocked_bindings(
    .msr_stop_unsafe_parallel_ram(inputs, fitted_mod = NULL, join_by = NULL,
                                  refit_fn = NULL, n_cores = 3L, PSOCK = TRUE),
    estimate_multiscale_ram = function(...) stop("RAM detection unavailable"),
    .package = "multiScaleR")
  expect_null(result)
})

test_that("over-budget requests stop before main setup, screening, or any worker", {
  fix <- make_core_fixture()
  for (strategy in c("single", "screen")) {
    started <- character()
    expect_error(with_mocked_bindings(
      multiScale_optim(fix$fitted_mod, fix$kernel_inputs, n_cores = 3,
                       start_strategy = strategy, verbose = FALSE),
      estimate_multiscale_ram = function(...) {
        list(max_cores_by_ram = 1L, backend = "PSOCK",
             component_bytes = data.frame(component = "peak_parallel_estimate",
                                           pretty = "233.9 Gb"))
      },
      build_opt_context = function(...) { started <<- c(started, "context") },
      .msr_prescreen_start = function(...) { started <<- c(started, "screen") },
      makeCluster = function(...) { started <<- c(started, "PSOCK") },
      makeForkCluster = function(...) { started <<- c(started, "fork") },
      cluster_prep = function(...) { started <<- c(started, "export") },
      optimParallel = function(...) { started <<- c(started, "parallel") },
      optim = function(...) { started <<- c(started, "serial") },
      .package = "multiScaleR"), "Parallel optimization was not started")
    expect_length(started, 0)
  }
})

test_that("RAM estimates include the same binned context used by optimization", {
  fix <- make_core_fixture()
  for (lean in c(FALSE, TRUE)) {
    inputs <- fix$kernel_inputs
    if (lean) {
      inputs$d_list <- inputs$raw_cov <- NULL
      inputs$cell_data_stored <- FALSE
    }
    context <- build_opt_context(fix$fitted_mod, cov_df = inputs$raw_cov,
      scale_vars = inputs$scale_vars, unit_conv = inputs$unit_conv,
      resolution = inputs$resolution, n_cols = inputs$n_cols, binned = inputs$binned)
    expect_false(is.null(context$binned))
    ram <- with_mocked_bindings(
      estimate_multiscale_ram(inputs, fitted_mod = fix$fitted_mod, n_cores = 2),
      .msr_detect_cores = function(logical = FALSE) 4L,
      .msr_total_ram_bytes = function() 8 * 1024^3,
      .package = "multiScaleR")
    expect_equal(ram$component_bytes$bytes[ram$component_bytes$component == "opt_context"],
                 as.numeric(utils::object.size(context)))
    expect_false(any(grepl("context could not be built", ram$notes, fixed = TRUE)))
  }
})

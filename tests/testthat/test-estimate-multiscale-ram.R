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

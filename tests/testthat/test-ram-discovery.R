test_that("Windows RAM discovery uses its direct physical-memory query", {
  calls <- list()
  local_mocked_bindings(.msr_system_stdout = function(command, args) {
    calls[[length(calls) + 1L]] <<- list(command = command, args = args)
    "8589934592"
  }, .package = "multiScaleR")
  expect_equal(.msr_windows_total_ram_bytes(), 8 * 1024^3)
  expect_length(calls, 1L)
  expect_identical(calls[[1L]]$command, "powershell")
  expect_match(paste(calls[[1L]]$args, collapse = " "),
    "Microsoft.VisualBasic.Devices.ComputerInfo", fixed = TRUE)
})

test_that("invalid direct RAM output falls back without inventing a budget", {
  for (invalid in list(character(), "unavailable", "0", "-1", "Inf",
      "diagnostic code 123", c("1", "2"))) {
    calls <- list()
    local_mocked_bindings(.msr_system_stdout = function(command, args) {
      calls[[length(calls) + 1L]] <<- args
      if (length(calls) == 1L) invalid else "8589934592"
    }, .package = "multiScaleR")
    expect_equal(.msr_windows_total_ram_bytes(), 8 * 1024^3)
    expect_length(calls, 2L)
    expect_match(paste(calls[[2L]], collapse = " "),
      "Win32_ComputerSystem", fixed = TRUE)
  }
})

test_that("unavailable RAM stays unknown when both Windows queries fail", {
  calls <- 0L
  local_mocked_bindings(.msr_system_stdout = function(...) {
    calls <<- calls + 1L
    character()
  }, .package = "multiScaleR")
  expect_true(is.na(.msr_windows_total_ram_bytes()))
  expect_equal(calls, 2L)
})

test_that("failed RAM discovery is retried and a successful value is cached", {
  cache <- new.env(parent = emptyenv())
  calls <- 0L
  local_mocked_bindings(.msr_hw_cache = cache,
    .msr_total_ram_bytes_uncached = function() {
      calls <<- calls + 1L
      if (calls == 1L) NA_real_ else 8 * 1024^3
    }, .package = "multiScaleR")
  expect_true(is.na(.msr_total_ram_bytes()))
  expect_false(exists("total_ram_bytes", cache, inherits = FALSE))
  expect_equal(.msr_total_ram_bytes(), 8 * 1024^3)
  expect_equal(.msr_total_ram_bytes(), 8 * 1024^3)
  expect_equal(calls, 2L)
})

test_that("old invalid RAM cache entries cannot prevent a later valid reading", {
  for (invalid in list(NA_real_, 0, -1, Inf, "invalid")) {
    cache <- new.env(parent = emptyenv())
    cache$total_ram_bytes <- invalid
    local_mocked_bindings(.msr_hw_cache = cache,
      .msr_total_ram_bytes_uncached = function() 8 * 1024^3,
      .package = "multiScaleR")
    expect_equal(.msr_total_ram_bytes(), 8 * 1024^3)
    expect_equal(cache$total_ram_bytes, 8 * 1024^3)
  }
})

test_that("unknown and insufficient memory still stop differentiation admission", {
  spec <- list(response_design = NULL,
    surface = list(buffers = list(habitat = list(list(v = 1, d = 1)))))
  cache <- new.env(parent = emptyenv())
  local_mocked_bindings(.msr_hw_cache = cache,
    .msr_total_ram_bytes_uncached = function() NA_real_,
    .package = "multiScaleR")
  expect_error(.msr_scale_surface_ram(spec), "Cannot detect RAM.*no tape was launched")
  local_mocked_bindings(.msr_total_ram_bytes_uncached = function() 1024^3,
    .package = "multiScaleR")
  expect_error(.msr_scale_surface_ram(spec), "exceed.*conservative differentiation RAM budget")
})

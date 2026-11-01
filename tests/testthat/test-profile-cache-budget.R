cache_budget_fixture <- function() structure(list(
  scale_est = data.frame(Mean = 2, SE = .2, row.names = "hab"),
  min_D = .1, max_D = 10, kernel_inputs = list(kernel = "gaussian", unit_conv = 1),
  optim_results = list(par = 2, value = 4),
  fitted_mod_original = list(y = 1:5)), class = "multiScaleR")

cache_budget_profile <- function(...) list(
  intervals = matrix(c(1, 3), 1, dimnames = list("hab", c("2.5%", "97.5%"))),
  method = c(hab = "profile"))

test_that("budgets preserve intervals, evict entries, and support disabling", {
  old <- options(multiScaleR.profile_cache_max_bytes = 1e6)
  on.exit(options(old), add = TRUE)
  on.exit(.msr_clear_profile_cache(), add = TRUE)
  .msr_clear_profile_cache()
  object <- cache_budget_fixture()
  count <- 0
  key <- function(object, ...) structure(paste0("key", object$max_D), payload = as.raw(1:100))
  profile <- function(...) { count <<- count + 1; cache_budget_profile() }
  with_mocked_bindings({
    first <- scale_ci_table(object, df = 20, profile = TRUE)
    expect_equal(scale_ci_table(object, df = 20, profile = TRUE), first)
    expect_equal(count, 1)
    bytes <- .msr_profile_cache_bytes()
    options(multiScaleR.profile_cache_max_bytes = bytes + 10)
    object$max_D <- 11
    expect_equal(scale_ci_table(object, df = 20, profile = TRUE), first)
    expect_equal(length(ls(.profile_scale_cache)), 1L)
    expect_false(exists("key10", .profile_scale_cache, inherits = FALSE))
    expect_lte(.msr_profile_cache_bytes(), bytes + 10)
    options(multiScaleR.profile_cache_max_bytes = 0)
    out <- scale_ci_table(object, df = 20, profile = TRUE)
    expect_equal(out[["2.5%"]], 1)
    expect_identical(attr(out, "profile_cache")$reason, "zero_budget")
    expect_equal(.msr_profile_cache_bytes(), 0)
    options(multiScaleR.profile_cache_max_bytes = 1)
    out <- scale_ci_table(object, df = 20, profile = TRUE)
    expect_identical(attr(out, "profile_cache")$reason, "payload_exceeds_budget")
    expect_equal(.msr_profile_cache_bytes(), 0)
    expect_equal(out[["97.5%"]], 3)
  }, profile_scale_cache_key = key, profile_sigma_intervals = profile)
})

test_that("key-storage failures bypass caching while numerical failures fail", {
  old <- options(multiScaleR.profile_cache_max_bytes = 1e6)
  on.exit(options(old), add = TRUE)
  object <- cache_budget_fixture()
  with_mocked_bindings({
    expect_warning(out <- scale_ci_table(object, 20, profile = TRUE),
      "computing intervals without caching: key storage denied")
    expect_identical(attr(out, "profile_cache")$status, "unavailable")
    expect_equal(out[["2.5%"]], 1)
  }, profile_scale_cache_key = function(...) stop("key storage denied"),
    profile_sigma_intervals = cache_budget_profile)
  object$opt_context <- list(refit_fn = function(...) NULL)
  with_mocked_bindings({
    out <- scale_ci_table(object, 20, profile = TRUE)
    expect_identical(attr(out, "profile_cache")$reason, "custom_refit")
    expect_equal(out[["97.5%"]], 3)
  }, profile_scale_cache_key = function(...) stop("Keys must not be built"),
    profile_sigma_intervals = cache_budget_profile)
  with_mocked_bindings(
    expect_error(scale_ci_table(object, 20, profile = TRUE), "numerical failure"),
    profile_sigma_intervals = function(...) stop("numerical failure"))
  for (invalid in list(-1, NA_real_, Inf, "64 MiB", c(1, 2))) {
    options(multiScaleR.profile_cache_max_bytes = invalid)
    expect_error(scale_ci_table(object, 20, profile = TRUE), "multiScaleR.profile_cache_max_bytes")
  }
})

test_that("digest collisions and corrupt entries never supply different intervals", {
  old <- options(multiScaleR.profile_cache_max_bytes = 1e6)
  on.exit(options(old), add = TRUE)
  on.exit(.msr_clear_profile_cache(), add = TRUE)
  .msr_clear_profile_cache()
  object <- cache_budget_fixture()
  count <- 0
  with_mocked_bindings({
    scale_ci_table(object, 20, profile = TRUE)
    object$max_D <- 11
    scale_ci_table(object, 20, profile = TRUE)
    expect_equal(count, 2)
    assign("collision", list(payload = as.raw(11), result = NULL), envir = .profile_scale_cache)
    out <- scale_ci_table(object, 20, profile = TRUE)
    expect_equal(count, 3)
    assign("collision", 42, envir = .profile_scale_cache)
    out <- scale_ci_table(object, 20, profile = TRUE)
    expect_equal(count, 4)
    expect_equal(out[["2.5%"]], 1)
  }, profile_scale_cache_key = function(object, ...)
    structure("collision", payload = as.raw(object$max_D)),
    profile_sigma_intervals = function(...) { count <<- count + 1; cache_budget_profile() })
})


test_that("real refit profile limits agree with caching disabled", {
  skip_if_not_installed("terra")
  old <- options(multiScaleR.profile_cache_max_bytes = 64 * 1024^2)
  on.exit(options(old), add = TRUE)
  on.exit(.msr_clear_profile_cache(), add = TRUE)
  .msr_clear_profile_cache()
  fit <- make_core_fixture()$opt
  cached <- scale_ci_table(fit, df = 100, profile = TRUE)
  expect_true(length(ls(.profile_scale_cache)) > 0)
  options(multiScaleR.profile_cache_max_bytes = 0)
  direct <- scale_ci_table(fit, df = 100, profile = TRUE)
  expect_equal(direct[c("2.5%", "97.5%")], cached[c("2.5%", "97.5%")])
  expect_identical(attr(direct, "interval_method"), attr(cached, "interval_method"))
  expect_equal(.msr_profile_cache_bytes(), 0)
})


test_that("oversized profile results are returned without retention", {
  old <- options(multiScaleR.profile_cache_max_bytes = 10000)
  on.exit(options(old), add = TRUE)
  on.exit(.msr_clear_profile_cache(), add = TRUE)
  .msr_clear_profile_cache()
  with_mocked_bindings({
    out <- scale_ci_table(cache_budget_fixture(), 20, profile = TRUE)
    expect_equal(out[["2.5%"]], 1)
    expect_identical(attr(out, "profile_cache")$reason, "result_exceeds_budget")
    expect_equal(.msr_profile_cache_bytes(), 0)
  }, profile_scale_cache_key = function(...) structure("key", payload = raw(1)),
    profile_sigma_intervals = function(...) {
      out <- cache_budget_profile()
      out$padding <- raw(20000)
      out
    })
})

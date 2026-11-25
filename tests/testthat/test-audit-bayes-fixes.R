# Regression tests for Bayesian findings of the 2026-11-24 integration audit
# (B-1, B-3, B-6). Where a model is needed it is built from sim_rast(),
# sim_dat(), and kernel_prep() inputs.

audit_bayes_inputs <- function(store_cell_data = TRUE) {
  rs <- sim_rast(dim = 60, resolution = 10, user_seed = 11)
  sim <- sim_dat(alpha = 0.5, beta = 0.8, kernel = "gaussian", sigma = 60,
    type = "count", n_points = 40, raster_stack = rs[["bin1"]], max_D = 250,
    user_seed = 11)
  inputs <- kernel_prep(pts = sim$pts, raster_stack = rs[["bin1"]], max_D = 250,
    store_cell_data = store_cell_data, verbose = FALSE)
  data <- data.frame(y = sim$obs, inputs$kernel_dat)
  folds <- grouped_cv_folds(rep(1:2, length.out = nrow(data)), k = 2, ids = rownames(data))
  list(inputs = inputs, data = data, folds = folds)
}

test_that("worker files that Windows briefly holds open are retried before failing (B-1)", {
  file <- tempfile(fileext = ".rds")
  stuck <- tempfile(fileext = ".rds")
  on.exit(unlink(c(file, paste0(file, ".ready"), stuck)), add = TRUE)
  saveRDS(list(answer = 42), file)
  file.create(paste0(file, ".ready"))
  # The first two removals leave the files in place, as a held handle would.
  calls <- 0L
  local_mocked_bindings(.msr_bayes_pool_unlink = function(files) {
    calls <<- calls + 1L
    if (calls <= 2L) return(1L)
    unlink(files)
  })
  expect_identical(.msr_bayes_pool_take(file), list(answer = 42))
  expect_identical(calls, 3L)
  expect_false(file.exists(file))
  expect_false(file.exists(paste0(file, ".ready")))

  # A file that never goes away still fails, after a bounded number of tries.
  saveRDS(1, stuck)
  calls <- 0L
  local_mocked_bindings(.msr_bayes_pool_unlink = function(files) {
    calls <<- calls + 1L
    1L
  })
  expect_false(.msr_bayes_pool_clear(stuck, attempts = 3L, wait = 0))
  expect_identical(calls, 3L)
  expect_true(file.exists(stuck))
  # A nonzero unlink() status is not fatal when the file is in fact gone.
  local_mocked_bindings(.msr_bayes_pool_unlink = function(files) {
    unlink(files)
    1L
  })
  expect_true(.msr_bayes_pool_clear(stuck, attempts = 1L, wait = 0))
})

test_that("a stuck worker file stops with the original message after about two seconds (B-1)", {
  skip_on_cran()
  stuck <- tempfile(fileext = ".rds")
  on.exit(unlink(stuck), add = TRUE)
  saveRDS(1, stuck)
  local_mocked_bindings(.msr_bayes_pool_unlink = function(files) 1L)
  started <- proc.time()[["elapsed"]]
  expect_error(.msr_bayes_pool_take(stuck), "Cannot clear Bayesian worker data at ", fixed = TRUE)
  expect_lt(proc.time()[["elapsed"]] - started, 10)
})

test_that("diagnostics of a prepared Bayesian object say it has no draws (B-3)", {
  x <- audit_bayes_inputs()
  prepared <- multiScale_bayes(y ~ bin1, x$inputs, x$data, family = "poisson",
    method = "prepare")
  expect_null(prepared$draws)
  reason <- "This object has no posterior draws. Fit with method = 'sample', 'slice', 'laplace', or 'marginal'."
  scales <- diagnose_bayes_scales(prepared)
  expect_false(scales$available)
  expect_identical(scales$reason, reason)
  expect_equal(nrow(scales$parameters), 0L)
  expect_true(is.na(scales$passed))
  probabilities <- diagnose_bayes_probabilities(prepared, list(bin1 = 100))
  expect_false(probabilities$available)
  expect_identical(probabilities$reason, reason)
  expect_equal(nrow(probabilities$probabilities), 0L)
  # The same wording as the other functions that need draws.
  expect_error(posterior_scales(prepared), "no posterior draws", fixed = TRUE)
})

test_that("Bayesian CV names bin_reference = 'fixed' for lean inputs (B-6)", {
  x <- audit_bayes_inputs(store_cell_data = FALSE)
  expect_null(x$inputs$raw_cov)
  expect_error(cross_validate_bayes(y ~ bin1, x$data, x$inputs, x$folds,
    family = "poisson", method = "laplace"),
    "or set bin_reference = 'fixed'", fixed = TRUE)
  skip_on_cran()
  skip_if_not_installed("RTMB")
  cv <- cross_validate_bayes(y ~ bin1, x$data, x$inputs, x$folds, family = "poisson",
    method = "laplace", bin_reference = "fixed", seed = 1)
  expect_identical(cv$settings$bin_reference, "fixed")
  expect_true(all(cv$fold_results$status %in% c("scored", "diagnostic_failed")))
})

test_that("Bayesian CV warns once when it returns unscored folds (B-6)", {
  skip_if_not_installed("RTMB")
  x <- audit_bayes_inputs()
  # With alternating folds, one analysis set holds a single site level (the
  # fit fails) and the other never sees level "b" (prediction fails), so both
  # folds fail with different reasons.
  x$data$site <- factor(ifelse(seq_len(nrow(x$data)) == 1L, "rare",
    rep(c("a", "b"), length.out = nrow(x$data))))
  warnings <- testthat::capture_warnings(cv <- cross_validate_bayes(y ~ bin1 + site,
    x$data, x$inputs, x$folds, family = "poisson", method = "laplace", seed = 1))
  expect_length(warnings, 1L)
  expect_match(warnings, "cross_validate_bayes() returned an incomplete result: 2 of 2 folds failed",
    fixed = TRUE)
  for (fold in names(cv$errors)) expect_match(warnings, cv$errors[[fold]]$message, fixed = TRUE)
  expect_match(warnings, "Reasons:", fixed = TRUE)
  expect_false(cv$complete)
  expect_true(is.na(cv$summary$estimate))
})

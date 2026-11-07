testthat::test_that("every repeated regression and variance coordinate moves with its block position", {
  fitted <- c(rep("beta_lambda", 2), "beta_alpha", rep("beta_phi", 2),
    rep("beta_dist", 2), "beta_scale", rep("beta_rem", 3),
    rep("lsigma_lambda", 2), "lsigma_rem")
  constructor <- c(rep("beta_lambda", 2), rep("beta_phi", 2), rep("beta_dist", 2),
    rep("beta_rem", 3), rep("lsigma_lambda", 2), "lsigma_rem", "beta_alpha", "beta_scale")
  values <- seq_along(fitted) / 10
  actual <- .msr_gdr_exact_start_map(constructor, fitted, values)
  testthat::expect_equal(actual$starts, values[c(1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 3, 8)])
  testthat::expect_identical(actual$map$block, constructor)
  testthat::expect_equal(actual$map$within_block, c(1, 2, 1, 2, 1, 2, 1, 2, 3, 1, 2, 1, 1, 1))
  restored <- numeric(length(values))
  restored[actual$map$fitted_position] <- actual$starts
  testthat::expect_identical(restored, values)
  testthat::expect_match(actual$qualification, "no_model_or_random_effect_admission")
})

testthat::test_that("incomplete, invalid and random-mode maps fail rather than discard parameters", {
  testthat::expect_error(.msr_gdr_exact_start_map(c("beta_lambda", "beta_alpha"), "beta_lambda", 0), "bijection")
  testthat::expect_error(.msr_gdr_exact_start_map("beta_lambda", "beta_lambda", Inf), "finite")
  testthat::expect_error(.msr_gdr_exact_start_map(c("beta_lambda", "lsigma_lambda"),
    c("beta_lambda", "lsigma_rem"), c(0, .1)), "bijection")
  testthat::expect_error(.msr_gdr_exact_start_map(c("beta_lambda", "b_lambda"),
    c("beta_lambda", "b_lambda"), c(0, .1)), "random modes")
  testthat::expect_error(.msr_gdr_exact_start_map(NA_character_, "beta_lambda", 0), "identities")
})

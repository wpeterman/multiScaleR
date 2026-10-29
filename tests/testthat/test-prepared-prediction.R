test_that("prepared GLM predictions freeze design and retain covariance safeguards", {
  f <- profile_glm_fixture(n = 180)
  fit <- profile_glm_fit(f, profile_control = list(n_starts = 3))
  prepared <- prepare_prediction(fit, f$data, f$inputs)
  for (type in c("link", "response")) {
    expect_equal(predict(prepared, type = type),
                 predict(fit, f$data, f$inputs, type = type), tolerance = 1e-12)
    expect_equal(predict(prepared, type = type, se.fit = TRUE, uncertainty = "conditional"),
                 predict(fit, f$data, f$inputs, type = type, se.fit = TRUE,
                         uncertainty = "conditional"), tolerance = 1e-12)
  }
  restored <- unserialize(serialize(prepared, NULL))
  expect_equal(predict(restored), predict(prepared))
  original <- predict(prepared)
  f$inputs$binned$vsum$hab[,] <- 0
  f$data$z <- 100
  fit$uncertainty$par[fit$uncertainty$beta_index[1]] <- 100
  expect_equal(predict(prepared), original)
  expect_error(predict(prepared, kernel_inputs = f$inputs), "fixed inputs")
  expect_error(prepare_prediction(fit, max_bytes = Inf), "max_bytes")
  expect_error(prepare_prediction(fit, max_bytes = 0), "max_bytes")
  weak <- prepared
  weak$uncertainty$available <- FALSE
  weak$uncertainty$reasons <- "weak scale identification"
  expect_error(predict(weak, se.fit = TRUE), "Joint Wald uncertainty unavailable")
  expect_equal(predict(weak, se.fit = TRUE, uncertainty = "conditional"),
               predict(prepared, se.fit = TRUE, uncertainty = "conditional"))
  if (prepared$uncertainty$available)
    expect_equal(predict(prepared, interval = "confidence"),
                 predict(restored, interval = "confidence"))
  expect_error(prepare_prediction(list(backend = "profile_glm")), "requires")
})

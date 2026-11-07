simulation_score_fixture <- function() data.frame(replicate = 1:3, target = "density",
  units = "birds/km2", truth = 0.4, estimate = c(0.2, 0.4, 0.6), eligible = TRUE)

test_that("simulation scores transform estimates and interval endpoints together", {
  x <- simulation_score_fixture(); x$lower <- x$estimate / 2; x$upper <- c(0.45, 0.8, 0.95)
  a <- score_simulation(x, "response")
  expect_equal(a$summary$bias, 0, tolerance = 1e-15)
  expect_equal(a$summary$relative_bias_pct, 0, tolerance = 1e-12)
  expect_equal(a$summary$cv_pct, 50)
  expect_equal(a$summary$rmse, sqrt(0.08 / 3))
  expect_equal(a$summary$coverage, 1)
  y <- x; y[c("estimate", "lower", "upper")] <- lapply(y[c("estimate", "lower", "upper")], log)
  b <- score_simulation(y, "log")
  expect_equal(b$summary, a$summary)
  expect_equal(b$attempts$response_lower, x$lower)
  expect_equal(b$attempts$response_upper, x$upper)
  y <- x; y[c("estimate", "lower", "upper")] <- lapply(y[c("estimate", "lower", "upper")], qlogis)
  expect_equal(score_simulation(y, "logit")$summary, a$summary)
  expect_identical(a$attempts$estimate, x$estimate)
})

test_that("truth near one is nonsingular and constant recovery has zero error", {
  x <- simulation_score_fixture(); x$truth <- x$estimate <- 1
  a <- score_simulation(x, "response")$summary
  expect_equal(a$relative_bias_pct, 0)
  expect_equal(a$rmse, 0)
  expect_equal(a$cv_pct, 0)
  x$truth <- 1 + 1e-10
  expect_true(is.finite(score_simulation(x, "response")$summary$relative_bias_pct))
})

test_that("failed attempts and unusable intervals remain visible", {
  x <- simulation_score_fixture(); x$eligible[1] <- FALSE
  x$reason <- c("native_hessian_failure", "", "")
  x$estimate[3] <- Inf; x$lower <- c(0.1, 0.5, 0.4); x$upper <- c(0.3, 0.3, 0.8)
  a <- score_simulation(x, "response")
  expect_equal(a$summary$n_attempted, 3)
  expect_equal(a$summary$n_eligible, 1)
  expect_equal(a$summary$n_failed, 2)
  expect_equal(a$summary$n_interval, 0)
  expect_equal(a$attempts$score_reason, c("native_hessian_failure", "", "nonfinite_response_estimate"))
  expect_equal(a$attempts$interval_status[2], "reversed_endpoints")
  expect_equal(a$summary$coverage_reason, "no_eligible_intervals")
  expect_equal(a$summary$cv_reason, "fewer_than_two_eligible_attempts")
  x$eligible[] <- FALSE
  expect_equal(score_simulation(x, "response")$summary$metric_reason, "no_eligible_attempts")
})

test_that("zero and unidentified targets retain explicit missing metric reasons", {
  x <- simulation_score_fixture(); x$truth <- 0; x$estimate <- c(-1, 0, 1)
  a <- score_simulation(x, "response")$summary
  expect_true(is.na(a$relative_bias_pct))
  expect_equal(a$relative_bias_reason, "zero_truth")
  expect_true(is.na(a$cv_pct))
  expect_equal(a$cv_reason, "zero_mean")
  x$target <- "ecological_scale_under_null_effect"; x$identified <- FALSE
  a <- score_simulation(x, "response")$summary
  expect_equal(a$n_attempted, 3)
  expect_equal(a$n_eligible, 3)
  expect_equal(a$metric_reason, "unidentified_target")
  expect_true(all(is.na(a[c("bias", "relative_bias_pct", "rmse", "cv_pct", "coverage")])))
})

test_that("simulation scoring groups named targets without coefficient-position assumptions", {
  x <- simulation_score_fixture(); y <- x; y$target <- "removal_probability"; y$units <- "probability"
  a <- score_simulation(rbind(x, y)[c(6, 1, 4, 3, 2, 5), ], "response")$summary
  expect_setequal(a$target, c("density", "removal_probability"))
  expect_identical(a$target, c("removal_probability", "density"))
  expect_equal(a$n_attempted, c(3, 3))
  expect_equal(a$cv_pct, c(50, 50))
  expect_error(score_simulation(rbind(x, x[1, ]), "response"), "only once")
  x$truth[2] <- 1
  expect_error(score_simulation(x, "response"), "one response-scale truth")
})

test_that("scoring requires explicit scale and valid interval and unit contracts", {
  x <- simulation_score_fixture()
  expect_error(score_simulation(x), "input_scale explicitly")
  expect_error(score_simulation(x, "bad"), "arg")
  expect_error(score_simulation(x, "response", level = 1), "level")
  x$lower <- 0
  expect_error(score_simulation(x, "response"), "both lower and upper")
  x$lower <- NULL; x$units[1] <- "birds/ha"
  expect_error(score_simulation(x, "response"), "one response-scale truth")
  x$units[] <- "birds/km2"; x$identified <- c(TRUE, TRUE, FALSE)
  expect_error(score_simulation(x, "response"), "identification status")
  x$identified <- NULL; x$eligible[1] <- NA
  expect_error(score_simulation(x, "response"), "eligible nonmissing")
  x <- simulation_score_fixture(); x$estimate <- matrix(1:6, 3, 2)
  expect_error(score_simulation(x, "response"), "atomic vector")
})

test_that("score printing retains all target and failure records while bounding output", {
  x <- simulation_score_fixture()[rep(1, 12), ]
  x$target <- sprintf("target_%02d", 1:12)
  a <- score_simulation(x, "response")
  expect_s3_class(a, "msr_simulation_score")
  text <- utils::capture.output(returned <- print(a))
  expect_match(paste(text, collapse = "\n"), "first 10 of 12")
  expect_false(any(grepl("target_12", text)))
  expect_identical(returned, a)
  expect_equal(nrow(a$attempts), 12)
})

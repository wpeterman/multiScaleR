# select_cv() must accept candidates whose default metric sets differ by
# family: Poisson records deviance by default, NB2 does not.

mixed_cv_fixture <- function() {
  set.seed(42)
  n <- 60
  dat <- data.frame(x = stats::rnorm(n))
  dat$y <- stats::rnbinom(n, size = 2, mu = exp(0.5 + 0.6 * dat$x))
  folds <- grouped_cv_folds(rep(seq_len(12), each = 5), k = 3)
  list(dat = dat, folds = folds)
}

test_that("select_cv scores mixed Poisson and NB2 GLM candidates on shared metrics", {
  skip_if_not_installed("MASS")
  f <- mixed_cv_fixture()
  pois <- cross_validate_glm(y ~ x, f$dat, f$folds, "poisson", keep_selection = TRUE)
  nb <- cross_validate_glm(y ~ x, f$dat, f$folds, "nbinom2", keep_selection = TRUE)
  expect_true("deviance" %in% pois$settings$metrics)
  expect_false("deviance" %in% nb$settings$metrics)

  shared <- c("rmse", "mae", "neg_log_score")
  for (models in list(list(pois = pois, nb = nb), list(nb = nb, pois = pois))) {
    selected <- select_cv(models, criterion = "AICc")
    expect_true(selected$complete)
    expect_identical(selected$selection$metrics, shared)
    expect_identical(selected$selection$metrics_dropped, "deviance")
    expect_identical(selected$settings$metrics, shared)
    expect_identical(selected$summary$metric, shared)
    expect_false("deviance" %in% names(selected$predictions))
    expect_true(all(is.finite(selected$summary$estimate)))
  }

  # Scores equal those from a Poisson run restricted to the shared metrics.
  pois_shared <- cross_validate_glm(y ~ x, f$dat, f$folds, "poisson",
    metrics = shared, keep_selection = TRUE)
  mixed <- select_cv(list(pois = pois, nb = nb), criterion = "AICc")
  matched <- select_cv(list(pois = pois_shared, nb = nb), criterion = "AICc")
  expect_identical(mixed$selection$chosen, matched$selection$chosen)
  expect_equal(mixed$summary, matched$summary)
  expect_equal(mixed$predictions, matched$predictions)
  expect_identical(matched$selection$metrics_dropped, character())
})

test_that("select_cv leaves unmixed candidate lists unchanged and rejects disjoint metrics", {
  f <- mixed_cv_fixture()
  null <- cross_validate_glm(y ~ 1, f$dat, f$folds, "poisson", keep_selection = TRUE)
  trend <- cross_validate_glm(y ~ x, f$dat, f$folds, "poisson", keep_selection = TRUE)
  selected <- select_cv(list(null = null, trend = trend), "AIC", ties = "first")
  expect_identical(selected$settings$metrics, null$settings$metrics)
  expect_identical(selected$selection$metrics_dropped, character())
  # Each fold's scores are the winning candidate's recorded scores.
  for (fold in names(f$folds)) {
    winner <- list(null = null, trend = trend)[[selected$selection$chosen[[fold]]]]
    expect_equal(selected$metrics[selected$metrics$fold == fold, ],
      winner$metrics[winner$metrics$fold == fold, ], ignore_attr = TRUE)
  }

  rmse_only <- cross_validate_glm(y ~ 1, f$dat, f$folds, "poisson",
    metrics = "rmse", keep_selection = TRUE)
  mae_only <- cross_validate_glm(y ~ x, f$dat, f$folds, "poisson",
    metrics = "mae", keep_selection = TRUE)
  expect_error(select_cv(list(a = rmse_only, b = mae_only), "AIC"), "share no scoring metric")
})

test_that("select_cv handles mixed Poisson and NB2 cross_validate() candidates", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  data("landscape_counts", package = "multiScaleR", envir = environment())
  data("surv_pts", package = "multiScaleR", envir = environment())
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  inputs <- kernel_prep(pts = surv_pts, raster_stack = land[["land1"]], max_D = 1700,
    kernel = "gaussian", bin = TRUE, verbose = FALSE)
  dat <- data.frame(landscape_counts, inputs$kernel_dat)
  folds <- grouped_cv_folds(rep(seq_len(10), length.out = nrow(dat)), k = 2)
  pois <- cross_validate(counts ~ land1, dat, inputs, folds, family = "poisson", keep_selection = TRUE)
  nb <- cross_validate(counts ~ land1, dat, inputs, folds, family = "nbinom2", keep_selection = TRUE)
  selected <- select_cv(list(pois = pois, nb = nb), criterion = "AICc")
  expect_true(selected$complete)
  expect_identical(selected$selection$metrics, c("rmse", "mae", "neg_log_score"))
  expect_identical(selected$selection$metrics_dropped, "deviance")
})

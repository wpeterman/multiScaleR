# Regression tests for cross-validation findings of the 2026-11-24
# integration audit (M-3, M-4, M-9, C-13, X-8, and the incomplete-result
# warning of B-6). Each test starts from sim_rast(), sim_dat(), and
# kernel_prep() inputs and real model fits.

# Simulated Poisson counts at 45 points with prepared kernel inputs. With
# `scale_vars`, every layer of the raster is available to the specifications.
audit_cv_inputs <- function(scale_vars = NULL, kernel = "gaussian") {
  rs <- sim_rast(dim = 60, resolution = 10, user_seed = 11)
  sim <- sim_dat(alpha = 0.5, beta = 0.8, kernel = "gaussian", sigma = 60,
    type = "count", n_points = 45, raster_stack = rs[["bin1"]], max_D = 250,
    user_seed = 11)
  inputs <- if (is.null(scale_vars))
    kernel_prep(pts = sim$pts, raster_stack = rs[["bin1"]], max_D = 250,
      kernel = kernel, verbose = FALSE) else
    kernel_prep(pts = sim$pts, raster_stack = rs, max_D = 250, kernel = kernel,
      scale_vars = scale_vars, verbose = FALSE)
  data <- data.frame(y = sim$obs, inputs$kernel_dat)
  folds <- grouped_cv_folds(rep(1:3, length.out = nrow(data)), k = 3,
    ids = rownames(data))
  list(inputs = inputs, data = data, folds = folds)
}

# A site factor whose level "rare" occurs on one row only. The fold that
# holds that row out cannot predict it, so exactly one of three folds fails.
audit_cv_rare_level <- function(data) {
  site <- rep(c("a", "b"), length.out = nrow(data))
  site[1L] <- "rare"
  data$site <- factor(site)
  data
}

test_that("profile_glm CV stops before fitting folds for predictors it cannot fit (M-3)", {
  skip_on_cran()  # about 3 s
  x <- audit_cv_inputs(msr_vars(kmean = kernel_var("bin1"),
    m = surface_var("cont2", metric = "sq")))
  expected <- "profile_glm supports optimized kernel means and weighted surface SQ only"
  # The full-data fit refuses the same model with the same text.
  full <- stats::glm(y ~ kmean + m, poisson, x$data)
  expect_error(multiScale_optim(full, x$inputs, backend = "profile_glm", verbose = FALSE),
    expected, fixed = TRUE)
  # No fold may start a scale fit before the refusal.
  with_mocked_bindings({
    expect_error(cross_validate(y ~ kmean + m, x$data, x$inputs, x$folds, family = "poisson"),
      expected, fixed = TRUE)
    # The NB2 route uses the same analytic contract.
    expect_error(cross_validate(y ~ kmean + m, x$data, x$inputs, x$folds, family = "nbinom2"),
      expected, fixed = TRUE)
  }, multiScale_optim = function(...) stop("a fold was fitted"))
  # A kernel-mean-only formula on the same inputs passes the admission check.
  expect_invisible(.msr_cv_analytic_admission(x$inputs, "kmean", "profile_glm"))
  expect_invisible(.msr_cv_analytic_admission(x$inputs, c("kmean", "m"), "refit"))

  # A uniform (fixed) kernel is refused for the same reason as in multiScale_optim().
  fixed <- audit_cv_inputs(kernel = "fixed")
  expect_error(cross_validate(y ~ bin1, fixed$data, fixed$inputs, fixed$folds,
    family = "poisson"), "profile_glm requires smooth kernels.", fixed = TRUE)
})

test_that("a fold that fails on an estimated hard radius names radius_grid (M-4)", {
  x <- audit_cv_inputs(msr_vars(kmean = kernel_var("bin1"),
    m = surface_var("cont2", metric = "sq")))
  base <- "Training GLM or scale optimizer did not converge."
  hard <- .msr_cv_convergence_message(base, x$inputs, c("kmean", "m"), list())
  expect_true(startsWith(hard, base))
  expect_match(hard, "hard radius of m", fixed = TRUE)
  expect_match(hard, "optim_args = list(radius_grid = ...)", fixed = TRUE)
  # No suggestion when the grid is already supplied or no hard radius is used.
  expect_identical(.msr_cv_convergence_message(base, x$inputs, c("kmean", "m"),
    list(radius_grid = list(m = c(50, 100)))), base)
  expect_identical(.msr_cv_convergence_message(base, x$inputs, "kmean", list()), base)

  # End to end: a real refit whose optimizer reports a nonzero code. The
  # convergence rule is unchanged; only the recorded message is more useful.
  skip_on_cran()
  real <- multiScale_optim
  local_mocked_bindings(multiScale_optim = function(...) {
    fit <- real(...)
    fit$optim_results$convergence <- 1L
    fit
  })
  cv <- suppressWarnings(cross_validate(y ~ kmean + m, x$data, x$inputs, x$folds[1L],
    family = "poisson", backend = "refit"))
  expect_false(cv$complete)
  expect_identical(cv$errors[[1L]]$stage, "scale_fit")
  expect_match(cv$errors[[1L]]$message, "radius_grid", fixed = TRUE)
})

test_that("an incomplete cross_validate() result gives one warning that names the reason", {
  skip_on_cran()  # about 3 s
  x <- audit_cv_inputs()
  data <- audit_cv_rare_level(x$data)
  warnings <- testthat::capture_warnings(
    cv <- cross_validate(y ~ bin1 + site, data, x$inputs, x$folds, family = "poisson"))
  expect_length(warnings, 1L)
  failed <- cv$fold_results$fold[cv$fold_results$status != "scored"]
  expect_length(failed, 1L)
  expect_match(warnings, "cross_validate() returned an incomplete result: 1 of 3 folds failed",
    fixed = TRUE)
  expect_match(warnings, cv$errors[[failed]]$message, fixed = TRUE)
  expect_match(warnings, "$errors", fixed = TRUE)
  # The returned record is the same as that of a silent run (timings aside).
  quiet <- suppressWarnings(cross_validate(y ~ bin1 + site, data, x$inputs, x$folds,
    family = "poisson"))
  for (item in c("complete", "summary", "metrics", "predictions", "fold_results", "errors"))
    expect_identical(cv[[item]], quiet[[item]])
  expect_false(cv$complete)
  expect_true(is.na(cv$summary$estimate[1L]))

  # A complete run stays silent.
  expect_no_warning(complete <- cross_validate(y ~ bin1, x$data, x$inputs, x$folds,
    family = "poisson"))
  expect_true(complete$complete)

  # The reference-model function warns the same way.
  expect_warning(cross_validate_glm(y ~ site, data, x$folds, family = "poisson"),
    "cross_validate_glm() returned an incomplete result: 1 of 3 folds failed", fixed = TRUE)
})

# tibble is installed with amt (Suggests) but is not a declared dependency,
# so fetch its constructor by name after checking that it is available.
audit_as_tibble <- function(x) getExportedValue("tibble", "as_tibble")(x)

test_that("cross_validate() treats a tibble like the equivalent data frame (C-13)", {
  skip_if_not_installed("tibble")
  x <- audit_cv_inputs()
  from_frame <- cross_validate(y ~ bin1, x$data, x$inputs, x$folds, family = "poisson")
  expect_no_warning(from_tibble <- cross_validate(y ~ bin1, audit_as_tibble(x$data),
    x$inputs, x$folds, family = "poisson"))
  expect_true(from_tibble$complete)
  expect_equal(from_tibble$summary, from_frame$summary)
  expect_equal(from_tibble$predictions, from_frame$predictions)
})

# Forty simulated choice sets of four alternatives that respond to cont2.
audit_choice_inputs <- function(kernel = "gaussian") {
  rs <- sim_rast(dim = 60, resolution = 10, user_seed = 21)
  terra::crs(rs) <- "EPSG:32617"
  set.seed(21)
  n_sets <- 40L; n_alt <- 4L
  sx <- rep(stats::runif(n_sets, 200, 400), each = n_alt)
  sy <- rep(stats::runif(n_sets, 200, 400), each = n_alt)
  len <- stats::rgamma(n_sets * n_alt, shape = 2, scale = 40)
  angle <- stats::runif(n_sets * n_alt, 0, 2 * pi)
  steps <- data.frame(step_id = rep(seq_len(n_sets), each = n_alt),
    x = sx + len * cos(angle), y = sy + len * sin(angle))
  pts <- terra::vect(steps, geom = c("x", "y"), crs = "EPSG:32617")
  gaussian <- kernel_prep(pts = pts, raster_stack = rs[["cont2"]], max_D = 150,
    kernel = "gaussian", verbose = FALSE)
  z <- gaussian$kernel_dat$cont2
  steps$case_ <- 0
  for (s in seq_len(n_sets)) {
    rows <- which(steps$step_id == s)
    steps$case_[sample(rows, 1L, prob = exp(1.2 * z[rows]))] <- 1
  }
  inputs <- if (kernel == "gaussian") gaussian else
    kernel_prep(pts = pts, raster_stack = rs[["cont2"]], max_D = 150, kernel = kernel,
      verbose = FALSE)
  data <- data.frame(steps, inputs$kernel_dat)
  list(inputs = inputs, data = data,
    folds = grouped_cv_folds(data$step_id, k = 3, ids = rownames(data)))
}

test_that("cross_validate_clogit() scores a tibble like the equivalent data frame (C-13)", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("tibble")
  x <- audit_choice_inputs()
  form <- case_ ~ cont2 + strata(step_id)
  # A model fitted to the tibble is accepted by multiScale_optim(), so CV
  # must accept the same tibble.
  steps <- audit_as_tibble(x$data)
  # clogit() rewrites the call to coxph() and evaluates strata() in the
  # formula environment, so make both visible here, as survival users would.
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  native <- survival::clogit(case_ ~ cont2 + strata(step_id), data = steps, model = TRUE)
  expect_s3_class(multiScale_optim(native, x$inputs, backend = "conditional_logit",
    verbose = FALSE), "multiScaleR")
  from_frame <- cross_validate_clogit(form, x$data, x$folds, x$inputs)
  expect_no_warning(from_tibble <- cross_validate_clogit(form, steps, x$folds, x$inputs))
  expect_true(from_tibble$complete)
  expect_equal(from_tibble$summary, from_frame$summary)
  expect_equal(from_tibble$predictions, from_frame$predictions)
})

test_that("cross_validate_clogit() refuses up front or warns for deterministic refusals (C-13)", {
  skip_on_cran()  # about 3 s
  skip_if_not_installed("survival")
  x <- audit_choice_inputs()
  form <- case_ ~ cont2 + strata(step_id)
  # An expow kernel cannot be fitted by conditional_logit in any fold.
  expow <- audit_choice_inputs(kernel = "expow")
  with_mocked_bindings(
    expect_error(cross_validate_clogit(form, expow$data, expow$folds, expow$inputs),
      "conditional_logit requires Gaussian or exponential kernels.", fixed = TRUE),
    multiScale_optim = function(...) stop("a fold was fitted"))
  # A poly() term is refused inside every fold; the call now warns about it.
  expect_warning(cv <- cross_validate_clogit(case_ ~ poly(cont2, 2) + strata(step_id),
    x$data, x$folds, x$inputs),
    "cross_validate_clogit() returned an incomplete result: 3 of 3 folds failed", fixed = TRUE)
  expect_false(cv$complete)
})

test_that("cross_validate_mixed() warns when it returns unscored folds", {
  skip_if_not_installed("glmmTMB")
  set.seed(10)
  site <- factor(rep(sprintf("site%02d", 1:12), each = 6))
  surveys <- data.frame(site, habitat = rnorm(length(site)))
  surveys$count <- rpois(nrow(surveys),
    exp(0.9 + 0.4 * surveys$habitat + rnorm(nlevels(site), sd = 0.7)[site]))
  # A cover class seen only at site01 cannot be predicted when that site is held out.
  surveys$cover <- factor(ifelse(surveys$site == "site01", "rare",
    rep(c("open", "closed"), length.out = nrow(surveys))))
  folds <- grouped_cv_folds(surveys$site, k = 3, ids = rownames(surveys))
  warnings <- testthat::capture_warnings(cv <- cross_validate_mixed(
    count ~ habitat + cover + (1 | site), surveys, folds, family = "poisson",
    prediction = "new_group"))
  expect_length(warnings, 1L)
  expect_match(warnings, "cross_validate_mixed() returned an incomplete result: 1 of 3 folds failed",
    fixed = TRUE)
  expect_match(warnings, "rare", fixed = TRUE)
  expect_false(cv$complete)
})

test_that("Gaussian field CV folds whose starts disagree are reported by a warning (X-8)", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")
  # Simulated counts plus Gaussian noise with no spatial field. The two
  # declared field starts disagree in each fold, a documented refusal.
  rs <- sim_rast(dim = 100, resolution = 10, user_seed = 555)
  rs <- terra::subset(rs, c(1, 4))
  sim <- sim_dat(alpha = 0.25, beta = c(0.75, -0.5), kernel = "gaussian",
    sigma = c(75, 150), type = "count", n_points = 120, raster_stack = rs,
    max_D = 400, user_seed = 555)
  set.seed(42)
  inputs <- kernel_prep(pts = sim$pts, raster_stack = rs, max_D = 400,
    kernel = "gaussian", verbose = FALSE)
  data <- data.frame(y = 0.25 + 0.75 * sim$df$bin1 - 0.5 * sim$df$cont2 +
    rnorm(120, sd = 0.5), inputs$kernel_dat)
  xy <- sf::st_coordinates(sim$pts); rownames(xy) <- rownames(data)
  folds <- spatial_cv_folds(xy, block_size = 250, k = 3, ids = rownames(data), crs = 32617)
  audit_dir <- tempfile("audit_field_cv_")
  on.exit(unlink(audit_dir, recursive = TRUE), add = TRUE)
  warnings <- testthat::capture_warnings(cv <- cross_validate_field(y ~ bin1 + cont2,
    data, xy, folds, audit_dir = audit_dir, crs = 32617,
    covariance = "exponential"))
  failed <- sum(cv$fold_results$status != "scored")
  if (failed == 0L) skip("Every field fold scored on this platform.")
  expect_length(warnings, 1L)
  expect_match(warnings, sprintf("cross_validate_field() returned an incomplete result: %d of 3 folds failed",
    failed), fixed = TRUE)
  expect_match(warnings, cv$errors[[cv$fold_results$fold[cv$fold_results$status != "scored"][1L]]]$message,
    fixed = TRUE)
})

test_that("cv_tab() and select_cv() say what was passed instead of a CV list (M-9)", {
  x <- audit_cv_inputs()
  cv <- cross_validate(y ~ bin1, x$data, x$inputs, x$folds, family = "poisson",
    keep_selection = TRUE)
  fit <- multiScale_optim(stats::glm(y ~ bin1, poisson, x$data), x$inputs, verbose = FALSE)
  native <- stats::glm(y ~ bin1, poisson, x$data)
  # One CV result not wrapped in list()
  expect_error(cv_tab(cv, metric = "rmse"),
    "mod_list is a single multiScaleR_cv result. Wrap the candidates in list()", fixed = TRUE)
  expect_error(select_cv(cv), "mod_list is a single multiScaleR_cv result", fixed = TRUE)
  # A fitted model in the list
  expect_error(cv_tab(list(cv = cv, fit = fit), metric = "rmse"),
    "Not a CV result: fit (multiScaleR)", fixed = TRUE)
  expect_error(cv_tab(list(cv = cv, native = native), metric = "rmse"),
    "Not a CV result: native (glm)", fixed = TRUE)
  expect_error(select_cv(list(cv = cv, fit = fit)), "Not a CV result: fit (multiScaleR)",
    fixed = TRUE)
  # A fitted model in place of the list
  expect_error(cv_tab(fit, metric = "rmse"), "mod_list is a multiScaleR object", fixed = TRUE)
  # A valid list still works, and a real incomplete result keeps its own message.
  expect_s3_class(cv_tab(list(cv = cv), metric = "rmse"), "data.frame")
  failed <- suppressWarnings(cross_validate(y ~ bin1 + site, audit_cv_rare_level(x$data),
    x$inputs, x$folds, family = "poisson"))
  expect_error(cv_tab(list(failed = failed), metric = "rmse"), "failed folds cannot be dropped",
    fixed = TRUE)
})

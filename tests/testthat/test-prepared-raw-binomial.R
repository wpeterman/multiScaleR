test_that("finite-MLE geometry distinguishes overlap, separation and unresolved fits", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  X <- cbind(intercept = 1, x = c(-2, -1, 1, 2))
  overlap <- .msr_binomial_mle_status(X, c(0, 1, 0, 1))
  expect_identical(overlap$status, "finite")
  expect_gt(overlap$witness$positivity_margin, 0)
  expect_lt(max(abs(crossprod(rbind(X[c(2, 4), ], -X[c(1, 3), ]), overlap$witness$weights))), 1e-8)
  expect_false(overlap$selection_eligible)
  expect_identical(overlap$solver_version, as.character(utils::packageVersion("lpSolve")))
  expect_match(overlap$design_id, "^[[:xdigit:]]{64}$")
  complete <- .msr_binomial_mle_status(X, c(0, 0, 1, 1))
  expect_identical(complete$status, "nonexistent")
  expect_identical(complete$separation, "complete")
  quasi <- .msr_binomial_mle_status(cbind(1, c(-1, 0, 0, 1)), c(0, 0, 1, 1))
  expect_identical(quasi$status, "nonexistent")
  expect_identical(quasi$separation, "quasi_complete")
  for (y in list(rep(0, 4), rep(1, 4)))
    expect_identical(.msr_binomial_mle_status(X, y)$status, "nonexistent")
  expect_identical(.msr_binomial_mle_status(cbind(X, X[, 2]), c(0, 1, 0, 1))$status, "non_identifiable")
  near <- cbind(X, X[, 2] + c(0, 1e-14, 0, 0))
  expect_identical(.msr_binomial_mle_status(near, c(0, 1, 0, 1))$status, "unknown")
  grouped <- .msr_binomial_mle_status(X, rep(1, 4), rep(3, 4))
  expect_identical(grouped$status, "finite")
  failed <- .msr_binomial_mle_status(X, c(0, 1, 0, 1), .solver = function(...) list(status = 2L))
  expect_identical(failed$status, "unknown")
  expect_equal(failed$dual_status, 2)
  expect_equal(failed$primal_status, 2)
  failed <- .msr_binomial_mle_status(X, c(0, 1, 0, 1), .solver = function(...) stop("retained LP failure"))
  expect_identical(failed$status, "unknown")
  expect_true(all(grepl("retained LP failure", unlist(failed$solver_errors))))
  expect_error(.msr_binomial_mle_status(X, c(0, 1, NA, 1)), "aligned integer")
  expect_error(.msr_binomial_mle_status(X, c(0, .5, 0, 1)), "aligned integer")
})


test_that("prepared raw ns2 values preserve native designs, likelihoods and nominal counts", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  for (kind in c("binary", "grouped", "factor")) {
    f <- raw_binomial_fixture(grouped = kind == "grouped", factor_response = kind == "factor")
    ctx <- f$context; before <- serialize(ctx, NULL)
    cases <- list(f$reference, f$reference^1.2, .1 + .8 * f$reference)
    boundary <- f$reference; boundary[c(1L, 80L), ] <- c(0, 1, 0, 1); cases[[4L]] <- boundary
    for (features in cases) {
      d <- f$data; d$forest <- features[, 1L]; d$developed <- features[, 2L]
      native <- glm(f$formula, data = d, weights = frequency, family = binomial(),
        control = glm.control(epsilon = 1e-11, maxit = 100L))
      result <- .msr_evaluate_raw_binomial(ctx, features, c(sigma_forest = .2, sigma_developed = .4))
      expect_identical(result$finite_mle$status, "finite")
      expect_true(result$coefficient_fit_accepted)
      expect_equal(.msr_response_design_matrix(ctx$design, features), model.matrix(native),
        tolerance = 1e-10, ignore_attr = TRUE)
      expect_equal(result$coefficients, coef(native), tolerance = 1e-8)
      expected <- -sum(ctx$target$frequency * dbinom(ctx$target$successes,
        ctx$target$trials, fitted(native), log = TRUE))
      expect_equal(result$value, expected, tolerance = 1e-10)
      expect_equal(result$nominal_parameter_count, 7)
      expect_equal(result$coefficient_count, 5)
      expect_false(result$selection_eligible)
      expect_false(result$uncertainty_available)
      expect_identical(result$target_id, ctx$target_id)
    }
    expect_identical(serialize(ctx, NULL), before)
    expect_lt(length(before), 128 * 1024)
    expect_error(.msr_evaluate_raw_binomial(ctx, f$reference[80:1, ], c(sigma_forest = .2, sigma_developed = .4)), "ordered point")
    changed <- ctx; changed$target$successes[1L] <- 1 - changed$target$successes[1L]
    expect_error(.msr_evaluate_raw_binomial(changed, f$reference, c(sigma_forest = .2, sigma_developed = .4)), "target was modified")
    changed <- ctx; changed$design$multipliers[1L, 1L] <- 2
    expect_error(.msr_evaluate_raw_binomial(changed, f$reference, c(sigma_forest = .2, sigma_developed = .4)), "design or recipe was modified")
    expect_error(.msr_evaluate_raw_binomial(ctx, f$reference, c(sigma_developed = .2, sigma_forest = .4)), "prepared order")
  }
})

test_that("IRLS convergence cannot promote an unresolved finite-MLE status", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  f <- raw_binomial_fixture()
  evaluate <- .msr_evaluate_raw_binomial
  original <- environment(evaluate)
  scope <- new.env(parent = original)
  scope$.msr_binomial_mle_status <- function(...) list(status = "unknown", selection_eligible = FALSE)
  environment(evaluate) <- scope
  point <- evaluate(f$context, f$reference, c(sigma_forest = .2, sigma_developed = .4))
  expect_true(point$diagnostics$converged)
  expect_true(is.finite(point$value))
  expect_identical(point$finite_mle$status, "unknown")
  expect_false(point$coefficient_fit_accepted)
  expect_false(point$selection_eligible)
  expect_false(point$uncertainty_available)
  expect_identical(environment(.msr_evaluate_raw_binomial), original)
})

test_that("prepared multiplication and factors use the original coding and retain guards", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  f <- raw_binomial_fixture(); d <- f$data
  contrasts(d$cover) <- contr.sum(2)
  fit <- glm(y ~ forest * cover + I(developed^2) + offset(log(exposure)), data = d, family = binomial())
  ctx <- .msr_prepare_raw_binomial(fit, colnames(f$reference), f$reference, paste(rep("b", 64), collapse = ""))
  old <- options(contrasts = c("contr.treatment", "contr.poly")); on.exit(options(old))
  result <- .msr_evaluate_raw_binomial(ctx, f$reference)
  expect_true(result$coefficient_fit_accepted)
  expect_equal(result$coefficients, coef(fit), tolerance = 1e-8)
  expect_equal(result$nominal_parameter_count, length(coef(fit)))
  expect_identical(ctx$transformation, "raw_identity")
  expect_identical(ctx$source_authentication, "caller_declared_geometry_unverified")
  expect_error(aic_tab(list(result), verbose = FALSE))
  expect_error(.msr_prepare_raw_binomial(fit, colnames(f$reference), f$reference[80:1, ], paste(rep("b", 64), collapse = "")), "reference rows")
  expect_error(.msr_prepare_raw_binomial(fit, colnames(f$reference), f$reference, "unverified"), "SHA256")
  expect_error(multiScale_optim(NULL, backend = "profile_glm", covariate_transform = "raw"), "require.*refit")
  cat("\nPREPARED_BINOMIAL raw_reference=1 geometry=exact_certificate public_admission=0 uncertainty=0\n")
})

test_that("three declared coordinates retain nominal seven-plus-three counts", {
  skip_if_not_installed("lpSolve")
  skip_if_not_installed("gmp")
  f <- raw_binomial_fixture(); d <- f$data
  d$metric <- sin(d$forest * 17)^2
  fit <- glm(y ~ splines::ns(forest, knots = .5, Boundary.knots = c(0, 1)) +
    splines::ns(developed, knots = .5, Boundary.knots = c(0, 1)) +
    splines::ns(metric, knots = .5, Boundary.knots = c(0, 1)), data = d, family = binomial())
  features <- as.matrix(d[c("forest", "developed", "metric")])
  ctx <- .msr_prepare_raw_binomial(fit, colnames(features), features, paste(rep("c", 64), collapse = ""),
    c("sigma_forest", "sigma_developed", "radius_metric"))
  point <- .msr_evaluate_raw_binomial(ctx, features,
    c(sigma_forest = .2, sigma_developed = .4, radius_metric = .3))
  expect_true(point$coefficient_fit_accepted)
  expect_equal(point$coefficient_count, 7)
  expect_equal(point$nominal_parameter_count, 10)
  expect_false(point$selection_eligible)
  expect_false(point$uncertainty_available)
  expect_equal(point$value, -as.numeric(logLik(fit)), tolerance = 1e-10)
  before <- serialize(ctx, NULL); restored <- unserialize(before)
  rm(ctx, fit, d, f)
  again <- .msr_evaluate_raw_binomial(restored, features,
    c(sigma_forest = .2, sigma_developed = .4, radius_metric = .3))
  expect_equal(again$value, point$value, tolerance = 1e-12)
  expect_identical(serialize(restored, NULL), before)
})

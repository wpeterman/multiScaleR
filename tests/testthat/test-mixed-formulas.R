# Existing native integration data and existing weighted-SQ unit inputs only.
mixed_formula_fixture <- function(family = "poisson", weighted = FALSE,
                                  kernel = "gaussian") {
  if (!weighted) {
    f <- readRDS(test_path("fixtures", "joint-mixed-existing.rds"))
    dat <- f$data
    dat$site <- factor(dat$site); dat$year <- factor(dat$year)
    f$inputs$kernel <- kernel
    engine <- .msr_prepare_kernel_engine(f$inputs$binned, "habitat", kernel = kernel)
    par <- log(.23)
    dat$habitat <- .msr_eval_kernel_engine(engine, par)$value[, 1]
    formula <- count ~ habitat * year + poly(temp, 2) + (1 | site) + offset(log(effort))
  } else {
    f <- formula_weighted_fixture(kernel, surface_only = TRUE)
    dat <- f$data
    dat$site <- factor(rep(seq_len(20), each = 6))
    # Identical data construction to the established weighted mixed regression.
    set.seed(322); effect <- rnorm(20, sd = .7)
    eta <- .2 + .7 * dat$rough + .2 * dat$z
    dat$y <- rpois(120, exp(eta + effect[dat$site]))
    engine <- .msr_prepare_surface_engine(f$inputs, "rough", seq_len(nrow(dat)))
    par <- log(.4)
    formula <- y ~ rough * group + z + (1 | site) + offset(log(exposure))
  }
  model <- glmmTMB::glmmTMB(formula, dat,
    family = if (family == "nbinom2") glmmTMB::nbinom2() else poisson())
  list(model = model, data = dat, inputs = f$inputs, engine = engine, par = par)
}

test_that("mixed products retain native likelihood values and AD gradients", {
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  cases <- list(list(family = "poisson"), list(family = "nbinom2"),
    list(weighted = TRUE), list(weighted = TRUE, kernel = "exp"))
  for (case in cases) {
    f <- do.call(mixed_formula_fixture, case)
    args <- .msr_native_mixed_design(f$model, f$engine, f$engine$covariates,
      f$par, rownames(f$data))
    expect_s3_class(args$response_design, "msr_response_design")
    native <- model.matrix(f$model, component = "cond")
    expect_identical(args$response_design$columns, colnames(native))
    expect_equal(length(args$parameters$beta), ncol(native))
    obj <- do.call(.msr_mixed_rtmb, c(args, list(fixed_scale = TRUE)))
    expect_lt(abs(obj$fn(obj$par) + as.numeric(logLik(f$model))), 1e-5)
    numerical <- vapply(seq_along(obj$par), function(j) {
      hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
      (obj$fn(hi) - obj$fn(lo)) / 2e-5
    }, 0)
    # At a native optimum the score is near zero. Use absolute error rather
    # than scaling finite-difference cancellation by that nearly zero score.
    expect_lt(max(abs(as.numeric(obj$gr(obj$par)) - numerical)), 1e-7)
    away <- obj$par + .02
    numerical_away <- vapply(seq_along(away), function(j) {
      hi <- lo <- away; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
      (obj$fn(hi) - obj$fn(lo)) / 2e-5
    }, 0)
    expect_lt(max(abs(as.numeric(obj$gr(away)) - numerical_away)), 2e-5)
  }
})

test_that("mixed factor products retain all cold starts and native uncertainty boundaries", {
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  f <- mixed_formula_fixture()
  starts <- matrix(c(.1, .23, .5), ncol = 1)
  fit <- multiScale_optim(f$model, f$inputs, backend = "joint_mixed", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = starts))
  expect_s3_class(fit, "multiScaleR")
  outcomes <- fit$mixed_fit$start_results
  expect_length(outcomes, 3)
  expect_equal(exp(fit$mixed_fit$starts), starts, ignore_attr = TRUE)
  eligible <- Filter(function(x) isTRUE(x$eligible), outcomes)
  expect_gte(length(eligible), 1)
  expect_true(all(vapply(eligible, function(x) x$native_delta < 1e-5 &&
    x$design_unchanged && isTRUE(x$joint_sdreport$pdHess), logical(1))))
  expect_true(all(is.na(fit$scale_est$SE)))
  expect_error(predict(fit, se.fit = TRUE), "conditional")
  expect_equal(predict(fit), predict(fit$opt_mod, type = "response"))
  expected <- predict(fit$opt_mod, type = "link", se.fit = TRUE)
  expect_equal(predict(fit, type = "link", se.fit = TRUE, uncertainty = "conditional"), expected)
  ids <- rev(seq_len(nrow(f$data)))
  expect_equal(unname(predict(fit, f$data[ids, ], f$inputs)), rev(unname(predict(fit))), tolerance = 1e-10)
  snapshot <- prepare_prediction(fit)
  expect_equal(predict(snapshot), predict(fit))
  expect_equal(predict(snapshot, type = "link", se.fit = TRUE, uncertainty = "conditional"), expected)
  expect_equal(.msr_selection_parameter_count(fit), as.numeric(attr(logLik(fit$opt_mod), "df")) + 1L)
  raw <- fit$mixed_fit$prediction_data
  expect_true(all(c("habitat", "temp", "year", "site", "effort") %in% names(raw)))
  current <- .msr_eval_kernel_engine(fit$mixed_fit$engine, fit$mixed_fit$log_scale)$value[, 1]
  expect_equal(raw$habitat, current, ignore_attr = TRUE)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- with_mocked_bindings(suppressMessages(plot_marginal_effects(fit, length.out = 7)),
    .model_data = function(...) stop("Unexpected raw-data environment lookup."), .package = "multiScaleR")
  expect_setequal(names(plots), c("habitat", "year", "temp"))
  nd <- raw[rep(1, 7), , drop = FALSE]
  nd$habitat <- seq(min(raw$habitat), max(raw$habitat), length.out = 7)
  nd$temp <- mean(raw$temp); nd$effort <- mean(raw$effort)
  native <- predict(fit$opt_mod, nd, type = "response", re.form = NA)
  expect_equal(unname(plots$habitat$data$fit), unname(native), tolerance = 1e-10)
})

test_that("weighted mixed products retain native predictions and fixed-input snapshots", {
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  f <- mixed_formula_fixture(weighted = TRUE)
  fit <- multiScale_optim(f$model, f$inputs, backend = "joint_mixed", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = matrix(c(.2, .4, .6), ncol = 1)))
  expect_length(fit$mixed_fit$start_results, 3)
  expect_equal(predict(fit), predict(fit$opt_mod, type = "response"))
  expect_equal(predict(prepare_prediction(fit)), predict(fit))
  ids <- c(20, 1, 19)
  expect_equal(unname(predict(fit, f$data[ids, ], f$inputs)), unname(predict(fit)[ids]), tolerance = 1e-10)
  expect_true(all(is.na(fit$scale_est$SE)))
  expect_error(predict(fit, se.fit = TRUE), "conditional")
})

test_that("expanded mixed designs keep target and additional-memory guards", {
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  f <- mixed_formula_fixture()
  args <- .msr_native_mixed_design(f$model, f$engine, "habitat", f$par, rownames(f$data))
  bad <- f$model; bad$modelInfo$priors <- data.frame(prior = "normal(0,1)")
  expect_error(.msr_native_mixed_design(bad, f$engine, "habitat", f$par, rownames(f$data)), "priors")
  bad <- f$model; bad$modelInfo$allForm$ziformula <- ~habitat
  expect_error(.msr_native_mixed_design(bad, f$engine, "habitat", f$par, rownames(f$data)), "zero inflation")
  local_mocked_bindings(.msr_total_ram_bytes = function() NA_real_, .package = "multiScaleR")
  expect_error(do.call(.msr_mixed_rtmb, args), "Cannot detect RAM")
})

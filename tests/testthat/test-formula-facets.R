test_that("weighted SQ products preserve streaming, family targets and scale scores", {
  for (kernel in c("gaussian", "exp")) for (only in c(FALSE, TRUE)) {
    f <- formula_weighted_fixture(kernel, only)
    for (family in list(gaussian(), poisson(), binomial(), MASS::negative.binomial(2))) {
      dat <- f$data
      if (family$family == "gaussian") dat$y <- log1p(dat$y)
      if (family$family == "binomial") dat$y <- dat$y %% 2L
      m <- glm(f$formula, dat, family = family)
      a <- .msr_profile_contract(m, f$inputs)$objective
      b <- .msr_profile_contract(m, f$inputs, "stream")$objective
      p <- c(f$parameters, if (grepl("Negative Binomial", family$family)) log(2))
      numeric <- vapply(seq_along(p), function(j) {
        hi <- lo <- p; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
        (a$fn(hi) - a$fn(lo)) / 2e-5
      }, 0)
      expect_equal(a$gr(p), numeric, tolerance = 2e-5, ignore_attr = TRUE)
      expect_equal(b$fn(p), a$fn(p), tolerance = 1e-9)
      expect_equal(b$gr(p), a$gr(p), tolerance = 1e-8, ignore_attr = TRUE)
      native <- if (family$family == "gaussian")
        -sum(dnorm(dat$y, fitted(m), sqrt(mean((dat$y - fitted(m))^2)), log = TRUE)) else
        -as.numeric(logLik(m))
      expect_equal(a$fn(p), native, tolerance = 1e-7)
    }
  }
})

test_that("weighted SQ conditional crossings retain independent native checks", {
  skip_if_not_installed("survival")
  for (kernel in c("gaussian", "exp")) {
    f <- formula_weighted_fixture(kernel)
    dat <- f$data
    dat$case_ <- rep(c(1, 0, 0), 40)
    dat$choice <- rep(seq_len(40), each = 3)
    form <- as.formula("case_ ~ hab * rough + rough * group + z + strata(choice)",
      env = asNamespace("survival"))
    coxph <- survival::coxph
    m <- survival::clogit(form, dat, model = TRUE, x = TRUE, y = TRUE)
    contract <- .msr_clogit_contract(m, f$inputs)
    obj <- .msr_clogit_objective(contract)
    p <- c(f$parameters, rep(.1, obj$nb))
    numeric <- vapply(seq_along(p), function(j) {
      hi <- lo <- p; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
      (obj$fn(hi) - obj$fn(lo)) / 2e-5
    }, 0)
    expect_equal(obj$gr(p), numeric, tolerance = 1e-6, ignore_attr = TRUE)
    native <- .msr_clogit_native(contract, f$parameters)
    expect_equal(obj$fn(c(f$parameters, unname(coef(native)[contract$original_columns]))),
      -as.numeric(logLik(native)), tolerance = 1e-10)
  }
})

test_that("nonstandard kernel names retain prediction, plotting and coefficient penalties", {
  f <- formula_integration_fixture()
  name <- "hab + space"
  f$data[[name]] <- f$data$hab; f$data$hab <- NULL
  f$inputs$binned$vsum[[name]] <- f$inputs$binned$vsum$hab
  f$inputs$binned$csum[[name]] <- f$inputs$binned$csum$hab
  f$inputs$binned$vsum$hab <- f$inputs$binned$csum$hab <- NULL
  m <- glm(y ~ `hab + space` * group + z + offset(log(exposure)), f$data, family = poisson())
  fit <- multiScale_optim(m, f$inputs, backend = "profile_glm", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1)))
  expect_identical(rownames(fit$scale_est), name)
  expect_equal(unname(predict(fit)), unname(predict(fit$opt_mod, type = "response")), tolerance = 1e-8)
  order <- rev(seq_len(nrow(f$data)))
  expect_equal(unname(predict(fit, f$data[order, ], f$inputs)), unname(predict(fit)[order]), tolerance = 1e-10)
  expect_identical(predict(prepare_prediction(fit)), predict(fit))
  expect_equal(.msr_selection_parameter_count(fit), length(coef(fit$opt_mod)) + 1L)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- suppressMessages(plot_marginal_effects(fit, length.out = 7))
  expect_true(all(c(name, "group", "z") %in% names(plots)))
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(vapply(plots, function(p) all(is.finite(p$data$fit)), logical(1))))
})

test_that("expanded NB2 designs preserve every serial and caller-worker start", {
  skip_on_cran()
  f <- formula_integration_fixture("nbinom2")
  control <- list(starts = matrix(c(.08, .25, .6), ncol = 1), coefficient_start = "warm")
  fit <- function(control) suppressWarnings(multiScale_optim(f$model, f$inputs,
    backend = "joint_nb2", n_cores = 1, verbose = FALSE, profile_control = control))
  serial <- fit(control)
  pool <- parallel::makePSOCKcluster(1)
  on.exit(parallel::stopCluster(pool), add = TRUE)
  control$worker_pool <- pool
  parallel <- fit(control)
  a <- serial$profile_fit$start_results; b <- parallel$profile_fit$start_results
  a$seconds <- b$seconds <- NULL
  expect_equal(b, a, tolerance = 1e-10)
  expect_identical(parallel$profile_fit$starts, serial$profile_fit$starts)
  expect_identical(parallel$profile_fit$winner, serial$profile_fit$winner)
  expect_equal(predict(parallel), predict(serial), tolerance = 1e-10)
  expect_equal(parallel$uncertainty, serial$uncertainty, tolerance = 1e-10)
  expect_equal(.msr_selection_parameter_count(parallel), length(coef(parallel$opt_mod)) + 2L)
  expect_length(parallel::clusterCall(pool, Sys.getpid), 1)
})

test_that("plot annotations distinguish quoted symbols from formula operators", {
  f <- formula_integration_fixture()
  f$data[["z[1]:+"]] <- f$data$z
  m <- glm(y ~ hab * group + I(`z[1]:+`^2) + offset(log(exposure)),
    f$data, family = poisson())
  object <- structure(list(opt_mod = m,
    scl_params = list(mean = c(hab = 0), sd = c(hab = 1))), class = "multiScaleR")
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- suppressMessages(plot_marginal_effects(object, length.out = 7))
  expect_true("z[1]:+" %in% names(plots))
  expect_match(plots[["z[1]:+"]]$labels$subtitle, "Includes:", fixed = TRUE)
  expect_false(grepl("At mean", plots[["z[1]:+"]]$labels$subtitle, fixed = TRUE))
})

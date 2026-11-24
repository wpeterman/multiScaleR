# Per-scale Wald SEs for backend = "joint_mixed". Unit tests mock the TMB
# sdreport of the selected start; integration tests use the frozen
# joint-mixed fixture and a few starts.

mock_mixed_sdreport <- function(p = 2L, pd = TRUE, scale = .01, seed = 31) {
  set.seed(seed)
  names <- c(rep("kernel_par", p), "beta", "beta", "log_sd")
  k <- length(names)
  A <- matrix(rnorm(k * k), k)
  V <- scale * (crossprod(A) / k + diag(k))
  dimnames(V) <- list(names, names)
  list(cov.fixed = V, pdHess = pd)
}

mock_mixed_outcomes <- function(log_scale, objective = 10) {
  list(list(eligible = TRUE, objective = objective, log_scale = log_scale))
}

test_that("joint_mixed scale SEs come from the kernel block of cov.fixed", {
  sdr <- mock_mixed_sdreport()
  kp <- log(c(.2, .3))
  lower <- rep(log(.05), 2); upper <- rep(log(.8), 2)
  out <- .msr_mixed_scale_se(sdr, kp, lower, upper, mock_mixed_outcomes(kp), c("a", "b"))
  expect_identical(out$status, "available")
  expect_true(out$pd_hessian)
  expect_equal(unname(out$log_se), unname(sqrt(diag(sdr$cov.fixed)[1:2])), tolerance = 1e-10)
  expect_identical(unname(out$reasons), rep("available", 2))
  expect_identical(out$conditional_on, character())
  expect_false(any(out$at_bound))
})

test_that("joint_mixed withholds only the scale at a search bound", {
  sdr <- mock_mixed_sdreport()
  V <- sdr$cov.fixed
  kp <- c(log(.05), log(.3))
  lower <- rep(log(.05), 2); upper <- rep(log(.8), 2)
  out <- .msr_mixed_scale_se(sdr, kp, lower, upper, mock_mixed_outcomes(kp), c("a", "b"))
  expect_identical(out$status, "partial")
  expect_true(is.na(out$log_se[["a"]]))
  expect_identical(out$reasons[["a"]], "parameter at an optimization bound")
  expect_identical(out$reasons[["b"]], "available")
  expect_identical(out$conditional_on, "sigma:a")
  # Holding the bounded scale fixed is the Schur complement of its block.
  f <- 2:5
  conditional <- V[f, f] - V[f, 1, drop = FALSE] %*% solve(V[1, 1]) %*% V[1, f, drop = FALSE]
  expect_equal(out$log_se[["b"]], sqrt(conditional[1, 1]), tolerance = 1e-8)
  expect_lte(out$log_se[["b"]], sqrt(V[2, 2]))
  # Every scale at a bound withholds every SE without a conditional label.
  both <- .msr_mixed_scale_se(sdr, c(log(.05), log(.8)), lower, upper,
    mock_mixed_outcomes(c(log(.05), log(.8))), c("a", "b"))
  expect_identical(both$status, "all_scales_at_bound")
  expect_true(all(is.na(both$log_se)))
  expect_identical(both$conditional_on, character())
})

test_that("joint_mixed withholds every scale SE without a positive-definite sdreport", {
  kp <- log(c(.2, .3))
  lower <- rep(log(.05), 2); upper <- rep(log(.8), 2)
  outcomes <- mock_mixed_outcomes(kp)
  not_pd <- .msr_mixed_scale_se(mock_mixed_sdreport(pd = FALSE), kp, lower, upper,
    outcomes, c("a", "b"))
  expect_identical(not_pd$status, "hessian_not_positive_definite")
  expect_false(not_pd$pd_hessian)
  expect_true(all(is.na(not_pd$log_se)))
  expect_identical(unname(not_pd$reasons), rep("joint Hessian is not positive definite", 2))
  missing <- .msr_mixed_scale_se(NULL, kp, lower, upper, outcomes, c("a", "b"))
  expect_identical(missing$status, "no_sdreport")
  expect_true(all(is.na(missing$log_se)))
  bad <- mock_mixed_sdreport(); bad$cov.fixed[1, 1] <- NaN
  expect_identical(.msr_mixed_scale_se(bad, kp, lower, upper, outcomes,
    c("a", "b"))$status, "invalid_covariance")
})

test_that("joint_mixed applies the analytic per-scale checks to interior scales", {
  kp <- log(c(.2, .3))
  lower <- rep(log(.05), 2); upper <- rep(log(.8), 2)
  # A log-scale SE above one is withheld for that scale only.
  wide <- mock_mixed_sdreport()
  wide$cov.fixed[1, ] <- wide$cov.fixed[1, ] * 30
  wide$cov.fixed[, 1] <- wide$cov.fixed[, 1] * 30
  out <- .msr_mixed_scale_se(wide, kp, lower, upper, mock_mixed_outcomes(kp), c("a", "b"))
  expect_gt(sqrt(wide$cov.fixed[1, 1]), 1)
  expect_identical(out$reasons[["a"]], "log-kernel standard error exceeds one")
  expect_true(is.na(out$log_se[["a"]]))
  expect_true(is.finite(out$log_se[["b"]]))
  expect_identical(out$status, "partial")
  # A competing start inside the 95% likelihood threshold that moves scale b
  # by more than 0.5 log units withholds b; one outside the threshold does not.
  outcomes <- c(mock_mixed_outcomes(kp),
    list(list(eligible = TRUE, objective = 10.5, log_scale = kp + c(0, .8)),
         list(eligible = TRUE, objective = 20, log_scale = kp + c(.8, 0)),
         list(eligible = FALSE, error = "failed start")))
  out <- .msr_mixed_scale_se(mock_mixed_sdreport(), kp, lower, upper, outcomes, c("a", "b"))
  expect_identical(out$reasons[["b"]], "competing scale modes within the 95% likelihood threshold")
  expect_identical(out$reasons[["a"]], "available")
  expect_true(is.finite(out$log_se[["a"]]))
})

test_that("joint_mixed fits report sdreport scale SEs and log-Wald limits", {
  skip_on_cran()
  skip_if_not_installed("RTMB")
  skip_if_not_installed("TMB")
  skip_if_not_installed("glmmTMB")
  f <- readRDS(test_path("fixtures", "joint-mixed-existing.rds"))
  d <- f$data; d$site <- factor(d$site); d$year <- factor(d$year)
  e <- .msr_prepare_kernel_engine(f$inputs$binned, "habitat")
  d$habitat <- .msr_eval_kernel_engine(e, log(.23))$value[, 1]
  m <- glmmTMB::glmmTMB(count ~ temp + year + habitat + (1 | site) + offset(log(effort)),
    data = d, family = glmmTMB::nbinom2())
  fit <- multiScale_optim(m, f$inputs, backend = "joint_mixed", n_cores = 1,
    profile_control = list(starts = matrix(c(.1, .23, .5), ncol = 1)), verbose = FALSE)
  jm <- fit$diagnostics$joint_mixed
  sdr <- fit$mixed_fit$start_results[[fit$mixed_fit$selected]]$joint_sdreport
  expect_true(isTRUE(sdr$pdHess))
  expect_identical(jm$scale_se_status, "available")
  expect_false(jm$uncertainty_available)
  expect_false(fit$uncertainty$available)
  log_se <- sqrt(sdr$cov.fixed["kernel_par", "kernel_par"])
  expect_equal(fit$scale_est$SE, fit$scale_est$Mean * log_se, tolerance = 1e-8)
  expect_equal(unname(jm$log_scale_se), log_se, tolerance = 1e-8)
  tab <- summary(fit)$opt_scale
  expect_equal(unname(unlist(tab[1, c("2.5%", "97.5%")])),
    fit$scale_est$Mean * exp(c(-1, 1) * qnorm(.975) * log_se), tolerance = 1e-8)
  expect_identical(unname(attr(tab, "interval_method")), "log_wald")
  expect_equal(fit$profile_scale_est[, c("2.5%", "97.5%")], tab[, c("2.5%", "97.5%")])
  expect_true(all(is.finite(unlist(kernel_dist(fit)[1, c("2.5%", "97.5%")]))))
  # Search-extent check and follow-up starts use the point estimate.
  expect_false(jm$triggered)
  expect_false(fit$diagnostics$max_distance$triggered)
  expect_identical(fit$warn_message, 0)
  expect_equal(unname(fit$next_run$start_sigma), fit$scale_est$Mean)
  expect_equal(fit$next_run$max_D, fit$max_D)
  expect_equal(unname(fit$next_run$start_par), fit$scale_est$Mean / fit$max_D)
  # Coefficient and prediction uncertainty still condition on the scales.
  expect_error(vcov(fit), "Scale SEs are in object\\$scale_est")
  expect_equal(vcov(fit, type = "conditional"), vcov(fit$opt_mod))
  expect_error(predict(fit, se.fit = TRUE), "conditional")
  expect_output(print(fit), "condition on the selected scales")
  expect_output(print(summary(fit)), "joint observed information on the log scale")

  # A fit at the upper search bound keeps its estimate and withholds the SE.
  inputs <- f$inputs; inputs$max_D <- 100
  d$habitat <- .msr_eval_kernel_engine(e, log(.09))$value[, 1]
  m <- glmmTMB::glmmTMB(count ~ temp + year + habitat + (1 | site) + offset(log(effort)),
    data = d, family = poisson())
  bounded <- multiScale_optim(m, inputs, backend = "joint_mixed", n_cores = 1,
    profile_control = list(starts = matrix(c(.06, .09), ncol = 1)), verbose = FALSE)
  jm <- bounded$diagnostics$joint_mixed
  expect_true(jm$boundary)
  expect_true(jm$triggered)
  expect_identical(jm$scale_se_status, "all_scales_at_bound")
  # An estimate at max_D = 100 warns and suggests a larger search extent.
  expect_true(bounded$diagnostics$max_distance$triggered)
  expect_true(1L %in% bounded$warn_message)
  expect_gt(bounded$next_run$max_D, 100)
  expect_true(all(c("joint_mixed", "max_distance") %in% bounded$next_run$flags))
  expect_output(print(bounded), "Consider increasing")
  expect_identical(unname(jm$scale_se_reasons), "parameter at an optimization bound")
  expect_equal(bounded$scale_est$Mean, 100)
  expect_true(is.na(bounded$scale_est$SE))
  tab <- summary(bounded)$opt_scale
  expect_true(all(is.na(tab[, c("2.5%", "97.5%")])))
  expect_identical(unname(attr(tab, "interval_method")), "unavailable")
  expect_output(print(bounded), "Scale SE withheld for habitat \\(parameter at an optimization bound\\)")
  expect_output(print(summary(bounded)), "diagnostics\\(x\\)\\$joint_mixed")
})

test_that("joint_mixed fits saved before per-scale SEs keep missing limits", {
  legacy <- structure(list(backend = "joint_mixed",
    scale_est = data.frame(Mean = 150, SE = NA_real_, row.names = "habitat")),
    class = "multiScaleR")
  tab <- scale_ci_table(legacy, df = 10)
  expect_true(all(is.na(tab[, c("2.5%", "97.5%")])))
  expect_identical(attr(tab, "interval_method"), c(habitat = "unavailable"))
  expect_silent(.msr_print_mixed_scale_se(NULL))
})

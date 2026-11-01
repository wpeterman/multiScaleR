test_that("expression predictions and coefficient curves freeze the original reference", {
  for (surface in c(FALSE, TRUE)) {
    f <- prepared_formula_fixture(surface = surface)
    form <- count ~ I(habitat^2) * year + I(habitat * other) + poly(temp, 2) + offset(log(effort))
    fit <- multiScale_bayes(form, f$inputs, f$data, f$scales, "poisson",
      row_id = "key", method = "prepare")
    original <- multiScale_bayes(f$formula, f$inputs, f$data, f$scales, "poisson",
      row_id = "key", method = "prepare")
    expect_identical(fit$specification[c("center", "spread")],
      original$specification[c("center", "spread")])
    fit$parameters$beta[] <- seq(-.2, .2, length.out = length(fit$parameters$beta))
    fit <- fake_bayesian_draws(fit, n = 4)
    normalization <- fit$specification[c("center", "spread")]
    expected <- predict(fit, summary = FALSE)
    ids <- c(7, 1, 13, 7)
    actual <- predict(fit, newdata = f$data[ids, ], kernel_inputs = f$inputs, summary = FALSE)
    expect_equal(as.numeric(actual), as.numeric(expected[, ids, drop = FALSE]), tolerance = 1e-12)
    expect_identical(dim(actual), dim(expected[, ids, drop = FALSE]))
    expect_equal(nrow(attr(actual, "new_groups")), 0)
    grDevices::pdf(tempfile(fileext = ".pdf"))
    plots <- tryCatch(plot_marginal_effects(fit, length.out = 5), finally = grDevices::dev.off())
    expect_setequal(names(plots), c("habitat", "other", "year", "temp"))
    expect_true(all(vapply(plots, function(x) all(is.finite(x$data$fit)), logical(1))))
    expect_identical(fit$specification[c("center", "spread")], normalization)
  }
})

test_that("expression Gaussian coefficient proposals retain their full conditional density", {
  skip_if_not_installed("RTMB")
  f <- prepared_formula_fixture(surface = TRUE)
  form <- count ~ I(habitat^2) * year + I(habitat * other) + poly(temp, 2) + offset(log(effort))
  fit <- multiScale_bayes(form, f$inputs, f$data, f$scales, "gaussian",
    row_id = "key", method = "prepare")
  s <- fit$specification; p <- s$parameters
  obj <- RTMB::MakeADFun(.msr_bayes_objective(s, TRUE), p, silent = TRUE)
  keep <- which(names(obj$par) != "beta")
  outer <- matrix(rep(obj$par[keep], each = 3), nrow = 3)
  set.seed(5811)
  original <- .msr_bayes_components; calls <- 0L
  candidate <- with_mocked_bindings(.msr_bayes_beta_candidates(outer, keep, s, p$beta),
    .msr_bayes_components = function(...) { calls <<- calls + 1L; original(...) },
    .package = "multiScaleR")
  expect_identical(calls, 3L)
  replay <- .msr_bayes_beta_candidates(outer, keep, s, p$beta, given = candidate$samples)
  expect_equal(candidate$logq, replay$logq, tolerance = 1e-10)
  difference <- vapply(seq_len(3), function(i) -obj$fn(candidate$samples[i, ]), 0) - candidate$logq
  expect_lt(diff(range(difference)), 1e-8)
  expect_equal(candidate$fallbacks, 0)
})

test_that("private expression group, zero and limiting targets keep native eta", {
  skip_if_not_installed("RTMB")
  f <- prepared_formula_fixture(two = FALSE)
  form <- count ~ I(habitat^2) * year + poly(temp, 2) + offset(log(effort))
  s <- .msr_marginal_prepare(form, f$inputs, f$data, f$scales, "gaussian",
    row_id = "key", reference_scale = 250)
  p <- s$parameters
  p$beta[] <- seq(-.2, .2, length.out = length(p$beta))
  p$scale_z[] <- seq(-.1, .1, length.out = length(p$scale_z))
  eta <- prepared_formula_native_eta(s, p, f$data)$eta
  expected <- -dnorm(s$y, eta, exp(p$log_dispersion), log = TRUE)
  expect_equal(.msr_factorized_response_nll(s, p), expected, tolerance = 1e-10, ignore_attr = TRUE)
  term <- s$scales[[1]]$random[[1]]
  for (j in seq_along(term$levels)) {
    g <- .msr_factorized_group(s, j); ids <- which(term$index == j)
    gp <- p; gp$scale_z <- p$scale_z[j]
    expect_equal(g$response_design$reference_kernel,
      s$response_design$reference_kernel[ids, , drop = FALSE], tolerance = 0)
    expect_equal(.msr_factorized_response_nll(g, gp), expected[ids],
      tolerance = 1e-10, ignore_attr = TRUE)
  }
  for (kind in c("zero", "broad", "nearest")) {
    endpoint <- if (kind == "zero") .msr_marginal_reduce(s, "scale") else
      .msr_factorized_kernel_limit(s, kind)
    q <- p; q$scale_z <- q$log_scale_sd <- numeric()
    if (kind != "zero") q$scale_beta <- numeric()
    reference <- prepared_formula_native_eta(endpoint, q, f$data)$eta
    expect_equal(sum(.msr_factorized_response_nll(endpoint, q)),
      -sum(dnorm(s$y, reference, exp(q$log_dispersion), log = TRUE)), tolerance = 1e-10)
  }
})

test_that("conditional expression cold frames agree with independent public native fits", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- formula_weighted_fixture(surface_only = TRUE)
  dat <- f$data; dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- as.formula("case_ ~ I(rough^2) * group + z + offset(log(exposure)) + strata(choice)",
    env = asNamespace("survival"))
  model <- survival::clogit(form, dat, model = TRUE, x = TRUE, y = TRUE)
  contract <- .msr_clogit_contract(model, f$inputs)
  kp <- log(.3); kd <- .msr_eval_kernel_engine(contract$engine, kp)
  dat$rough <- kd$value[, "rough"]
  native <- survival::clogit(form, dat, model = TRUE, x = TRUE, y = TRUE)
  cold <- .msr_clogit_native(contract, kp)
  expect_lt(abs(as.numeric(logLik(cold)) - as.numeric(logLik(native))), 1e-8)
  obj <- .msr_clogit_objective(contract)
  beta <- coef(native)[contract$response_design$columns]
  expect_lt(abs(obj$fn(c(kp, beta)) + as.numeric(logLik(native))), 1e-8)
  p <- c(kp, rep(.1, obj$nb))
  score <- vapply(seq_along(p), function(j) {
    hi <- lo <- p; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
    (obj$fn(hi) - obj$fn(lo)) / 2e-5
  }, 0)
  expect_lt(max(abs(as.numeric(obj$gr(p)) - score)), 2e-5)
})


test_that("quoted Cox symbols are not confused with special calls", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- profile_glm_fixture(); dat <- f$data
  name <- "hab:itat [cluster(]"
  names(f$inputs$binned$vsum) <- names(f$inputs$binned$csum) <- name
  names(dat)[names(dat) == "hab"] <- name
  dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
  dat[["choice:id"]] <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  form <- as.formula(paste0("case_ ~ ", deparse1(as.name(name), backtick = TRUE),
    " * group + z + offset(log(exposure)) + strata(`choice:id`)"),
    env = asNamespace("survival"))
  model <- survival::clogit(form, dat, model = TRUE, x = TRUE, y = TRUE)
  contract <- .msr_clogit_contract(model, f$inputs)
  expect_identical(contract$strata_name, "choice:id")
  expect_identical(contract$covariates, name)
  native <- .msr_clogit_native(contract, log(.3))
  obj <- .msr_clogit_objective(contract)
  expect_lt(abs(obj$fn(c(log(.3), coef(native)[contract$response_design$columns])) +
    as.numeric(logLik(native))), 1e-8)
  # Genuine special interactions remain outside the admitted choice model.
  bad <- model
  bad$terms <- terms(as.formula("Surv(rep(1, 240), case_) ~ z * strata(`choice:id`)",
    env = asNamespace("survival")), specials = "strata")
  expect_error(.msr_clogit_contract(bad, f$inputs), "one bare strata")
})

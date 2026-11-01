test_that("prepared products preserve native designs, fixed references and normalized priors", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) for (surface in c(FALSE, TRUE))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
      f <- prepared_formula_fixture(kernel, surface)
      if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
      form <- update(f$formula, . ~ . + (1 | site))
      fit <- multiScale_bayes(form, f$inputs, f$data, f$scales, family,
        row_id = "key", method = "prepare", priors = list(scale_median = 250))
      s <- fit$specification; p <- s$parameters
      expect_s3_class(s$response_design, "msr_response_design")
      expect_identical(s$priors, .msr_bayes_priors(list(scale_median = 250), s$covariates, f$inputs))
      p$beta[] <- seq(-.2, .2, length.out = length(p$beta))
      p$scale_z[] <- seq(-.1, .2, length.out = length(p$scale_z))
      p$response_z[] <- seq(-.1, .1, length.out = length(p$response_z))
      native <- prepared_formula_native_eta(s, p, f$data)
      z <- .msr_bayes_components(s, p)
      expect_identical(s$beta_names, colnames(native$X))
      site_effect <- setNames(p$response_z, levels(f$data$site))[as.character(f$data$site)]
      expected_eta <- native$eta + exp(p$log_response_sd) * site_effect
      expect_equal(z$eta, expected_eta, tolerance = 1e-12, ignore_attr = TRUE)
      ll <- switch(family, gaussian = dnorm(s$y, expected_eta, exp(p$log_dispersion), log = TRUE),
        poisson = dpois(s$y, exp(expected_eta), log = TRUE),
        binomial = dbinom(s$y, s$trials, plogis(expected_eta), log = TRUE),
        nbinom2 = dnbinom(s$y, mu = exp(expected_eta), size = exp(p$log_dispersion), log = TRUE))
      expect_equal(.msr_bayes_objective(s)(p), -sum(ll) - prepared_formula_prior(s, p), tolerance = 1e-10)
      obj <- RTMB::MakeADFun(.msr_bayes_objective(s, TRUE), p, silent = TRUE)
      expect_equal(as.numeric(obj$fn(obj$par)), .msr_bayes_objective(s)(p), tolerance = 1e-10)
      score <- vapply(seq_along(obj$par), function(j) {
        hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
        (obj$fn(hi) - obj$fn(lo)) / 2e-5
      }, 0)
      expect_equal(as.numeric(obj$gr(obj$par)), score, tolerance = 2e-5)
      expect_true(all(is.finite(obj$he(obj$par))))
      blocks <- .msr_bayes_parameter_blocks(rbind(unlist(p), unlist(p)), p)
      expect_equal(.msr_bayes_batch_eta(s, blocks), rbind(native$eta, native$eta),
        tolerance = 1e-12, ignore_attr = TRUE)
      if (family == "gaussian") {
        engine <- if (surface) .msr_prepare_surface_engine(f$inputs, s$covariates,
          s$rows, standardize = FALSE) else .msr_prepare_kernel_engine(f$inputs$binned,
            s$covariates, s$rows, kernel, standardize = FALSE)
        raw <- .msr_eval_kernel_engine(engine, rep(log(.25), length(s$covariates)))$raw
        expect_equal(s$center, colMeans(raw))
        expect_equal(s$spread, apply(raw, 2L, sd))
      }
      if (surface) expect_error(.msr_scale_check_capability(s, "marginal_likelihood", "factorized"), "Weighted SQ")
    }
})

test_that("expanded coefficient proposals use the complete native Gaussian geometry", {
  skip_if_not_installed("RTMB")
  for (surface in c(FALSE, TRUE)) {
    f <- prepared_formula_fixture(surface = surface)
    fit <- multiScale_bayes(update(f$formula, . ~ . + (1 | site)), f$inputs, f$data, f$scales, "gaussian",
      row_id = "key", method = "prepare")
    s <- fit$specification; p <- s$parameters
    obj <- RTMB::MakeADFun(.msr_bayes_objective(s, TRUE), p, silent = TRUE)
    keep <- which(names(obj$par) != "beta")
    outer <- matrix(rep(obj$par[keep], each = 3), nrow = 3)
    outer[, names(obj$par)[keep] == "response_z"] <- .1
    set.seed(5811)
    original <- .msr_bayes_components; calls <- 0L
    candidate <- with_mocked_bindings(.msr_bayes_beta_candidates(outer, keep, s, p$beta),
      .msr_bayes_components = function(...) {
        calls <<- calls + 1L; original(...)
      }, .package = "multiScaleR")
    expect_identical(calls, 3L)
    replay <- .msr_bayes_beta_candidates(outer, keep, s, p$beta, given = candidate$samples)
    expect_equal(candidate$logq, replay$logq, tolerance = 1e-10)
    difference <- vapply(seq_len(3), function(i) -obj$fn(candidate$samples[i, ]), 0) - candidate$logq
    expect_lt(diff(range(difference)), 1e-8)
    expect_equal(candidate$fallbacks, 0)
    expect_equal(ncol(candidate$samples), length(obj$par))
  }
})

test_that("prepared interaction predictions retain basis, row identity and uncertainty contracts", {
  old <- options("contrasts")
  on.exit(options(old), add = TRUE)
  for (surface in c(FALSE, TRUE)) {
    options(old)
    f <- prepared_formula_fixture(surface = surface)
    fit <- multiScale_bayes(f$formula, f$inputs, f$data, f$scales, "poisson",
      row_id = "key", method = "prepare")
    fit$parameters$beta[] <- seq(-.2, .2, length.out = length(fit$parameters$beta))
    fit <- fake_bayesian_draws(fit, n = 4)
    expected <- predict(fit, summary = FALSE)
    normalization <- fit$specification[c("center", "spread")]
    ids <- c(7, 1, 13, 7)
    options(contrasts = c("contr.sum", "contr.poly"))
    actual <- predict(fit, newdata = f$data[ids, ], kernel_inputs = f$inputs, summary = FALSE)
    expect_equal(as.numeric(actual), as.numeric(expected[, ids, drop = FALSE]), tolerance = 1e-12)
    expect_identical(dim(actual), dim(expected[, ids, drop = FALSE]))
    expect_equal(nrow(attr(actual, "new_groups")), 0)
    expect_equal(nrow(summary(fit)$coefficients), length(fit$parameters$beta))
    grDevices::pdf(tempfile(fileext = ".pdf"))
    plots <- tryCatch(plot_marginal_effects(fit, length.out = 5), finally = grDevices::dev.off())
    expect_true(all(c("habitat", "other", "year", "temp") %in% names(plots)))
    expect_true(all(vapply(plots, function(p) all(is.finite(p$data$fit)), logical(1))))
    contract <- attr(plots, "prediction_contract")
    expect_identical(contract$inference_target, "posterior")
    expect_identical(fit$specification[c("center", "spread")], normalization)
  }
})

test_that("private group, zero and fixed-kernel paths use the expanded coefficient layout", {
  skip_if_not_installed("RTMB")
  f <- prepared_formula_fixture(two = FALSE)
  s <- .msr_marginal_prepare(f$formula, f$inputs, f$data, f$scales, "gaussian",
    row_id = "key", reference_scale = 250)
  expect_null(s$priors)
  p <- s$parameters; p$beta[] <- seq(-.15, .15, length.out = length(p$beta))
  p$scale_z[] <- seq(-.1, .1, length.out = length(p$scale_z))
  ref <- prepared_formula_native_eta(s, p, f$data)$eta
  expected <- -dnorm(s$y, ref, exp(p$log_dispersion), log = TRUE)
  expect_equal(.msr_factorized_response_nll(s, p), expected, tolerance = 1e-10, ignore_attr = TRUE)
  term <- s$scales[[1]]$random[[1]]
  for (group in seq_along(term$levels)) {
    g <- .msr_factorized_group(s, group)
    ids <- which(term$index == group)
    expect_equal(g$response_design$n, length(ids))
    gp <- p; gp$scale_z <- p$scale_z[group]
    expect_equal(.msr_factorized_response_nll(g, gp), expected[ids], tolerance = 1e-10, ignore_attr = TRUE)
  }
  zero <- .msr_marginal_reduce(s, "scale")
  q <- p; q$scale_z <- q$log_scale_sd <- numeric()
  zero_eta <- prepared_formula_native_eta(zero, q, f$data)$eta
  expect_equal(sum(.msr_factorized_response_nll(zero, q)),
    -sum(dnorm(s$y, zero_eta, exp(q$log_dispersion), log = TRUE)), tolerance = 1e-10)
  for (limit in c("broad", "nearest")) {
    endpoint <- .msr_factorized_kernel_limit(s, limit)
    q <- p; q$scale_beta <- q$scale_z <- q$log_scale_sd <- numeric()
    reference <- prepared_formula_native_eta(endpoint, q, f$data)$eta
    expect_equal(sum(.msr_factorized_response_nll(endpoint, q)),
      -sum(dnorm(s$y, reference, exp(q$log_dispersion), log = TRUE)), tolerance = 1e-10)
  }
  obj <- RTMB::MakeADFun(function(p) sum(.msr_factorized_response_nll(s, p, TRUE)), p, silent = TRUE)
  expect_equal(as.numeric(obj$fn(obj$par)), sum(expected), tolerance = 1e-10)
  score <- vapply(seq_along(obj$par), function(j) {
    hi <- lo <- obj$par; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
    (obj$fn(hi) - obj$fn(lo)) / 2e-5
  }, 0)
  expect_equal(as.numeric(obj$gr(obj$par)), score, tolerance = 2e-5)
})

test_that("expanded response AD dimensions are checked before a tape launch", {
  f <- prepared_formula_fixture()
  fit <- multiScale_bayes(f$formula, f$inputs, f$data, f$scales, "poisson",
    row_id = "key", method = "prepare")
  s <- fit$specification
  expect_gt(.msr_scale_surface_ram(s), 1024^3)
  broken <- s; broken$response_design$p <- broken$response_design$p + 1L
  expect_error(.msr_scale_surface_ram(broken), "dimensions")
  local_mocked_bindings(.msr_total_ram_bytes = function() NA_real_, .package = "multiScaleR")
  expect_error(.msr_scale_surface_ram(s), "Cannot detect RAM")
})

test_that("quoted Bayesian kernel names retain predictions and shared scale coordinates", {
  f <- prepared_formula_fixture(two = FALSE)
  name <- "hab \u03bb[]:+"
  names(f$inputs$binned$vsum) <- names(f$inputs$binned$csum) <- name
  form <- as.formula(paste0("count ~ ", deparse1(as.name(name), backtick = TRUE),
    " * year + poly(temp, 2) + offset(log(effort))"))
  fit <- multiScale_bayes(form, f$inputs, f$data,
    scales = setNames(list(~1), name), family = "poisson", row_id = "key", method = "prepare")
  expect_identical(fit$specification$covariates, name)
  expect_length(fit$parameters$scale_z, 0)
  expect_length(fit$parameters$log_scale_sd, 0)
  fit$parameters$beta[] <- seq(-.1, .2, length.out = length(fit$parameters$beta))
  fit <- fake_bayesian_draws(fit, n = 3)
  ids <- c(12, 2, 1)
  actual <- predict(fit, newdata = f$data[ids, ], kernel_inputs = f$inputs, summary = FALSE)
  expected <- predict(fit, summary = FALSE)[, ids, drop = FALSE]
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
  expect_identical(dim(actual), dim(expected))
  expect_equal(nrow(attr(actual, "new_groups")), 0)
  # Unicode text rendering needs a device that supports Unicode fonts. The
  # prediction assertions above also run on platforms without Cairo.
  if (!capabilities("cairo")) skip("Cairo is required to render Unicode labels.")
  grDevices::cairo_pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(fit, length.out = 5)
  expect_true(name %in% names(plots))
  expect_true(all(is.finite(plots[[name]]$data$fit)))
})

test_that("prepared dot and term removal optimize only active kernel variables", {
  f <- prepared_formula_fixture()
  form <- count ~ . - key - site - effort - other + offset(log(effort))
  fit <- multiScale_bayes(form, f$inputs, f$data, family = "poisson",
    row_id = "key", method = "prepare")
  expect_identical(fit$specification$covariates, "habitat")
  expect_length(fit$specification$scales, 1)
  expect_false(any(grepl("other", fit$specification$beta_names, fixed = TRUE)))
  fit <- fake_bayesian_draws(fit, n = 3)
  ids <- c(7, 1, 13)
  actual <- predict(fit, newdata = f$data[ids, ], kernel_inputs = f$inputs, summary = FALSE)
  expected <- predict(fit, summary = FALSE)[, ids, drop = FALSE]
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
  expect_identical(dim(actual), dim(expected))
  expect_equal(nrow(attr(actual, "new_groups")), 0)
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  plots <- plot_marginal_effects(fit, length.out = 5)
  expect_setequal(names(plots), c("habitat", "year", "temp"))
  bad_data <- f$data; bad_data$habitat <- bad_data$temp
  expect_error(multiScale_bayes(count ~ habitat + temp, f$inputs, bad_data,
    scales = list(habitat = ~. - count - key - site - year - temp - effort),
    row_id = "key", method = "prepare"), "Scale formulas cannot use kernel")
  expect_s3_class(multiScale_bayes(count ~ I(habitat^2), f$inputs, f$data,
    row_id = "key", method = "prepare")$specification$response_design, "msr_response_design")
  expect_error(multiScale_bayes(count ~ habitat + offset(habitat), f$inputs, f$data,
    row_id = "key", method = "prepare"), "offsets")
})


test_that("scale dot uses native active predictors and the existing proper priors", {
  f <- prepared_formula_fixture(two = FALSE)
  formula <- count ~ habitat + temp + offset(log(effort))
  dot <- ~. - count - key - site - year - effort
  fit <- multiScale_bayes(formula, f$inputs, f$data, scales = list(habitat = dot),
    family = "poisson", row_id = "key", method = "prepare")
  reference <- multiScale_bayes(formula, f$inputs, f$data, scales = list(habitat = ~temp),
    family = "poisson", row_id = "key", method = "prepare")
  expect_identical(fit$parameters, reference$parameters)
  expect_identical(fit$specification$priors, reference$specification$priors)
  expect_equal(fit$specification$scales$habitat$X, reference$specification$scales$habitat$X)
  expect_identical(fit$specification$scales$habitat$labels, reference$specification$scales$habitat$labels)
  fit$parameters$beta[3] <- 1
  fit$parameters$scale_beta[2] <- .2
  fit <- fake_bayesian_draws(fit, n = 3)
  ids <- c(7, 1, 13)
  nd <- f$data[ids, ]; nd$count <- NULL
  actual <- predict(fit, newdata = nd, kernel_inputs = f$inputs, summary = FALSE)
  expect_equal(as.numeric(actual), as.numeric(predict(fit, summary = FALSE)[, ids, drop = FALSE]),
    tolerance = 1e-12)
})

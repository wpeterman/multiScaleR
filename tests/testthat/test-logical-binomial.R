# One failure-retaining fixed-scale cache. No draws, bandwidth search or refits
# are triggered by repeated extraction/comparison assertions.
logical_binomial_cases <- local({
  value <- NULL
  function() {
    if (is.null(value)) value <<- tryCatch({
      n <- 64L; i <- seq_len(n); d <- seq(.02, 1, length.out = 12L)
      C <- 1 + outer(i, seq_along(d), function(a, b) (a + b) %% 5)
      radial <- outer(sin(i * .7), sin(d * 8)) + outer(cos(i * .9), cos(d * 13)) +
        outer(seq(-1, 1, length.out = n), d)
      ids <- paste0("b", i)
      out <- list()
      for (kernel in c("gaussian", "exp", "expow")) {
        inputs <- list(binned = list(bin_rep = d, vsum = list(hab = radial * C),
          csum = list(hab = C), point_ids = ids, n_points = n), min_D = 50,
          max_D = 800, unit_conv = 1000, kernel = kernel)
        engine <- .msr_prepare_kernel_engine(inputs$binned, kernel = kernel)
        par <- c(log(.2), if (kernel == "expow") log(1.6))
        dat <- data.frame(hab = .msr_eval_kernel_engine(engine, par)$value[, 1L],
          z = seq(-2, 2, length.out = n), group = factor(rep(c("a", "b"), n / 2L)),
          exposure = 1 + i / (4 * n), key = ids, row.names = ids)
        eta <- -.2 + .3 * dat$hab + .2 * dat$z + .1 * (dat$group == "b") +
          .1 * dat$hab * dat$z + .05 * dat$z^2 + log(dat$exposure)
        binary <- ((i * 17) %% 67 + .5) / 67 < plogis(eta)
        encodings <- list(logical = binary, numeric = as.integer(binary),
          factor = factor(binary, levels = c(FALSE, TRUE), labels = c("absent", "present")),
          counts = cbind(as.integer(binary), as.integer(!binary)),
          logical_counts = cbind(binary, !binary))
        if (kernel == "gaussian") {
          grouped <- cbind(binary, !binary | i %% 2 == 0)
          encodings$grouped_logical <- grouped
          encodings$grouped_numeric <- matrix(as.integer(grouped), ncol = 2L)
        }
        for (scheme in c("unit", "frequency")) for (name in names(encodings)) {
          rows <- dat; rows$response <- encodings[[name]]
          w <- if (scheme == "unit") rep(1, n) else rep(c(1, 2, 3, 1), length.out = n)
          warnings <- character()
          model <- withCallingHandlers(glm(response ~ hab * z + group + I(z^2) +
            offset(log(exposure)), data = rows, family = binomial(), weights = w,
            control = glm.control(epsilon = 1e-11, maxit = 100L)), warning = function(w) {
              warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
            })
          contract <- .msr_profile_contract(model, inputs)
          fitted <- contract$objective$evaluate(par)
          out[[paste(kernel, scheme, name, sep = "_")]] <- list(inputs = inputs,
            data = rows, model = model, contract = contract, point = fitted,
            par = par, weights = w, encoding = name, warnings = warnings)
        }
      }
      out
    }, error = identity)
    if (inherits(value, "error")) stop(value)
    value
  }
})

logical_glm_numbers <- function(x) x[c("coefficients", "fitted.values",
  "linear.predictors", "residuals", "deviance", "null.deviance", "aic",
  "df.residual", "df.null", "prior.weights", "weights", "rank", "converged",
  "boundary", "iter", "y")]

test_that("logical profile responses preserve fixed-scale native likelihood and rows", {
  cases <- logical_binomial_cases()
  for (entry in cases) {
    model <- entry$model; before <- serialize(logical_glm_numbers(model), NULL)
    response_before <- model.response(model$model)
    obj <- entry$contract$objective; point <- entry$point
    trials <- if (is.matrix(response_before)) rowSums(response_before) else rep(1, nrow(entry$data))
    y <- if (is.matrix(response_before)) as.numeric(response_before[, 1L]) else
      if (is.factor(response_before)) as.numeric(response_before) - 1 else as.numeric(response_before)
    expect_true(model$converged)
    expect_length(entry$warnings, 0)
    expect_equal(max(abs(crossprod(model.matrix(model), entry$weights *
      (y - trials * fitted(model))))), 0, tolerance = 1e-4)
    expect_true(all(is.finite(vcov(model))))
    expect_gt(min(eigen(vcov(model), symmetric = TRUE)$values), 0)
    expect_identical(obj$y, as.numeric(y))
    expect_equal(obj$trials, trials)
    expect_equal(obj$weights, entry$weights)
    expect_identical(entry$contract$rows, seq_len(nrow(entry$data)))
    expect_equal(point$value, -as.numeric(logLik(model)), tolerance = 1e-7)
    expect_equal(point$value, -sum(entry$weights * dbinom(y, trials,
      point$fit$fitted.values, log = TRUE)), tolerance = 1e-10)
    expect_equal(point$beta, coef(model), tolerance = 1e-8)
    expect_equal(unname(point$fit$fitted.values), unname(fitted(model)), tolerance = 1e-8)
    expect_identical(serialize(logical_glm_numbers(model), NULL), before)
    expect_identical(model.response(model$model), response_before)
    expect_equal(attr(logLik(model), "df"), length(coef(model)))
  }
  # Match the native representation: matrix initialization can retain a
  # different last IRLS covariance from a vector at the same optimum.
  for (kernel in c("gaussian", "exp", "expow")) for (scheme in c("unit", "frequency")) {
    key <- paste(kernel, scheme, sep = "_")
    reference <- cases[[paste0(key, "_numeric")]]
    for (name in c("logical", "factor")) {
      entry <- cases[[paste0(key, "_", name)]]
      expect_equal(vcov(entry$model), vcov(reference$model), tolerance = 1e-10)
      expect_identical(entry$point$fit$weights, reference$point$fit$weights)
      expect_equal(entry$point$gradient, reference$point$gradient, tolerance = 1e-10)
    }
    a <- cases[[paste0(key, "_logical_counts")]]
    b <- cases[[paste0(key, "_counts")]]
    expect_equal(vcov(a$model), vcov(b$model), tolerance = 1e-10)
    expect_identical(a$point$fit$weights, b$point$fit$weights)
  }
})

test_that("shared logical preparation preserves priors, normalized targets and fixed references", {
  cases <- logical_binomial_cases()
  fields <- c("y", "trials", "rows", "observation_ids", "V", "C", "delta",
    "center", "spread", "X_fixed", "beta_names", "parameters")
  for (kernel in c("gaussian", "exp")) {
    entry <- cases[[paste0(kernel, "_unit_numeric")]]
    formula <- response ~ hab * z + group + I(z^2) + offset(log(exposure))
    numeric <- .msr_scale_prepare(formula, entry$inputs, entry$data,
      list(hab = ~1), "binomial", "key", 200)
    posterior <- multiScale_bayes(formula, entry$inputs, entry$data,
      scales = list(hab = ~1), row_id = "key", family = "binomial",
      priors = list(scale_median = 200), method = "prepare")$specification
    for (name in c("logical", "factor", "counts", "logical_counts")) {
      rows <- cases[[paste0(kernel, "_unit_", name)]]$data
      original <- serialize(rows, NULL)
      neutral <- .msr_scale_prepare(formula, entry$inputs, rows,
        list(hab = ~1), "binomial", "key", 200)
      expect_identical(neutral[fields], numeric[fields])
      expect_null(neutral$priors)
      prepared <- multiScale_bayes(formula, entry$inputs, rows,
        scales = list(hab = ~1), row_id = "key", family = "binomial",
        priors = list(scale_median = 200), method = "prepare")$specification
      expect_identical(prepared[fields], posterior[fields])
      expect_identical(prepared$priors, posterior$priors)
      expect_equal(.msr_bayes_objective(prepared)(prepared$parameters),
        .msr_bayes_objective(posterior)(posterior$parameters), tolerance = 1e-10)
      expect_identical(serialize(rows, NULL), original)
    }
    a <- .msr_marginal_prepare(formula, entry$inputs, entry$data,
      list(hab = ~1), "binomial", "key", 200)
    rows <- entry$data; rows$response <- as.logical(rows$response)
    b <- .msr_marginal_prepare(formula, entry$inputs, rows,
      list(hab = ~1), "binomial", "key", 200)
    expect_identical(a[fields], b[fields])
    expect_null(b$priors)
    expect_equal(.msr_marginal_objective(a)(a$parameters),
      .msr_marginal_objective(b)(b$parameters), tolerance = 1e-10)
    for (outcome in c(FALSE, TRUE)) {
      # Only Gaussian input is needed for the constant-outcome MLE guard.
      if (kernel != "gaussian") next
      logical <- rows; logical$response <- rep(outcome, nrow(logical))
      binary <- logical; binary$response <- as.integer(binary$response)
      x <- .msr_marginal_prepare(formula, entry$inputs, logical,
        list(hab = ~1), "binomial", "key", 200)
      y <- .msr_marginal_prepare(formula, entry$inputs, binary,
        list(hab = ~1), "binomial", "key", 200)
      expect_match(x$identification, "No finite binomial-intercept MLE")
      expect_identical(x$identification, y$identification)
    }
  }
})

test_that("logical weighted-SQ preparation retains frozen targets and private admission guards", {
  base <- logical_binomial_cases()$gaussian_unit_numeric
  for (kernel in c("gaussian", "exp")) {
    inputs <- base$inputs; inputs$kernel <- kernel
    inputs$scale_vars <- msr_vars(hab = kernel_var("hab"),
      rough = surface_var("elev", "sq", weighted = TRUE))
    inputs$d_list <- rep(list(inputs$binned$bin_rep), nrow(base$data))
    inputs$raw_cov <- lapply(seq_len(nrow(base$data)), function(i) {
      raw <- inputs$binned$vsum$hab[i, ] / inputs$binned$csum$hab[i, ]
      data.frame(hab = raw, elev = raw^2 + sin(i * inputs$binned$bin_rep))
    })
    inputs$kernel_dat <- data.frame(row.names = row.names(base$data))
    formula <- response ~ hab + rough + z * group + offset(log(exposure))
    numeric <- base$data; logical <- numeric; logical$response <- as.logical(logical$response)
    prepare <- function(data) multiScale_bayes(formula, inputs, data,
      scales = list(hab = ~1, rough = ~1), family = "binomial", row_id = "key",
      priors = list(scale_median = 200), method = "prepare")$specification
    a <- prepare(numeric); b <- prepare(logical)
    for (name in c("y", "trials", "center", "spread", "parameters", "priors", "surface"))
      expect_identical(a[[name]], b[[name]])
    expect_equal(.msr_bayes_objective(a)(a$parameters),
      .msr_bayes_objective(b)(b$parameters), tolerance = 1e-10)
    for (data in list(numeric, logical)) expect_error(.msr_marginal_prepare(
      formula, inputs, data, list(hab = ~1, rough = ~1), "binomial", "key", 200),
      "Weighted SQ is not yet supported")
  }
})

test_that("logical formulas and invalid response boundaries keep explicit contracts", {
  entry <- logical_binomial_cases()$gaussian_unit_numeric
  response_formula <- response ~ hab + z
  prepare <- function(data, family = "binomial", formula = response_formula)
    .msr_scale_prepare(formula, entry$inputs, data, list(hab = ~1), family, "key", 200)
  dat <- entry$data; dat$response <- as.logical(dat$response)
  expect_identical(prepare(dat)$y, as.numeric(dat$response))
  expressed <- prepare(dat, formula = I(response) ~ hab + z)
  expect_identical(expressed$y, as.numeric(dat$response))
  # Native model.response removes the dimension of a one-column response.
  column <- dat; column$response <- matrix(dat$response, ncol = 1L)
  expect_identical(prepare(column)$y, as.numeric(dat$response))
  expect_identical(prepare(column)$trials, rep(1, nrow(dat)))
  for (bad in list(matrix(rep(dat$response, 3), ncol = 3L),
                  as.character(dat$response), structure(dat$response, class = "custom_boolean"))) {
    data <- dat; data$response <- bad
    expect_error(prepare(data), "binary|appropriate")
  }
  missing <- dat; missing$response[1L] <- NA
  expect_error(prepare(missing), "missing")
  expect_error(prepare(dat, "poisson"), "numeric vector")
  expect_error(prepare(dat, "gaussian"), "numeric vector")
  empty <- dat; empty$response <- cbind(rep(FALSE, nrow(dat)), rep(FALSE, nrow(dat)))
  expect_error(prepare(empty), "binary")
  for (w in c(0, .5)) {
    model <- entry$model; model$model$`(weights)` <- rep(w, nrow(model$model))
    model$model$response <- as.logical(model$model$response)
    expect_error(.msr_profile_contract(model, entry$inputs), "frequency weights")
  }
})

test_that("logical posterior preparation retains ordinary and AD point targets", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp")) {
    entry <- logical_binomial_cases()[[paste0(kernel, "_unit_numeric")]]
    formula <- response ~ hab * z + group + I(z^2) + offset(log(exposure))
    rows <- entry$data
    prepare <- function(data) multiScale_bayes(formula, entry$inputs, data,
      scales = list(hab = ~1), family = "binomial", row_id = "key",
      priors = list(scale_median = 200), method = "prepare")$specification
    a <- prepare(rows); rows$response <- as.logical(rows$response); b <- prepare(rows)
    a$parameters$beta <- b$parameters$beta <- seq(-.2, .3, length.out = length(a$parameters$beta))
    left <- RTMB::MakeADFun(.msr_bayes_objective(a, TRUE), a$parameters, silent = TRUE)
    right <- RTMB::MakeADFun(.msr_bayes_objective(b, TRUE), b$parameters, silent = TRUE)
    expect_identical(left$par, right$par)
    expect_equal(left$fn(left$par), right$fn(right$par), tolerance = 1e-10)
    expect_equal(left$gr(left$par), right$gr(right$par), tolerance = 1e-10)
    expect_equal(right$fn(right$par), .msr_bayes_objective(b)(b$parameters), tolerance = 1e-9)
  }
})

test_that("a logical public profile fit retains full starts and prediction uncertainty", {
  f <- profile_glm_fixture(n = 400, seed = 6221)
  set.seed(6222)
  f$data$presence <- rbinom(400, 1, plogis(-.4 + 1.2 * f$data$hab)) == 1L
  f$model <- glm(presence ~ hab, family = binomial(), data = f$data)
  fit <- profile_glm_fit(f)
  independent <- glm(presence ~ hab, family = binomial(), data = fit$opt_mod$model,
    control = glm.control(epsilon = 1e-11, maxit = 100L))
  expect_equal(-fit$optim_results$value, as.numeric(logLik(independent)), tolerance = 1e-7)
  expect_equal(unname(predict(fit)), unname(fitted(independent)), tolerance = 1e-8)
  expect_true(fit$uncertainty$available)
  expect_identical(fit$profile_fit$control$n_starts, 6L)
  expect_equal(nrow(fit$profile_fit$starts), 6)
  expect_true(is.logical(fit$opt_mod$model$presence))
  expect_identical(fit$profile_fit$contract$objective$y, as.numeric(f$data$presence))
  prepared <- prepare_prediction(fit)
  for (uncertainty in c("conditional", "joint")) {
    value <- predict(fit, interval = "confidence", uncertainty = uncertainty)
    expect_equal(predict(prepared, interval = "confidence", uncertainty = uncertainty), value)
    expect_true(all(is.finite(value$se.fit)))
    expect_true(all(value$lower >= 0 & value$upper <= 1))
  }
  expect_equal(aic_tab(list(fit), mod_names = "logical", verbose = FALSE)$K,
    length(coef(fit$opt_mod)) + nrow(fit$scale_est))
  expect_equal(bic_tab(list(fit), mod_names = "logical", verbose = FALSE)$K,
    length(coef(fit$opt_mod)) + nrow(fit$scale_est))
})

testthat::skip_if_not_installed("RTMB")

# Reuse a declared fixture, restoring RNG state. These tests fit no models.
likelihood_prediction_fixture <- function(...) {
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, .GlobalEnv) else
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  bayesian_scale_fixture(n_site = 4, n_year = 3, ...)
}

# Independent bin-weight reduction and manually trained polynomial design.
likelihood_prediction_reference <- function(s, p, d, inputs, training) {
  group <- match(as.character(d$year), s$scales$habitat$random[[1]]$levels)
  sigma <- exp(p$scale_beta + exp(p$log_scale_sd)*p$scale_z[group])
  rows <- match(d$key, inputs$binned$point_ids)
  raw <- vapply(seq_len(nrow(d)), function(i) {
    r <- rows[i]; distance <- inputs$binned$bin_rep*inputs$unit_conv/s$unit_conv
    w <- if (s$kernel == "gaussian") exp(-distance^2/(2*sigma[i]^2)) else exp(-distance/sigma[i])
    sum(w*inputs$binned$vsum$habitat[r, ])/sum(w*inputs$binned$csum$habitat[r, ])
  }, 0)
  d$habitat <- (raw-s$center[["habitat"]])/s$spread[["habitat"]]
  basis <- predict(poly(training$temp, 2), d$temp)
  d$b1 <- basis[, 1]; d$b2 <- basis[, 2]
  X <- model.matrix(~habitat*region+b1+b2, d)
  eta <- drop(X %*% p$beta) + log(d$effort)
  list(eta = eta, scale = sigma*s$unit_conv, metric = d$habitat)
}

test_that("private snapshots match independent kernel means and frozen formula designs", {
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- likelihood_prediction_fixture(kernel = kernel)
    f$data$region <- factor(rep(c("a", "b"), 6))
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    form <- count ~ habitat*region+poly(temp, 2)+offset(log(effort))
    s <- .msr_marginal_prepare(form, f$inputs, f$data,
      list(habitat = ~1+(1|year)), family, "key", 200)
    p <- s$parameters; p$beta <- seq(-.2, .4, length.out = length(p$beta))
    p$scale_z <- c(-.5, .1, .7); p$log_scale_sd <- log(.3)
    for (method in c("laplace", "factorized")) {
      before <- get0(".Random.seed", .GlobalEnv); original <- serialize(s, NULL)
      nd <- f$data[c(8, 2, 8, 5), ]; nd$temp <- nd$temp+.7
      # New spatial IDs and bin values, preserving original training inputs.
      prediction_inputs <- f$inputs
      prediction_inputs$binned$point_ids <- paste0(prediction_inputs$binned$point_ids, "_new")
      prediction_inputs$binned$vsum$habitat <- prediction_inputs$binned$vsum$habitat +
        .15*prediction_inputs$binned$csum$habitat
      nd$key <- paste0(nd$key, "_new")
      x <- .msr_likelihood_point_snapshot(s, p, nd, prediction_inputs, f$inputs, method = method)
      ref <- likelihood_prediction_reference(s, p, nd, prediction_inputs, f$data)
      expected <- switch(family, gaussian = ref$eta, binomial = plogis(ref$eta), exp(ref$eta))
      expect_equal(unname(x$linear_predictor), unname(ref$eta), tolerance = 1e-10)
      expect_equal(unname(x$expected_response), unname(expected), tolerance = 1e-10)
      expect_equal(unname(x$kernel_scale[, 1]), ref$scale, tolerance = 1e-10)
      expect_equal(unname(x$standardized_metric[, 1]), ref$metric, tolerance = 1e-10)
      expect_identical(x$point_ids, as.character(nd$key))
      expect_identical(names(x$linear_predictor), rownames(nd))
      expect_identical(names(x$expected_response), rownames(nd))
      expect_identical(rownames(x$kernel_scale), rownames(nd))
      expect_identical(get0(".Random.seed", .GlobalEnv), before); expect_identical(serialize(s, NULL), original)
      # Assessment outcomes and attempted replacement normalization cannot enter.
      changed <- nd; changed$count <- rep(1e8, nrow(nd))
      ki <- prediction_inputs; ki$scl_params <- list(mean = c(habitat = 1e8), sd = c(habitat = .001))
      same <- .msr_likelihood_point_snapshot(s, p, changed, ki, f$inputs, method = method)
      expect_identical(same$expected_response, x$expected_response)
      physical <- prediction_inputs; physical$unit_conv <- 500; physical$binned$bin_rep <- physical$binned$bin_rep*2
      converted <- .msr_likelihood_point_snapshot(s, p, nd, physical, f$inputs, method = method)
      expect_equal(converted$expected_response, x$expected_response, tolerance = 1e-12)
      contract <- attr(.msr_likelihood_predict(x), "prediction_contract")
      expect_identical(contract$inference_target, "conditional_likelihood")
      expect_identical(contract$diagnostics_status, "unavailable")
      expect_false(contract$metadata_establishes_eligibility)
      expect_false(contract$uncertainty$observation_noise)
      expect_identical(contract$uncertainty$intervals, "not_provided")
      expect_setequal(contract$uncertainty$omitted, c("population_parameter_uncertainty",
        "conditional_latent_variation", "new_group_variation", "observation_noise"))
      # Repeated outputs do not rebuild an engine or retain its fitted closure.
      reuse <- function(x) {
        testthat::local_mocked_bindings(.msr_bayes_components = function(...) stop("rebuild"),
          .msr_prepare_kernel_engine = function(...) stop("rebuild"), .package = "multiScaleR")
        .msr_likelihood_predict(x, "linear_predictor")
      }
      expect_equal(as.numeric(reuse(x)), unname(ref$eta), tolerance = 1e-10)
      restored <- unserialize(serialize(x, NULL))
      expect_identical(.msr_likelihood_predict(restored), .msr_likelihood_predict(x))
      forbidden <- c("specification", "V", "C", "delta", "surface", "model", "fit", "objective")
      expect_false(any(names(x) %in% forbidden))
      functions <- function(value) if (is.function(value)) TRUE else
        if (is.list(value)) any(vapply(value, functions, TRUE)) else FALSE
      expect_false(functions(x))
      expect_error(.msr_likelihood_predict(x, se.fit = TRUE), "no parameter/latent SEs")
      expect_error(.msr_likelihood_predict(x, interval = "confidence"), "no parameter/latent SEs")
      expect_error(.msr_likelihood_predict(x, newdata = nd), "prepare again")
    }
  }
})

test_that("source, group, parameter and private admission failures remain explicit", {
  f <- likelihood_prediction_fixture()
  s <- .msr_marginal_prepare(count ~ habitat+temp+(1|site), f$inputs, f$data,
    list(habitat = ~1+(1|year)), "gaussian", "key", 200)
  p <- s$parameters; p$response_z[] <- c(-.5, .2, .4, .7)
  p$scale_z[] <- c(-.4, .2, .5)
  nd <- f$data[1:3, ]; nd$year <- as.character(nd$year); nd$site <- as.character(nd$site)
  bad <- f$inputs; bad$binned$vsum$habitat[1, 3] <- bad$binned$vsum$habitat[1, 3]+1
  expect_error(.msr_likelihood_point_snapshot(s, p, nd, f$inputs, bad), "do not reproduce")
  expect_error(.msr_likelihood_point_snapshot(s, p, nd, f$inputs), "original prepared")
  bad <- p; bad$scale_z <- numeric()
  expect_error(.msr_likelihood_point_snapshot(s, bad), "parameter lengths")
  bad <- p; bad$beta <- matrix(bad$beta, nrow = 1L)
  expect_error(.msr_likelihood_point_snapshot(s, bad), "parameter lengths")
  bad <- p; bad$beta[1] <- Inf
  expect_error(.msr_likelihood_point_snapshot(s, bad), "finite parameter")
  expect_error(.msr_likelihood_point_snapshot(s, p, method = "factorized"), "coupling")
  nd$year[1] <- "future"; nd$site[1] <- "newsite"
  expect_error(.msr_likelihood_point_snapshot(s, p, nd, f$inputs, f$inputs), "New scale group")
  x <- .msr_likelihood_point_snapshot(s, p, nd, f$inputs, f$inputs,
    new_scale_levels = "zero", new_response_levels = "zero")
  expect_setequal(x$new_groups$action, "zero")
  expect_equal(x$kernel_scale[1, 1], exp(p$scale_beta)*s$unit_conv)
  zero <- .msr_likelihood_point_snapshot(s, p, scale_re.form = NA, re.form = NA)
  expect_identical(zero$prediction_contract$expected_response$effects$scale, "zero_effect_reference")
  posterior <- s; posterior$priors <- list(beta_sd = 5)
  expect_error(.msr_likelihood_point_snapshot(posterior, p), "must not carry")
  spatial <- s; spatial$scales$habitat$random[[1]]$gaussian_block <- list()
  expect_error(.msr_likelihood_point_snapshot(spatial, p), "private joint spatial")
  sq <- s; sq$surface$buffers <- list(habitat = list())
  expect_error(.msr_likelihood_point_snapshot(sq, p), "Weighted SQ")
  failed <- structure(list(specification = s, diagnostics = list(passed = FALSE)), class = "msr_marginal_prototype")
  expect_error(.msr_likelihood_prepare_prediction(failed), "original fit must pass")
  # Synthetic incomplete checkpoint tests rejection only, not fit eligibility.
  independent <- .msr_marginal_reduce(s, "response")
  integrated <- independent$parameters; integrated$scale_z <- numeric()
  partial <- structure(list(specification = independent, parameters = integrated,
    diagnostics = list(passed = TRUE), selected = 1L, zero_sd = independent$zero_sd,
    attempts = list(list(eligible = TRUE, parameters = integrated, zero_sd = independent$zero_sd))),
    class = "msr_factorized_prototype")
  for (selected in c(Inf, 1e100)) {
    invalid <- partial; invalid$selected <- selected
    expect_error(.msr_likelihood_prepare_prediction(invalid), "selected eligible start")
  }
  expect_error(.msr_likelihood_prepare_prediction(partial), "Integrated scale effects")
  expect_error(.msr_likelihood_prepare_prediction(partial, scale_effects = c(-.2, .1, .3)), "named realized scale")
})

test_that("source declarations and exact zero-SD shared-scale limits are retained", {
  f <- likelihood_prediction_fixture()
  s <- .msr_marginal_prepare(count ~ habitat+temp, f$inputs, f$data,
    list(habitat = ~1), "gaussian", "key", 200)
  # No year effect is requested or needed by the prediction scale formula.
  p <- s$parameters; nd <- f$data[1:4, ]; nd$year <- NULL
  train <- f$inputs; train$scale_vars <- msr_vars(habitat = kernel_var("terrain"))
  fresh <- train
  x <- .msr_likelihood_point_snapshot(s, p, nd, fresh, train)
  expect_identical(x$prediction_contract$kernel_scale$effects$scale, "no_active_scale_effects")
  expect_equal(unname(x$kernel_scale[, 1]), rep(exp(p$scale_beta)*s$unit_conv, 4))
  fresh$scale_vars$source <- "different"
  expect_error(.msr_likelihood_point_snapshot(s, p, nd, fresh, train), "source and metric")
  reduced <- .msr_marginal_reduce(s, character())
  expect_identical(.msr_likelihood_point_snapshot(reduced, p)$expected_response,
    .msr_likelihood_point_snapshot(s, p)$expected_response)
})

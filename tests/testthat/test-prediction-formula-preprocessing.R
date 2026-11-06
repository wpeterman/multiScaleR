# Software design/prediction references only, not a new inference study.
prediction_formula_fixture <- function(...) {
  had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, .GlobalEnv) else
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
  bayesian_scale_fixture(n_site = 4, n_year = 3, ...)
}

# Recreate original preparation without the new prediction-only learner walk.
prediction_formula_legacy <- function(code) {
  testthat::local_mocked_bindings(.msr_cv_prediction_terms = function(terms, data, qualify = FALSE) terms,
    .package = "multiScaleR")
  force(code)
}

# Match the production preparer's compact formula environments before testing
# immutability. Serializing a raw caller formula would capture the test itself.
prediction_formula_freeze_model <- function(model) {
  environment(model$formula) <- baseenv()
  environment(model$fixed) <- baseenv()
  attr(model$terms, ".Environment") <- baseenv()
  attr(model$frame, "terms") <- model$terms
  model
}

prediction_formula_basis <- function(name, training, newdata) switch(name,
  scale = (newdata-mean(training))/stats::sd(training),
  poly = predict(poly(training, 2), newdata),
  ns = predict(splines::ns(training, df = 3), newdata),
  bs = predict(splines::bs(training, df = 3), newdata))

test_that("nested fixed learners preserve training definitions and chunk invariance", {
  f <- prediction_formula_fixture(); train <- f$data
  nd <- train[c(8, 2, 8, 5), ]; nd$temp <- nd$temp+.7
  expressions <- c(scale = "scale(temp)", poly = "poly(temp, 2)",
    ns = "splines::ns(temp, df = 3)", bs = "splines::bs(temp, df = 3)")
  for (response in c(FALSE, TRUE)) for (wrapper in c("bare", "I")) for (name in names(expressions)) {
    ex <- if (wrapper == "I") paste0("I(", expressions[[name]], ")") else expressions[[name]]
    form <- as.formula(paste(if (response) "count ~" else "~", ex))
    current <- .msr_bayes_formula(form, train, response)
    legacy <- prediction_formula_legacy(.msr_bayes_formula(form, train, response))
    X <- model.matrix(current$terms, current$frame)
    expect_identical(X, model.matrix(legacy$terms, legacy$frame))
    current$prediction <- .msr_bayes_design_info(current$frame, X, response)
    legacy$prediction <- current$prediction
    expected <- cbind(1, prediction_formula_basis(name, train$temp, nd$temp))
    for (model in list(current, legacy)) {
      model <- prediction_formula_freeze_model(model)
      saved <- serialize(model, NULL)
      actual <- .msr_bayes_prediction_matrix(model, nd, X, response, train)
      expect_equal(unname(actual$X), unname(expected), tolerance = 1e-12)
      singles <- lapply(seq_len(nrow(nd)), function(i)
        .msr_bayes_prediction_matrix(model, nd[i, , drop = FALSE], X, response, train)$X)
      expect_equal(do.call(rbind, singles), actual$X, tolerance = 1e-12)
      expect_identical(serialize(model, NULL), saved)
      if (response) {
        changed <- nd; changed$count <- 1e8
        expect_identical(.msr_bayes_prediction_matrix(model, changed, X, TRUE, train)$X, actual$X)
      }
    }
  }
})

test_that("nested offsets, learner aliases and trusted spline bindings stay fixed", {
  f <- prediction_formula_fixture(); train <- f$data
  nd <- train[c(8, 2, 8, 5), ]; nd$temp <- nd$temp+.7
  form <- count ~ I(base::scale(temp, cen = TRUE, sc = FALSE)) +
    offset(log1p(abs(scale(effort))))
  for (legacy in c(FALSE, TRUE)) {
    model <- if (legacy) prediction_formula_legacy(.msr_bayes_formula(form, train, TRUE)) else
      .msr_bayes_formula(form, train, TRUE)
    X <- model.matrix(model$terms, model$frame)
    model$prediction <- .msr_bayes_design_info(model$frame, X, TRUE)
    actual <- .msr_bayes_prediction_matrix(model, nd, X, TRUE, train)
    expect_equal(unname(actual$X[, 2]), nd$temp-mean(train$temp), tolerance = 1e-12)
    expected <- log1p(abs((nd$effort-mean(train$effort))/sd(train$effort)))
    offset <- model.offset(actual$frame)
    expect_equal(as.numeric(offset), expected, tolerance = 1e-12)
    expect_identical(dim(offset), c(nrow(nd), 1L))
    expect_equal(attr(offset, "scaled:center"), mean(train$effort), tolerance = 0)
    expect_equal(attr(offset, "scaled:scale"), sd(train$effort), tolerance = 1e-12)
  }
  ns <- splines::ns
  model <- .msr_bayes_formula(count ~ I(ns(temp, df = 3)), train, TRUE)
  X <- model.matrix(model$terms, model$frame)
  model$prediction <- .msr_bayes_design_info(model$frame, X, TRUE)
  environment(model$formula) <- baseenv(); attr(model$terms, ".Environment") <- baseenv()
  actual <- .msr_bayes_prediction_matrix(model, nd, X, TRUE)
  expect_equal(unname(actual$X), unname(cbind(1, prediction_formula_basis("ns", train$temp, nd$temp))), tolerance = 1e-12)
  expect_match(paste(deparse(attr(model$terms, "predvars")), collapse = ""), "splines::ns", fixed = TRUE)
  scale <- function(x) x+4
  expect_error(.msr_bayes_formula(count ~ I(scale(temp)), train, TRUE), "Unqualified custom learner")
})

test_that("legacy reconstruction rejects missing, altered and misaligned originals", {
  f <- prediction_formula_fixture(); train <- f$data; nd <- train[1:3, ]
  model <- prediction_formula_legacy(.msr_bayes_formula(count ~ I(scale(temp)), train, TRUE))
  X <- model.matrix(model$terms, model$frame)
  model$prediction <- .msr_bayes_design_info(model$frame, X, TRUE)
  expect_error(.msr_bayes_prediction_matrix(model, nd, X, TRUE), "original training")
  altered <- train; altered$temp[1] <- altered$temp[1]+.1
  expect_error(.msr_bayes_prediction_matrix(model, nd, X, TRUE, altered), "retained training frame")
  expect_error(.msr_bayes_prediction_matrix(model, nd, X, TRUE, train[12:1, ]), "retained training frame")
  missing <- train; missing$temp <- NULL
  expect_error(.msr_bayes_prediction_matrix(model, nd, X, TRUE, missing), "original training")
  no_frame <- model; no_frame$frame <- NULL
  expect_error(.msr_bayes_prediction_matrix(no_frame, nd, X, TRUE, train), "original aligned")
  # Top-level native metadata needs no reconstruction or original data.
  model <- prediction_formula_legacy(.msr_bayes_formula(count ~ scale(temp), train, TRUE))
  X <- model.matrix(model$terms, model$frame)
  model$prediction <- .msr_bayes_design_info(model$frame, X, TRUE)
  expect_equal(unname(.msr_bayes_prediction_matrix(model, nd, X, TRUE)$X[, 2]),
    (nd$temp-mean(train$temp))/sd(train$temp), tolerance = 1e-12)
})

# The reference computes kernel means directly from bins and a training-frozen
# polynomial scale design. It does not use the prediction or fitting engine.
prediction_formula_reference <- function(s, p, training, d, inputs) {
  basis <- predict(poly(training$temp, 2), d$temp)
  scale <- exp(drop(cbind(1, basis) %*% p$scale_beta))
  rows <- match(d$key, inputs$binned$point_ids)
  habitat <- vapply(seq_len(nrow(d)), function(i) {
    distance <- inputs$binned$bin_rep*inputs$unit_conv/s$unit_conv
    w <- if (s$kernel == "gaussian") exp(-distance^2/(2*scale[i]^2)) else exp(-distance/scale[i])
    sum(w*inputs$binned$vsum$habitat[rows[i], ])/sum(w*inputs$binned$csum$habitat[rows[i], ])
  }, 0)
  d$habitat <- (habitat-s$center[["habitat"]])/s$spread[["habitat"]]
  d$fixed <- (d$temp-mean(training$temp))/sd(training$temp)
  eta <- drop(model.matrix(~habitat*region+fixed, d) %*% p$beta) +
    log1p(abs((d$effort-mean(training$effort))/sd(training$effort)))
  list(eta = eta, scale = scale*s$unit_conv)
}

test_that("Bayesian and private likelihood predictions reuse all trained formula components", {
  for (kernel in c("gaussian", "exp")) for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    f <- prediction_formula_fixture(kernel = kernel)
    f$data$region <- factor(rep(c("a", "b"), 6))
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    form <- count ~ habitat*region+I(scale(temp))+offset(log1p(abs(scale(effort))))
    scales <- list(habitat = ~1+I(poly(temp, 2)))
    fit <- multiScale_bayes(form, f$inputs, f$data, scales, family, row_id = "key",
      method = "prepare", priors = list(scale_median = 200))
    legacy <- prediction_formula_legacy(multiScale_bayes(form, f$inputs, f$data, scales, family,
      row_id = "key", method = "prepare", priors = list(scale_median = 200)))
    neutral <- .msr_marginal_prepare(form, f$inputs, f$data, scales, family, "key", 200)
    old <- prediction_formula_legacy(.msr_marginal_prepare(form, f$inputs, f$data, scales, family, "key", 200))
    expect_identical(fit$parameters, legacy$parameters)
    expect_identical(fit$specification$priors, legacy$specification$priors)
    expect_identical(neutral$parameters, old$parameters)
    expect_identical(neutral$X_fixed, old$X_fixed)
    expect_identical(neutral$response_design$multipliers, old$response_design$multipliers)
    expect_identical(neutral$offset, old$offset)
    expect_identical(neutral$scales$habitat$X, old$scales$habitat$X)
    expect_identical(lengths(neutral$parameters), lengths(old$parameters))
    p <- neutral$parameters
    p$beta <- seq(-.2, .4, length.out = length(p$beta))
    p$scale_beta <- c(log(.2), .35, -.2)
    nd <- f$data[c(8, 2, 8, 5), ]; nd$temp <- nd$temp+.7; nd$effort <- nd$effort+.1
    expected <- prediction_formula_reference(neutral, p, f$data, nd, f$inputs)
    response <- switch(family, gaussian = expected$eta, binomial = plogis(expected$eta), exp(expected$eta))
    for (s in list(neutral, old)) {
      x <- .msr_likelihood_point_snapshot(s, p, nd, f$inputs, f$inputs)
      expect_null(dim(x$linear_predictor))
      expect_null(dim(x$expected_response))
      expect_equal(unname(x$linear_predictor), unname(expected$eta), tolerance = 1e-10)
      expect_equal(unname(x$kernel_scale[, 1]), expected$scale, tolerance = 1e-10)
      expect_equal(unname(x$expected_response), unname(response), tolerance = 1e-10)
      singles <- vapply(seq_len(nrow(nd)), function(i)
        unname(.msr_likelihood_point_snapshot(s, p, nd[i, , drop = FALSE], f$inputs, f$inputs)$expected_response), 0)
      expect_equal(singles, unname(response), tolerance = 1e-10)
    }
    for (b in list(fit, legacy)) {
      b <- fake_bayesian_draws(b, 3L)
      b$draws <- lapply(p, function(value) matrix(rep(value, each = 3L), nrow = 3L))
      saved <- serialize(b, NULL)
      actual <- predict(b, nd, kernel_inputs = f$inputs, summary = FALSE)
      expect_equal(as.numeric(actual), rep(unname(response), each = 3L), tolerance = 1e-10)
      sigma <- predict_scales(b, nd, summary = FALSE)$habitat
      expect_equal(as.numeric(sigma), rep(expected$scale, each = 3L), tolerance = 1e-10)
      expect_identical(serialize(b, NULL), saved)
      if (kernel == "gaussian" && family == "gaussian") {
        get_curves <- function(object) {
          grDevices::pdf(tempfile(fileext = ".pdf"))
          on.exit(grDevices::dev.off())
          .msr_bayes_marginal_effects(object, "Expected response", 7L, 3L)
        }
        curves <- get_curves(b)
        grid <- curves$temp$data$x
        d <- data.frame(habitat = rep(mean(b$specification$response$frame$habitat), length(grid)),
          region = factor(rep(as.character(f$data$region[1]), length(grid)), levels = levels(f$data$region)),
          fixed = (grid-mean(f$data$temp))/sd(f$data$temp))
        expected_curve <- drop(model.matrix(~habitat*region+fixed, d) %*% p$beta)
        expect_equal(curves$temp$data$fit, unname(expected_curve), tolerance = 1e-10)
        expect_equal(curves$temp$data$lwr, curves$temp$data$upr, tolerance = 0)
        expect_identical(serialize(b, NULL), saved)
      }
    }
    if (family == "gaussian" && requireNamespace("RTMB", quietly = TRUE)) {
      a <- RTMB::MakeADFun(.msr_bayes_objective(fit$specification, TRUE), p, silent = TRUE)
      b <- RTMB::MakeADFun(.msr_bayes_objective(legacy$specification, TRUE), p, silent = TRUE)
      expect_equal(a$fn(a$par), b$fn(b$par), tolerance = 0)
      expect_equal(a$gr(a$par), b$gr(b$par), tolerance = 0)
      expect_equal(a$he(a$par), b$he(b$par), tolerance = 0)
    }
  }
})


test_that("dynamic learned kernel bases retain the original derivative guard", {
  f <- prediction_formula_fixture()
  expect_error(multiScale_bayes(count ~ I(scale(habitat)), f$inputs, f$data,
    row_id = "key", method = "prepare", family = "gaussian"), "Unsupported dynamic function")
})

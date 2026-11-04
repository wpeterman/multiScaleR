# The native fits are small, deterministic software cases, shared across tests.
# Retain a failed initialization so later tests cannot launch replacement fits.
logical_selection_cache <- new.env(parent = emptyenv())
logical_selection_fixture <- function() {
  if (exists("failure", logical_selection_cache, inherits = FALSE))
    stop(logical_selection_cache$failure)
  if (!exists("value", logical_selection_cache, inherits = FALSE)) {
    logical_selection_cache$value <- tryCatch({
      n <- 64L; i <- seq_len(n); d <- seq(.02, 1, length.out = 12L)
      counts <- 1 + outer(i, seq_along(d), function(a, b) (a + b) %% 5)
      radial <- outer(sin(i * .7), sin(d * 8)) + outer(cos(i * .9), cos(d * 13)) +
        outer(seq(-1, 1, length.out = n), d)
      ids <- paste0("b", i)
      bins <- list(bin_rep = d, vsum = list(hab = radial * counts),
        csum = list(hab = counts), point_ids = ids, n_points = n)
      engine <- .msr_prepare_kernel_engine(bins, kernel = "gaussian")
      dat <- data.frame(hab = .msr_eval_kernel_engine(engine, log(.2))$value[, 1L],
        z = seq(-2, 2, length.out = n), group = factor(rep(c("a", "b"), n / 2L)),
        exposure = 1 + i / (4 * n), row.names = ids)
      eta <- -.2 + .3 * dat$hab + .2 * dat$z + .1 * (dat$group == "b") +
        .1 * dat$hab * dat$z + .05 * dat$z^2 + log(dat$exposure)
      binary <- ((i * 17) %% 67 + .5) / 67 < plogis(eta)
      encodings <- list(logical = binary, numeric = as.integer(binary),
        factor = factor(binary, levels = c(FALSE, TRUE), labels = c("absent", "present")),
        counts = cbind(as.integer(binary), as.integer(!binary)),
        logical_counts = cbind(binary, !binary),
        grouped_logical = cbind(binary, !binary | i %% 2 == 0))
      encodings$grouped_numeric <- matrix(as.integer(encodings$grouped_logical), ncol = 2L)
      frequency <- rep(c(1, 2, 3, 1), length.out = n)
      calls <- 0L; fits <- list(); warnings <- list()
      one <- function(key, response, weights, reverse = FALSE, as_is = FALSE) {
        rows <- dat; rows$response <- response; rows$case_weights <- weights
        if (reverse) rows <- rows[n:1L, , drop = FALSE]
        calls <<- calls + 1L
        stopifnot(calls <= 21L)
        messages <- character()
        formula <- if (as_is) I(response) ~ hab * z + group + I(z^2) + offset(log(exposure)) else
          response ~ hab * z + group + I(z^2) + offset(log(exposure))
        fit <- withCallingHandlers(glm(formula, data = rows, family = binomial(),
          weights = case_weights, control = glm.control(epsilon = 1e-11, maxit = 100L)),
          warning = function(w) {
            messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
          })
        warnings[[key]] <<- messages
        fit
      }
      for (scheme in c("unit", "frequency")) for (name in names(encodings)) {
        key <- paste(scheme, name, sep = "_")
        fits[[key]] <- one(key, encodings[[name]],
          if (scheme == "unit") rep(1, n) else frequency)
      }
      changed <- as.integer(binary); changed[3L] <- 1L - changed[3L]
      fits$changed_response <- one("changed_response", changed, frequency)
      changed <- frequency; changed[3L] <- changed[3L] + 1
      fits$changed_weights <- one("changed_weights", as.integer(binary), changed)
      changed <- encodings$grouped_numeric; changed[3L, 2L] <- changed[3L, 2L] + 1
      fits$changed_trials <- one("changed_trials", changed, frequency)
      fits$reversed_logical <- one("reversed_logical", binary, frequency, reverse = TRUE)
      zero <- frequency; zero[4L] <- 0
      fits$zero_logical <- one("zero_logical", binary, zero)
      fits$zero_numeric <- one("zero_numeric", as.integer(binary), zero)
      fits$as_is_logical <- one("as_is_logical", binary, frequency, as_is = TRUE)
      list(fits = fits, warnings = warnings, calls = calls)
    }, error = function(e) {
      logical_selection_cache$failure <- e
      stop(e)
    })
  }
  logical_selection_cache$value
}

logical_selection_pairs <- function() {
  pairs <- list()
  for (scheme in c("unit", "frequency")) {
    for (name in c("logical", "factor")) pairs[[paste(scheme, name, sep = "_")]] <-
      c(paste(scheme, name, sep = "_"), paste0(scheme, "_numeric"))
    pairs[[paste0(scheme, "_logical_counts")]] <- c(paste0(scheme, "_logical_counts"),
      paste0(scheme, "_counts"))
    pairs[[paste0(scheme, "_grouped")]] <- c(paste0(scheme, "_grouped_logical"),
      paste0(scheme, "_grouped_numeric"))
  }
  c(pairs, list(zero_weights = c("zero_logical", "zero_numeric"),
    as_is = c("as_is_logical", "frequency_numeric")))
}

test_that("logical selection fixtures preserve verified native binomial targets", {
  f <- logical_selection_fixture()
  expect_identical(f$calls, 21L)
  for (name in names(f$fits)) {
    fit <- f$fits[[name]]
    raw <- model.response(fit$model)
    y <- if (is.matrix(raw)) as.numeric(raw[, 1L]) else if (is.factor(raw))
      as.integer(raw != levels(raw)[1L]) else as.integer(raw)
    trials <- if (is.matrix(raw)) rowSums(raw) else rep(1, NROW(raw))
    w <- model.weights(fit$model); X <- model.matrix(fit); mu <- fitted(fit)
    information <- crossprod(X, X * (w * trials * mu * (1 - mu)))
    expect_true(fit$converged, info = name)
    expect_identical(f$warnings[[name]], character(), info = name)
    expect_equal(as.numeric(logLik(fit)), sum(w * dbinom(y, trials, mu, log = TRUE)),
      tolerance = 1e-8, info = name)
    expect_lt(max(abs(crossprod(X, w * (y - trials * mu)))), 1e-4)
    expect_gt(min(eigen(information, symmetric = TRUE)$values), 0)
    expect_gt(min(eigen(vcov(fit), symmetric = TRUE)$values), 0)
    expect_equal(attr(logLik(fit), "df"), ncol(X), tolerance = 0, info = name)
    expect_equal(nobs(fit), sum(w > 0), tolerance = 0, info = name)
  }
})

test_that("equivalent logical and numeric binomial outcomes share selection signatures", {
  fits <- logical_selection_fixture()$fits
  before <- fits
  for (pair in c(logical_selection_pairs(), list(reversed = c("reversed_logical", "frequency_numeric")))) {
    a <- fits[[pair[1L]]]; b <- fits[[pair[2L]]]
    expect_equal(.msr_selection_signature(a), .msr_selection_signature(b), tolerance = 0,
      info = paste(pair, collapse = "/"))
    expect_no_error(.msr_assert_selection_targets(list(a, b)))
    expect_no_error(.msr_validate_comparable_models(list(a, b)))
    expect_equal(coef(a), coef(b), tolerance = 1e-10)
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-8)
  }
  expect_identical(fits, before)
  expect_type(model.response(fits$frequency_logical$model), "logical")
  expect_type(model.response(fits$frequency_logical_counts$model), "logical")
  expect_s3_class(fits$as_is_logical$model[[1L]], "AsIs")
  expect_type(model.response(fits$as_is_logical$model), "logical")
})

test_that("AIC and BIC retain native penalties and expected redundancy diagnostics", {
  fits <- logical_selection_fixture()$fits
  before <- fits
  for (pair in logical_selection_pairs()) {
    models <- fits[pair]
    # Equivalent encodings of the same model retain the provider's redundancy warning.
    expect_warning(aic <- aic_tab(models, AICc = FALSE, mod_names = pair), "redundant")
    expect_warning(bic <- bic_tab(models, mod_names = pair), "redundant")
    expect_s3_class(aic, "aictab")
    expect_s3_class(bic, "bictab")
    K <- vapply(models, function(x) attr(logLik(x), "df"), numeric(1))
    ll <- vapply(models, function(x) as.numeric(logLik(x)), numeric(1))
    n <- vapply(models, nobs, numeric(1))
    expect_equal(sort(aic$K), sort(unname(K)), tolerance = 0)
    expect_equal(sort(bic$K), sort(unname(K)), tolerance = 0)
    expect_equal(sort(aic$AIC), sort(unname(2 * K - 2 * ll)), tolerance = 1e-10)
    expect_equal(sort(bic$BIC), sort(unname(log(n) * K - 2 * ll)), tolerance = 1e-10)
  }
  expect_identical(fits, before)
})

test_that("logical normalization preserves target and private-admission rejections", {
  fits <- logical_selection_fixture()$fits
  pairs <- list(response = c("frequency_logical", "changed_response"),
    weights = c("frequency_logical", "changed_weights"),
    trials = c("frequency_grouped_logical", "changed_trials"))
  messages <- c(response = "response values", weights = "trial counts", trials = "trial counts")
  for (name in names(pairs)) {
    models <- fits[pairs[[name]]]
    expect_error(aic_tab(models, AICc = FALSE), messages[[name]])
    expect_error(bic_tab(models), messages[[name]])
  }
  # Repeated one-trial groups and multiplied trial totals keep their existing distinction.
  expect_error(aic_tab(fits[c("frequency_logical_counts", "frequency_numeric")], AICc = FALSE),
    "response values|likelihood weights|trial counts")
  posterior <- structure(list(), class = "multiScaleR_bayes")
  spatial <- structure(list(), class = "msr_spatial_scale_prototype")
  expect_error(aic_tab(list(posterior)), "Bayesian posterior fits")
  expect_error(bic_tab(list(spatial)), "Private spatial prototypes")
})

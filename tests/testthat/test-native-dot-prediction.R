# Six fixed native software fits, never an ecological recovery/interval study.
# A failed fixture is cached; no retry, changed seed, or altered optimizer.
native_dot_cache <- new.env(parent = emptyenv())
native_dot_fixture <- function() {
  skip_if_not_installed("survival")
  if (exists("failure", native_dot_cache, inherits = FALSE)) stop(native_dot_cache$failure)
  if (!exists("value", native_dot_cache, inherits = FALSE)) {
    native_dot_cache$value <- tryCatch({
      had <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      if (had) seed <- get(".Random.seed", .GlobalEnv)
      on.exit(if (had) assign(".Random.seed", seed, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv))
      set.seed(10171)
      n <- 180L
      raw <- data.frame(x = rnorm(n), z = runif(n, -.8, .8),
        f = factor(rep(letters[1:3], length.out = n)), omit = cos(seq_len(n)/9))
      eta <- .2 + .25 * raw$x - .2 * raw$z + .1 * (as.integer(raw$f) - 2)
      values <- list(lm = eta + rnorm(n, sd = .5),
        gaussian = eta + rnorm(n, sd = .7), binomial = rbinom(n, 1L, plogis(eta)),
        poisson = rpois(n, exp(eta)), gamma = rgamma(n, shape = 4, scale = exp(eta)/4))
      event <- rexp(n, rate = exp(.2 * raw$x)); censor <- rexp(n, rate = .25)
      form <- y ~ . - omit + x:f + I(z^2) + offset(.1 * x)
      fits <- data <- snapshots <- list(); calls <- 0L
      for (name in c(names(values), "cox")) {
        rows <- raw
        if (name == "cox") {
          # A single native Surv response avoids R's varlist warnings when dot
          # and explicit subtraction both mention multi-input response columns.
          rows$y <- survival::Surv(pmin(event, censor), as.integer(event <= censor))
          current <- form
        } else {rows$y <- values[[name]]; current <- form}
        calls <- calls + 1L; stopifnot(calls <= 6L)
        fit <- if (name == "lm") stats::lm(current, rows) else if (name == "cox")
          survival::coxph(current, rows, model = TRUE, x = TRUE, y = TRUE) else
          stats::glm(current, rows, family = switch(name, gaussian = gaussian(),
            binomial = binomial(), poisson = poisson(), gamma = Gamma(link = "log")))
        stopifnot(all(is.finite(coef(fit))), all(is.finite(vcov(fit))))
        stopifnot(min(eigen(vcov(fit), symmetric = TRUE, only.values = TRUE)$values) > 0)
        if (inherits(fit, "glm")) stopifnot(isTRUE(fit$converged), !isTRUE(fit$boundary))
        if (inherits(fit, "coxph")) stopifnot(is.null(fit$fail))
        fits[[name]] <- fit; data[[name]] <- rows
        snapshots[[name]] <- list(coefficients = coef(fit), covariance = vcov(fit),
          formula = formula(fit), terms = terms(fit), likelihood = logLik(fit),
          K = .msr_parameter_count(fit))
      }
      list(fits = fits, data = data, snapshots = snapshots, calls = calls)
    }, error = function(e) {native_dot_cache$failure <- e; stop(e)})
  }
  native_dot_cache$value
}

# The original native frame still evaluates inactive terms. Give its oracle
# explicit unused context; the wrapper's checked local terms must not need it.
native_dot_original_rows <- function(rows, name) {
  rows$omit <- 999
  rows
}

# All references call the ORIGINAL provider on the ORIGINAL fitted model.
# The six models have no data-learning terms, so single/batched rows agree.
test_that("native dot formulas validate the fitted variable expansion", {
  f <- native_dot_fixture()
  expect_identical(f$calls, 6L)
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; object <- structure(list(opt_mod = model,
      scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
    expect_setequal(.model_predictors(model), c("x", "z", "f"))
    rows <- f$data[[name]][c(2L, 17L, 53L, 102L, 179L), c("x", "z", "f"), drop = FALSE]
    attr(rows$f, "contrasts") <- NULL
    types <- if (name == "cox") c("lp", "risk") else c("response", "link")
    for (indices in list(seq_len(nrow(rows)), 3L, c(3L, 3L, 1L))) for (type in types) {
      nd <- rows[indices, , drop = FALSE]
      native_type <- if (name == "lm" && type == "link") "response" else type
      controls <- if (name == "cox") list(reference = "zero") else list()
      expected <- do.call(stats::predict, c(list(object = model, newdata = native_dot_original_rows(nd, name),
        type = native_type, se.fit = TRUE), controls))
      actual <- .msr_native_prediction_rows(object, nd, type, TRUE, controls)
      contract <- attr(actual, "prediction_contract"); attr(actual, "prediction_contract") <- NULL
      expect_equal(actual, expected, tolerance = 1e-12)
      expect_identical(contract$uncertainty$kernel_parameters, "fixed_at_estimate")
      extra <- nd; extra$omit <- 999; extra$new_prediction_column <- letters[seq_len(nrow(nd))]
      augmented <- .msr_native_prediction_rows(object, extra, type, TRUE, controls)
      attr(augmented, "prediction_contract") <- NULL
      expect_equal(augmented, expected, tolerance = 1e-12)
      thin <- list(object = object, type = type, se.fit = TRUE,
        native_args = controls, logical_names = character())
      raster <- .msr_native_raster_eval(thin, nd)
      expect_equal(as.numeric(raster$fit), as.numeric(expected$fit), tolerance = 1e-12)
      expect_equal(as.numeric(raster$se.fit), as.numeric(expected$se.fit), tolerance = 1e-12)
    }
    for (variable in c("x", "z", "f")) {
      absent <- rows; absent[[variable]] <- NULL
      # Missing prediction columns must be rejected before native evaluation.
      expect_error(.msr_native_prediction_rows(object, absent, types[1L], FALSE),
        paste0("explicitly contain.*", variable))
    }
    before <- f$snapshots[[name]]
    expect_identical(coef(model), before$coefficients)
    expect_identical(vcov(model), before$covariance)
    expect_identical(formula(model), before$formula)
    expect_identical(terms(model), before$terms)
    expect_identical(logLik(model), before$likelihood)
    expect_identical(.msr_parameter_count(model), before$K)
  }
})

test_that("native dot snapshots retain original fitted outputs and parameter counts", {
  f <- native_dot_fixture()
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; object <- structure(list(opt_mod = model,
      scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
    types <- if (name == "cox") c("lp", "risk") else c("response", "link")
    for (type in types) {
      controls <- if (name == "cox") list(reference = "zero") else list()
      snapshot <- prepare_prediction(object, component = type, native_args = controls)
      restored <- unserialize(serialize(snapshot, NULL))
      actual <- predict(restored, type = type, se.fit = TRUE, uncertainty = "conditional")
      expected <- .msr_native_prediction_rows(object, NULL, type, TRUE, controls)
      expect_equal(actual, expected, tolerance = 1e-12)
      expect_null(restored$opt_mod)
      expect_lt(as.numeric(object.size(restored)), 128 * 1024)
      expect_identical(.msr_parameter_count(model), f$snapshots[[name]]$K)
    }
  }
})

test_that("dot curves retain active offset covariates and interaction annotations", {
  f <- native_dot_fixture()
  open_null_device(); on.exit(grDevices::dev.off(), add = TRUE)
  for (name in names(f$fits)) {
    model <- f$fits[[name]]; object <- structure(list(opt_mod = model,
      scl_params = list(mean = numeric(), sd = numeric())), class = "multiScaleR")
    original <- f$data[[name]][c("x", "z", "f")]
    defaults <- lapply(original, function(x) if (is.numeric(x)) mean(x) else x[1L])
    expect_message(plots <- plot_marginal_effects(object, length.out = 7L, link = TRUE), "Interaction term")
    expect_setequal(names(plots), c("x", "z", "f"))
    expect_match(plots$x$labels$subtitle, "At mean of: f")
    expect_match(plots$z$labels$subtitle, "I\\(z\\^2\\)")
    for (variable in names(plots)) {
      count <- if (variable == "f") 3L else 7L
      rows <- as.data.frame(lapply(defaults, rep, times = count))
      rows[[variable]] <- if (variable == "f") levels(original$f) else
        seq(min(original[[variable]]), max(original[[variable]]), length.out = count)
      controls <- if (name == "cox") list(reference = "zero") else list()
      expected <- do.call(stats::predict, c(list(object = model, newdata = native_dot_original_rows(rows, name),
        type = if (name == "cox") "lp" else if (name == "lm") "response" else "link",
        se.fit = TRUE), controls))
      inverse <- if (name == "cox") exp else if (name == "lm") identity else family(model)$linkinv
      expect_equal(plots[[variable]]$data$fit, unname(inverse(expected$fit)), tolerance = 1e-12)
      expect_equal(plots[[variable]]$data$lwr,
        unname(inverse(expected$fit + qnorm(.025) * expected$se.fit)), tolerance = 1e-12)
      expect_equal(plots[[variable]]$data$upr,
        unname(inverse(expected$fit + qnorm(.975) * expected$se.fit)), tolerance = 1e-12)
    }
    expect_identical(.msr_parameter_count(model), f$snapshots[[name]]$K)
  }
})

test_that("upstream unmarked dot fitting stays outside the supported contract", {
  skip_if_not_installed("unmarked")
  frame <- unmarked::unmarkedFramePCount(matrix(0, 12, 3),
    siteCovs = data.frame(x = seq_len(12)))
  forms <- list(det = ~1, state = ~.)
  # No optimization or fabricated fitted model: this is the provider's boundary.
  expect_error(unmarked:::getDesign(frame, forms), "Variable.*not found in covariates")
})


# Metadata-only models below have no coefficients/covariance/criterion and are
# never used for inference. All evaluated designs remain local, not checkpoints.
test_that("inactive native frame pruning retains learned calls and extra classes", {
  raw <- data.frame(y = 1:24, x = sin((1:24)/4), z = cos((1:24)/7),
    f = factor(rep(letters[1:3], 8)), omit = (1:24)/10, effort = 1 + (1:24)/20)
  form <- y ~ f * I(base::scale(x)) + splines::ns(z, df = 3) + omit - omit + offset(log(effort))
  frame <- model.frame(form, raw)
  model <- structure(list(terms = terms(frame), model = frame), class = "lm")
  model <- .msr_native_fixed_metadata(model, nested = TRUE)
  classes <- attr(model$terms, "dataClasses")
  attr(model$terms, "dataClasses") <- c(classes, "(weights)" = "numeric")
  original <- model
  active <- .msr_native_active_metadata(model)
  expect_identical(model, original)
  expect_identical(.msr_native_active_metadata(active), active)
  expect_identical(attr(active$terms, "term.labels"), attr(original$terms, "term.labels"))
  expect_identical(attr(active$terms, "order"), attr(original$terms, "order"))
  expect_identical(tail(attr(active$terms, "dataClasses"), 1L), c("(weights)" = "numeric"))
  variables <- as.list(attr(original$terms, "variables"))[-1L]
  new_variables <- as.list(attr(active$terms, "variables"))[-1L]
  keep <- vapply(new_variables, function(expr) which(vapply(variables, identical, logical(1), expr)), integer(1))
  expect_identical(as.list(attr(active$terms, "predvars"))[-1L],
    as.list(attr(original$terms, "predvars"))[-1L][keep])
  expect_false("omit" %in% .msr_active_formula_variables(active$terms))
  for (indices in list(c(2, 17, 8, 17), 8L)) {
    full <- raw[indices, , drop = FALSE]; reduced <- full; reduced$omit <- NULL
    old_terms <- delete.response(original$terms); new_terms <- delete.response(active$terms)
    expected <- model.matrix(old_terms, model.frame(old_terms, full))
    actual <- model.matrix(new_terms, model.frame(new_terms, reduced))
    expect_equal(actual, expected, tolerance = 1e-12)
    expect_equal(as.numeric(model.offset(model.frame(new_terms, reduced))), log(full$effort), tolerance = 1e-12)
    reduced$effort <- NULL
    expect_error(model.frame(new_terms, reduced), "effort")
  }
  bad <- original; attr(bad$terms, "predvars") <- quote(list(y))
  expect_error(.msr_native_active_metadata(bad), "Cannot align original native")
  bad <- original; attr(bad$terms, "factors") <- attr(bad$terms, "factors")[-1L, , drop = FALSE]
  expect_error(.msr_native_active_metadata(bad), "Cannot align original native")
  for (kind in list(c("gam", "glm", "lm"), "glmmTMB", c("coxph.penal", "coxph"))) {
    untouched <- original; class(untouched) <- kind
    expect_identical(.msr_native_active_metadata(untouched), untouched)
  }
})


test_that("inactive Cox metadata preserves native specials and their indices", {
  skip_if_not_installed("survival")
  strata <- survival::strata
  raw <- data.frame(y = 1:12, x = sin(1:12), omit = cos(1:12),
    group = factor(rep(1:3, 4)), effort = 1 + (1:12)/20)
  tt <- terms(y ~ x + omit - omit + strata(group) + offset(log(effort)),
    specials = c("strata", "tt"))
  frame <- model.frame(tt, raw)
  # A terms-only prototype, without a fitted Cox coefficient or covariance.
  model <- structure(list(terms = terms(frame), model = frame), class = "coxph")
  original <- model
  active <- .msr_native_active_metadata(model)
  expect_identical(model, original)
  variables <- as.list(attr(original$terms, "variables"))[-1L]
  keep <- which(!vapply(variables, identical, logical(1), as.name("omit")))
  expect_identical(attr(active$terms, "specials")$strata,
    match(attr(original$terms, "specials")$strata, keep))
  expect_identical(attr(active$terms, "offset"), match(attr(original$terms, "offset"), keep))
  reduced <- raw; reduced$omit <- NULL
  expect_equal(model.matrix(delete.response(active$terms), model.frame(delete.response(active$terms), reduced)),
    model.matrix(delete.response(original$terms), model.frame(delete.response(original$terms), raw)), tolerance = 1e-12)
  untouched <- original; attr(untouched$terms, "specials")$tt <- 2L
  expect_identical(.msr_native_active_metadata(untouched), untouched)
  bad <- original; attr(bad$terms, "dataClasses") <- "numeric"
  expect_error(.msr_native_active_metadata(bad), "Cannot align original native")
})

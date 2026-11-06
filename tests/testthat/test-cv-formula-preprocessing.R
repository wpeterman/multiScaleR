# Deterministic contract checks, not a recovery/calibration study. Fifteen cases
# use one CV fit and two native references each. The literal-coefficient case
# adds three, for48 ordinary fits in total.
cv_formula_data <- function() {
  i <- 1:96
  d <- data.frame(x = sin(i*.73), z = cos(i*.29), g = factor(rep(letters[1:3], 32)),
    exposure = 1 + (i %% 5)/5, row.names = paste0("formulaAudit", i))
  d$continuous <- .3 + .4*d$x + .2*d$z + .15*(d$g == "b") + .1*cos(i*.91)
  d$count <- rep(c(0, 1, 0, 2, 9, 0, 1, 5, 0, 2, 14, 1), 8)
  d$binary <- as.numeric(i %% 5 < 2)
  d$positive <- exp(.2 + .3*d$x + .1*d$z)*rep(c(.3, .6, 1.1, 2.4, .8, 1.8, .4, 1.6), 12)
  d[65:96, "x"] <- d[65:96, "x"] + 1.7
  d[65:96, "z"] <- d[65:96, "z"] + 5
  d
}

test_that("ordinary CV freezes nested learned expressions and offsets for five families", {
  skip_if_not_installed("MASS")
  d <- cv_formula_data(); ids <- rownames(d)
  folds <- list(audit = list(analysis = ids[1:64], assessment = ids[65:96]))
  forms <- list(square = y ~ I(base::scale(z)^2) + x*g + offset(log(exposure)),
    offset = y ~ x*g + offset(base::scale(z)) + offset(log(exposure)),
    compound = y ~ log1p(I(base::scale(z)^2)) + stats::poly(x, 2) + g + offset(log(exposure)))
  manual <- list(square = y ~ squared + x*g + offset(log(exposure)),
    offset = y ~ x*g + offset(scaled_z) + offset(log(exposure)),
    compound = y ~ compound_z + b1 + b2 + g + offset(log(exposure)))
  for (family in c("gaussian", "poisson", "binomial", "nbinom2", "Gamma")) for (label in names(forms)) {
    d$y <- d[[switch(family, gaussian = "continuous", poisson = "count", binomial = "binary",
      nbinom2 = "count", Gamma = "positive")]]
    controls <- switch(family, nbinom2 = list(theta = 3), Gamma = list(dispersion = .4), list())
    cv <- do.call(cross_validate_glm, c(list(formula = forms[[label]], data = d, folds = folds,
      family = family, keep_fits = TRUE, keep_selection = TRUE), controls))
    expect_true(cv$complete, info = paste(family, label, paste(unlist(cv$errors), collapse = "; ")))
    if (!cv$complete) next
    fit <- cv$fits[[1]]; train <- d[ids[1:64], ]; assess <- d[ids[65:96], ]
    original <- glm(forms[[label]], train, family = fit$family, na.action = na.fail)
    basis <- poly(train$x, 2); held <- predict(basis, assess$x)
    recipe <- function(data, b) {
      data$scaled_z <- (data$z - mean(train$z))/sd(train$z)
      data$squared <- data$scaled_z^2; data$compound_z <- log1p(data$scaled_z^2)
      data$b1 <- b[, 1]; data$b2 <- b[, 2]; data
    }
    ref <- glm(manual[[label]], recipe(train, basis), family = fit$family, na.action = na.fail)
    expected <- as.numeric(predict(ref, recipe(assess, held), type = "response"))
    expect_equal(cv$predictions$prediction, expected, tolerance = 1e-8, info = paste(family, label))
    expect_identical(coef(fit), coef(original))
    expect_identical(vcov(fit), vcov(original))
    expect_identical(logLik(fit), logLik(original))
    expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
    expect_equal(unname(vcov(fit)), unname(vcov(ref)), tolerance = 1e-8)
    expect_identical(formula(fit), formula(original))
    expect_equal(cv$selection_records[[1]]$data$native_K, attr(logLik(original), "df"))
    # Predicting a row alone must not relearn its normalization or become NA.
    one <- as.numeric(predict(fit, assess[1, , drop = FALSE], type = "response"))
    expect_equal(one, expected[1], tolerance = 1e-8)
    expect_identical(cv$settings$parameter_uncertainty, "not_integrated")
    expect_identical(cv$settings$observation_intervals, "not_provided")
    compact <- cv; compact$fits <- NULL
    expect_identical(select_cv(list(model = cv), "AIC")$selection,
      select_cv(list(model = compact), "AIC")$selection)
  }
})

test_that("prediction terms retain trusted bare/namespaced splines and nested learners", {
  skip_if_not_installed("splines")
  d <- data.frame(z = seq(-2, 2, length.out = 31))
  nd <- data.frame(z = c(-1.7, .2, 1.3))
  ns <- splines::ns; bs <- splines::bs
  forms <- list(~I(scale(z)^2), ~I(base::scale(z)^2), ~I(stats::poly(z, 2)^2),
    ~I(splines::ns(z, df = 4)^2), ~I(splines::bs(z, df = 4)^2),
    ~I(ns(z, df = 4)^2), ~I(bs(z, df = 4)^2), ~I(scale(poly(z, 2))^2))
  learners <- list(scale(d$z), scale(d$z), poly(d$z, 2), ns(d$z, df = 4), bs(d$z, df = 4),
    ns(d$z, df = 4), bs(d$z, df = 4))
  for (j in seq_along(forms)) {
    frame <- model.frame(forms[[j]], d); tt <- attr(frame, "terms")
    fixed <- .msr_cv_prediction_terms(tt, d)
    actual <- model.frame(fixed, nd)[[1]]
    expected <- if (j %in% 1:2) ((nd$z - mean(d$z))/sd(d$z))^2 else if (j == 8L) {
      p <- poly(d$z, 2); new <- predict(p, nd$z)
      sweep(sweep(new, 2, colMeans(p)), 2, apply(p, 2, sd), "/")^2
    } else predict(learners[[j]], nd$z)^2
    expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), ignore_attr = TRUE, tolerance = 1e-10)
    expect_identical(attr(fixed, "variables"), attr(tt, "variables"))
    expect_identical(attr(fixed, "term.labels"), attr(tt, "term.labels"))
    expect_identical(.msr_cv_prediction_terms(fixed, d), fixed)
  }
})

test_that("metadata recovery does not evaluate unknown user functions", {
  d <- data.frame(z = seq(-2, 2, length.out = 15)); calls <- 0L
  wrapper <- function(x) { calls <<- calls + 1L; x }
  frame <- model.frame(~I(scale(wrapper(z))^2), d)
  before <- calls
  expect_error(.msr_cv_prediction_terms(attr(frame, "terms"), d), "supported pointwise")
  expect_identical(calls, before)
  # A user-defined function named scale keeps its own native semantics.
  scale <- wrapper
  frame <- model.frame(~I(scale(z)^2), d); before <- calls
  expect_identical(.msr_cv_prediction_terms(attr(frame, "terms"), d), attr(frame, "terms"))
  expect_identical(calls, before)
})

test_that("native scale abbreviations and numeric polynomial constants retain their meanings", {
  d <- data.frame(z = seq(-1, 3, length.out = 31)); nd <- data.frame(z = c(-.7, 1.2, 2.3))
  coefs <- attr(poly(d$z, 2), "coefs")
  forms <- list(~I(scale(z, cent = FALSE, sc = TRUE)^2),
    ~I(base::scale(z, cent = FALSE, sc = TRUE)^2),
    ~I(scale(z, cent = TRUE, sc = 2)^2), ~I(base::scale(z, cent = TRUE, sc = 2)^2),
    ~scale(z, cent = FALSE, sc = TRUE), ~base::scale(z, cent = FALSE, sc = TRUE),
    ~I(stats::poly(z, 2, coefs = coefs)^2),
    as.formula(paste0("~I(stats::poly(z, 2, coefs = ", paste(deparse(coefs), collapse = " "), ")^2)")),
    ~I(splines::ns(z, d = 4)^2))
  skip_if_not_installed("splines")
  for (j in seq_along(forms)) {
    tt <- attr(model.frame(forms[[j]], d), "terms")
    fixed <- .msr_cv_prediction_terms(tt, d)
    actual <- model.frame(fixed, nd)[[1]]
    rms <- sqrt(sum(d$z^2)/(nrow(d)-1))
    expected <- if (j <= 2) (nd$z/rms)^2 else if (j <= 4) ((nd$z-mean(d$z))/2)^2 else
      if (j <= 6) nd$z/rms else if (j <= 8) predict(poly(d$z, 2), nd$z)^2 else
        predict(splines::ns(d$z, df = 4), nd$z)^2
    expect_equal(unname(as.matrix(actual)), unname(as.matrix(expected)), ignore_attr = TRUE,
      tolerance = 1e-10, info = paste("argument case", j))
    expect_identical(attr(fixed, "variables"), attr(tt, "variables"))
    expect_identical(attr(fixed, "term.labels"), attr(tt, "term.labels"))
    expect_identical(.msr_cv_prediction_terms(fixed, d), fixed)
  }
})

test_that("literal polynomial coefficients work through the guarded public CV path", {
  d <- cv_formula_data(); d$y <- d$continuous; ids <- rownames(d)
  folds <- list(audit = list(analysis = ids[1:64], assessment = ids[65:96]))
  train <- d[ids[1:64], ]; assess <- d[ids[65:96], ]
  basis <- poly(train$z, 2); held <- predict(basis, assess$z)
  coefs <- attr(basis, "coefs")
  form <- as.formula(paste0("y ~ I(stats::poly(z, 2, coefs = ",
    paste(deparse(coefs), collapse = " "), ")^2) + x*g"))
  expect_setequal(all.vars(form), c("y", "z", "x", "g"))
  cv <- cross_validate_glm(form, d, folds, keep_fits = TRUE, keep_selection = TRUE)
  expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = "; "))
  if (!cv$complete) return(invisible(NULL))
  fit <- cv$fits[[1]]
  train$b1 <- basis[, 1]^2; train$b2 <- basis[, 2]^2
  assess$b1 <- held[, 1]^2; assess$b2 <- held[, 2]^2
  original <- glm(form, data = train)
  ref <- glm(y ~ b1 + b2 + x*g, data = train)
  expect_equal(cv$predictions$prediction, as.numeric(predict(ref, assess)), tolerance = 1e-8)
  expect_identical(coef(fit), coef(original)); expect_identical(vcov(fit), vcov(original))
  expect_identical(logLik(fit), logLik(original))
  expect_identical(formula(fit), formula(original))
  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = 1e-8)
  expect_equal(unname(vcov(fit)), unname(vcov(ref)), tolerance = 1e-8)
  expect_equal(cv$selection_records[[1]]$data$K, attr(logLik(original), "df"))
})

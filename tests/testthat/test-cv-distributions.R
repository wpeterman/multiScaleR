testthat::skip_if_not_installed("MASS")

# Native software references, not a response-recovery or calibration study.
cv_distribution_data <- function() {
  i <- 1:96
  d <- data.frame(x = sin(i * .73), z = cos(i * .29),
    g = factor(rep(letters[1:3], 32)), exposure = 1 + (i %% 5)/5,
    row.names = paste0("distribution", i))
  d$nb <- rep(c(0, 1, 0, 2, 9, 0, 1, 5, 0, 2, 14, 1), 8)
  d$positive <- exp(.2 + .35*d$x + .25*(d$g == "b") - .15*(d$g == "c") + .1*d$z) *
    rep(c(.3, .6, 1.1, 2.4, .8, 1.8, .4, 1.6), 12) * d$exposure
  list(data = d, folds = grouped_cv_folds(rep(letters[1:12], each = 8), 3, rownames(d)))
}

cv_distribution_cache <- new.env(parent = emptyenv())
cv_distribution_fixture <- function() {
  if (exists("failure", cv_distribution_cache, inherits = FALSE)) stop(cv_distribution_cache$failure)
  if (!exists("value", cv_distribution_cache, inherits = FALSE)) {
    cv_distribution_cache$value <- tryCatch({
      z <- cv_distribution_data(); values <- list(); planned <- 0L
      for (family in c("nbinom2", "Gamma")) for (known in c(FALSE, TRUE)) for (weighted in c(FALSE, TRUE)) for (rich in c(FALSE, TRUE)) {
        d <- z$data[c("x", "z", "g", "exposure")]
        d$y <- z$data[[if (family == "nbinom2") "nb" else "positive"]]
        form <- if (rich) y ~ x*g + poly(z, 2) + offset(log(exposure)) else
          y ~ . - exposure + I(z^2) + offset(log(exposure))
        w <- if (weighted) rep(c(0, 2, 1, 1), 24) else rep(1, 96)
        control <- list(theta = if (family == "nbinom2" && known) 3 else NULL,
          dispersion = if (family == "Gamma" && known) .4 else NULL)
        cv <- do.call(cross_validate_glm, c(list(form, d, z$folds, family, w,
          keep_fits = TRUE, keep_selection = TRUE), control))
        native_warnings <- list()
        native <- lapply(names(z$folds), function(name) {
          f <- z$folds[[name]]; messages <- character()
          train <- d[f$analysis, ]; weights <- w[match(f$analysis, rownames(d))]
          result <- withCallingHandlers(if (family == "nbinom2" && !known) do.call(MASS::glm.nb, list(formula = form, data = train, weights = weights)) else {
            fam <- if (family == "nbinom2") MASS::negative.binomial(3) else stats::Gamma("log")
            if (family == "Gamma" && known) fam$dispersion <- .4
            do.call(glm, list(formula = form, data = train, weights = weights, family = fam))
          }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") })
          native_warnings[name] <<- list(messages)
          result
        })
        names(native) <- names(z$folds)
        planned <- planned + 2L*length(z$folds)
        values[[paste(family, known, weighted, rich, sep = "_")]] <-
          list(cv = cv, native = native, native_warnings = native_warnings, data = d, w = w, known = known, family = family)
      }
      stopifnot(planned == 96L)
      list(values = values, planned_native_fits = planned)
    }, error = function(e) { cv_distribution_cache$failure <- e; stop(e) })
  }
  cv_distribution_cache$value
}

test_that("NB/Gamma CV matches native means, coefficient uncertainty, proper densities and nuisance K", {
  skip_if_not_installed("MASS")
  z <- cv_distribution_fixture()
  expect_identical(z$planned_native_fits, 96L)
  for (case in z$values) {
    cv <- case$cv
    expected <- all(vapply(case$native, function(m) isTRUE(m$converged) && is.null(m$th.warn), logical(1)))
    expect_identical(cv$complete, expected, info = paste(unlist(cv$errors), collapse = "; "))
    expect_s3_class(cv, "multiScaleR_distribution_cv")
    expect_identical(cv$settings$parameter_uncertainty, "not_integrated")
    expect_identical(cv$settings$observation_intervals, "not_provided")
    for (name in names(cv$folds)) {
      m <- cv$fits[[name]]; ref <- case$native[[name]]; d <- cv$distribution[[name]]
      p <- cv$predictions[cv$predictions$fold == name, ]; r <- cv$selection_records[[name]]
      expect_identical(cv$warnings[[name]], case$native_warnings[[name]])
      if (!isTRUE(ref$converged) || !is.null(ref$th.warn)) {
        expect_identical(cv$fold_results$status[cv$fold_results$fold == name], "failed")
        expect_match(cv$errors[[name]]$message, "theta warnings|converged")
        expect_true(all(is.na(p$prediction)))
        expect_identical(m$th.warn, ref$th.warn)
        next
      }
      expect_identical(cv$fold_results$status[cv$fold_results$fold == name], "scored")
      expect_equal(coef(m), coef(ref), tolerance = 1e-10)
      # Native Gamma summary warns about excluded zero weights. Both coefficient
      # covariance oracles retain that native convention, separately from density.
      expect_equal(suppressWarnings(vcov(m)), suppressWarnings(vcov(ref)), tolerance = 1e-10)
      expect_equal(p$prediction, as.numeric(predict(ref, case$data[p$row_id, ], type = "response")), tolerance = 1e-10)
      expect_identical(d$estimated, !case$known)
      expect_equal(r$data$K, ref$rank + !case$known)
      expect_equal(r$data$nobs, sum(case$w[match(cv$folds[[name]]$analysis, rownames(case$data))] > 0))
      expect_identical(r$version, 3L)
      mu <- fitted(ref); y <- model.response(ref$model); w <- model.weights(ref$model); used <- w > 0
      if (case$family == "nbinom2") {
        theta <- if (case$known) 3 else ref$theta
        expect_identical(d$theta, as.numeric(theta))
        train_ll <- sum(w[used]*dnbinom(y[used], size = theta, mu = mu[used], log = TRUE))
        loss <- -dnbinom(p$observed, size = theta, mu = p$prediction, log = TRUE)
      } else {
        phi <- if (case$known) .4 else ref$deviance/sum(w)
        expect_identical(d$dispersion, phi)
        train_ll <- sum(w[used]*dgamma(y[used], shape = 1/phi, scale = mu[used]*phi, log = TRUE))
        loss <- -dgamma(p$observed, shape = 1/phi, scale = p$prediction*phi, log = TRUE)
      }
      expect_equal(as.numeric(logLik(m)), train_ll, tolerance = 1e-10)
      expect_equal(r$data$logLik, train_ll, tolerance = 1e-10)
      expect_equal(AIC(m), -2*train_ll + 2*(ref$rank + !case$known), tolerance = 1e-10)
      loss[p$weight == 0] <- NA_real_
      expect_equal(p$neg_log_score, loss, tolerance = 1e-12)
      expect_equal(p$eligible, p$weight > 0)
      expect_equal(p$rmse[p$eligible], (p$observed[p$eligible] - p$prediction[p$eligible])^2)
      expect_identical(rownames(m$model), cv$folds[[name]]$analysis)
      expect_true(is.null(r$error))
    }
    if (cv$complete) expect_equal(cv$summary$estimate[cv$summary$metric == "neg_log_score"],
      weighted.mean(cv$predictions$neg_log_score[cv$predictions$eligible], cv$predictions$weight[cv$predictions$eligible])) else
      expect_true(all(is.na(cv$summary$estimate)))
  }
})

test_that("NB/Gamma compact and full native model selection agree for every nuisance convention", {
  z <- cv_distribution_fixture()
  for (case in z$values) {
    full <- case$cv; compact <- full; compact$fits <- NULL
    if (!full$complete) {
      expect_error(select_cv(list(candidate = full), "AIC"), "failed folds")
      expect_error(select_cv(list(candidate = compact), "AIC"), "failed folds")
      next
    }
    for (criterion in c("AIC", "AICc", "BIC")) {
      a <- select_cv(list(candidate = full), criterion)
      b <- select_cv(list(candidate = compact), criterion)
      expect_true(a$complete, info = paste(unlist(a$errors), collapse = "; "))
      expect_true(b$complete, info = paste(unlist(b$errors), collapse = "; "))
      for (field in c("predictions", "summary", "distribution", "diagnostics", "selection")) expect_identical(a[[field]], b[[field]])
      for (name in names(a$folds)) expect_identical(a$settings$fold_settings[[name]], full$settings)
      expect_identical(cv_tab(list(candidate = compact, selected = b), "neg_log_score")$Delta, c(0, 0))
    }
    small <- function(x) {
      if (is.environment(x) || typeof(x) == "externalptr" || inherits(x, c("glm", "multiScaleR", "formula"))) return(FALSE)
      if (is.list(x)) return(all(vapply(x, small, logical(1))))
      TRUE
    }
    expect_true(small(compact$selection_records))
  }
})

test_that("nuisance parameters and training transformations never learn from assessment rows", {
  z <- cv_distribution_data(); fold <- z$folds[1L]; held <- fold[[1]]$assessment
  for (family in c("nbinom2", "Gamma")) {
    d <- z$data; d$y <- d[[if (family == "nbinom2") "nb" else "positive"]]
    form <- y ~ x*g + base::scale(z) + offset(log(exposure))
    a <- cross_validate_glm(form, d, fold, family, keep_fits = TRUE)
    changed <- d; changed[held, "y"] <- changed[held, "y"]*2 + 1; changed[held, "z"] <- changed[held, "z"] + 7
    b <- cross_validate_glm(form, changed, fold, family, keep_fits = TRUE)
    expect_true(a$complete, info = paste(unlist(a$errors), collapse = "; "))
    expect_true(b$complete, info = paste(unlist(b$errors), collapse = "; "))
    expect_identical(a$distribution, b$distribution)
    expect_identical(coef(a$fits[[1]]), coef(b$fits[[1]]))
    # Independent design: fixed training mean/SD and native fitted coefficient order.
    m <- a$fits[[1]]; train <- d[fold[[1]]$analysis, ]; nd <- d[held, ]
    design <- model.matrix(~ x*g + scaled, transform(nd, scaled = (z - mean(train$z))/sd(train$z)))
    expected <- exp(as.vector(design %*% coef(m)) + log(nd$exposure))
    expect_equal(a$predictions$prediction, expected, tolerance = 1e-10)
  }
})

test_that("NB/Gamma support, backend, likelihood and record failures remain explicit", {
  z <- cv_distribution_data(); d <- z$data; folds <- z$folds[1L]
  expect_error(cross_validate_glm(nb ~ x, d, folds, "nbinom2", theta = 0), "finite positive")
  expect_error(cross_validate_glm(positive ~ x, d, folds, "Gamma", dispersion = Inf), "finite positive")
  expect_error(cross_validate_glm(nb ~ x, d, folds, "poisson", theta = 3), "theta")
  expect_error(cross_validate_glm(positive ~ x, d, folds, "Gamma", theta = 3), "theta")
  expect_error(cross_validate_glm(nb ~ x, d, folds, "nbinom2", metrics = "deviance"), "nuisance-dependent")
  expect_error(cross_validate_glm(positive ~ x, d, folds, "Gamma", weights = rep(.5, 96)), "integer")
  for (family in c("nbinom2", "Gamma")) {
    changed <- d; changed$y <- changed[[if (family == "nbinom2") "nb" else "positive"]]
    changed$y[1] <- if (family == "nbinom2") .5 else 0
    bad <- cross_validate_glm(y ~ x, changed, z$folds, family, weights = c(0, rep(1, 95)))
    expect_false(bad$complete); expect_true(all(is.na(bad$summary$estimate)))
    expect_equal(nrow(bad$predictions), 96)
    expect_true(any(grepl("requires", unlist(bad$errors))))
  }
  good <- cv_distribution_fixture()$values$nbinom2_FALSE_FALSE_FALSE$cv
  compact <- good; compact$fits <- NULL
  bad <- compact; bad$predictions$theta[1] <- bad$predictions$theta[1]*2
  expect_error(cv_tab(list(bad = bad), "neg_log_score"), "nuisance")
  bad <- compact; bad$selection_records[[1]]$data$native_K <- 99
  expect_error(select_cv(list(bad = bad), "AIC"), "count")
  bad <- compact; bad$distribution[[1]]$poisson_limit <- TRUE
  expect_error(cv_tab(list(bad = bad), "rmse"), "nuisance")
  bad <- compact; bad$selection_records[[1]]$data$nuisance$theta <- NULL
  expect_error(select_cv(list(bad = bad), "AIC"), "nuisance")
  broken <- cv_distribution_fixture()$values$Gamma_FALSE_FALSE_FALSE$cv$fits[[1]]
  broken$family$dispersion <- .4
  expect_error(.msr_cv_distribution_parameters(broken), "unverified fixed Gamma")
  m <- good$fits[[1]]; m$th.warn <- "kept native theta warning"
  expect_error(.msr_cv_distribution_parameters(m), "theta warnings")
  # Underflow-safe likelihoods retain an infinite loss rather than clipping it.
  p <- data.frame(observed = c(0, 1), prediction = c(0, 0), weight = c(0, 1), trials = 1,
    distribution_family = "nbinom2", theta = 3, dispersion = NA_real_, poisson_limit = FALSE)
  scored <- .msr_cv_distribution_score(p, "neg_log_score")
  expect_true(is.na(scored$neg_log_score[1])); expect_identical(scored$neg_log_score[2], Inf)
})

test_that("mixed-family selection retains actual fold distributions and comparable score measures", {
  z <- cv_distribution_fixture(); nb <- z$values$nbinom2_FALSE_FALSE_FALSE$cv
  case <- z$values$nbinom2_FALSE_FALSE_FALSE
  po <- cross_validate_glm(y ~ . - exposure + I(z^2) + offset(log(exposure)), case$data, nb$folds,
    "poisson", metrics = nb$settings$metrics, keep_fits = TRUE, keep_selection = TRUE)
  candidates <- list(poisson = po, nb = nb)
  full <- select_cv(candidates, "AIC", ties = "first", keep_fits = TRUE)
  compact <- lapply(candidates, function(x) { x$fits <- NULL; x })
  small <- select_cv(compact, "AIC", ties = "first")
  expect_true(full$complete); expect_true(small$complete)
  expect_identical(full$predictions, small$predictions)
  expect_identical(full$selection, small$selection)
  expect_equal(nrow(cv_tab(c(candidates, list(selected = small)), "neg_log_score")), 3)
  for (name in names(full$folds)) {
    winner <- candidates[[full$selection$chosen[[name]]]]
    expect_identical(full$settings$fold_settings[[name]], winner$settings)
    expect_true(all(full$predictions$distribution_family[full$predictions$fold == name] == winner$settings$family))
  }
  # A selected result remains a valid comparison candidate with fold-specific controls.
  expect_equal(nrow(cv_tab(list(selected = small, nb = nb), "rmse")), 2)
  ga <- z$values$Gamma_FALSE_FALSE_FALSE$cv
  expect_error(cv_tab(list(nb = nb, gamma = ga), "neg_log_score"), "identical splits|observations")
  same <- ga; same$targets$observed <- nb$targets$observed
  expect_error(cv_tab(list(nb = nb, gamma = same), "neg_log_score"), "disagree|positive")
})

test_that("continuous mean selection and density-measure guards apply across qualified families", {
  z <- cv_distribution_data(); d <- z$data; folds <- z$folds
  d$y <- d$positive
  candidates <- list(
    gamma = cross_validate_glm(y ~ x + g, d, folds, "gamma", metrics = c("rmse", "mae"), keep_selection = TRUE),
    gaussian = cross_validate_glm(y ~ x + g, d, folds, "gaussian", metrics = c("rmse", "mae"), keep_selection = TRUE))
  expect_true(all(vapply(candidates, `[[`, logical(1), "complete")))
  selected <- select_cv(candidates, "BIC", ties = "first")
  expect_true(selected$complete)
  expect_equal(nrow(cv_tab(c(candidates, list(selected = selected)), "rmse")), 3)
  for (name in names(folds)) expect_identical(selected$settings$fold_settings[[name]],
    candidates[[selected$selection$chosen[[name]]]]$settings)
  d$y <- d$nb + 1
  measures <- list(nb = cross_validate_glm(y ~ x, d, folds, "nbinom2", keep_selection = TRUE),
    gamma = cross_validate_glm(y ~ x, d, folds, "Gamma", keep_selection = TRUE))
  expect_true(all(vapply(measures, `[[`, logical(1), "complete")))
  expect_error(cv_tab(measures, "neg_log_score"), "continuous/discrete")
  expect_error(select_cv(measures, "AIC"), "continuous/discrete")
  expect_equal(nrow(cv_tab(measures, "mae")), 2)
})

test_that("known nuisance parameters preserve aliased regression rank and conditional SEs", {
  z <- cv_distribution_data(); d <- z$data; d$alias <- d$x
  for (family in c("nbinom2", "Gamma")) {
    d$y <- d[[if (family == "nbinom2") "nb" else "positive"]]
    controls <- if (family == "nbinom2") list(theta = 3) else list(dispersion = .4)
    cv <- do.call(cross_validate_glm, c(list(y ~ x + alias + g, d, z$folds[1L], family,
      keep_fits = TRUE, keep_selection = TRUE), controls))
    expect_true(cv$complete)
    m <- cv$fits[[1]]
    expect_true(is.na(coef(m)["alias"]))
    expect_equal(cv$selection_records[[1]]$data$K, m$rank)
    expect_false(cv$distribution[[1]]$estimated)
    expect_true(all(is.finite(vcov(m, complete = FALSE))))
  }
})

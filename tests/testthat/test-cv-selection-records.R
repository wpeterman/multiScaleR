# Deterministic native-GLM oracles; no scale fits or random-number generation.
cv_record_data <- function() {
  i <- 1:24
  d <- data.frame(x = i/24, group = factor(rep(c("a", "b"), 12)),
    exposure = 1 + (i %% 3)/10, y = rep(c(1, 3, 0, 4, 2, 1), 4),
    yes = rep(0:3, 6), no = 4 - rep(0:3, 6), binary = as.logical(i %% 2),
    continuous = 2 + i/30 + sin(i)/10, row.names = paste0("record", i))
  list(data = d, folds = grouped_cv_folds(rep(letters[1:6], each = 4), 3, rownames(d)))
}

test_that("compact selection matches native tables, choices and held-out scores", {
  z <- cv_record_data(); d <- z$data
  models <- list(null = cross_validate_glm(y ~ 1, d, z$folds, "poisson", keep_fits = TRUE, keep_selection = TRUE),
    trend = cross_validate_glm(y ~ poly(x, 2) * group + offset(log(exposure)), d, z$folds,
      "poisson", keep_fits = TRUE, keep_selection = TRUE))
  compact <- lapply(models, function(x) { x$fits <- NULL; x })
  mixed <- compact; mixed$null <- models$null
  for (criterion in c("AIC", "AICc", "BIC")) {
    full <- select_cv(models, criterion, ties = "first")
    for (inputs in list(compact, mixed)) {
      selected <- select_cv(inputs, criterion, ties = "first")
      for (item in c("predictions", "summary", "diagnostics", "warnings", "errors", "training", "prediction_contracts"))
        expect_identical(selected[[item]], full[[item]])
      expect_identical(selected$selection, full$selection)
    }
  }
  expect_lt(as.numeric(object.size(compact)), as.numeric(object.size(models)))
  plain <- cross_validate_glm(y ~ x, d, z$folds, "poisson", keep_selection = TRUE)
  expect_null(plain$fits); expect_true(select_cv(list(plain = plain), "AIC")$complete)
  expect_error(select_cv(compact, "AIC", keep_fits = TRUE), "cannot reconstruct")
  legacy <- lapply(models, function(x) { x$selection_records <- NULL; x })
  expect_identical(select_cv(legacy, "AIC")$selection$chosen, select_cv(compact, "AIC")$selection$chosen)
  # Records have no model/design/formula environments or external pointers.
  small <- function(x) {
    if (is.environment(x) || typeof(x) == "externalptr" || inherits(x, c("glm", "multiScaleR", "formula"))) return(FALSE)
    if (is.list(x)) return(all(vapply(x, small, logical(1))))
    TRUE
  }
  expect_true(small(plain$selection_records))
  restored <- unserialize(serialize(plain, NULL))
  a <- select_cv(list(plain = plain), "AIC"); b <- select_cv(list(plain = restored), "AIC")
  expect_identical(a[setdiff(names(a), "call")], b[setdiff(names(b), "call")])
  expect_error(cross_validate_glm(y ~ x, d, z$folds, "poisson", keep_selection = NA), "keep_selection")
  defaults <- cross_validate_glm(y ~ x, d, z$folds, "poisson")
  expect_null(defaults$selection_records)
  expect_identical(defaults$predictions, plain$predictions)
})

test_that("native response representations, zero rows and rank counts survive compact capture", {
  z <- cv_record_data(); d <- z$data; d$alias <- d$x
  d$yes[c(1, 5)] <- d$no[c(1, 5)] <- 0
  d$factor_y <- factor(ifelse(d$binary, "yes", "no"), levels = c("no", "yes"))
  w <- rep(c(0, 1, 2), 8)
  cases <- list(gaussian = list(continuous ~ x + alias + group, "gaussian", NULL),
    poisson = list(y ~ x + offset(log(exposure)), "poisson", w),
    counts = list(cbind(yes, no) ~ x + group, "binomial", w),
    logical = list(I(binary) ~ x, "binomial", w),
    factor = list(factor_y ~ x, "binomial", w))
  for (spec in cases) {
    cv <- cross_validate_glm(spec[[1]], d, z$folds, spec[[2]], spec[[3]], keep_fits = TRUE, keep_selection = TRUE)
    expect_true(cv$complete)
    compact <- cv; compact$fits <- NULL
    for (criterion in c("AIC", "AICc", "BIC")) {
      a <- select_cv(list(candidate = cv), criterion)
      b <- select_cv(list(candidate = compact), criterion)
      expect_identical(a$selection$tables, b$selection$tables)
      expect_identical(a$summary, b$summary)
    }
    for (name in names(z$folds)) {
      record <- cv$selection_records[[name]]$data; fit <- cv$fits[[name]]
      expect_equal(record$K, attr(logLik(fit), "df"))
      expect_equal(record$logLik, as.numeric(logLik(fit)))
      expect_identical(record$nobs, as.integer(nobs(fit)))
      expect_identical(record$formula, paste(deparse(formula(fit)), collapse = " "))
    }
  }
})

test_that("compact selection preserves failures and rejects stale or substituted records", {
  z <- cv_record_data(); d <- z$data
  cv <- cross_validate_glm(y ~ x, d, z$folds, "poisson", keep_fits = TRUE, keep_selection = TRUE)
  compact <- cv; compact$fits <- NULL
  tie <- select_cv(list(first = compact, second = compact), "AIC")
  expect_false(tie$complete); expect_match(tie$errors[[1]]$message, "tie")
  expect_true(all(is.na(tie$summary$estimate)))
  expect_true(all(select_cv(list(first = compact, second = compact), "AIC", "first")$selection$chosen == "first"))
  small <- cross_validate_glm(y ~ poly(x, 2), d[1:8, ], grouped_cv_folds(rep(letters[1:4], each = 2), 2,
    rownames(d)[1:8]), "poisson", keep_selection = TRUE)
  result <- select_cv(list(small = small), "AICc")
  expect_false(result$complete); expect_match(result$errors[[1]]$message, "n > K")
  changed <- compact; changed$predictions$prediction[] <- 100
  expect_identical(select_cv(list(changed = changed), "AIC")$selection$chosen,
    select_cv(list(changed = compact), "AIC")$selection$chosen)
  bad <- compact; bad$selection_records[[1]]$data$targets$observed[1] <- 999
  expect_error(select_cv(list(bad = bad), "AIC"), "training responses")
  bad <- compact; bad$selection_records[[1]]$data$signature$weights[1] <- 999
  expect_error(select_cv(list(bad = bad), "AIC"), "signature")
  bad <- compact; bad$selection_records[[1]]$data$family <- "gaussian"
  expect_error(select_cv(list(bad = bad), "AIC"), "family and link")
  bad <- compact; bad$selection_records[[1]]$data$observation_ids <- rev(bad$folds[[1]]$analysis)
  expect_error(select_cv(list(bad = bad), "AIC"), "likelihood row IDs")
  bad <- compact; bad$selection_records[[1]]$data$K <- 999
  expect_error(select_cv(list(bad = bad), "AIC"), "parameter counts")
  bad <- compact; bad$settings$glm_control$maxit <- 999
  expect_error(select_cv(list(bad = bad), "AIC"), "controls")
  bad <- compact; bad$diagnostics[[1]]$fitted_glm$converged <- FALSE
  expect_error(select_cv(list(bad = bad), "AIC"), "diagnostics")
  bad <- cv; bad$selection_records[[1]]$data$logLik <- -999
  expect_error(select_cv(list(bad = bad), "AIC"), "disagree with their retained")
  bad <- cv; bad$selection_records <- bad$selection_records[-1L]
  expect_error(select_cv(list(bad = bad), "AIC"), "every training selection record")
  bad <- cv; bad$fits <- bad$fits[-1L]
  expect_error(select_cv(list(bad = bad), "AIC"), "every training fit")
  broken <- cv$fits[[1]]; broken$aic <- Inf
  metadata <- .msr_cv_selection_record(broken, cv$diagnostics[[1]], .msr_cv_selection_controls(cv$settings))
  expect_null(metadata$data)
  expect_match(metadata$error, "finite scalar log-likelihood")
  # Metadata extraction errors do not erase otherwise scored predictions.
  bad <- compact
  bad$selection_records[[1]]$data <- NULL
  bad$selection_records[[1]]$error <- "Native likelihood metadata unavailable"
  bad$warnings[[1]] <- "original candidate warning"
  bad$selection_records[[1]]$warnings <- "native metadata warning"
  selected <- select_cv(list(bad = bad), "AIC")
  expect_false(selected$complete); expect_match(selected$errors[[1]]$message, "metadata unavailable")
  expect_identical(selected$selection$candidate_warnings$bad, bad$warnings)
  expect_identical(selected$warnings[[1]], "native metadata warning")
  expect_true(compact$complete); expect_true(all(is.finite(compact$predictions$prediction)))
})

# Small deterministic native-GLM checks; no scale optimization or simulation.
cv_reference_data <- function() {
  n <- 24L; i <- seq_len(n)
  d <- data.frame(x = seq(-1, 1, length.out = n), g = factor(rep(c("a", "b"), n / 2)),
    effort = 1 + (i %% 3) / 10, count = rep(0:5, 4), yes = rep(0:3, 6), no = 4 - rep(0:3, 6),
    binary = as.logical(i %% 2))
  d$continuous <- 2 + .6 * d$x + .1 * d$x^2 + sin(i * 1.3) / 10
  rownames(d) <- paste0("r", i)
  list(data = d, folds = grouped_cv_folds(rep(letters[1:6], each = 4), 3, rownames(d)))
}

test_that("ordinary Gaussian equations and training transformations match native GLMs", {
  z <- cv_reference_data(); d <- z$data
  formula <- continuous ~ poly(x, 2) * g + offset(log(effort))
  cv <- cross_validate_glm(formula, d, z$folds, keep_fits = TRUE)
  expect_true(cv$complete)
  for (name in names(z$folds)) {
    split <- z$folds[[name]]
    fit <- glm(formula, data = d[split$analysis, ], family = gaussian())
    expect_equal(coef(cv$fits[[name]]), coef(fit))
    p <- predict(fit, d[split$assessment, ], type = "response")
    expect_equal(cv$predictions$prediction[cv$predictions$fold == name], unname(p))
    expect_identical(cv$prediction_contracts[[name]]$uncertainty$intervals, "not_provided")
    expect_false(cv$prediction_contracts[[name]]$uncertainty$observation_noise)
  }
  p <- cv$predictions
  expect_equal(cv_tab(list(reference = cv), "rmse")$Score, sqrt(mean((p$observed - p$prediction)^2)))
  expect_true(select_cv(list(reference = cv), "AIC")$complete)
  expect_error(cross_validate_glm(formula, d, z$folds, weights = rep(2, nrow(d))), "Gaussian weights")
  expect_error(cross_validate_glm(formula, d, z$folds, metrics = "brier"), "supported")
  expect_error(cross_validate_glm(formula, d, z$folds, glm_control = list(invented = 1)), "glm_control")
  # No assessment rows enter learned poly metadata or the native training fit.
  changed <- d; changed[z$folds[[1]]$assessment, "continuous"] <- 999
  changed[z$folds[[1]]$assessment, "x"] <- 50
  other <- cross_validate_glm(formula, changed, z$folds[1], keep_fits = TRUE)
  expect_identical(coef(other$fits[[1]]), coef(cv$fits[[1]]))
  expect_identical(other$fits[[1]]$model, cv$fits[[1]]$model)
  expect_identical(other$fits[[1]]$terms, cv$fits[[1]]$terms)
})

test_that("ordinary GLM scoring retains frequency, zero-trial and failure cases", {
  z <- cv_reference_data(); d <- z$data; weights <- rep(c(0L, 1L, 2L), length.out = nrow(d))
  d$yes[c(1, 5)] <- d$no[c(1, 5)] <- 0
  cv <- cross_validate_glm(cbind(yes, no) ~ x + g, d, z$folds, "binomial", weights, keep_fits = TRUE)
  expect_true(cv$complete)
  p <- cv$predictions; keep <- p$weight > 0 & p$trials > 0
  expect_identical(p$eligible, keep)
  expect_equal(cv_tab(list(current = cv), "neg_log_score")$Score,
    weighted.mean(-dbinom(p$observed[keep], p$trials[keep], p$prediction[keep], log = TRUE), p$weight[keep]))
  expect_true(select_cv(list(current = cv), "AIC")$complete)
  legacy <- cv; legacy$targets <- NULL
  expect_equal(cv_tab(list(current = cv, legacy = legacy), "brier")$Score, rep(cv$summary$estimate[2], 2))
  expect_true(select_cv(list(legacy = legacy), "AIC")$complete)
  # Equivalent factor and logical responses have the same recorded binary target.
  binary <- cross_validate_glm(binary ~ x, d, z$folds, "binomial")
  d$factor_response <- factor(ifelse(d$binary, "yes", "no"), levels = c("no", "yes"))
  factor_fit <- cross_validate_glm(factor_response ~ x, d, z$folds, "binomial")
  expect_equal(cv_tab(list(logical = binary, factor = factor_fit), "brier")$Score,
    rep(binary$summary$estimate[2], 2))
  # A native prediction error is retained, with no retries or dropped fold.
  d$g <- factor(ifelse(rownames(d) %in% z$folds[[1]]$assessment, "new", as.character(d$g)))
  failed <- cross_validate_glm(count ~ g + x, d, z$folds, "poisson")
  expect_false(failed$complete)
  expect_identical(failed$fold_results$status[1], "failed")
  expect_match(failed$errors[[1]]$message, "new level")
  expect_true(all(is.na(failed$summary$estimate)))
  expect_error(cv_tab(list(failed = failed), "rmse"), "complete")
})

test_that("training selection retains criterion failures and validates its actual family", {
  d <- data.frame(y = c(1, 2, 1, 3, 2, 4, 1, 2), x = 1:8)
  f <- grouped_cv_folds(rep(letters[1:4], each = 2), 2)
  full <- cross_validate_glm(y ~ poly(x, 2), d, f, "poisson", keep_fits = TRUE)
  expect_true(full$complete)
  selected <- select_cv(list(full = full), "AICc")
  expect_false(selected$complete)
  expect_true(all(is.na(selected$summary$estimate)))
  expect_match(selected$errors[[1]]$message, "n > K")
  expect_true(select_cv(list(full = full), "AIC")$complete)
  altered <- full; altered$fits[[1]]$family <- gaussian()
  expect_error(select_cv(list(altered = altered), "AIC"), "family and link")
  compact <- cross_validate_glm(y ~ 1, d, f, "poisson")
  compact$predictions$weight[1] <- -1
  expect_error(cv_tab(list(altered = compact), "rmse"), "frequency weights")
})

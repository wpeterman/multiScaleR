# Deterministic choice-set data, not a recovery or simulation panel.
choice_cv_data <- function() {
  i <- 1:72; set <- (i - 1L) %/% 3L
  d <- data.frame(case_ = as.numeric((i - 1L) %% 3L == set %% 3L), choice = set + 1L,
    x = sin(i * 1.7), z = cos(i * .8), g = factor(rep(c("a", "b"), 36)),
    exposure = 1 + (i %% 4)/10, row.names = paste0("alt", i))
  list(data = d, folds = grouped_cv_folds(d$choice, 3, rownames(d)))
}

test_that("native choice references, null likelihoods and event penalties agree", {
  skip_if_not_installed("survival")
  z <- choice_cv_data(); d <- z$data; coxph <- survival::coxph; strata <- survival::strata
  for (method in c("exact", "efron", "breslow")) {
    forms <- list(null = case_ ~ strata(choice), offset = case_ ~ offset(log(exposure)) + strata(choice),
      trend = case_ ~ poly(x, 2) * g + offset(log(exposure)) + strata(choice))
    models <- lapply(forms, function(f) cross_validate_clogit(f, d, z$folds, method = method,
      keep_fits = TRUE, keep_selection = TRUE))
    expect_true(all(vapply(models, function(x) x$complete, logical(1))),
      info = paste(unlist(lapply(models, `[[`, "errors")), collapse = "; "))
    for (nm in names(models)) {
      cv <- models[[nm]]
      expect_s3_class(cv, "multiScaleR_choice_cv")
      for (fold in names(z$folds)) {
        p <- cv$predictions[cv$predictions$fold == fold, ]
        assessment <- d[p$row_id, ]; fit <- cv$fits[[fold]]
        if (length(coef(fit))) {
          # LP reference='zero' has no fitted stratum effects. The provider still
          # validates its stratum levels, so use an existing label solely for
          # this independent native LP oracle; normalize over the actual sets.
          assessment$choice <- d[z$folds[[fold]]$analysis[1L], "choice"]
          eta <- as.numeric(predict(fit, assessment, type = "lp", reference = "zero"))
        } else eta <- if (nm == "offset") log(assessment$exposure) else rep(0, nrow(assessment))
        expected <- numeric(nrow(p))
        for (g in split(seq_len(nrow(p)), p$choice)) {
          m <- max(eta[g]); expected[g] <- eta[g] - m - log(sum(exp(eta[g] - m)))
        }
        expect_equal(p$log_probability, expected, tolerance = 1e-10)
        expect_equal(as.numeric(rowsum(p$prediction, p$choice)), rep(1, 8), tolerance = 1e-12)
        assessment <- d[p$row_id, ]; assessment$fixed_lp <- eta
        # No parameters are fitted on held-out data. The public survival
        # fitter supplies an independent complete-set likelihood at fixed LP.
        oracle <- survival::clogit(case_ ~ offset(fixed_lp) + strata(choice), assessment,
          method = method, model = TRUE, x = TRUE, y = TRUE)
        expect_equal(sum(p$log_probability[p$observed == 1]), unname(oracle$loglik[1L]), tolerance = 1e-10)
        record <- cv$selection_records[[fold]]$data
        expect_identical(record$nobs, 16L)
        expect_equal(record$K, length(coef(fit)))
        expect_equal(record$logLik, as.numeric(logLik(fit)))
      }
      expect_identical(cv$summary$unit, "choice_set")
      expect_equal(cv$summary$n_assessment, 24L)
      expect_equal(cv$summary$n_alternatives, 72L)
      expect_equal(cv$summary$estimate, mean(cv$predictions$neg_log_score[cv$predictions$eligible]))
    }
    compact <- lapply(models, function(x) { x$fits <- NULL; x })
    expect_equal(cv_tab(compact, "neg_log_score")$Unit, rep("choice_set", 3))
    for (criterion in c("AIC", "AICc", "BIC")) {
      a <- select_cv(models, criterion, ties = "first")
      b <- select_cv(compact, criterion, ties = "first")
      expect_true(a$complete); expect_true(b$complete, info = paste(unlist(b$errors), collapse = "; "))
      for (field in c("predictions", "summary", "diagnostics", "warnings", "prediction_contracts", "selection"))
        expect_identical(a[[field]], b[[field]])
      tab <- a$selection$tables[[1L]]
      k <- tab$K[tab$Modnames == "trend"]; ll <- tab$LL[tab$Modnames == "trend"]
      expected <- switch(criterion, AIC = -2*ll + 2*k, AICc = -2*ll + 2*k + 2*k*(k+1)/(16-k-1), BIC = -2*ll + log(16)*k)
      expect_equal(tab[[criterion]][tab$Modnames == "trend"], expected)
    }
  }
})

test_that("choice scoring remains stable and repeated assessment sets stay separate", {
  p <- data.frame(row_id = letters[1:4], fold = c("a", "a", "b", "b"), choice = "same",
    observed = c(0, 1, 0, 1), linear_predictor = c(1000, -1000, 1001, 1000))
  p <- .msr_choice_score(p, "neg_log_score")
  expect_equal(p$neg_log_score[p$eligible], c(2000, log1p(exp(1))), tolerance = 1e-12)
  expect_identical(p$prediction[2L], 0)
  expect_true(all(is.finite(p$log_probability)))
  expect_equal(.msr_choice_metric(p, TRUE)$n_assessment, 2)
  z <- choice_cv_data(); folds <- c(z$folds, list(repeat_a = z$folds[[1L]]))
  cv <- cross_validate_clogit(case_ ~ x + strata(choice), z$data, folds, keep_selection = TRUE)
  expect_true(cv$complete)
  expect_equal(cv$summary$n_assessment, 32)
  expect_true(is.finite(cv_tab(list(repeated = cv), "neg_log_score")$Score))
})

test_that("choice CV retains native failures and rejects split or target substitutions", {
  z <- choice_cv_data(); d <- z$data
  partial <- z$folds; partial[[1L]]$assessment <- partial[[1L]]$assessment[-1L]
  expect_error(cross_validate_clogit(case_ ~ x + strata(choice), d, partial), "complete")
  broken <- d; broken$case_[1:3] <- 0
  expect_error(cross_validate_clogit(case_ ~ x + strata(choice), broken, z$folds), "exactly one")
  expect_error(cross_validate_clogit(case_ ~ x + survival::strata(choice), d, z$folds), "bare strata")
  expect_error(cross_validate_clogit(case_ ~ x + strata(choice) + cluster(g), d, z$folds), "special")
  expect_error(cross_validate_clogit(case_ ~ x + strata(choice), d, z$folds, weights = rep(1, 72)), "unused argument")
  cv <- cross_validate_clogit(case_ ~ x + strata(choice), d, z$folds, keep_fits = TRUE, keep_selection = TRUE)
  expect_true(cv$complete)
  bad <- cv; bad$targets$choice[1] <- "changed"
  expect_error(cv_tab(list(bad = bad), "neg_log_score"), "exactly one|two alternatives|complete|membership")
  bad <- cv; bad$predictions$prediction[1] <- .5
  expect_error(cv_tab(list(bad = bad), "neg_log_score"), "disagree")
  bad <- cv; bad$selection_records[[1L]]$data$nobs <- 48L
  expect_error(select_cv(list(bad = bad)), "event count")
  bad <- cv; bad$selection_records[[1L]]$data$signature$strata <- rep(1L, 48)
  expect_error(select_cv(list(bad = bad)), "signature")
  expect_error(cv_tab(list(cv = cv), "brier"), "neg_log_score")
  glm <- cross_validate_glm(case_ ~ x, d, z$folds, "binomial")
  expect_error(cv_tab(list(choice = cv, independent = glm), "neg_log_score"), "conditional-choice")
  tie <- select_cv(list(first = cv, second = cv), "AIC")
  expect_false(tie$complete); expect_true(all(is.na(tie$predictions$linear_predictor)))
  expect_true(all(is.na(tie$summary$estimate)))
  failed_data <- d; failed_data$g <- factor(ifelse(rownames(d) %in% z$folds[[1L]]$assessment, "new", as.character(d$g)))
  failed <- cross_validate_clogit(case_ ~ x + g + strata(choice), failed_data, z$folds)
  # survival retains an absent factor level as an aliased coefficient, so
  # its native training guard fires before the prediction new-level check.
  expect_false(failed$complete); expect_match(failed$errors[[1L]]$message, "nonfinite coefficients")
  expect_identical(failed$errors[[1L]]$stage, "training_clogit")
  expect_equal(nrow(failed$predictions), nrow(d)); expect_true(all(is.na(failed$summary$estimate)))
  expect_error(cv_tab(list(failed = failed), "neg_log_score"), "complete")
  # A null adapter must not make inconsistent native metadata selectable.
  null <- cross_validate_clogit(case_ ~ strata(choice), d, z$folds, keep_fits = TRUE)
  null$fits[[1L]]$loglik <- 999
  expect_error(logLik(null$fits[[1L]]), "disagrees")
})

test_that("namespaced fixed training transforms retain native metadata", {
  z <- choice_cv_data(); d <- z$data
  form <- case_ ~ stats::poly(x, 2) + base::scale(z) + g + strata(choice)
  cv <- cross_validate_clogit(form, d, z$folds, keep_fits = TRUE)
  expect_true(cv$complete)
  for (name in names(z$folds)) {
    p <- cv$predictions[cv$predictions$fold == name, ]; nd <- d[p$row_id, ]
    nd$choice <- d[z$folds[[name]]$analysis[1L], "choice"]
    native <- .msr_native_fixed_metadata(cv$fits[[name]])
    eta <- as.numeric(predict(native, nd, type = "lp", reference = "zero"))
    expected <- numeric(nrow(p))
    for (g in split(seq_len(nrow(p)), p$choice)) {
      m <- max(eta[g]); expected[g] <- eta[g]-m-log(sum(exp(eta[g]-m)))
    }
    expect_equal(p$log_probability, expected, tolerance = 1e-10)
  }
})

test_that("native dot formulas and unrepresentably large losses keep their targets", {
  z <- choice_cv_data(); d <- z$data[c("case_", "choice", "x", "z")]
  a <- cross_validate_clogit(case_ ~ . - choice + strata(choice), d, z$folds, keep_fits = TRUE)
  b <- cross_validate_clogit(case_ ~ x + z + strata(choice), d, z$folds, keep_fits = TRUE)
  expect_true(a$complete); expect_true(b$complete)
  expect_identical(a$predictions, b$predictions)
  for (name in names(z$folds)) expect_identical(coef(a$fits[[name]]), coef(b$fits[[name]]))
  a$predictions$linear_predictor <- ifelse(a$predictions$observed == 1, -1e308, 1e308)
  a$predictions <- .msr_choice_score(a$predictions, "neg_log_score")
  tab <- cv_tab(list(first = a, second = a), "neg_log_score")
  expect_true(all(is.infinite(tab$Score))); expect_true(all(is.na(tab$Delta)))
  expect_equal(tab$Rank, c(1, 1))
})

test_that("native choice admission does not depend on retained representations", {
  z <- choice_cv_data(); d <- z$data
  strata <- survival::strata; ridge <- survival::ridge; coxph <- survival::coxph
  for (method in c("efron", "breslow")) {
    native <- survival::clogit(case_ ~ ridge(x, theta = 1) + strata(choice),
      data = d[z$folds[[1L]]$analysis, ], method = method, model = TRUE, x = TRUE, y = TRUE)
    expect_true(inherits(native, "coxph") && (!is.null(native$penalty) ||
      any(vapply(native$model, inherits, logical(1), "coxph.penalty"))))
    expect_true(all(is.finite(coef(native))) && is.finite(as.numeric(logLik(native))))
    cv <- cross_validate_clogit(case_ ~ ridge(x, theta = 1) + strata(choice), d,
      z$folds[1L], method = method, keep_fits = FALSE, keep_selection = FALSE)
    expect_false(cv$complete)
    expect_identical(cv$errors[[1L]]$stage, "training_identity")
    expect_match(cv$errors[[1L]]$message, "unpenalized")
    expect_null(cv$fits); expect_null(cv$selection_records)
    expect_true(all(is.na(cv$predictions$prediction)))
    expect_true(all(is.na(cv$summary$estimate)))
  }
})

test_that("prepared fixed formulas freeze qualified scaling and shield fixed logic", {
  dat <- data.frame(y = seq_len(12), z = c(-2, -1, -.5, 0, .2, .7, 1, 1.3, 2, 3, 4, 5),
    a = rep(c(TRUE, FALSE, FALSE), 4), b = rep(c(FALSE, TRUE), 6),
    site = factor(rep(seq_len(4), each = 3)))
  bare <- .msr_bayes_formula(y ~ scale(z) + I(a | b) + (1 | site), dat,
    response = TRUE, active_frame = TRUE)
  qualified <- .msr_bayes_formula(y ~ base::scale(z) + I(a | b) + (1 | site), dat,
    response = TRUE, active_frame = TRUE)
  expect_equal(model.matrix(bare$terms, bare$frame),
    model.matrix(qualified$terms, qualified$frame), ignore_attr = TRUE)
  expect_length(qualified$random, 1L)
  expect_identical(qualified$random[[1]]$group, "site")
  nd <- dat[c(9, 2, 5, 9), , drop = FALSE]
  expected <- (nd$z - mean(dat$z)) / sd(dat$z)
  tt <- delete.response(qualified$terms)
  frame <- model.frame(tt, nd, na.action = na.fail)
  X <- model.matrix(tt, frame)
  expect_equal(as.numeric(X[, "base::scale(z)"]), expected, tolerance = 1e-12)
  expect_identical(frame[["I(a | b)"]], I(nd$a | nd$b))
  expect_identical(.msr_fixed_prediction_terms(qualified$terms, qualified$frame), qualified$terms)
  scale <- .msr_bayes_formula(~ base::scale(z) + I(a | b), dat, active_frame = TRUE)
  expect_equal(as.numeric(model.matrix(scale$terms, model.frame(scale$terms, nd))[, "base::scale(z)"]),
    expected, tolerance = 1e-12)
  expect_error(.msr_bayes_formula(y ~ z + (z | site), dat, response = TRUE), "additive random intercepts")
  expect_error(.msr_bayes_formula(y ~ z + (1 || site), dat, response = TRUE), "additive random intercepts")
})

test_that("profile fixed metadata retains native family targets and new-row bases", {
  for (family in c("poisson", "gaussian", "binomial", "nbinom2")) {
    f <- formula_integration_fixture(family)
    make <- function(qualified) {
      form <- update.formula(formula(f$model), if (qualified)
        . ~ . - z + base::scale(z) else . ~ . - z + scale(z))
      model <- glm(form, f$data, family = f$model$family)
      contract <- .msr_profile_contract(model, f$inputs)
      kp <- log(.3)
      solved <- contract$objective$evaluate(c(kp, if (family == "nbinom2") log(2)))
      native <- .msr_profile_conditional_model(contract, solved, model, kp)
      object <- structure(list(backend = "profile_glm", opt_mod = native,
        kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = kp),
        uncertainty = list(par = c(kp, solved$beta), beta_index = 1L + seq_along(solved$beta))),
        class = "multiScaleR")
      list(object = object, solved = solved)
    }
    b <- make(FALSE); q <- make(TRUE)
    expect_equal(b$solved$value, q$solved$value, tolerance = 1e-10)
    expect_equal(unname(b$solved$beta), unname(q$solved$beta), tolerance = 1e-10)
    for (rows in list(seq_len(20), rev(seq_len(30)), c(4, 1, 8))) {
      nd <- f$data[rows, , drop = FALSE]
      bd <- .msr_profile_prediction_design(b$object, nd, f$inputs)
      qd <- .msr_profile_prediction_design(q$object, nd, f$inputs)
      expect_equal(qd$X, bd$X, tolerance = 1e-12, ignore_attr = TRUE)
      expect_equal(drop(qd$X %*% q$solved$beta) + qd$offset,
        drop(bd$X %*% b$solved$beta) + bd$offset, tolerance = 1e-10, ignore_attr = TRUE)
    }
    expect_equal(predict(q$object$opt_mod, type = "link", se.fit = TRUE)$fit,
      predict(b$object$opt_mod, type = "link", se.fit = TRUE)$fit, tolerance = 1e-10)
  }
})


test_that("conditional fixed metadata retains original native expression identities", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  for (contrast in c("contr.treatment", "contr.sum")) {
    for (basis in c("scale", "poly")) {
      f <- profile_glm_fixture(); dat <- f$data
      dat$case_ <- rep(c(1, 0, 0), length.out = nrow(dat))
      dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
      contrasts(dat$group) <- get(contrast)(nlevels(dat$group))
      make <- function(qualified) {
        expr <- if (basis == "scale") if (qualified) "base::scale(z)" else "scale(z)" else
          if (qualified) "stats::poly(z,2)" else "poly(z,2)"
        form <- as.formula(paste("case_ ~ hab * group +", expr,
          "+ offset(log(exposure)) + strata(choice)"), env = asNamespace("survival"))
        model <- survival::clogit(form, data = dat, model = TRUE, x = TRUE, y = TRUE)
        contract <- .msr_clogit_contract(model, f$inputs)
        contract$objective <- list(engine = contract$engine, X_fixed = contract$X_fixed,
          offset = contract$offset, response_design = contract$response_design)
        object <- structure(list(backend = "conditional_logit", opt_mod = model,
          kernel_inputs = f$inputs, profile_fit = list(contract = contract, par = log(.3))),
          class = "multiScaleR")
        list(model = model, object = object, contract = contract)
      }
      b <- make(FALSE); q <- make(TRUE)
      expect_equal(as.numeric(logLik(q$model)), as.numeric(logLik(b$model)), tolerance = 1e-10)
      expect_equal(unname(coef(q$model)), unname(coef(b$model)), tolerance = 1e-10)
      ob <- .msr_clogit_objective(b$contract); oq <- .msr_clogit_objective(q$contract)
      columns <- if (is.null(b$contract$response_design))
        c(b$contract$fixed_columns, b$contract$covariates) else b$contract$response_design$columns
      beta <- coef(b$model)[columns]
      par <- c(log(.3), beta)
      expect_equal(oq$fn(par), ob$fn(par), tolerance = 1e-10)
      expect_equal(oq$gr(par), ob$gr(par), tolerance = 1e-10, ignore_attr = TRUE)
      for (rows in list(seq_len(12), rev(seq_len(12)), c(4, 1, 8))) {
        nd <- dat[rows, , drop = FALSE]
        # New rows carry factor levels. The fitted contrast is supplied by
        # the prediction contract, not an embedded training-data attribute.
        attr(nd$group, "contrasts") <- NULL
        bd <- .msr_profile_prediction_design(b$object, nd, f$inputs)
        qd <- .msr_profile_prediction_design(q$object, nd, f$inputs)
        expect_equal(qd$X, bd$X, tolerance = 1e-12, ignore_attr = TRUE)
        expect_equal(qd$offset, bd$offset, tolerance = 1e-12)
        if (basis == "scale")
          expect_equal(as.numeric(qd$X[, "base::scale(z)"]),
            (nd$z - mean(dat$z)) / sd(dat$z), tolerance = 1e-12)
      }
    }
  }
})

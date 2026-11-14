# Fixed-knot splines with backend = "conditional_logit", fitted the way a user
# does it: a plain clogit(..., model = TRUE) call with library(survival) and
# library(splines) bindings. survival::coxph() passes its model frame through
# na.omit(), which strips the ns()/bs() class and knot attributes, so the
# basis must be checked against the saved prediction call instead.

# ssf_fixture() and ssf_native_loglik() are in helper-step-selection.R.

test_that("fixed-knot ns() fits through conditional_logit from a plain clogit() call", {
  skip_on_cran()
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  # The bindings a user gets from library(survival) and library(splines).
  # clogit() evaluates its coxph() call in the caller's frame.
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  ns <- splines::ns
  f_ns <- case_ ~ ns(forest, knots = 0, Boundary.knots = c(-3, 3)) + sl_km + strata(step_id)
  m_ns <- survival::clogit(f_ns, data = fx$df, model = TRUE)

  fit <- multiScale_optim(m_ns, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE)
  sig <- fit$scale_est["forest", "Mean"]
  expect_true(is.finite(sig) && sig > 0)
  expect_equal(fit$optim_results$value, -as.numeric(logLik(fit$opt_mod)),
               tolerance = 1e-5)
  expect_length(grep("^ns\\(forest", names(coef(fit$opt_mod))), 2L)

  # Independent check: refit clogit() on kernel values prepared at a given
  # sigma. The joint fit must match the native likelihood at its estimate
  # (up to distance binning) and beat nearby scales.
  native_ll <- function(s) ssf_native_loglik(fx, f_ns, s)
  ll_hat <- native_ll(sig)
  expect_equal(ll_hat, as.numeric(logLik(fit$opt_mod)), tolerance = 1e-4)
  expect_gt(ll_hat, native_ll(0.75 * sig))
  expect_gt(ll_hat, native_ll(1.33 * sig))

  p <- predict(fit)
  expect_true(all(is.finite(p)))
})

test_that("fixed-knot bs() fits through conditional_logit from a plain clogit() call", {
  skip_on_cran()
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  bs <- splines::bs
  m_bs <- survival::clogit(case_ ~ bs(forest, knots = 0, Boundary.knots = c(-4, 4)) +
                             sl_km + strata(step_id), data = fx$df, model = TRUE)
  fit <- multiScale_optim(m_bs, fx$inputs, backend = "conditional_logit",
                          n_cores = 1, verbose = FALSE)
  expect_true(is.finite(fit$scale_est["forest", "Mean"]))
  expect_equal(fit$optim_results$value, -as.numeric(logLik(fit$opt_mod)),
               tolerance = 1e-5)
})

test_that("a saved spline call that differs from the fixed recipe is still rejected", {
  skip_on_cran()
  skip_if_not_installed("survival")
  fx <- ssf_fixture()
  coxph <- survival::coxph; Surv <- survival::Surv; strata <- survival::strata
  ns <- splines::ns
  m_ns <- survival::clogit(case_ ~ ns(forest, knots = 0, Boundary.knots = c(-3, 3)) +
                             sl_km + strata(step_id), data = fx$df, model = TRUE)
  # Alter the knot stored in the prediction call; the recipe comes from the
  # formula, so the two no longer agree and the fallback must refuse.
  tt <- m_ns$terms
  pv <- attr(tt, "predvars")
  pv[[3L]]$knots <- 0.5
  attr(tt, "predvars") <- pv
  m_ns$terms <- tt
  attr(m_ns$model, "terms") <- tt
  # Some survival versions keep the basis attributes in the model frame, and a
  # frame that carries them is checked directly. Strip them so the check must
  # use the saved prediction call, whatever the installed survival version.
  nm <- grep("^ns[(]", names(m_ns$model), value = TRUE)
  basis <- m_ns$model[[nm]]
  m_ns$model[[nm]] <- matrix(as.numeric(basis), nrow = nrow(basis),
                             dimnames = list(NULL, colnames(basis)))
  expect_error(
    multiScale_optim(m_ns, fx$inputs, backend = "conditional_logit",
                     n_cores = 1, verbose = FALSE),
    "does not match its fixed basis recipe"
  )
})

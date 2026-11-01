# Deterministic choice sets expose alternative-row versus likelihood sample size.
cox_selection_fixture <- function() {
  coxph <- survival::coxph
  i <- seq_len(160L); group <- rep(seq_len(40L), each = 4L)
  d <- data.frame(choice = group, x = sin(i * 1.7) + group / 20,
    z = cos(i * .61), case_ = 0L)
  d$case_[(seq_len(40L) - 1L) * 4L + (seq_len(40L) - 1L) %% 4L + 1L] <- 1L
  formula <- as.formula("case_ ~ x + z + strata(choice)", env = asNamespace("survival"))
  list(data = d, model = survival::clogit(formula, d, model = TRUE, x = TRUE, y = TRUE),
    formula = formula)
}

test_that("Cox information criteria use native event count without losing alternative rows", {
  skip_if_not_installed("survival")
  f <- cox_selection_fixture(); m <- f$model
  before <- serialize(m, NULL)
  ll <- as.numeric(logLik(m)); k <- as.numeric(attr(logLik(m), "df"))
  n <- attr(logLik(m), "nobs")
  expect_equal(n, 40)
  expect_equal(.msr_model_nobs(m), 160L)
  expect_equal(bic_tab(list(m), mod_names = "native")$BIC, as.numeric(BIC(m)), tolerance = 1e-12)
  expect_equal(aic_tab(list(m), mod_names = "native")$AICc,
    -2 * ll + 2 * k + 2 * k * (k + 1) / (n - k - 1), tolerance = 1e-12)
  expect_equal(aic_tab(list(m), AICc = FALSE, mod_names = "native")$AIC,
    as.numeric(AIC(m)), tolerance = 1e-12)
  expect_identical(serialize(m, NULL), before)
})

test_that("estimated scale penalties remain separate from Cox sample-size metadata", {
  skip_if_not_installed("survival")
  m <- cox_selection_fixture()$model
  # Count-only fixture: no new optimization or scientific result is asserted.
  wrapped <- structure(list(model = m), class = c("fit_clogit", "list"))
  fitted <- structure(list(opt_mod = wrapped,
    scale_est = data.frame(Mean = 20, row.names = "habitat"), shape_est = NULL,
    kernel_inputs = list(kernel = "gaussian")), class = "multiScaleR")
  k <- as.numeric(attr(logLik(m), "df")) + 1
  ll <- as.numeric(logLik(m)); n <- attr(logLik(m), "nobs")
  expect_equal(.msr_selection_parameter_count(fitted), k)
  expect_equal(bic_tab(list(fitted), mod_names = "count fixture")$BIC,
    -2 * ll + log(n) * k, tolerance = 1e-12)
  expect_equal(aic_tab(list(fitted), mod_names = "count fixture")$AICc,
    -2 * ll + 2 * k + 2 * k * (k + 1) / (n - k - 1), tolerance = 1e-12)
})

test_that("actual amt wrappers retain native likelihood criterion sample size", {
  skip_if_not_installed("amt")
  f <- cox_selection_fixture()
  wrapped <- amt::fit_clogit(f$data, f$formula, model = TRUE, x = TRUE, y = TRUE,
    more = list(marker = "original metadata"))
  expect_s3_class(wrapped, "fit_clogit")
  expect_equal(bic_tab(list(wrapped), mod_names = "amt")$BIC,
    as.numeric(BIC(wrapped$model)), tolerance = 1e-12)
  expect_equal(.msr_model_nobs(wrapped), 160L)
  expect_identical(wrapped$more$marker, "original metadata")
})

test_that("Cox event metadata is required and alternative rows remain comparable", {
  skip_if_not_installed("survival")
  m <- cox_selection_fixture()$model
  for (bad in list(NULL, 0, Inf, NA_real_, 2.5, c(40, 41), "40")) {
    broken <- m; broken$nevent <- bad
    expect_error(bic_tab(list(broken), mod_names = "invalid"),
      "native logLik event count")
  }
  other <- m
  rownames(other$model)[1] <- "different alternative"
  expect_error(bic_tab(list(m, other), mod_names = c("a", "b")),
    "different observation sets")
})

test_that("AICc rejects too few events even with many alternatives", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  d <- data.frame(choice = rep(1:2, each = 3), x = rep(-1:1, 2),
    case_ = c(1, 0, 0, 0, 0, 1))
  formula <- as.formula("case_ ~ x + strata(choice)", env = asNamespace("survival"))
  m <- survival::clogit(formula, d, model = TRUE)
  expect_equal(.msr_model_nobs(m), 6L)
  expect_equal(.msr_selection_nobs(m), 2L)
  expect_error(aic_tab(list(m), mod_names = "small"), "n > K \\+ 1")
  expect_equal(aic_tab(list(m), AICc = FALSE, mod_names = "small")$AIC,
    as.numeric(AIC(m)), tolerance = 1e-12)
})

test_that("non-Cox likelihood sample sizes and penalties are unchanged", {
  m <- lm(y ~ x, data = data.frame(y = sin(1:20), x = 1:20))
  expect_identical(.msr_selection_nobs(m), .msr_model_nobs(m))
  expect_equal(bic_tab(list(m), mod_names = "lm")$BIC, as.numeric(BIC(m)), tolerance = 1e-12)
})

test_that("finite Poisson bin scores retain the native cutoff without renormalization", {
  cases <- data.frame(lambda = c(0, 8, 8, 8, 100, 400, 800, 4, 30, 5, 12, 50, 0, 40),
    q = c(0, .2, .2, .9, .1, .7, .01, 1, .95, 1e-8, .5, .99, .5, .01),
    n = c(0, 0, 3, 3, 4, 10, 0, 2, 5, 0, 3, 10, 1, 0),
    K = c(0, 40, 3, 40, 6, 15, 0, 20, 30, 0, 3, 10, 3, 2))
  actual <- direct <- numeric(nrow(cases))
  for (i in seq_len(nrow(cases))) {
    z <- cases[i, ]; y <- matrix(z$n, 1L); pd <- matrix(z$q, 1L); pr <- matrix(1, 1L)
    value <- .msr_gdr_finite_poisson_score_bins(z$lambda, pd, pr, y, y, z$K)
    actual[i] <- value$logLik
    count <- seq.int(z$n, z$K)
    logs <- dpois(count, z$lambda, log = TRUE) + dbinom(z$n, count, z$q, log = TRUE)
    high <- max(logs)
    direct[i] <- if (is.finite(high)) high + log(sum(exp(logs - high))) else high
    expect_identical(value$finite_count_target, "native_finite_K_not_renormalized")
    expect_identical(value$uncertainty, "none_no_fitting_or_backend_admission")
  }
  expect_equal(actual, direct, tolerance = 1e-12)
  expect_equal(actual[1L], 0)
  expect_equal(actual[7L], -800, tolerance = 1e-12)
  expect_identical(actual[13L], -Inf)
  # With no detection and K=0 the native sum still contains only P(N=0).
  zero <- .msr_gdr_finite_poisson_score_bins(8, matrix(0), matrix(1), matrix(0), matrix(0), 0)
  expect_equal(zero$logLik, -8)
  expect_equal(zero$per_event$expected_observed, 0)
})

test_that("finite identity and the retained infinite scorer have distinct target metadata", {
  yd <- rbind(c(1, 0), c(0, 2)); yr <- rbind(c(0, 1), c(1, 1))
  pd <- rbind(c(.1, .2), c(.2, .3)); pr <- rbind(c(.3, .1), c(.1, .2))
  args <- list(expected_abundance = c(8, 15), pd = pd, pr = pr,
    y_distance = yd, y_removal = yr, probability_method = "declared_numeric_bins")
  before <- args
  full <- do.call(.msr_gdr_exact_score_bins, args)
  finite <- do.call(.msr_gdr_finite_poisson_score_bins, c(args, list(cutoff = 2)))
  correction <- ppois(2 - rowSums(yd), args$expected_abundance *
    (1 - rowSums(pd) * rowSums(pr)), log.p = TRUE)
  expect_equal(finite$per_event$logLik, full$per_event$logLik + correction, tolerance = 1e-12)
  expect_equal(finite$per_event$log_abundance_cdf, correction)
  expect_identical(finite$expected_distance, full$expected_distance)
  expect_identical(finite$expected_removal, full$expected_removal)
  expect_identical(finite$per_event$expected_abundance, full$per_event$expected_abundance)
  expect_identical(finite$probability_method, full$probability_method)
  expect_gt(full$logLik, finite$logLik)
  expect_identical(args, before)
  expect_false("abundance_cutoff" %in% names(full))
  for (bad in list(-1, Inf, NA_real_, .5, TRUE, matrix(2), 1:2, .Machine$integer.max + 1))
    expect_error(do.call(.msr_gdr_finite_poisson_score_bins, c(args, list(cutoff = bad))), "cutoff")
  expect_error(do.call(.msr_gdr_finite_poisson_score_bins, c(args, list(cutoff = 1))), "cover")
  above <- args; above$pd[,] <- .5 + 1e-13; above$pr[,] <- .5
  expect_error(do.call(.msr_gdr_finite_poisson_score_bins, c(above, list(cutoff = 2))), "exceed one")
})

test_that("stable log-CDF program preserves four derivatives and vector cross terms", {
  skip_if_not_installed("RTMB")
  reference <- function(x, m) {
    mu <- exp(x); y <- ppois(m, mu, log.p = TRUE)
    r <- exp((m + 1) * x - mu - lgamma(m + 1) - y)
    A <- m + 1 - mu + r; rp <- r * A; Ap <- -mu + rp
    rpp <- rp * A + r * Ap; App <- -mu + rpp
    c(y, -r, -r * A, -rp * A - r * Ap, -rpp * A - 2 * rp * Ap - r * App)
  }
  cases <- data.frame(m = c(rep(0L, 4), rep(3L, 4), rep(40L, 4), rep(400L, 3)),
    mu = c(.001, 1, 10, 800, .001, 3, 50, 800, 10, 40, 100, 1000, 100, 400, 1000))
  gaps <- matrix(0, nrow(cases), 5L)
  for (i in seq_len(nrow(cases))) {
    z <- cases[i, ]; x <- log(z$mu)
    program <- .msr_gdr_poisson_logcdf_program(z$m)
    tape <- RTMB::MakeTape(function(p) sum(program$evaluate(p)), x)
    actual <- vapply(0:4, function(order) {
      if (order == 0L) return(as.numeric(tape(x)))
      derivative <- tape
      for (j in seq_len(order - 1L)) derivative <- derivative$jacfun()
      as.numeric(derivative$jacobian(x))
    }, 0)
    gaps[i, ] <- abs(actual - reference(x, z$m))
  }
  expect_lt(max(gaps[, 1L]), 1e-12)
  expect_lt(max(gaps[, 2L]), 2e-5)
  expect_lt(max(gaps[, 3:5]), 1e-4)
  m <- c(0L, 3L, 40L, 400L); x <- log(c(800, 50, 1000, 400))
  program <- .msr_gdr_poisson_logcdf_program(m)
  tape <- RTMB::MakeTape(function(p) sum(program$evaluate(p)), x)
  expected <- vapply(seq_along(x), function(i) reference(x[i], m[i])[2:3], numeric(2L))
  expect_equal(as.numeric(tape$jacobian(x)), expected[1L, ], tolerance = 2e-5)
  expect_equal(tape$jacfun()$jacobian(x), diag(expected[2L, ]), tolerance = 1e-4)
  expect_true(program$requires_serial_R_callbacks)
  expect_false(program$portable_objective)
  expect_false(program$fit_admission_established)
  # The old nonlog route underflows. Its finite log target is retained.
  expect_identical(log(ppois(0, 800)), -Inf)
  expect_equal(.msr_gdr_poisson_logcdf_program(0)$evaluate(log(800)), -800, tolerance = 1e-12)
  expect_equal(program$evaluate(rep(-1000, 4)), rep(0, 4))
  expect_error(program$evaluate(1), "log-rate")
  expect_error(program$evaluate(c(Inf, 0, 0, 0)), "log-rate")
  expect_error(program$evaluate(rep(1000, 4)), "representable")
  for (bad in list(numeric(), -1, NA_real_, Inf, .2, TRUE, matrix(2), .Machine$integer.max + 1))
    expect_error(.msr_gdr_poisson_logcdf_program(bad), "cutoffs")
})

test_that("finite Poisson audit adds an exact reference while retaining native gates", {
  model <- gdr_guard_fixture()$P
  before <- gdr_guard_state(model)
  audit <- audit_gdr_integration(model, c(5, 400))
  expect_identical(audit$finite_poisson$status, "evaluated_closed_form_reference")
  expect_identical(audit$finite_poisson$comparison$cutoff, audit$cutoffs$cutoff)
  expect_equal(audit$finite_poisson$comparison$nll, audit$cutoffs$nll, tolerance = 1e-12)
  expect_true(all(audit$finite_poisson$comparison$agrees))
  expect_false(audit$finite_poisson$fitting_admission_established)
  expect_false(audit$cutoffs$truncation_qualified[2L])
  expect_identical(audit$parameter_count, .msr_parameter_count(model))
  expect_identical(gdr_guard_state(model), before)
  nb <- audit_gdr_integration(gdr_guard_fixture()$NB)
  expect_identical(nb$finite_poisson$status, "not_qualified_for_NB")
  expect_null(nb$finite_poisson$comparison)
  for (key in c("exp", "uniform")) {
    native <- unmarked::gdistremoval(data = model@data, keyfun = key, mixture = "P",
      output = "abund", engine = "C", K = 40, se = FALSE, threads = 1,
      starts = c(log(8), if (key != "uniform") log(40), qlogis(.3)), control = list(maxit = 0))
    result <- audit_gdr_integration(native, 5)
    expect_true(all(result$finite_poisson$comparison$agrees))
    expect_equal(result$finite_poisson$comparison$nll, result$cutoffs$nll, tolerance = 1e-12)
    if (key == "exp") expect_false(result$distance$analytic_equivalence_qualified)
  }
  cat("GDR_FINITE_K_REFERENCE native_C_only=1 AD_derivatives=4 fitting_admission=0\n")
})

test_that("failed native probability records remain visible in the closed-form audit", {
  local_mocked_bindings(.msr_gdr_audit_native_bins = function(...) stop("retained native-bin failure"))
  result <- audit_gdr_integration(gdr_guard_fixture()$P, 400)
  expect_identical(result$references$native_bins$error, "retained native-bin failure")
  expect_true(all(is.finite(result$cutoffs$nll)))
  expect_true(all(is.na(result$finite_poisson$comparison$nll)))
  expect_false(any(result$finite_poisson$comparison$agrees))
  expect_match(result$finite_poisson$evaluations[[1L]]$error, "Native-bin reference failed")
  expect_false(result$fitting_admission_established)
})

test_that("a cutoff ladder constructs each invariant score reference once", {
  original <- .msr_gdr_exact_score_bins
  calls <- 0L
  local_mocked_bindings(.msr_gdr_exact_score_bins = function(...) {
    calls <<- calls + 1L
    original(...)
  })
  audit_gdr_integration(gdr_guard_fixture()$P, c(5, 400))
  expect_identical(calls, 2L) # Stable and native bins, independent of ladder length.
})

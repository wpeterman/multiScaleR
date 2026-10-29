test_that("stream kernel arithmetic preserves values, derivatives and defaults", {
  # A fixed arithmetic fixture, without responses or model fitting.
  counts <- matrix(rep(1:4, 18), 8, 9)
  values <- matrix(sin(seq_len(72)), 8, 9) * counts
  b <- list(bin_rep = seq(0, 1, length.out = 9),
    vsum = list(a = values, b = values[8:1, ]),
    csum = list(a = counts, b = counts[8:1, ]))
  for (kernel in c("gaussian", "exp", "expow")) {
    for (standardize in c(FALSE, TRUE)) {
      e <- .msr_prepare_kernel_engine(b, rows = c(8, 2, 4, 1, 6),
        kernel = kernel, standardize = standardize)
      expect_identical(e$implementation, "reference")
      par <- c(log(c(.13, .4)), if (kernel == "expow") log(c(.8, 3)))
      old <- .msr_eval_kernel_engine(e, par)
      e$implementation <- "stream"
      new <- .msr_eval_kernel_engine(e, par)
      expect_equal(new, old, tolerance = 1e-12)
      for (j in seq_along(par)) {
        plus <- minus <- par
        plus[j] <- plus[j] + 1e-5; minus[j] <- minus[j] - 1e-5
        numerical <- (.msr_eval_kernel_engine(e, plus)$value -
          .msr_eval_kernel_engine(e, minus)$value) / 2e-5
        analytic <- matrix(0, nrow(numerical), ncol(numerical))
        k <- (j - 1L) %% e$p + 1L
        analytic[, k] <- if (j <= e$p) new$d_scale[, k] else new$d_shape[, k]
        expect_equal(unname(numerical), analytic, tolerance = 1e-7)
      }
      e$implementation <- NULL
      expect_equal(.msr_eval_kernel_engine(e, par), old)
    }
  }
  expect_error(.msr_prepare_kernel_engine(b, implementation = "invalid"), "arg")
  expect_identical(.msr_profile_control(list())$kernel_engine, "reference")
  expect_identical(.msr_profile_control(list(kernel_engine = "stream"))$kernel_engine, "stream")
  expect_error(.msr_profile_control(list(kernel_engine = NA)), "kernel_engine")
})

test_that("cached NB2 terms agree across dispersion changes and independent responses", {
  y <- c(0, 1, 4, 12, 3); w <- c(1, 2, 1, 4, 1); eta <- log(c(.2, 2, 3, 9, 4))
  reference <- .msr_nb2_normalizer(y, w)
  cached <- .msr_nb2_normalizer(y, w, cache = TRUE)
  for (alpha in c(0, 1e-10, .5, .5, 100, 0, .5)) {
    expect_identical(cached(alpha), reference(alpha))
    expect_equal(.msr_nb2_terms(y, eta, alpha, w, cached),
                 .msr_nb2_terms(y, eta, alpha, w, reference), tolerance = 1e-14)
  }
  other <- .msr_nb2_normalizer(y + 1, w, cache = TRUE)
  expect_false(identical(other(.5), cached(.5)))
  expect_identical(cached(.5), reference(.5))
})

test_that("stream validation retains finite cell and selected-row contracts", {
  b <- list(bin_rep = c(0, 1), vsum = list(a = matrix(1:6, 3)),
            csum = list(a = matrix(1, 3, 2)))
  for (mode in c("reference", "stream")) {
    bad <- b; bad$vsum$a[1] <- Inf
    expect_error(.msr_prepare_kernel_engine(bad, implementation = mode), "sum/count")
    bad <- b; bad$csum$a[1] <- -1
    expect_error(.msr_prepare_kernel_engine(bad, implementation = mode), "sum/count")
    bad <- b; bad$csum$a[1] <- 0
    expect_error(.msr_prepare_kernel_engine(bad, implementation = mode), "sum/count")
    bad$vsum$a[1, ] <- bad$csum$a[1, ] <- 0
    expect_error(.msr_prepare_kernel_engine(bad, implementation = mode), "support")
    expect_s3_class(.msr_prepare_kernel_engine(bad, rows = 2:3,
      implementation = mode), "msr_kernel_engine")
  }
  original <- b
  e <- .msr_prepare_kernel_engine(b, implementation = "stream")
  e$V$a[1] <- 10
  expect_identical(b, original)
})

test_that("binned all-row reuse preserves subset and duplicate row semantics", {
  b <- list(bin_rep = c(0, .2, .4), n_points = 4L,
    vsum = list(a = matrix(1:12, 4)), csum = list(a = matrix(2, 4, 3)))
  for (kernel in c("gaussian", "exp", "expow", "fixed")) {
    evaluate <- function(rows) .msr_binned_kernel_means(b, "a", "a", .3, 2,
                                                       kernel, rows)
    all <- evaluate(NULL)
    expect_identical(evaluate(1:4), all)
    expect_identical(evaluate(4:1), all[4:1, , drop = FALSE])
    expect_identical(evaluate(c(2L, 2L, 1L)), all[c(2L, 2L, 1L), , drop = FALSE])
  }
})

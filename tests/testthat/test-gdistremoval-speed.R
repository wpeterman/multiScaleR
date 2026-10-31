test_that("closed distance integrals agree with independent quadrature", {
  breaks <- c(0, 0.1, 5, 25, 60, 99.999, 100)
  scales <- c(0.05, 1, 10, 40, 100, 1000, 1e6)
  for (key in c("halfnorm", "exp")) {
    got <- .msr_gdr_distance(key, scales, breaks, length(scales))
    expected <- t(vapply(scales, function(scale) {
      vapply(seq_len(length(breaks) - 1L), function(j) {
        stats::integrate(function(r) 2 * r / 100^2 * if (key == "halfnorm")
          exp(-r^2 / (2 * scale^2)) else exp(-r / scale),
          breaks[j], breaks[j + 1L], rel.tol = 1e-10, abs.tol = 1e-14)$value
      }, numeric(1))
    }, numeric(length(breaks) - 1L)))
    expect_equal(got, expected, tolerance = 1e-10)
    # Repeated scales, missing rows, and one-bin frames preserve their layouts.
    repeated <- .msr_gdr_distance(key, c(40, NA, 40, 0, Inf), breaks, 5)
    expect_equal(repeated[c(1, 3), ], got[rep(4, 2), ])
    expect_true(all(is.na(repeated[2, ])))
    expect_equal(repeated[4, ], rep(0, length(breaks) - 1L))
    expect_equal(repeated[5, ], diff((breaks / 100)^2))
    expect_identical(dim(.msr_gdr_distance(key, 40, c(0, 100), 1)), c(1L, 1L))
    expect_true(all(is.na(.msr_gdr_distance(key, c(NA, NA), c(0, 100), 2))))
    expect_error(.msr_gdr_distance(key, -1, c(0, 100), 1), "nonnegative")
  }
})

test_that("exponential integrals retain tiny and nearly uniform probabilities", {
  scales <- c(1e-100, 0.01, 100, 1e100, 1e300, 0, Inf, NA)
  got <- .msr_gdr_distance("exp", scales, c(0, 100), length(scales))[, 1]
  # When the scale is tiny relative to the circle, the untruncated integral
  # 2*(scale/radius)^2 is an accurate independent limiting expression.
  expect_equal(got[1] / (2 * (scales[1] / 100)^2), 1, tolerance = 1e-12)
  expect_equal(got[2], 2 * (scales[2] / 100)^2, tolerance = 1e-12)
  expect_equal(got[3], 2 * (1 - 2 / exp(1)), tolerance = 1e-12)
  expect_equal(got[4:5], c(1, 1), tolerance = 1e-12)
  expect_equal(got[6:7], c(0, 1))
  expect_true(is.na(got[8]))
  expect_true(all(got[!is.na(got)] >= 0 & got[!is.na(got)] <= 1 + 1e-12))
})

test_that("the GDR refitter skips only objective Hessians", {
  f <- gdr_fixture()
  frame <- f$frame[, 1]
  model <- unmarked::gdistremoval(~hab, ~1, ~tod, ~1, data = frame,
    keyfun = "uniform", K = 35, se = TRUE)
  context <- build_opt_context(model, f$inputs$raw_cov,
    unit_conv = f$inputs$unit_conv, refit_fn = refit_gdistremoval)
  run <- function(return_model = NULL, ctx = context) kernel_scale_fn(
    30 / f$inputs$unit_conv, f$inputs$d_list, f$inputs$raw_cov,
    f$inputs$kernel, model, mod_return = return_model, opt_context = ctx)
  scalar <- run()
  returned <- run(TRUE)$mod
  reference_context <- context
  reference_context$refit_fn <- NULL
  reference <- run(TRUE, reference_context)$mod
  expect_equal(scalar, -as.numeric(unmarked::logLik(reference)), tolerance = 1e-10)
  expect_equal(unmarked::coef(returned), unmarked::coef(reference), tolerance = 1e-10)
  expect_equal(unmarked::vcov(returned), unmarked::vcov(reference), tolerance = 1e-10)
  expect_true(all(is.finite(unmarked::vcov(returned))))
  expect_identical(returned@formlist, model@formlist)
  expect_identical(returned@data@yDistance, frame@yDistance)
  expect_identical(returned@data@yRemoval, frame@yRemoval)
  expect_identical(model@call$se, TRUE)
  lean <- refit_gdistremoval(model, frame, list(return_model = FALSE))
  expect_error(unmarked::vcov(lean))
  expect_equal(unmarked::coef(lean), unmarked::coef(model), tolerance = 1e-10)
  lean_final <- refit_gdistremoval(lean, frame)
  expect_error(unmarked::vcov(lean_final))
  restored <- unserialize(serialize(refit_gdistremoval, NULL))
  expect_equal(unmarked::coef(restored(model, frame)), unmarked::coef(model))
  expect_error(refit_gdistremoval(stats::lm(1:3 ~ 1), frame), "unmarkedFitGDR")
  expect_error(refit_gdistremoval(model, frame, context = FALSE), "context")
})

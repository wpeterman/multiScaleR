# Skipped on CRAN to keep the package check within CRAN's time limit. CI sets
# NOT_CRAN and runs this file.
skip_on_cran()

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
    # Repeated scales, missing rows, and one-bin frames keep their layouts.
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

test_that("exponential integrals keep tiny and nearly uniform probabilities", {
  scales <- c(1e-100, 0.01, 100, 1e100, 1e300, 0, Inf, NA)
  got <- .msr_gdr_distance("exp", scales, c(0, 100), length(scales))[, 1]
  # For a scale that is tiny relative to the circle, the untruncated integral
  # 2 * (scale / radius)^2 is an accurate independent limit.
  expect_equal(got[1] / (2 * (scales[1] / 100)^2), 1, tolerance = 1e-12)
  expect_equal(got[2], 2 * (scales[2] / 100)^2, tolerance = 1e-12)
  expect_equal(got[3], 2 * (1 - 2 / exp(1)), tolerance = 1e-12)
  expect_equal(got[4:5], c(1, 1), tolerance = 1e-12)
  expect_equal(got[6:7], c(0, 1))
  expect_true(is.na(got[8]))
  expect_true(all(got[!is.na(got)] >= 0 & got[!is.na(got)] <= 1 + 1e-12))
})

test_that("combined probabilities match unmarked's components visit by visit", {
  model <- gdr_survey_fixture()$model
  out <- predict_gdistremoval(model)
  p <- unmarked::getP(model)
  for (t in 1:2) {
    rows <- out$period == t
    expect_equal(out$availability[rows], rowSums(p$rem[, , t]), tolerance = 1e-10)
    expect_equal(out$phi[rows], p$phi[, t], tolerance = 1e-10)
  }
  # unmarked integrates the first visit with its own distance scale.
  expect_equal(out$perceptibility[out$period == 1], rowSums(p$dist[, , 1]),
    tolerance = 1e-6)
  # Every visit uses its own scale; here the two visits differ.
  sigma <- unmarked::predict(model, "dist", level = NULL)$Predicted
  expect_gt(abs(sigma[2] - sigma[1]), 1)
  quadrature <- vapply(sigma, function(s) stats::integrate(function(r)
    2 * r / 100^2 * exp(-r^2 / (2 * s^2)), 0, 100, rel.tol = 1e-10)$value, numeric(1))
  expect_equal(out$perceptibility, quadrature, tolerance = 1e-8)
  expect_equal(out$detection, out$phi * out$availability * out$perceptibility)
  expect_equal(out$expected_count, out$abundance * out$detection)
})

test_that("combined means agree with independent analytical probabilities", {
  model <- gdr_survey_oracle()
  out <- predict_gdistremoval(model)
  sigma <- rep(c(50, 100), nrow(model@data@yDistance))
  pd <- 2 * sigma^2 / 100^2 * (1 - exp(-100^2 / (2 * sigma^2)))
  pa <- 1 - (1 - 0.3)^5
  expect_equal(out$phi, rep(0.7, nrow(out)))
  expect_equal(out$availability, rep(pa, nrow(out)))
  expect_equal(out$perceptibility, pd, tolerance = 1e-12)
  expect_equal(out$abundance, rep(8, nrow(out)))
  expect_equal(out$density, rep(8 / pi, nrow(out)))
  expect_equal(out$detection, 0.7 * pa * pd)
  expect_equal(out$expected_count, 8 * 0.7 * pa * pd)
  expect_identical(attr(out, "density_units"), "ha")
  expect_equal(attr(out, "survey_area"), pi)
  expect_identical(attr(out, "prediction_contract")$uncertainty$intervals, "not_provided")
  expect_gt(out$perceptibility[2], out$perceptibility[1])
})

test_that("density units and the ZIP mean follow the native likelihood", {
  model <- gdr_survey_oracle("ZIP", "density")
  out <- predict_gdistremoval(model)
  expect_equal(out$lambda, rep(8, nrow(out)))
  expect_equal(out$zero_probability, rep(0.4, nrow(out)))
  expect_equal(out$density, rep(8 * 0.6, nrow(out)))
  expect_equal(out$abundance, rep(8 * 0.6 * pi, nrow(out)))
  model@unitsOut <- "kmsq"
  out <- predict_gdistremoval(model)
  expect_equal(out$abundance, rep(8 * 0.6 * pi / 100, nrow(out)))
  km <- model
  km@data@dist.breaks <- km@data@dist.breaks / 1000
  km@data@unitsIn <- "km"
  km@estimates@estimates$dist@estimates[1] <- log(0.05)
  converted <- predict_gdistremoval(km)
  expect_equal(converted$perceptibility, out$perceptibility)
  expect_equal(converted$abundance, out$abundance)
})

test_that("distance keys and zero removal probabilities stay valid", {
  for (key in c("uniform", "exp", "hazard")) {
    model <- gdr_survey_oracle(keyfun = key)
    out <- predict_gdistremoval(model)
    expected <- if (key == "uniform") 1 else vapply(rep(c(50, 100), 24), function(scale) {
      stats::integrate(function(r) 2 * r / 100^2 * if (key == "exp")
        exp(-r / scale) else (1 - exp(-(r / scale)^-2)), 0, 100)$value
    }, numeric(1))
    expect_equal(out$perceptibility, rep(expected, length.out = nrow(out)), tolerance = 1e-7)
  }
  # Zero probabilities would divide by zero in unmarked's removal cells.
  p <- rbind(c(0, 0.5, 0), c(1, 0, 1), c(0, 0, 0))
  expect_equal(.msr_gdr_removal(p, c(1, 2, 2)), rbind(c(0, 0.75, 0), c(1, 0, 0), c(0, 0, 0)))
  expect_equal(.msr_gdr_distance("halfnorm", c(0, Inf, NA), c(0, 100), 3),
    matrix(c(0, 1, NA), ncol = 1))
})

test_that("new survey protocols keep native rows, fitted normalization and inputs", {
  f <- gdr_survey_fixture()
  frame <- f$object$opt_mod@data
  original <- frame
  frame@period.lengths <- c(2, 2, 3)
  keys <- data.frame(site = seq_len(nrow(frame@yDistance)))
  manual <- frame
  manual@siteCovs <- prepare_unmarked_data(f$object, frame@siteCovs, f$inputs, join_by = keys)
  got <- predict_gdistremoval(f$object, frame, f$inputs, join_by = keys)
  expected <- predict_gdistremoval(f$object$opt_mod, manual)
  expect_equal(as.matrix(got), as.matrix(expected))
  expect_identical(attr(got, "prediction_contract")$uncertainty$kernel_parameters, "fixed_at_estimate")
  expect_equal(frame@siteCovs, original@siteCovs)
  expect_equal(frame@period.lengths, c(2, 2, 3))
  # A longer count raises availability.
  expect_gt(got$availability[1], predict_gdistremoval(f$object)$availability[1])
  expect_error(predict_gdistremoval(f$object, frame, f$inputs), "join_by")
  expect_error(predict_gdistremoval(f$object, kernel_inputs = f$inputs), "explicit prediction frame")
  expect_error(predict_gdistremoval(f$model, newdata = data.frame(x = 1)), "unmarkedFrameGDR")
  expect_error(predict_gdistremoval(stats::lm(1:3 ~ 1)), "unmarkedFitGDR")
  changed <- frame
  changed@dist.breaks <- c(0, 20, 60, 100)
  expect_error(predict_gdistremoval(f$model, changed), "survey footprint")
  changed <- frame
  changed@siteCovs$hab[2] <- NA_real_
  predicted <- predict_gdistremoval(f$model, changed)
  expect_true(all(is.na(predicted$abundance[3:4])))
  expect_equal(predicted$perceptibility[-(3:4)], got$perceptibility[-(3:4)])
})

test_that("single-visit and uniform fits need no inactive component", {
  f <- gdr_survey_fixture()
  frame <- f$frame[, 1]
  model <- unmarked::gdistremoval(~hab, ~1, ~tod, ~1, data = frame,
    keyfun = "uniform", K = 35, se = FALSE)
  out <- predict_gdistremoval(model)
  expect_equal(out$phi, rep(1, nrow(out)))
  expect_equal(out$perceptibility, rep(1, nrow(out)))
  expect_identical(.msr_gdr_components(model), c("lambda", "rem"))
})

test_that("native random effects keep explicit prediction policies", {
  skip_if_not_installed("TMB")
  f <- gdr_survey_fixture()
  frame <- f$frame
  frame@siteCovs$route <- factor(rep(seq_len(6), each = 4))
  model <- suppressWarnings(unmarked::update(f$model, data = frame,
    lambdaformula = ~hab + (1 | route), engine = "TMB", se = FALSE))
  expect_false(is.null(model@TMB))
  for (policy in list(NULL, NA)) {
    expected <- unmarked::predict(model, "lambda", newdata = frame@siteCovs,
      re.form = policy, level = NULL)$Predicted
    got <- predict_gdistremoval(model, re.form = policy)
    expect_equal(got$lambda, rep(expected, each = frame@numPrimary))
  }
  # A zero group deviation is an explicit scenario for unsampled groups.
  unseen <- frame
  unseen@siteCovs$route <- factor(rep("unsampled", nrow(frame@yDistance)))
  expected <- unmarked::predict(model, "lambda", newdata = unseen@siteCovs,
    re.form = NA, level = NULL)$Predicted
  got <- predict_gdistremoval(model, unseen, re.form = NA)
  expect_equal(got$lambda, rep(expected, each = frame@numPrimary))
  expect_identical(attr(got, "prediction_contract")$effects$response,
    "zero_deviation_reference")
})

test_that("simulation reproduces paired totals and conditional means", {
  model <- gdr_survey_oracle("ZIP")
  # Four sites suffice to check Monte Carlo moments without a slow refit.
  model@data <- model@data[1:4, ]
  set.seed(81)
  simulated <- simulate_gdistremoval(model, nsim = 600)
  set.seed(81)
  expect_identical(simulate_gdistremoval(model, nsim = 2),
    structure(simulated[1:2], prediction_contract = attr(simulated, "prediction_contract")))
  expect_true(all(vapply(simulated, function(frame) {
    totals <- function(x) cbind(rowSums(x[, 1:3, drop = FALSE]), rowSums(x[, 4:6, drop = FALSE]))
    identical(totals(frame@yDistance), totals(frame@yRemoval)) && identical(frame@y, frame@yRemoval)
  }, logical(1))))
  empirical <- Reduce(`+`, lapply(simulated, function(frame) cbind(
    rowSums(frame@yRemoval[, 1:3, drop = FALSE]), rowSums(frame@yRemoval[, 4:6, drop = FALSE])))) / length(simulated)
  expected <- matrix(predict_gdistremoval(model)$expected_count, ncol = 2, byrow = TRUE)
  expect_lt(max(abs(empirical - expected)), 0.45)
  expect_true(attr(simulated, "prediction_contract")$uncertainty$observation_noise)
  expect_identical(attr(simulated, "prediction_contract")$quantity, "replicate_response")
  expect_error(simulate_gdistremoval(model, nsim = 0), "nsim")
  missing <- model@data
  missing@yDistance[2, 4] <- NA
  simulated <- simulate_gdistremoval(model, newdata = missing)[[1]]
  expect_true(all(is.na(simulated@yDistance[2, 4:6])))
  expect_true(all(is.na(simulated@yRemoval[2, 4:6])))
})

test_that("NB simulations keep overdispersion and shared site abundance", {
  model <- gdr_survey_oracle("NB")
  model@data <- model@data[1, ]
  expected <- predict_gdistremoval(model)
  set.seed(482)
  simulated <- simulate_gdistremoval(model, nsim = 1500)
  counts <- t(vapply(simulated, function(frame)
    c(sum(frame@yRemoval[1, 1:3]), sum(frame@yRemoval[1, 4:6])), numeric(2)))
  mean <- expected$expected_count
  # NB size 2: variance of N is 8 + 8^2 / 2, thinned by each visit's detection.
  theoretical_variance <- mean + mean^2 / 2
  theoretical_covariance <- expected$detection[1] * expected$detection[2] * (8 + 8^2 / 2)
  expect_lt(max(abs(colMeans(counts) - mean)), 0.35)
  expect_lt(max(abs(apply(counts, 2, var) / theoretical_variance - 1)), 0.25)
  expect_lt(abs(cov(counts)[1, 2] / theoretical_covariance - 1), 0.3)
})

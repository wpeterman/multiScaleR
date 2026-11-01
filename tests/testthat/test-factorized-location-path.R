location_path_fixture <- function(kernel, family, response, surface, design) {
  if (surface) f <- private_sq_fixture(kernel, family, response = response) else {
    f <- bayesian_scale_fixture(n_site = 4, n_year = 3, annual = FALSE, kernel = kernel)
    if (family == "binomial") f$data$count <- as.integer(f$data$count > 2)
    f$data$batch <- interaction(f$data$year, as.integer(f$data$site) <= 2)
    f$formula <- if (response) count ~ habitat + temp + (1 | batch) + offset(log(effort)) else
      count ~ habitat + temp + offset(log(effort))
  }
  f$data$band <- factor(as.integer(f$data$site) <= 2)
  scale <- if (design == "contrast") ~1 + temp + (1 | year) else ~0 + band + (1 | year)
  .msr_scale_prepare(f$formula, f$inputs, f$data, stats::setNames(list(scale), if (surface) f$spec$covariates else "habitat"), family,
    "key", reference_scale = 250)
}

test_that("common location paths preserve normalized group limits across supported designs", {
  cases <- 0L; maximum <- 0; pointwise <- 0
  ratios <- c(-.75, 0, .75)
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (response_zero in if (response) c(FALSE, TRUE) else FALSE)
          for (surface in c(FALSE, TRUE)) for (design in c("contrast", "cell_means")) {
            cases <- cases + 1L
            s <- location_path_fixture(kernel, family, response, surface, design)
            if (family == "binomial") s$trials[] <- 4
            if (response_zero) s <- .msr_marginal_reduce(s, "response")
            p <- private_sq_parameters(s); p$beta[3L] <- .12
            p$scale_beta[] <- c(.2, -.3); p$log_scale_sd[] <- log(8)
            p$scale_z <- p$response_z <- numeric()
            before <- list(s, p, .Random.seed)
            direction <- .msr_factorized_location_direction(s$scales[[1L]], s$parameters)
            expect_equal(unname(drop(s$scales[[1L]]$X %*% direction)), rep(1, length(s$y)), tolerance = 0)
            broad <- mixture_native_nll(s, p, "broad")
            nearest <- mixture_native_nll(s, p, "nearest")
            expected <- vapply(pnorm(ratios), function(w) mixture_native_total(broad, nearest, w), 0)
            actual <- vapply(ratios, function(ratio) sum(vapply(1:3, function(i)
              .msr_factorized_group_mixture_limit(.msr_factorized_group(s, i), p, ratio)$nll, 0)), 0)
            maximum <- max(maximum, abs(actual - expected))
            label <- paste(kernel, family, response, response_zero, surface, design)
            expect_equal(unname(actual), unname(expected), tolerance = 1e-10, info = label)
            if (!surface) {
              full <- .msr_factorized_path_limits(s, p)
              zero <- p; zero$log_scale_sd <- numeric()
              boundary <- .msr_factorized_path_limits(s, zero)
              expect_equal(unname(vapply(full$references[1:3], `[[`, 0, "objective")),
                unname(expected), tolerance = 1e-10, info = label)
              expect_identical(full$references, boundary$references)
            } else expect_error(.msr_factorized_path_limits(s, p), "Weighted SQ")
            # Large finite paths verify the predictor direction away from each
            # switch boundary; they do not certify a finite optimum or infinity.
            g <- .msr_factorized_group(s, 2L)
            finite <- matrix(0, 2L, length(ratios))
            for (j in seq_along(ratios)) for (i in 1:2) {
              q <- p; q$scale_beta <- p$scale_beta + 1e6 * ratios[j] * direction
              q$log_scale_sd <- log(1e6); q$scale_z <- c(-2, 2)[i]
              finite[i, j] <- sum(.msr_factorized_response_nll(g, q))
            }
            reference <- matrix(rep(c(nearest[2L], broad[2L]), length(ratios)), 2L)
            pointwise <- max(pointwise, abs(finite - reference))
            expect_equal(unname(finite), unname(reference), tolerance = 1e-10, info = label)
            expect_identical(list(s, p, .Random.seed), before)
          }
  expect_equal(cases, 48L)
  cat(sprintf("\nLOCATION_PATH cases=48 native_error=%.12g pointwise_error=%.12g\n", maximum, pointwise))
})

test_that("location directions retain strict layout and unsupported-design failures", {
  s <- location_path_fixture("gaussian", "gaussian", FALSE, FALSE, "contrast")
  sm <- s$scales[[1L]]
  direction <- .msr_factorized_location_direction(sm, s$parameters)
  expect_identical(direction, c(1, 0))
  bad <- sm; bad$X[, 1L] <- seq(-1, 1, length.out = nrow(bad$X))
  expect_null(.msr_factorized_location_direction(bad, s$parameters))
  bad <- sm; bad$X[1L, 1L] <- NA_real_
  expect_null(.msr_factorized_location_direction(bad, s$parameters))
  bad <- sm; bad$X <- NULL
  expect_null(.msr_factorized_location_direction(bad, s$parameters))
  bad <- sm; bad$beta_index <- rev(bad$beta_index)
  expect_null(.msr_factorized_location_direction(bad, s$parameters))
  bad <- sm; bad$random[[1L]]$sd_index <- 2L
  expect_null(.msr_factorized_location_direction(bad, s$parameters))
  bad <- sm; bad$X <- cbind(bad$X, 1)
  expect_null(.msr_factorized_location_direction(bad, s$parameters))
  cell <- location_path_fixture("gaussian", "gaussian", FALSE, FALSE, "cell_means")
  expect_identical(.msr_factorized_location_direction(cell$scales[[1L]], cell$parameters), c(1, 1))
  # Prepared groups may retain unobserved factor columns without full local rank.
  reduced <- sm; reduced$X <- cbind(rep(0, 2), rep(1, 2))
  expect_equal(drop(reduced$X %*% .msr_factorized_location_direction(reduced, s$parameters)), c(1, 1))
})

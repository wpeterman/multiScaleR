test_that("prepared group paths retain native endpoint and group-mixture likelihoods", {
  ratios <- c(-.75, 0, .75)
  mass <- vapply(ratios, function(b) integrate(dnorm, -b, Inf, rel.tol = 1e-12)$value, 0)
  expect_equal(mass, pnorm(ratios), tolerance = 1e-12)
  cases <- 0L; maximum <- 0
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (kind in c("mean", "sq", "unique"))
          for (response_zero in if (response) c(FALSE, TRUE) else FALSE) {
            cases <- cases + 1L
            f <- mixture_fixture(kernel, family, kind != "mean", response)
            s <- f$spec
            if (kind == "unique") {
              d <- f$inputs$d_list[[2L]]; nearest <- min(d)
              d[2L] <- (nearest + min(d[d > nearest])) / 2
              f$inputs$d_list[[2L]] <- d
              s <- .msr_scale_prepare(f$formula, f$inputs, f$data, f$scales,
                family, "key", reference_scale = 250)
            }
            if (family == "binomial") s$trials[] <- 4
            if (response_zero) s <- .msr_marginal_reduce(s, "response")
            p <- private_sq_parameters(s); p$beta[3L] <- .12
            p$scale_beta[] <- 0; p$log_scale_sd[] <- log(8)
            p$scale_z <- p$response_z <- numeric()
            before <- s; parameters_before <- p; rng <- .Random.seed
            broad <- as.numeric(mixture_native_nll(s, p, "broad"))
            nearest <- as.numeric(mixture_native_nll(s, p, "nearest"))
            expected <- vapply(mass, function(w) vapply(seq_along(broad),
              function(i) mixture_native_total(broad[i], nearest[i], w), 0), numeric(3L))
            actual <- matrix(0, 3L, length(ratios)); endpoints <- matrix(0, 3L, 2L)
            for (j in seq_along(ratios)) for (i in 1:3) {
              result <- .msr_factorized_group_mixture_limit(.msr_factorized_group(s, i), p, ratios[j])
              actual[i, j] <- result$nll
              if (j == 1L) endpoints[i, ] <- result$component_nll
            }
            label <- paste(kernel, family, response, kind, response_zero)
            # Aggregate groups/ratios so repeated cases retain useful failure
            # locations without inflating metadata assertions per evaluation.
            expect_equal(endpoints, unname(cbind(broad, nearest)), tolerance = 1e-10, info = label)
            expect_equal(actual, expected, tolerance = 1e-10, info = label)
            expect_equal(colSums(actual), vapply(mass, function(w)
              mixture_native_total(broad, nearest, w), 0), tolerance = 1e-10, info = label)
            expect_equal(colSums(actual), vapply(pnorm(ratios), function(w)
              .msr_factorized_mixture_nll(broad, nearest, w), 0), tolerance = 1e-10, info = label)
            expect_identical(list(specification = s, parameters = p, rng = .Random.seed),
              list(specification = before, parameters = parameters_before, rng = rng))
            maximum <- max(maximum, abs(actual - expected), abs(endpoints - cbind(broad, nearest)))
          }
  expect_equal(cases, 36L)
  cat(sprintf("\nSQ_GROUP_LIMIT cases=36 groups=3 ratios=3 native_error=%.12g\n", maximum))
})

test_that("group path metadata and restrictions keep likelihood limits separate from fitting", {
  s <- private_sq_fixture(response = TRUE)$spec; g <- .msr_factorized_group(s, 2L)
  p <- private_sq_parameters(g); p$scale_z <- p$response_z <- numeric()
  result <- .msr_factorized_group_mixture_limit(g, p, 0)
  expect_identical(result$target, "marginal_likelihood_limit")
  expect_false(result$finite_scale_estimate)
  expect_equal(result$broad_weight, .5)
  shifted <- p; shifted$scale_beta <- shifted$scale_beta + 10
  shifted$log_scale_sd <- shifted$log_scale_sd + 1
  expect_identical(.msr_factorized_group_mixture_limit(g, shifted, 0), result)
  expect_error(.msr_factorized_group_mixture_limit(s, p, 0), "one recoded")
  expect_error(.msr_factorized_group_mixture_limit(.msr_marginal_reduce(g, "scale"), p, 0), "one scale grouping")
  bad <- g; bad$priors <- list(beta_sd = 2.5)
  expect_error(.msr_factorized_group_mixture_limit(bad, p, 0), "neutral specification")
  bad <- g; bad$scales[[1L]]$X[,] <- 2
  expect_error(.msr_factorized_group_mixture_limit(bad, p, 0), "canonical scale intercept")
  bad <- p; bad$scale_z <- 0
  expect_error(.msr_factorized_group_mixture_limit(g, bad, 0), "population layout")
  bad <- p; bad$beta[1L] <- Inf
  expect_error(.msr_factorized_group_mixture_limit(g, bad, 0), "finite population")
  bad <- p; names(bad)[1L] <- "wrong"
  expect_error(.msr_factorized_group_mixture_limit(g, bad, 0), "population layout")
  for (ratio in list(NA_real_, Inf, c(0, 1), "0"))
    expect_error(.msr_factorized_group_mixture_limit(g, p, ratio), "finite real")
  with_mocked_bindings({
    expect_error(.msr_factorized_group_mixture_limit(g, p, 0), "Finite normalized endpoint")
  }, .msr_factorized_response_nll = function(...) Inf, .package = "multiScaleR")
  expect_error(.msr_factorized_limit_spec(s, "broad"), "Weighted SQ")
  expect_error(.msr_factorized_limit_audit(s, p, result$nll), "Weighted SQ")
})

test_that("prepared paths preserve raw physical units and exact Poisson reductions", {
  for (kernel in c("gaussian", "exp")) {
    f <- private_sq_fixture(kernel, response = TRUE); s <- f$spec
    moved <- f$inputs
    moved$d_list <- lapply(moved$d_list, function(d) d * moved$unit_conv)
    moved$unit_conv <- 1
    other <- .msr_scale_prepare(f$formula, moved, f$data, f$scales,
      "gaussian", "key", reference_scale = 250)
    p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
    q <- p; q$scale_beta <- p$scale_beta + log(s$unit_conv)
    expect_equal(.msr_factorized_group_mixture_limit(.msr_factorized_group(s, 2L), p, .4),
      .msr_factorized_group_mixture_limit(.msr_factorized_group(other, 2L), q, .4), tolerance = 1e-10)
    nb <- private_sq_fixture(kernel, "nbinom2")$spec
    poisson <- nb; poisson$family <- "poisson"; poisson$parameters$log_dispersion <- numeric()
    q <- private_sq_parameters(nb); q$scale_z <- q$response_z <- q$log_dispersion <- numeric()
    broad <- as.numeric(mixture_native_nll(poisson, q, "broad"))
    nearest <- as.numeric(mixture_native_nll(poisson, q, "nearest"))
    for (ratio in c(-.75, 0, .75)) {
      actual <- .msr_factorized_group_mixture_limit(.msr_factorized_group(poisson, 2L), q, ratio)
      expect_equal(actual$nll, mixture_native_total(broad[2L], nearest[2L], pnorm(ratio)), tolerance = 1e-10)
    }
  }
})

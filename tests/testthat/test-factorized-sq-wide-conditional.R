test_that("large-SD conditional caches retain native probabilities and uncertainty scope", {
  skip_if_not_installed("RTMB")
  maximum <- c(normalizer = 0, cdf = 0, quantile = 0, inverse = 0)
  cases <- 0L
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (unique_nearest in c(FALSE, TRUE))
          for (response_zero in if (response) c(FALSE, TRUE) else FALSE) {
            cases <- cases + 1L
            f <- private_sq_fixture(kernel, family, response = response); s <- f$spec
            if (unique_nearest) {
              distance <- f$inputs$d_list[[2L]]; nearest <- min(distance)
              distance[2L] <- (nearest + min(distance[distance > nearest])) / 2
              f$inputs$d_list[[2L]] <- distance
              s <- .msr_scale_prepare(f$formula, f$inputs, f$data, f$scales,
                family, "key", reference_scale = 250)
            }
            expect_identical(s$rows, f$spec$rows)
            expect_identical(s$observation_ids, f$spec$observation_ids)
            expect_identical(s$offset, f$spec$offset)
            if (family == "binomial") s$trials[] <- 4
            if (response_zero) s <- .msr_marginal_reduce(s, "response")
            before <- s; g <- .msr_factorized_group(s, 2L)
            p <- private_sq_parameters(g); p$beta[3] <- .12
            p$scale_beta[] <- 0; p$log_scale_sd[] <- log(8)
            p$scale_z <- p$response_z <- numeric()
            rng <- .Random.seed; configuration <- RTMB::TapeConfig()
            template <- p; template$scale_z <- 0
            initial <- unlist(template, use.names = FALSE)
            unpack <- local({ fixed <- template; function(x) private_sq_unpack(x, fixed) })
            tape <- .msr_factorized_group_log_tape(g, initial, unpack)
            cached <- function(z) {
              q <- p; q$scale_z <- z
              as.numeric(tape(unlist(q, use.names = FALSE)))
            }
            ordinary <- function(z) { q <- p; q$scale_z <- z; -.msr_factorized_group_joint_nll(g, q) }
            bound <- if (family == "gaussian")
              -length(g$y) * (p$log_dispersion + .5 * log(2 * pi)) else 0
            locations <- unique(drop(g$scales[[1L]]$X %*% p$scale_beta[g$scales[[1L]]$beta_index]))
            knots <- as.vector(outer(seq(-30, 30, by = 1), locations, "-") / exp(p$log_scale_sd))
            a <- .msr_factorized_distribution(ordinary, bound, knots)
            b <- .msr_factorized_distribution(cached, bound, knots)
            log_normalizer <- -wide_sq_integral(g, p)
            native_density <- function(z) vapply(z, function(t) {
              if (exp(bound + dnorm(t, log = TRUE) - log_normalizer) == 0) return(0)
              q <- p; q$scale_z <- t
              exp(-wide_sq_nll(g, q) + dnorm(t, log = TRUE) - log_normalizer)
            }, 0)
            native_cdf <- function(q, lower = TRUE) vapply(q, function(t) {
              integrate(native_density, if (lower) -Inf else t, if (lower) t else Inf,
                rel.tol = 1e-11, abs.tol = 0, subdivisions = 1000L, stop.on.error = TRUE)$value
            }, 0)
            nodes <- c(-3, 0, 3); probabilities <- c(.025, .5, .975)
            lower <- native_cdf(nodes); upper <- native_cdf(nodes, FALSE)
            quantiles <- b$quantile(probabilities)
            ordinary_quantiles <- a$quantile(probabilities)
            inverse <- native_cdf(quantiles)
            expect_equal(b$log_normalizer, a$log_normalizer, tolerance = 1e-8)
            expect_equal(b$log_normalizer, log_normalizer, tolerance = 1e-8)
            expect_equal(b$cdf(nodes), a$cdf(nodes), tolerance = 1e-9)
            expect_equal(b$cdf(nodes, FALSE), a$cdf(nodes, FALSE), tolerance = 1e-9)
            expect_equal(b$cdf(nodes), lower, tolerance = 1e-9)
            expect_equal(b$cdf(nodes, FALSE), upper, tolerance = 1e-9)
            expect_equal(quantiles, ordinary_quantiles, tolerance = 1e-7)
            expect_equal(inverse, probabilities, tolerance = 1e-7)
            expect_equal(b$cdf(c(-Inf, Inf)), c(0, 1))
            expect_equal(b$quantile(c(0, 1)), c(-Inf, Inf))
            expect_true(b$infinite_support)
            expect_lte(b$relative_error, 1e-10)
            expect_identical(RTMB::TapeConfig(), configuration)
            expect_identical(.Random.seed, rng)
            expect_identical(g$center, s$center)
            expect_identical(g$spread, s$spread)
            expect_identical(s, before)
            expect_error(.msr_factorized_conditional_cache(s), "Weighted SQ")
            no_annual <- .msr_marginal_reduce(g, "scale")
            expect_error(.msr_factorized_group_log_tape(no_annual, initial, unpack), "one scale grouping")
            maximum <- pmax(maximum, c(normalizer = max(abs(b$log_normalizer - a$log_normalizer),
              abs(b$log_normalizer - log_normalizer)),
              cdf = max(abs(b$cdf(nodes) - lower), abs(b$cdf(nodes, FALSE) - upper)),
              quantile = max(abs(quantiles - ordinary_quantiles)), inverse = max(abs(inverse - probabilities))))
          }
  expect_equal(cases, 24L)
  # Fixed-population latent intervals include no population SE or observation noise.
  cat(sprintf("\nSQ_WIDE_CONDITIONAL SD=8 cases=24 normalizer_error=%.12g cdf_error=%.12g quantile_error=%.12g inverse_error=%.12g\n",
    maximum[1L], maximum[2L], maximum[3L], maximum[4L]))
})

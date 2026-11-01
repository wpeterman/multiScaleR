test_that("large-SD signed-score moments retain full normalized nuisance identities", {
  skip_if_not_installed("RTMB")
  maximum <- c(mass = 0, value = 0, gradient = 0, hessian = 0)
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (unique in c(FALSE, TRUE)) {
          f <- private_sq_fixture(kernel, family, response = response)
          s <- f$spec
          if (unique) {
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
          before <- s; g <- .msr_factorized_group(s, 2L)
          p <- private_sq_parameters(g); p$beta[3] <- .12
          p$scale_beta[] <- 0; p$log_scale_sd[] <- log(8)
          p$scale_z <- p$response_z <- numeric()
          rng <- .Random.seed; configuration <- RTMB::TapeConfig()
          result <- wide_sq_score_moments(g, p)
          fn <- .msr_factorized_group_integral(g, TRUE)
          obj <- RTMB::MakeADFun(function(q) {
            part <- fn(q); -log(part$value) - part$shift
          }, p, silent = TRUE)
          native <- function(x) wide_sq_integral(g, private_sq_unpack(x, p))
          gradient <- private_sq_fd_gradient(native, obj$par)
          hessian <- private_sq_fd_hessian(native, obj$par)
          value <- as.numeric(obj$fn(obj$par))
          ag <- as.numeric(obj$gr(obj$par)); ah <- unname(obj$he(obj$par))
          expect_true(all(is.finite(c(result$mass, result$normalizer, result$score,
            result$hessian, value, gradient, hessian, ag, ah))))
          expect_equal(result$mass, 1, tolerance = 1e-8)
          expect_identical(names(obj$par), result$parameter_names)
          expect_equal(value, result$normalizer, tolerance = 1e-6)
          expect_equal(result$score, ag, tolerance = 2e-5)
          expect_equal(unname(result$hessian), ah, tolerance = 2e-5)
          expect_equal(result$score, gradient, tolerance = 2e-5)
          expect_equal(unname(result$hessian), hessian, tolerance = 2e-5)
          expect_identical(RTMB::TapeConfig(), configuration)
          expect_identical(.Random.seed, rng)
          expect_identical(g$center, s$center)
          expect_identical(g$spread, s$spread)
          expect_identical(s, before)
          expect_error(.msr_factorized_score_adfun(s), "Weighted SQ")
          maximum <- pmax(maximum, c(mass = abs(result$mass - 1),
            value = abs(value - result$normalizer),
            gradient = max(abs(result$score - ag), abs(result$score - gradient)),
            hessian = max(abs(result$hessian - ah), abs(result$hessian - hessian))))
        }
  cat(sprintf("\nSQ_WIDE_SCORE SD=8 cases=20 mass_error=%.12g value_error=%.12g gradient_error=%.12g hessian_error=%.12g\n",
    maximum[1L], maximum[2L], maximum[3L], maximum[4L]))
})

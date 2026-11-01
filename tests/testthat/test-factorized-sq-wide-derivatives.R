test_that("large-SD AD coordinates retain full native nuisance gradients and Hessians", {
  skip_if_not_installed("RTMB")
  maximum <- c(value = 0, gradient = 0, hessian = 0)
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        s <- private_sq_fixture(kernel, family, response = response)$spec
        if (family == "binomial") s$trials[] <- 4
        before <- s; g <- .msr_factorized_group(s, 2L)
        p <- private_sq_parameters(g); p$beta[3] <- .12
        p$scale_beta[] <- 0; p$log_scale_sd[] <- log(8)
        p$scale_z <- p$response_z <- numeric()
        rng <- .Random.seed; configuration <- RTMB::TapeConfig()
        physical <- .msr_factorized_group_scale_integral(g, TRUE)
        normal <- .msr_factorized_group_integral(g, TRUE)
        make <- function(fn) RTMB::MakeADFun(function(q) {
          part <- fn(q); -log(part$value) - part$shift
        }, p, silent = TRUE)
        a <- make(physical); b <- make(normal)
        native <- function(x) wide_sq_integral(g, private_sq_unpack(x, p))
        value <- native(a$par)
        gradient <- private_sq_fd_gradient(native, a$par)
        hessian <- private_sq_fd_hessian(native, a$par)
        av <- as.numeric(a$fn(a$par)); bv <- as.numeric(b$fn(b$par))
        ag <- as.numeric(a$gr(a$par)); bg <- as.numeric(b$gr(b$par))
        ah <- unname(a$he(a$par)); bh <- unname(b$he(b$par))
        expect_true(all(is.finite(c(value, gradient, hessian, av, bv, ag, bg, ah, bh))))
        expect_equal(av, value, tolerance = 1e-6)
        expect_equal(bv, value, tolerance = 1e-6)
        expect_equal(ag, gradient, tolerance = 2e-5)
        expect_equal(bg, gradient, tolerance = 2e-5)
        expect_equal(ah, hessian, tolerance = 2e-5)
        expect_equal(bh, hessian, tolerance = 2e-5)
        expect_equal(ag, bg, tolerance = 2e-5)
        expect_equal(ah, bh, tolerance = 2e-5)
        expect_identical(RTMB::TapeConfig(), configuration)
        expect_identical(.Random.seed, rng)
        expect_identical(g$center, s$center)
        expect_identical(g$spread, s$spread)
        expect_identical(s, before)
        expect_error(.msr_factorized_scale_objective(s), "Weighted SQ")
        maximum <- pmax(maximum, c(value = max(abs(c(av, bv) - value)),
          gradient = max(abs(ag - gradient), abs(bg - gradient), abs(ag - bg)),
          hessian = max(abs(ah - hessian), abs(bh - hessian), abs(ah - bh))))
      }
  cat(sprintf("\nSQ_WIDE_AD SD=8 cases=10 value_error=%.12g gradient_error=%.12g hessian_error=%.12g\n",
    maximum[1L], maximum[2L], maximum[3L]))
})

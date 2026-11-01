test_that("independent log weights retain finite references and exact SQ endpoints", {
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        s <- private_sq_fixture(kernel, family, response = response)$spec
        before <- s; p <- private_sq_parameters(s); p$beta[3] <- .12
        p$scale_z[] <- 0
        buffers <- s$surface$buffers[[1L]]
        for (t in c(-8, -2, 0, 2, 8)) {
          p$scale_beta[] <- t
          native <- wide_sq_reference(s, p)
          old <- private_sq_reference(s, p, response = FALSE)
          expect_equal(native$raw, old$raw, tolerance = 1e-12)
          expect_equal(native$kernel, old$kernel, tolerance = 1e-12)
          expect_equal(native$eta, old$eta, tolerance = 1e-12)
          expect_equal(native$raw, .msr_scale_surface_raw(buffers, rep(t, length(s$y)), kernel), tolerance = 1e-12)
        }
        for (t in c(-1000, 1000, -Inf, Inf)) {
          expected <- vapply(buffers, function(b) {
            v <- if (t < 0) b$v[b$d == min(b$d)] else b$v
            sqrt(mean((v - mean(v))^2))
          }, 0)
          actual <- wide_sq_raw(buffers, rep(t, length(s$y)), kernel)
          expect_true(all(is.finite(actual)))
          expect_equal(actual, unname(expected), tolerance = 1e-12)
          # Infinite log scales are mathematical endpoints, not valid fitted
          # parameters. Production arithmetic is compared at finite extremes.
          if (is.finite(t)) expect_equal(actual,
            .msr_scale_surface_raw(buffers, rep(t, length(s$y)), kernel), tolerance = 1e-12)
        }
        expect_identical(s, before)
      }
})

test_that("large annual SD preserves independent normalized coordinate values", {
  errors <- setNames(rep(0, 3), c("2", "8", "32"))
  ranges <- matrix(c(rep(Inf, 3), rep(-Inf, 3)), ncol = 2L)
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        s <- private_sq_fixture(kernel, family, response = response)$spec
        if (family == "binomial") s$trials[] <- 4
        before <- s; g <- .msr_factorized_group(s, 2L)
        p <- private_sq_parameters(g); p$beta[3] <- .12
        p$scale_z <- p$response_z <- numeric()
        physical <- .msr_factorized_group_scale_integral(g)
        normal <- .msr_factorized_group_integral(g)
        endpoint <- function(t) {
          q <- p; q$scale_beta[] <- t; q$scale_z <- 0
          wide_sq_nll(g, q)
        }
        broad <- endpoint(Inf); nearest <- endpoint(-Inf)
        for (j in 1:3) for (b in c(-.75, 0, .75)) {
          sd <- c(2, 8, 32)[j]
          p$log_scale_sd[] <- log(sd); p$scale_beta[] <- b * sd
          native <- wide_sq_integral(g, p)
          pv <- physical(p); pv <- -log(pv$value) - pv$shift
          nv <- normal(p); nv <- -log(nv$value) - nv$shift
          expect_true(all(is.finite(c(native, pv, nv))))
          expect_equal(pv, native, tolerance = 1e-6)
          expect_equal(nv, native, tolerance = 1e-6)
          weight <- pnorm(b)
          limit <- .msr_factorized_mixture_nll(broad, nearest, weight)
          errors[j] <- max(errors[j], abs(pv - native), abs(nv - native))
          ranges[j, ] <- c(min(ranges[j, 1L], native - limit), max(ranges[j, 2L], native - limit))
        }
        expect_identical(g$center, s$center)
        expect_identical(g$spread, s$spread)
        expect_identical(s, before)
        expect_error(.msr_factorized_scale_objective(s), "Weighted SQ")
      }
  # Compact diagnostics only; finite SD is not a limit or identification gate.
  for (j in 1:3) cat(sprintf("\nSQ_WIDE SD=%s coordinate_error=%.12g finite_minus_limit=[%.12g,%.12g]\n",
    names(errors)[j], errors[j], ranges[j, 1L], ranges[j, 2L]))
})

test_that("large-SD paths give the correct normalized group-mixture endpoints", {
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        s <- private_sq_fixture(kernel, family, response = response)$spec
        if (family == "binomial") s$trials[] <- 4
        g <- .msr_factorized_group(s, 2L); p <- private_sq_parameters(g)
        p$beta[3] <- .12; p$log_scale_sd[] <- log(1e6); p$response_z <- numeric()
        for (b in c(-.75, 0, .75)) {
          p$scale_beta[] <- b * 1e6
          nll <- numeric(2L)
          for (i in 1:2) {
            q <- p; q$scale_z <- -b + c(-.25, .25)[i]
            actual <- wide_sq_reference(g, q)
            expected <- vapply(g$surface$buffers[[1L]], function(cell) {
              v <- if (i == 1L) cell$v[cell$d == min(cell$d)] else cell$v
              sqrt(mean((v - mean(v))^2))
            }, 0)
            expect_equal(actual$raw, unname(expected), tolerance = 1e-12)
            nll[i] <- wide_sq_nll(g, q)
            expect_equal(.msr_factorized_group_joint_nll(g, q),
              nll[i] - dnorm(q$scale_z, log = TRUE), tolerance = 1e-8)
          }
          mass <- integrate(dnorm, -b, Inf, rel.tol = 1e-12)$value
          expect_equal(mass, pnorm(b), tolerance = 1e-12)
          native <- -log(sum(c(mass, 1 - mass) * exp(-rev(nll))))
          expect_equal(.msr_factorized_mixture_nll(nll[2L], nll[1L], pnorm(b)),
            native, tolerance = 1e-12)
        }
        expect_error(.msr_factorized_conditional_cache(s), "Weighted SQ")
      }
})

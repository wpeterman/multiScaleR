# Diagnostic arithmetic only. Full-model raw-SQ constructors stay guarded.
zero_sq_adfun <- function(g, p) {
  .msr_scale_surface_ram(g)
  if (length(g$scales[[1L]]$random)) {
    fn <- .msr_factorized_group_integral(g, TRUE)
    return(RTMB::MakeADFun(function(q) {
      part <- fn(q); -log(part$value) - part$shift
    }, p, silent = TRUE))
  }
  RTMB::MakeADFun(function(q) sum(.msr_factorized_response_nll(g, q, TRUE)), p, silent = TRUE)
}

zero_sq_native <- function(g, p) {
  if (length(g$scales[[1L]]$random)) wide_sq_integral(g, p) else wide_sq_nll(g, p)
}

test_that("finite-bandwidth SQ zero-variance faces retain normalized targets and population counts", {
  skip_if_not_installed("RTMB")
  population_count <- function(p) sum(lengths(p[!names(p) %in% c("scale_z", "response_z")]))
  maximum <- c(value = 0, gradient = 0, hessian = 0, near_value = 0,
    near_gradient = 0, near_hessian = 0)
  faces <- 0L
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (unique in c(FALSE, TRUE)) {
          f <- private_sq_fixture(kernel, family, response = response); s <- f$spec
          if (unique) {
            distance <- f$inputs$d_list[[2L]]; nearest <- min(distance)
            distance[2L] <- (nearest + min(distance[distance > nearest])) / 2
            f$inputs$d_list[[2L]] <- distance
            s <- .msr_scale_prepare(f$formula, f$inputs, f$data, f$scales,
              family, "key", reference_scale = 250)
          }
          if (family == "binomial") s$trials[] <- 4
          before <- s; g <- .msr_factorized_group(s, 2L)
          selections <- if (response) list("scale", "response", c("scale", "response")) else list("scale")
          for (face in selections) {
            faces <- faces + 1L
            p <- private_sq_parameters(g); p$beta[3] <- .12
            p$scale_beta[] <- 0; p$log_scale_sd[] <- log(8)
            p$scale_z <- p$response_z <- numeric()
            if ("scale" %in% face) p$log_scale_sd[] <- log(1e-7)
            if ("response" %in% face) p$log_response_sd[] <- log(1e-7)
            reduced <- .msr_marginal_reduce(g, face); q <- p
            for (block in face) {
              q[[paste0("log_", block, "_sd")]] <- numeric()
              q[[paste0(block, "_z")]] <- numeric()
              expect_true(isTRUE(reduced$zero_sd[[block]]))
              expect_length(reduced$parameters[[paste0("log_", block, "_sd")]], 0)
              expect_length(reduced$parameters[[paste0(block, "_z")]], 0)
            }
            expect_equal(population_count(q), population_count(p) - length(face))
            expect_equal(population_count(reduced$parameters), population_count(q))
            expect_identical(q$scale_beta, p$scale_beta)
            expect_identical(reduced$parameters$scale_beta, g$parameters$scale_beta)
            expect_identical(.msr_marginal_reduce(reduced, face), reduced)
            for (nm in c("rows", "observation_ids", "X_fixed", "y", "offset", "trials", "center", "spread", "surface"))
              expect_identical(reduced[[nm]], g[[nm]])
            if (!"response" %in% face) expect_identical(reduced$response, g$response)
            if (!"scale" %in% face) expect_identical(reduced$scales, g$scales)
            expect_null(reduced$priors)
            rng <- .Random.seed; configuration <- RTMB::TapeConfig()
            exact <- zero_sq_adfun(reduced, q); near <- zero_sq_adfun(g, p)
            expect_length(exact$par, population_count(q))
            expect_length(near$par, population_count(p))
            native <- function(x) zero_sq_native(reduced, private_sq_unpack(x, q))
            value <- native(exact$par)
            gradient <- private_sq_fd_gradient(native, exact$par)
            hessian <- private_sq_fd_hessian(native, exact$par)
            ev <- as.numeric(exact$fn(exact$par)); eg <- as.numeric(exact$gr(exact$par))
            eh <- unname(exact$he(exact$par))
            remove <- paste0("log_", face, "_sd")
            keep <- which(!names(near$par) %in% remove)
            nv <- as.numeric(near$fn(near$par)); ng <- as.numeric(near$gr(near$par))[keep]
            nh <- unname(near$he(near$par)[keep, keep, drop = FALSE])
            expect_identical(names(near$par)[keep], names(exact$par))
            expect_true(all(is.finite(c(value, gradient, hessian, ev, eg, eh, nv, ng, nh))))
            expect_equal(ev, value, tolerance = 1e-6)
            expect_equal(eg, gradient, tolerance = 2e-5)
            expect_equal(eh, hessian, tolerance = 2e-5)
            expect_lt(abs(nv - ev), 1e-7)
            expect_equal(ng, eg, tolerance = 2e-5)
            expect_equal(nh, eh, tolerance = 2e-5)
            expect_identical(RTMB::TapeConfig(), configuration)
            expect_identical(.Random.seed, rng)
            expect_error(.msr_factorized_adfun(reduced, q), "Weighted SQ")
            maximum <- pmax(maximum, c(value = abs(ev - value),
              gradient = max(abs(eg - gradient)), hessian = max(abs(eh - hessian)),
              near_value = abs(nv - ev), near_gradient = max(abs(ng - eg)), near_hessian = max(abs(nh - eh))))
          }
          if (!response) expect_error(.msr_marginal_reduce(g, "response"), "absent SD")
          expect_identical(s, before)
        }
  expect_equal(faces, 28L)
  cat(sprintf("\nSQ_ZERO faces=28 value_error=%.12g gradient_error=%.12g hessian_error=%.12g near_value_error=%.12g near_gradient_error=%.12g near_hessian_error=%.12g\n",
    maximum[1L], maximum[2L], maximum[3L], maximum[4L], maximum[5L], maximum[6L]))
})

test_that("ordinary density components match native full integrals in both coordinates", {
  preparations <- 0L; points <- 0L; maximum <- 0
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (response_zero in if (response) c(FALSE, TRUE) else FALSE) {
          preparations <- preparations + 1L
          design <- if (response) "cell_means" else if (family %in% c("poisson", "binomial")) "contrast" else "canonical"
          s <- path_component_fixture(kernel, family, response,
            if (kernel == "gaussian") "tied" else "unique", design)
          if (family == "binomial") s$trials[] <- 4
          if (response_zero) s <- .msr_marginal_reduce(s, "response")
          p <- private_sq_parameters(s); p$beta[3L] <- .12
          p$scale_beta[] <- seq_along(p$scale_beta) * .1
          p$scale_z <- p$response_z <- numeric()
          before <- list(s, p, .Random.seed)
          standard <- .msr_factorized_density_component(s)
          physical <- .msr_factorized_density_component(s, "physical")
          metadata <- attr(standard, "component")
          expect_identical(metadata[c("target", "coordinate", "rel.tol", "step", "finite_scale_estimate", "uncertainty_available")],
            list(target = "normalized_full_marginal_likelihood", coordinate = "standard", rel.tol = 1e-8,
              step = 1, finite_scale_estimate = FALSE, uncertainty_available = FALSE))
          for (sd in c(.3, 8)) {
            points <- points + 1L
            q <- p; q$log_scale_sd <- log(sd)
            reference <- sum(vapply(1:3, function(i)
              wide_sq_integral(.msr_factorized_group(s, i), q), 0))
            actual <- c(standard(q), physical(q))
            maximum <- max(maximum, abs(actual - reference))
            expect_equal(actual, rep(reference, 2L), tolerance = 1e-6,
              info = paste(kernel, family, response, response_zero, design, sd))
          }
          fixed <- .msr_marginal_reduce(s, "scale")
          zero <- p; zero$log_scale_sd <- numeric()
          reference <- wide_sq_nll(fixed, zero)
          expect_equal(c(.msr_factorized_density_component(fixed)(zero),
            .msr_factorized_density_component(fixed, "physical")(zero)), rep(reference, 2L), tolerance = 1e-6)
          expect_identical(list(s, p, .Random.seed), before)
        }
  expect_equal(c(preparations, points), c(12L, 24L))
  cat(sprintf("\nORDINARY_DENSITY_COMPONENT preparations=12 points=24 native_error=%.12g\n", maximum))
})

test_that("mean wrapper parity and strict population layouts preserve evaluation scope", {
  s <- mixture_fixture("gaussian", "gaussian", FALSE, TRUE)$spec
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  ordinary <- .msr_factorized_density_component(s, rel.tol = 1e-10, step = .5)
  expect_equal(ordinary(p), .msr_factorized_objective(s, FALSE, 1e-10, .5)(p), tolerance = 1e-10)
  expect_equal(.msr_factorized_density_component(s, "physical")(p),
    .msr_factorized_scale_objective(s, FALSE)(p), tolerance = 1e-10)
  bad <- p; bad$scale_z <- 0
  expect_error(ordinary(bad), "without realized latent")
  bad <- p; bad$response_z <- 0
  expect_error(ordinary(bad), "without realized latent")
  bad <- p; bad$beta[1L] <- Inf
  expect_error(ordinary(bad), "finite numeric population")
  bad <- p; bad$beta <- matrix(bad$beta, ncol = 1L)
  expect_error(ordinary(bad), "finite numeric population")
  bad <- p; bad$beta <- as.complex(bad$beta)
  expect_error(ordinary(bad), "finite numeric population")
  bad <- s; bad$priors <- list(beta = 1)
  expect_error(.msr_factorized_density_component(bad), "neutral")
  bad <- s; bad$response$random[[1L]]$index[] <- 1L
  expect_error(.msr_factorized_density_component(bad), "coupling")
  expect_error(.msr_factorized_density_component(s, "invalid"), "arg")
  expect_error(.msr_factorized_density_component(s, step = 0), "> 0")
  expect_error(.msr_factorized_density_component(s, rel.tol = NA_real_), "finite")
})

test_that("ordinary raw evaluations reach neither tape nor optimizer boundaries", {
  skip_if_not_installed("RTMB")
  s <- .msr_marginal_reduce(private_sq_fixture()$spec, "scale")
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  launched <- 0L
  testthat::local_mocked_bindings(MakeADFun = function(...) {
    launched <<- launched + 1L; stop("Tape reached")
  }, MakeTape = function(...) {
    launched <<- launched + 1L; stop("Tape reached")
  }, .package = "RTMB")
  testthat::local_mocked_bindings(nlminb = function(...) {
    launched <<- launched + 1L; stop("Optimizer reached")
  }, .package = "stats")
  expect_true(is.finite(.msr_factorized_density_component(s)(p)))
  expect_error(.msr_factorized_objective(s, FALSE), "Weighted SQ")
  expect_error(.msr_factorized_fit(s), "Weighted SQ")
  expect_error(.msr_factorized_adfun(s, p), "Weighted SQ")
  expect_identical(launched, 0L)
})

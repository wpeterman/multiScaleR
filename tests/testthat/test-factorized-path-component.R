test_that("explicit raw exact-limit components retain normalized native whole-model densities", {
  cases <- 0L; maximum <- 0
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (response_zero in if (response) c(FALSE, TRUE) else FALSE)
          for (support in c("tied", "unique"))
            for (design in c("canonical", "contrast", "cell_means")) {
              cases <- cases + 1L
              s <- path_component_fixture(kernel, family, response, support, design)
              if (family == "binomial") s$trials[] <- 4
              if (response_zero) s <- .msr_marginal_reduce(s, "response")
              p <- private_sq_parameters(s); p$beta[3L] <- .12
              p$scale_beta[] <- if (design == "canonical") .2 else c(.2, -.3)
              p$log_scale_sd[] <- log(8); p$scale_z <- p$response_z <- numeric()
              before <- list(s, p, .Random.seed)
              component <- .msr_factorized_path_component(s, p)
              expected <- numeric()
              for (poisson in if (family == "nbinom2") c(FALSE, TRUE) else FALSE) {
                g <- s; q <- p
                if (poisson) {g$family <- "poisson"; q$log_dispersion <- numeric()}
                broad <- mixture_native_nll(g, q, "broad")
                nearest <- mixture_native_nll(g, q, "nearest")
                expected <- c(expected, vapply(pnorm(c(-.75, 0, .75)), function(w)
                  mixture_native_total(broad, nearest, w), 0))
              }
              actual <- unname(vapply(component$references, `[[`, 0, "objective"))
              maximum <- max(maximum, abs(actual - expected))
              label <- paste(kernel, family, response, response_zero, support, design)
              expect_equal(actual, unname(expected), tolerance = 1e-10, info = label)
              zero <- p; zero$log_scale_sd <- numeric()
              expect_identical(component$references, .msr_factorized_path_component(s, zero)$references)
              expect_true(all(vapply(component$references, function(x) x$completed && !x$finite_scale_estimate, TRUE)))
              expect_identical(list(s, p, .Random.seed), before)
            }
  expect_equal(cases, 72L)
  cat(sprintf("\nEXACT_PATH_COMPONENT cases=72 native_error=%.12g\n", maximum))
})

test_that("component structure preserves dependence, neutral target and full observation provenance", {
  s <- path_component_fixture("gaussian", "gaussian", TRUE, "tied", "canonical")
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  expect_identical(.msr_factorized_limit_component_validate(s), s)
  bad <- s; bad$priors <- list(beta = 1)
  expect_error(.msr_factorized_path_component(bad, p), "neutral")
  bad <- s; bad$response$random[[1L]]$index[] <- 1L
  expect_error(.msr_factorized_path_component(bad, p), "coupling")
  bad <- s; bad$scales[[1L]]$random[[1L]]$gaussian_block <- 1L
  expect_error(.msr_factorized_path_component(bad, p), "Gaussian projection")
  bad <- s; bad$covariates <- c(bad$covariates, "other")
  expect_error(.msr_factorized_path_component(bad, p), "one kernel")
  bad <- s; bad$family <- "poisson"
  expect_error(.msr_factorized_path_component(bad, p), "Gaussian response")
  bad <- .msr_marginal_reduce(s, "response"); bad$parameters$response_z <- 0
  expect_error(.msr_factorized_path_component(bad, p), "no corresponding grouping term")
  bad <- s; bad$rows <- bad$rows[-1L]
  expect_error(.msr_factorized_path_component(bad, p), "observation metadata")
  bad <- s; bad$observation_ids[1L] <- NA_character_
  expect_error(.msr_factorized_path_component(bad, p), "observation metadata")
  bad <- s; bad$data <- bad$data[-1L, , drop = FALSE]
  expect_error(.msr_factorized_path_component(bad, p), "observation metadata")
  fixed <- .msr_marginal_reduce(s, "scale")
  expect_false(.msr_factorized_path_component(fixed, p)$applicable)
  mean <- mixture_fixture("gaussian", "gaussian", FALSE, TRUE)$spec
  q <- private_sq_parameters(mean); q$scale_z <- q$response_z <- numeric()
  expect_identical(.msr_factorized_path_component(mean, q), .msr_factorized_path_limits(mean, q))
})

test_that("exact components cannot admit raw fitting or reach tape constructors", {
  skip_if_not_installed("RTMB")
  s <- path_component_fixture("gaussian", "gaussian", FALSE, "tied", "canonical")
  p <- private_sq_parameters(s); p$scale_z <- p$response_z <- numeric()
  launched <- 0L
  testthat::local_mocked_bindings(
    MakeADFun = function(...) { launched <<- launched + 1L; stop("Tape reached") },
    MakeTape = function(...) { launched <<- launched + 1L; stop("Tape reached") }, .package = "RTMB")
  expect_true(.msr_factorized_path_component(s, p)$applicable)
  constructors <- list(function() .msr_factorized_fit(s), function() .msr_marginal_fit(s),
    function() .msr_factorized_objective(s, FALSE), function() .msr_factorized_adfun(s),
    function() .msr_factorized_score_adfun(s), function() .msr_factorized_conditional_cache(s),
    function() .msr_factorized_limit_audit(s, p, 0), function() .msr_factorized_path_limits(s, p))
  for (construct in constructors) expect_error(construct(), "Weighted SQ")
  expect_identical(launched, 0L)
  group_limit <- .msr_factorized_group_mixture_limit
  testthat::local_mocked_bindings(.msr_factorized_group_mixture_limit = function(group, parameters, ratio) {
    if (ratio == 0) stop("mock component path failure")
    group_limit(group, parameters, ratio)
  }, .package = "multiScaleR")
  partial <- .msr_factorized_path_component(s, p)$references
  expect_match(partial$`annual-sd-infinite-gaussian-zero`$error, "mock component path failure")
  expect_true(partial$`annual-sd-infinite-gaussian-negative`$completed)
})

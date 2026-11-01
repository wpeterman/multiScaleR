test_that("exact kernel components preserve native targets and existing constructor dimensions", {
  skip_if_not_installed("RTMB")
  preparations <- 0L; endpoints <- 0L; maximum <- 0
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (response_zero in if (response) c(FALSE, TRUE) else FALSE)
          for (support in c("tied", "unique"))
            for (design in c("canonical", "contrast", "cell_means")) {
              preparations <- preparations + 1L
              s <- path_component_fixture(kernel, family, response, support, design)
              if (family == "binomial") s$trials[] <- 4
              if (response_zero) s <- .msr_marginal_reduce(s, "response")
              p <- private_sq_parameters(s); p$beta[3L] <- .12
              p$scale_z <- p$response_z <- numeric()
              before <- list(s, p); pure <- TRUE
              for (limit in c("broad", "nearest"))
                for (poisson in if (family == "nbinom2") c(FALSE, TRUE) else FALSE) {
                  endpoints <- endpoints + 1L
                  rng <- .Random.seed
                  endpoint <- .msr_factorized_kernel_component(s, limit, poisson)
                  pure <- pure && identical(.Random.seed, rng)
                  q <- p; q$scale_beta <- q$scale_z <- q$log_scale_sd <- numeric()
                  native_spec <- s
                  if (poisson) {native_spec$family <- "poisson"; q$log_dispersion <- numeric()}
                  expected <- sum(mixture_native_nll(native_spec, q, limit))
                  ordinary <- .msr_factorized_objective(endpoint, FALSE)(q)
                  obj <- .msr_factorized_adfun(endpoint, q)
                  value <- as.numeric(obj$fn(obj$par))
                  maximum <- max(maximum, abs(c(ordinary, value) - expected))
                  label <- paste(kernel, family, response, response_zero, support, design, limit, poisson)
                  expect_equal(c(ordinary, value), rep(expected, 2L), tolerance = 1e-9, info = label)
                  expect_true(all(is.finite(c(obj$gr(obj$par), obj$he(obj$par)))), info = label)
                  expected_parameters <- sum(lengths(endpoint$parameters)[
                    !names(endpoint$parameters) %in% c("scale_z", "response_z")])
                  fields <- c("center", "spread", "offset", "trials", "rows", "observation_ids", "data", "response")
                  # One structured endpoint contract retains all metadata fields
                  # and useful failure paths without multiplying their counts.
                  metadata <- list(population_dimension = length(obj$par),
                    scale_blocks = lengths(endpoint$parameters[c("scale_beta", "scale_z", "log_scale_sd")]),
                    retained = endpoint[fields], surface_specs = endpoint$surface$specs,
                    operation = endpoint$exact_limit_component,
                    raw_buffers = length(endpoint$surface$buffers), raw_RAM_estimate = .msr_scale_surface_ram(endpoint))
                  expected_metadata <- list(population_dimension = as.integer(expected_parameters),
                    scale_blocks = c(scale_beta = 0L, scale_z = 0L, log_scale_sd = 0L),
                    retained = s[fields], surface_specs = s$surface$specs,
                    operation = list(operation = "exact_kernel_limit", target = "marginal_likelihood_limit",
                      source_family = s$family, kernel = limit, poisson = poisson, finite_scale_estimate = FALSE),
                    raw_buffers = 0L, raw_RAM_estimate = NULL)
                  expect_identical(metadata, expected_metadata, info = label)
                  rm(obj)
                }
              expect_identical(list(s, p), before)
              expect_true(pure)
            }
  expect_equal(c(preparations, endpoints), c(72L, 168L))
  cat(sprintf("\nEXACT_KERNEL_COMPONENT preparations=72 endpoints=168 native_error=%.12g\n", maximum))
})

test_that("kernel component validation retains original-model guards and explicit scope", {
  s <- path_component_fixture("gaussian", "gaussian", TRUE, "tied", "contrast")
  expect_error(.msr_factorized_kernel_component(s, "invalid"), "kernel limit")
  expect_error(.msr_factorized_kernel_component(s, "broad", TRUE), "Poisson limit requires NB2")
  bad <- s; bad$priors <- list(beta = 1)
  expect_error(.msr_factorized_kernel_component(bad, "broad"), "neutral")
  bad <- s; bad$response$random[[1L]]$index[] <- 1L
  expect_error(.msr_factorized_kernel_component(bad, "nearest"), "coupling")
  bad <- s; bad$rows <- bad$rows[-1L]
  expect_error(.msr_factorized_kernel_component(bad, "broad"), "observation metadata")
  fixed <- .msr_marginal_reduce(s, "scale")
  expect_length(.msr_factorized_kernel_component(fixed, "broad")$parameters$log_scale_sd, 0L)
  expect_error(.msr_factorized_limit_spec(s, "broad"), "Weighted SQ")
  expect_error(.msr_factorized_fit(s), "Weighted SQ")
})

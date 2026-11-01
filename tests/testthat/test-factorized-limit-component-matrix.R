test_that("exact audit components preserve every supported family and zero-face model universe", {
  # Mock optimization/tapes only: retain actual normalized point/path values.
  # The deliberate nonconvergence/singular curvature is not a fitted result.
  constructed <- list()
  testthat::local_mocked_bindings(.msr_factorized_adfun = function(spec, parameters, ...) {
    constructed[[length(constructed) + 1L]] <<- c(raw = length(spec$surface$buffers),
      scale_beta = length(spec$parameters$scale_beta), annual_sd = length(spec$parameters$log_scale_sd))
    par <- unlist(parameters, use.names = FALSE)
    env <- new.env(); env$last.par <- par
    env$parList <- function(x, par) private_sq_unpack(x, parameters)
    ordinary <- .msr_factorized_density_component(spec)
    list(par = par, env = env, fn = function(x) ordinary(private_sq_unpack(x, parameters)),
      gr = function(x) rep(0, length(x)), he = function(x) matrix(0, length(x), length(x)))
  }, .msr_marginal_polish = function(...) list(passed = FALSE), .package = "multiScaleR")
  testthat::local_mocked_bindings(nlminb = function(start, objective, gradient, control) {
    stopifnot(identical(control, list(iter.max = 250L, eval.max = 750L, rel.tol = 1e-12)))
    list(par = start, objective = objective(start), convergence = 1L,
      message = "deliberate compatibility mock")
  }, .package = "stats")
  actual <- native <- numeric(); checks <- list(); preparations <- references <- 0L
  for (kernel in c("gaussian", "exp"))
    for (role in c("gaussian", "gaussian_response", "gaussian_response_zero", "poisson", "binomial", "nbinom2")) {
      family <- if (grepl("^gaussian", role)) "gaussian" else role
      response <- grepl("response", role)
      original <- path_component_fixture(kernel, family, response,
        if (kernel == "gaussian") "tied" else "unique",
        if (response) "cell_means" else if (family %in% c("binomial", "poisson")) "contrast" else "canonical")
      if (family == "binomial") original$trials <- rep(c(3, 5, 4), length.out = length(original$y))
      preparations <- preparations + 1L
      for (mode in c("positive", "selected_zero", "fixed")) {
        s <- if (mode == "fixed") .msr_marginal_reduce(original, "scale") else original
        p <- private_sq_parameters(s); p$beta[3L] <- .12
        p$scale_beta[] <- seq_along(p$scale_beta) * .1
        p$scale_z <- p$response_z <- numeric()
        p$log_scale_sd <- if (mode == "positive") log(9) else numeric()
        if (role == "gaussian_response_zero") p$log_response_sd <- numeric()
        before <- list(s, p, .Random.seed)
        audit <- .msr_factorized_limit_component_audit(s, p, 1e5)
        key <- paste(kernel, role, mode, sep = "/")
        oracle <- if (role == "gaussian_response_zero")
          .msr_marginal_reduce(original, "response") else original
        labels <- character(); attempts <- list(); dimensions <- logical()
        for (endpoint in c("broad", "nearest"))
          for (poisson in if (family == "nbinom2") c(FALSE, TRUE) else FALSE) {
            label <- paste(endpoint, if (poisson) "poisson" else family, sep = "-")
            labels <- c(labels, label); references <- references + 1L
            a <- audit$references[[label]]; attempts <- c(attempts, a$attempts)
            g <- oracle
            if (poisson) g$family <- "poisson"
            expected_size <- 3L + as.integer(g$family %in% c("gaussian", "nbinom2")) +
              as.integer(role == "gaussian_response")
            dimensions <- c(dimensions, length(a$attempts) == 3L,
              vapply(a$attempts, function(x) length(unlist(x$parameters, use.names = FALSE)) == expected_size &&
                all(lengths(x$parameters)[c("scale_beta", "scale_z", "log_scale_sd", "response_z")] == 0L), TRUE))
            for (j in seq_along(a$attempts)) {
              n <- paste(key, label, j, sep = "/")
              actual[n] <- a$attempts[[j]]$objective
              native[n] <- sum(mixture_native_nll(g, a$attempts[[j]]$parameters, endpoint))
            }
          }
        path_names <- grep("^annual-sd-infinite-", names(audit$references), value = TRUE)
        expected_path <- if (mode == "fixed") 0L else if (family == "nbinom2") 6L else 3L
        checks[[key]] <- c(completed = isTRUE(audit$completed),
          source = identical(list(s, p, .Random.seed), before), dimensions = all(dimensions),
          diagnostics = all(vapply(attempts, function(x) isTRUE(x$audited) &&
            identical(x$optimizer$convergence, 1L) && !x$curvature$passed && !x$refinement$passed, TRUE)),
          model_universe = identical(audit$infinite_sd_path$applicable, mode != "fixed") &&
            (mode == "fixed" || identical(audit$infinite_sd_path$candidate_zero_sd, mode == "selected_zero")),
          path_count = length(path_names) == expected_path,
          poisson = family != "nbinom2" || isTRUE(audit$references$`fitted-scales-poisson`$completed),
          scope = isTRUE(audit$component_only) && identical(audit$target, "marginal_likelihood_limit") &&
            identical(audit$finite_scale_estimate, FALSE) && identical(audit$uncertainty_available, FALSE),
          witnesses = all(labels %in% audit$witnesses))
      }
    }
  # Aggregate related contracts without multiplying identical metadata checks.
  contract <- do.call(rbind, checks)
  for (name in colnames(contract))
    expect_true(all(contract[, name]), info = paste(name, paste(rownames(contract)[!contract[, name]], collapse = ", ")))
  expect_identical(c(preparations, length(checks), references, length(actual)), c(12L, 36L, 84L, 252L))
  expect_identical(names(actual), names(native))
  expect_lt(max(abs(actual - native)), 1e-6)
  expect_identical(length(constructed), 84L)
  expect_true(all(vapply(constructed, function(x) all(x == 0L), TRUE)))
  cat(sprintf("\nLIMIT_AUDIT_COMPATIBILITY preparations=12 cases=36 endpoints=84 mock_attempts=252 native_error=%.12g\n",
    max(abs(actual - native))))
})

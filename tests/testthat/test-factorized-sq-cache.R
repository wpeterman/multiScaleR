sq_numeric_group <- function(group, population) {
  template <- population; template$scale_z <- 0; template$response_z <- numeric()
  force(template)
  unpack <- function(x) private_sq_unpack(x, template)
  initial <- unlist(template, use.names = FALSE)
  list(tape = .msr_factorized_group_log_tape(group, initial, unpack),
    initial = initial, unpack = unpack, template = template)
}

sq_numeric_log_density <- function(record, population) {
  force(record); force(population)
  function(z) {
    p <- population; p$scale_z <- z
    as.numeric(record$tape(unlist(p, use.names = FALSE)))
  }
}

test_that("numeric group tapes preserve native joint targets across groups and parameter changes", {
  skip_if_not_installed("RTMB")
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE)
        for (sd in c(.3, 2)) {
          s <- private_sq_fixture(kernel, family, response = response)$spec
          if (family == "binomial") s$trials[] <- 4
          before <- s
          for (i in 1:3) {
            g <- .msr_factorized_group(s, i); p <- private_sq_parameters(g)
            p$beta[3] <- .12; p$log_scale_sd[] <- log(sd)
            p$scale_z <- p$response_z <- numeric()
            rng <- .Random.seed; configuration <- RTMB::TapeConfig()
            record <- sq_numeric_group(g, p)
            expect_identical(RTMB::TapeConfig(), configuration)
            expect_identical(.Random.seed, rng)
            score <- .msr_factorized_group_score_tape(g, p)
            for (change in c(FALSE, TRUE)) {
              q <- p
              if (change) {
                q$beta <- q$beta + c(.07, -.03, .02)
                q$scale_beta <- q$scale_beta + .1
                q$log_scale_sd <- q$log_scale_sd + .02
                q$log_response_sd <- q$log_response_sd + .03
                q$log_dispersion <- q$log_dispersion + .05
              }
              log_density <- sq_numeric_log_density(record, q)
              for (z in c(-8, -2, 0, 2, 8)) {
                node <- q; node$scale_z <- z
                expected <- -sq_group_native_response(g, node) + dnorm(z, log = TRUE)
                v <- score$obj$par; v[score$outer] <- unlist(q, use.names = FALSE); v[score$z] <- z
                expect_equal(log_density(z), expected, tolerance = 1e-8)
                expect_equal(log_density(z), -as.numeric(score$obj$fn(v)), tolerance = 1e-8)
                expect_equal(log_density(z), -.msr_factorized_group_joint_nll(g, node), tolerance = 1e-8)
              }
            }
            expect_identical(RTMB::TapeConfig(), configuration)
            expect_identical(.Random.seed, rng)
          }
          expect_identical(s, before)
          expect_error(.msr_factorized_conditional_cache(s), "Weighted SQ")
        }
})

test_that("cached conditional latent distributions match ordinary and independent native integrals", {
  skip_if_not_installed("RTMB")
  cases <- expand.grid(kernel = c("gaussian", "exp"), kind = c("gaussian", "poisson", "response"))
  for (i in seq_len(nrow(cases))) {
    kind <- as.character(cases$kind[i]); family <- if (kind == "response") "gaussian" else kind
    s <- private_sq_fixture(as.character(cases$kernel[i]), family, response = kind == "response")$spec
    g <- .msr_factorized_group(s, 2L); p <- private_sq_parameters(g)
    p$beta[3] <- .12; p$log_scale_sd[] <- log(.3); p$scale_z <- p$response_z <- numeric()
    rng <- .Random.seed; configuration <- RTMB::TapeConfig()
    record <- sq_numeric_group(g, p)
    cached <- sq_numeric_log_density(record, p)
    ordinary <- function(z) { q <- p; q$scale_z <- z; -.msr_factorized_group_joint_nll(g, q) }
    bound <- if (family == "gaussian") -length(g$y) * (p$log_dispersion + .5 * log(2 * pi)) else 0
    tau <- exp(p$log_scale_sd)
    locations <- unique(drop(g$scales[[1L]]$X %*% p$scale_beta[g$scales[[1L]]$beta_index]))
    knots <- as.vector(outer(seq(-30, 30, by = 1), locations, "-") / tau)
    a <- .msr_factorized_distribution(ordinary, bound, knots)
    b <- .msr_factorized_distribution(cached, bound, knots)
    expect_equal(b$log_normalizer, a$log_normalizer, tolerance = 1e-8)
    expect_equal(b$log_normalizer, -sq_group_native_integral(g, p), tolerance = 1e-8)
    expect_equal(b$cdf(c(-3, 0, 3)), a$cdf(c(-3, 0, 3)), tolerance = 1e-9)
    expect_equal(b$cdf(c(-3, 0, 3), FALSE), a$cdf(c(-3, 0, 3), FALSE), tolerance = 1e-9)
    probabilities <- c(.025, .5, .975)
    expect_equal(b$quantile(probabilities), a$quantile(probabilities), tolerance = 1e-7)
    log_normalizer <- -sq_group_native_integral(g, p)
    native_density <- function(z) vapply(z, function(t) {
      if (exp(bound + dnorm(t, log = TRUE) - log_normalizer) == 0) return(0)
      q <- p; q$scale_z <- t
      exp(-sq_group_native_response(g, q) + dnorm(t, log = TRUE) - log_normalizer)
    }, 0)
    native_cdf <- function(q, lower = TRUE) vapply(q, function(t) {
      integrate(native_density, if (lower) -Inf else t, if (lower) t else Inf,
        rel.tol = 1e-11, abs.tol = 0, subdivisions = 1000L, stop.on.error = TRUE)$value
    }, 0)
    expect_equal(b$cdf(c(-3, 0, 3)), native_cdf(c(-3, 0, 3)), tolerance = 1e-9)
    expect_equal(b$cdf(c(-3, 0, 3), FALSE), native_cdf(c(-3, 0, 3), FALSE), tolerance = 1e-9)
    expect_equal(native_cdf(b$quantile(probabilities)), probabilities, tolerance = 1e-7)
    expect_equal(b$cdf(c(-Inf, Inf)), c(0, 1))
    expect_true(b$infinite_support)
    expect_lte(b$relative_error, 1e-10)
    expect_identical(RTMB::TapeConfig(), configuration)
    expect_identical(.Random.seed, rng)
    expect_error(.msr_factorized_conditional_cache(s), "Weighted SQ")
  }
})

test_that("group numeric tape freezes source data and rejects malformed layouts before recording", {
  skip_if_not_installed("RTMB")
  s <- private_sq_fixture()$spec; g <- .msr_factorized_group(s, 2L)
  p <- private_sq_parameters(g); p$scale_z <- p$response_z <- numeric()
  record <- sq_numeric_group(g, p); density <- sq_numeric_log_density(record, p)
  original <- density(0)
  altered <- g; altered$y[1] <- altered$y[1] + 1
  changed <- sq_numeric_log_density(sq_numeric_group(altered, p), p)
  expect_identical(density(0), original)
  expect_gt(abs(changed(0) - original), 1e-4)
  expect_equal(density(0), -sq_group_native_response(g, record$template) + dnorm(0, log = TRUE), tolerance = 1e-8)
  expect_equal(changed(0), -sq_group_native_response(altered, record$template) + dnorm(0, log = TRUE), tolerance = 1e-8)
  expect_error(.msr_factorized_group_log_tape(s, record$initial, record$unpack), "one recoded")
  bad <- g; bad$observation_ids <- s$observation_ids
  expect_error(.msr_factorized_group_log_tape(bad, record$initial, record$unpack), "observation metadata")
  expect_error(.msr_factorized_group_log_tape(g, record$initial, NULL), "unpack function")
  expect_error(.msr_factorized_group_log_tape(g, c(record$initial, Inf), record$unpack), "finite numeric")
  expect_error(.msr_factorized_group_log_tape(g, record$initial, function(x) unname(record$template)), "matching")
  wrong <- record$template; names(wrong)[names(wrong) == "scale_z"] <- "scale_z_extra"
  expect_error(.msr_factorized_group_log_tape(g, record$initial, function(x) wrong), "matching")
  wrong <- record$template; wrong$response_z <- 0
  expect_error(.msr_factorized_group_log_tape(g, record$initial, function(x) wrong), "matching")
  expect_error(.msr_factorized_group_joint_nll(g, wrong), "integrated response")
  # A matching total length must not hide incompatible coefficient/SD blocks.
  wrong <- record$template
  wrong$scale_beta <- c(wrong$scale_beta, tail(wrong$beta, 1L))
  wrong$beta <- head(wrong$beta, -1L)
  malformed <- unlist(wrong, use.names = FALSE)
  with_mocked_bindings({
    expect_error(.msr_factorized_group_log_tape(g, malformed, function(x) wrong), "matching")
  }, MakeTape = function(...) stop("Tape reached"), .package = "RTMB")
  launched <- 0L
  for (ram in c(NA_real_, 1024^3)) with_mocked_bindings({
    with_mocked_bindings({
      expect_error(.msr_factorized_group_log_tape(g, record$initial, record$unpack),
        if (is.na(ram)) "Cannot detect RAM" else "RAM budget")
    }, MakeTape = function(...) { launched <<- launched + 1L; stop("Tape reached") }, .package = "RTMB")
  }, .msr_total_ram_bytes = function() ram, .package = "multiScaleR")
  expect_equal(launched, 0L)
  expect_error(.msr_factorized_conditional_cache(s), "Weighted SQ")
})

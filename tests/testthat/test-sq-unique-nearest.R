test_that("unique-nearest SQ retains scaled analytic values and two log-scale derivatives", {
  skip_if_not_installed("RTMB")
  buffers <- list(list(d = c(0, 1), v = c(0, 1)))
  for (kernel in c("gaussian", "exp")) for (rate in c(3, 30, 300, 600, 700)) {
    power <- if (kernel == "gaussian") 2 else 1
    t <- if (power == 2) -log(2 * rate) / 2 else -log(rate)
    multiplier <- exp(rate / 2)
    w <- exp(-rate); value <- 1 / (1 + w)
    first <- power * rate * (1 - w) / (2 * (1 + w))
    second <- -power^2 * rate * (1 - w) / (2 * (1 + w)) -
      power^2 * rate^2 * w / (1 + w)^2
    expected <- c(value, value * first, value * (first^2 + second))
    .msr_scale_surface_ram(list(surface = list(buffers = list(x = buffers))))
    configuration <- RTMB::TapeConfig()
    obj <- RTMB::MakeADFun(function(p)
      multiplier * .msr_scale_surface_raw(buffers, p$t, kernel, TRUE), list(t = t), silent = TRUE)
    actual <- c(as.numeric(obj$fn(obj$par)), as.numeric(obj$gr(obj$par)), as.numeric(obj$he(obj$par)))
    expect_true(all(is.finite(actual)))
    expect_equal(multiplier * .msr_scale_surface_raw(buffers, t, kernel), value, tolerance = 1e-6)
    expect_equal(actual[1L], expected[1L], tolerance = 1e-6)
    expect_equal(actual[2L], expected[2L], tolerance = 2e-5)
    expect_equal(actual[3L], expected[3L], tolerance = 2e-5)
    expect_identical(RTMB::TapeConfig(), configuration)
  }
})

test_that("unique-nearest finite saturation keeps explicit finite AD endpoints", {
  skip_if_not_installed("RTMB")
  buffers <- list(list(d = c(0, 1), v = c(0, 1)))
  for (kernel in c("gaussian", "exp")) for (t in c(-1000, 1000)) {
    .msr_scale_surface_ram(list(surface = list(buffers = list(x = buffers))))
    obj <- RTMB::MakeADFun(function(p)
      .msr_scale_surface_raw(buffers, p$t, kernel, TRUE), list(t = t), silent = TRUE)
    actual <- c(as.numeric(obj$fn(obj$par)), as.numeric(obj$gr(obj$par)), as.numeric(obj$he(obj$par)))
    expect_true(all(is.finite(actual)))
    expect_equal(actual, c(if (t < 0) 0 else .5, 0, 0), tolerance = 1e-12)
    expect_equal(actual[1L], .msr_scale_surface_raw(buffers, t, kernel), tolerance = 1e-12)
  }
})

test_that("unique-nearest prepared SQ retains full large-SD nuisance derivatives", {
  skip_if_not_installed("RTMB")
  maximum <- c(value = 0, gradient = 0, hessian = 0)
  for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        f <- private_sq_fixture(kernel, family, response = response)
        distance <- f$inputs$d_list[[2L]]
        nearest <- min(distance)
        distance[2L] <- (nearest + min(distance[distance > nearest])) / 2
        f$inputs$d_list[[2L]] <- distance
        # Re-prepare the changed synthetic support once. Its original recipe
        # supplies a new fixed reference, never a pooled or moving transform.
        s <- .msr_scale_prepare(f$formula, f$inputs, f$data, f$scales,
          family, "key", reference_scale = 250)
        expect_identical(s$rows, f$spec$rows)
        expect_identical(s$observation_ids, f$spec$observation_ids)
        expect_identical(s$offset, f$spec$offset)
        if (family == "binomial") s$trials[] <- 4
        before <- s; g <- .msr_factorized_group(s, 2L)
        cell <- g$surface$buffers[[1L]][[which(g$rows == 2L)[1L]]]
        expect_equal(sum(cell$d == min(cell$d)), 1)
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
  cat(sprintf("\nSQ_UNIQUE_AD SD=8 cases=10 value_error=%.12g gradient_error=%.12g hessian_error=%.12g\n",
    maximum[1L], maximum[2L], maximum[3L]))
})

test_that("analytic SQ engines retain independent unique-nearest raw values and slopes", {
  for (kernel in c("gaussian", "exp")) for (rate in c(3, 30, 300, 600, 700)) {
    power <- if (kernel == "gaussian") 2 else 1
    t <- if (power == 2) -log(2 * rate) / 2 else -log(rate)
    multiplier <- exp(rate / 2); w <- exp(-rate)
    value <- 1 / (1 + w)
    slope <- value * power * rate * (1 - w) / (2 * (1 + w))
    inputs <- list(kernel = kernel,
      scale_vars = msr_vars(rough = surface_var("elev", "sq", weighted = TRUE)),
      d_list = list(c(0, 1)), raw_cov = list(data.frame(elev = c(0, 1))))
    for (implementation in c("reference", "stream")) {
      # Raw comparison only. Each inference method keeps its own reference
      # and support policy; no new-data normalization is inferred here.
      engine <- .msr_prepare_surface_engine(inputs, "rough", 1L,
        standardize = FALSE, implementation = implementation)
      before <- engine
      z <- .msr_eval_surface_engine(engine, t)
      actual <- multiplier * c(z$raw, z$d_scale)
      expect_true(all(is.finite(actual)))
      expect_equal(actual[1L], value, tolerance = 1e-6)
      expect_equal(actual[2L], slope, tolerance = 2e-5)
      expect_identical(engine, before)
    }
  }
})

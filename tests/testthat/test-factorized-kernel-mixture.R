test_that("group mixtures retain exact endpoints and stable zero-density arithmetic", {
  fn <- .msr_factorized_mixture_nll
  a <- c(.5, 700, Inf); b <- c(900, .2, 3)
  expect_identical(fn(a, b, 0), sum(b))
  expect_identical(fn(a, b, 1), sum(a))
  expect_equal(fn(a, b, .3), mixture_native_total(a, b, .3), tolerance = 1e-12)
  expect_identical(fn(Inf, Inf, .5), Inf)
  expect_equal(fn(Inf, 2, .25), 2 - log(.75), tolerance = 1e-12)
  expect_equal(fn(2, Inf, .25), 2 - log(.25), tolerance = 1e-12)
  expect_equal(fn(c(1e300, -1e300), c(1e300, -1e300), .7), 0)
  expect_equal(fn(c(1000, 1500), c(1000, 1500), .7), 2500)
  for (w in list(NA_real_, Inf, NaN, -.1, 1.1, numeric(), c(.2, .5), "0.5"))
    expect_error(fn(1, 2, w), "weight")
  expect_error(fn(numeric(), numeric(), .5), "nonempty")
  expect_error(fn(1:2, 1, .5), "matching")
  expect_error(fn(matrix(1), 2, .5), "vectors")
  expect_error(fn(c(1, NA), 1:2, .5), "missing")
  expect_error(fn(-Inf, 2, .5), "negative-infinite")
})

test_that("mixtures are over annual group joint densities with preserved native normalizers", {
  skip_if_not_installed("RTMB")
  for (surface in c(FALSE, TRUE)) for (kernel in c("gaussian", "exp"))
    for (family in c("gaussian", "poisson", "binomial", "nbinom2"))
      for (response in if (family == "gaussian") c(FALSE, TRUE) else FALSE) {
        f <- mixture_fixture(kernel, family, surface, response)
        s <- f$spec; if (family == "binomial") s$trials[] <- 4
        before <- s; p <- private_sq_parameters(s); p$beta[3] <- .12
        parts <- lapply(seq_len(3), function(i) {
          g <- .msr_factorized_group(s, i)
          lapply(c("broad", "nearest"), function(k) .msr_factorized_kernel_limit(g, k))
        })
        q <- p; q$scale_beta <- q$scale_z <- q$log_scale_sd <- q$response_z <- numeric()
        values <- function(parameters, ad = FALSE) lapply(1:2, function(k)
          do.call(c, lapply(parts, function(pair) sum(.msr_factorized_response_nll(pair[[k]], parameters, ad)))))
        a <- mixture_native_nll(s, q, "broad"); b <- mixture_native_nll(s, q, "nearest")
        got <- values(q)
        expect_equal(as.numeric(got[[1]]), unname(a), tolerance = 1e-10)
        expect_equal(as.numeric(got[[2]]), unname(b), tolerance = 1e-10)
        for (w in c(0, .001, .25, .5, .75, .999, 1))
          expect_equal(.msr_factorized_mixture_nll(got[[1]], got[[2]], w),
            mixture_native_total(a, b, w), tolerance = 1e-10)
        for (pair in parts) for (lim in pair) {
          expect_identical(lim$center, s$center)
          expect_identical(lim$spread, s$spread)
          expect_length(lim$surface$buffers, 0)
          expect_length(lim$parameters$log_scale_sd, 0)
          expect_identical(lim$parameters$log_response_sd, s$parameters$log_response_sd)
          .msr_scale_surface_ram(lim)
        }
        obj <- RTMB::MakeADFun(function(parameters) {
          v <- values(parameters, TRUE)
          .msr_factorized_mixture_nll(v[[1]], v[[2]], .35, TRUE)
        }, q, silent = TRUE)
        native <- function(x) {
          v <- private_sq_unpack(x, q)
          mixture_native_total(mixture_native_nll(s, v, "broad"), mixture_native_nll(s, v, "nearest"), .35)
        }
        expect_equal(as.numeric(obj$fn(obj$par)), native(obj$par), tolerance = 1e-9)
        expect_equal(as.numeric(obj$gr(obj$par)), private_sq_fd_gradient(native, obj$par), tolerance = 2e-5)
        expect_equal(unname(obj$he(obj$par)), private_sq_fd_hessian(native, obj$par), tolerance = 2e-5)
        expect_identical(s, before)
        if (surface) expect_error(.msr_factorized_validate(s), "Weighted SQ")
      }
})

test_that("shared group mixing differs from mixing observation likelihoods separately", {
  # Two perfectly informative observations within one group share its kernel.
  # Both group components have joint density.09; observation-wise mixing gives.25.
  broad <- -log(c(.9, .1)); near <- -log(c(.1, .9))
  grouped <- .msr_factorized_mixture_nll(sum(broad), sum(near), .5)
  wrong <- .msr_factorized_mixture_nll(broad, near, .5)
  expect_equal(grouped, -log(.09), tolerance = 1e-12)
  expect_equal(wrong, -log(.25), tolerance = 1e-12)
  expect_gt(abs(grouped - wrong), .5)
  # Sum the utility's implied probabilities over two exhaustive histories.
  for (w in c(0, .001, .35, .999, 1)) {
    probability <- vapply(1:2, function(i) exp(-.msr_factorized_mixture_nll(
      -log(c(.1, .9))[i], -log(c(.7, .3))[i], w)), 0)
    expect_equal(sum(probability), 1, tolerance = 1e-12)
  }
})

test_that("exact NB2-to-Poisson limits commute with group mixture arithmetic", {
  for (surface in c(FALSE, TRUE)) for (kernel in c("gaussian", "exp")) {
    s <- mixture_fixture(kernel, "nbinom2", surface)$spec
    before <- s; p <- private_sq_parameters(s); p$beta[3] <- .12
    p$scale_beta <- p$scale_z <- p$log_scale_sd <- p$log_dispersion <- p$response_z <- numeric()
    limits <- lapply(1:3, function(i) {
      g <- .msr_factorized_group(s, i)
      lapply(c("broad", "nearest"), function(k) .msr_factorized_kernel_limit(g, k, poisson = TRUE))
    })
    for (pair in limits) for (lim in pair) {
      expect_identical(lim$family, "poisson")
      expect_length(lim$parameters$log_dispersion, 0)
      expect_identical(lim$center, s$center)
      expect_identical(lim$spread, s$spread)
    }
    a <- vapply(limits, function(pair) sum(.msr_factorized_response_nll(pair[[1]], p)), 0)
    b <- vapply(limits, function(pair) sum(.msr_factorized_response_nll(pair[[2]], p)), 0)
    native <- s; native$family <- "poisson"
    ra <- mixture_native_nll(native, p, "broad"); rb <- mixture_native_nll(native, p, "nearest")
    for (w in c(0, .35, 1)) expect_equal(.msr_factorized_mixture_nll(a, b, w),
      mixture_native_total(ra, rb, w), tolerance = 1e-10)
    expect_identical(s, before)
    if (surface) expect_error(.msr_factorized_validate(s), "Weighted SQ")
  }
})

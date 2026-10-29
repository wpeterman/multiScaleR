smooth_kernel_fixture <- function(n = 60L, p = 2L, nb = 24L, seed = 19L) {
  set.seed(seed)
  d <- seq(0, 1, length.out = nb)
  nm <- paste0("hab", seq_len(p))
  V <- C <- stats::setNames(vector("list", p), nm)
  for (j in seq_len(p)) {
    C[[j]] <- matrix(sample(0:5, n * nb, TRUE), n, nb)
    C[[j]][, 1L] <- 1
    radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(13 * d)) +
      outer(rnorm(n), d) + matrix(rnorm(n * nb, sd = 0.2), n, nb)
    V[[j]] <- radial * C[[j]]
  }
  list(bin_rep = d, vsum = V, csum = C, n_points = n, covariates = nm)
}

smooth_numeric_gradient <- function(fn, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    a <- b <- x
    a[j] <- a[j] + h
    b[j] <- b[j] - h
    (fn(a) - fn(b)) / (2 * h)
  }, numeric(1))
}

test_that("shared kernels preserve finite-cell means and row/parameter order", {
  b <- smooth_kernel_fixture()
  rows <- c(10:20, 1:9)
  for (kernel in c("gaussian", "exp", "expow")) {
    e <- .msr_prepare_kernel_engine(b, c("hab2", "hab1"), rows, kernel, FALSE)
    par <- c(log(c(0.18, 0.4)), if (kernel == "expow") log(c(1.3, 2.4)))
    actual <- .msr_eval_kernel_engine(e, par)
    expected <- .msr_binned_kernel_means(
      b, e$covariates, e$covariates, c(0.18, 0.4), c(1.3, 2.4), kernel, rows)
    expect_equal(actual$raw, expected, tolerance = 1e-12)
    expect_equal(actual$value, expected, tolerance = 1e-12)
  }
})

test_that("analytic kernel derivatives include quotient and sample scaling", {
  b <- smooth_kernel_fixture()
  for (kernel in c("gaussian", "exp", "expow")) {
    for (standardize in c(FALSE, TRUE)) {
      e <- .msr_prepare_kernel_engine(b, kernel = kernel, standardize = standardize)
      par <- c(log(c(0.13, 0.35)), if (kernel == "expow") log(c(0.8, 3)))
      ans <- .msr_eval_kernel_engine(e, par)
      for (j in seq_along(par)) {
        plus <- minus <- par
        plus[j] <- plus[j] + 1e-5
        minus[j] <- minus[j] - 1e-5
        numerical <- (.msr_eval_kernel_engine(e, plus)$value -
                        .msr_eval_kernel_engine(e, minus)$value) / 2e-5
        analytical <- matrix(0, e$n, e$p)
        col <- (j - 1L) %% e$p + 1L
        analytical[, col] <- if (j <= e$p) ans$d_scale[, col] else ans$d_shape[, col]
        expect_equal(unname(numerical), analytical, tolerance = 2e-7)
      }
      if (standardize) {
        expect_equal(unname(colMeans(ans$value)), rep(0, e$p), tolerance = 1e-12)
        expect_equal(unname(apply(ans$value, 2L, sd)), rep(1, e$p), tolerance = 1e-12)
      }
    }
  }
})

test_that("invalid support, constant predictors and unsupported kernels fail clearly", {
  b <- smooth_kernel_fixture()
  expect_error(.msr_prepare_kernel_engine(b, kernel = "fixed"), "smooth")
  expect_error(.msr_prepare_kernel_engine(b, rows = 0), "indices")
  bad <- b
  bad$csum[[1]][1, 1] <- -1
  expect_error(.msr_prepare_kernel_engine(bad), "sum/count")
  bad <- b
  bad$vsum[[1]] <- bad$vsum[[1]][, -1]
  expect_error(.msr_prepare_kernel_engine(bad), "sum/count")
  e <- .msr_prepare_kernel_engine(b)
  expect_error(.msr_eval_kernel_engine(e, c(a = -1, b = -2)), "names")
  b$vsum[[1]] <- b$csum[[1]]
  expect_error(.msr_eval_kernel_engine(.msr_prepare_kernel_engine(b), c(-1, -2)),
               "Undefined")
  b$vsum[[1]][1, ] <- b$csum[[1]][1, ] <- 0
  expect_error(.msr_prepare_kernel_engine(b), "support")
})

test_that("profile likelihood scores match differences for all four families", {
  skip_if_not_installed("MASS")
  b <- smooth_kernel_fixture(n = 100L)
  for (kernel in c("gaussian", "exp", "expow")) {
    e <- .msr_prepare_kernel_engine(b, kernel = kernel)
    kp <- c(log(c(0.15, 0.4)), if (kernel == "expow") log(c(1.2, 2)))
    x <- .msr_eval_kernel_engine(e, kp)$value
    offset <- seq(-0.2, 0.2, length.out = e$n)
    eta <- as.vector(0.3 + x %*% c(0.6, -0.4) + offset)
    set.seed(17)
    ys <- list(gaussian = rnorm(e$n, eta, 0.7), poisson = rpois(e$n, exp(eta)),
               binomial = rbinom(e$n, 5, plogis(eta)),
               nbinom2 = rnbinom(e$n, mu = exp(eta), size = 1.5))
    for (family in names(ys)) {
      obj <- .msr_profile_glm_objective(e, ys[[family]], family,
        offset = offset, weights = rep(c(1, 2), length.out = e$n),
        trials = rep(5, e$n))
      par <- c(kp, if (family == "nbinom2") log(1.5))
      expect_equal(obj$gr(par), smooth_numeric_gradient(obj$fn, par), tolerance = 3e-5)
      fit <- obj$evaluate(par)
      expect_true(is.finite(fit$value))
      expect_equal(obj$evaluations(), 1L + 2L * length(par) + 1L)
    }
  }
})

test_that("NB2 estimated dispersion matches an independent MASS fit at fixed scales", {
  skip_if_not_installed("MASS")
  e <- .msr_prepare_kernel_engine(smooth_kernel_fixture(n = 200L))
  kp <- log(c(0.15, 0.4))
  x <- .msr_eval_kernel_engine(e, kp)$value
  offset <- seq(-0.3, 0.3, length.out = e$n)
  weights <- rep(c(1, 2), length.out = e$n)
  set.seed(941)
  y <- rnbinom(e$n, mu = exp(0.4 + x %*% c(0.7, -0.5) + offset), size = 0.7)
  obj <- .msr_profile_glm_objective(e, y, "nbinom2", offset = offset, weights = weights)
  fit <- nlminb(log(1), function(z) obj$fn(c(kp, z)),
                function(z) tail(obj$gr(c(kp, z)), 1L),
                lower = log(0.01), upper = log(1000),
                control = list(rel.tol = 1e-12))
  ref <- MASS::glm.nb(y ~ x + offset(offset), weights = weights,
                      control = glm.control(epsilon = 1e-11, maxit = 100))
  expect_equal(fit$objective, -as.numeric(logLik(ref)), tolerance = 1e-7)
  expect_equal(exp(fit$par), ref$theta, tolerance = 1e-5)
  expect_equal(unname(obj$evaluate(c(kp, fit$par))$beta), unname(coef(ref)),
               tolerance = 1e-6)
})

test_that("profile fitting rejects invalid responses and singular designs", {
  e <- .msr_prepare_kernel_engine(smooth_kernel_fixture())
  expect_error(.msr_profile_glm_objective(e, rep(0.5, e$n), "nbinom2"), "integers")
  expect_error(.msr_profile_glm_objective(e, rep(2, e$n), "binomial"), "successes")
  expect_error(.msr_profile_glm_objective(e, rep(0, e$n), weights = rep(0, e$n)),
               "positive weights")
  obj <- .msr_profile_glm_objective(e, seq_len(e$n),
                                   X_fixed = matrix(1, e$n, 2))
  expect_error(obj$fn(log(c(0.1, 0.3))), "full-rank")
})

test_that("multistart selects the best eligible fit and retains failed starts", {
  obj <- list(fn = function(x) {
    if (x > 2.9) stop("intentional failure")
    (x^2 - 1)^2 + 0.1 * x
  }, gr = function(x) 4 * x * (x^2 - 1) + 0.1)
  fit <- .msr_multistart(obj, matrix(c(0.8, -0.8, 3), ncol = 1), -4, 4)
  expect_equal(fit$winner, 2L)
  expect_lt(fit$best$par, 0)
  expect_false(fit$records$eligible[3])
  expect_match(fit$records$message[3], "intentional failure")
  expect_length(fit$fits, 3L)
  none <- .msr_multistart(obj, matrix(3, ncol = 1), -4, 4)
  expect_null(none$best)
  expect_true(is.na(none$winner))
  bound <- .msr_multistart(list(fn = function(x) x, gr = function(x) 1),
                           matrix(0.5, ncol = 1), 0, 1)
  expect_true(bound$records$eligible)
  expect_equal(bound$best$par, 0)
  polished <- .msr_multistart(
    list(fn = function(x) (x - 0.3)^2, gr = function(x) 2 * (x - 0.3)),
    matrix(2, ncol = 1), -4, 4, control = list(iter.max = 1L, eval.max = 5L))
  expect_true(polished$records$eligible)
  expect_true(polished$records$polished)
  expect_equal(polished$records$primary_convergence, 1L)
  expect_equal(polished$best$par, 0.3, tolerance = 1e-6)
})

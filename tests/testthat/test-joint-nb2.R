testthat::skip_if_not_installed("MASS")

joint_nb2_fit <- function(f, ...) suppressWarnings(multiScale_optim(
  f$model, f$inputs, backend = "joint_nb2", n_cores = 1, verbose = FALSE, ...))

test_that("stable NB2 likelihood includes the exact boundary and large counts", {
  y <- c(0, 1, 2, 19, 10001, 20000)
  w <- c(1, 2, 3, 1, 2, 1); eta <- log(c(.2, 1, 3, 18, 8000, 25000))
  norm <- .msr_nb2_normalizer(y, w)
  for (a in c(0, 1e-10, 1e-6, .009, .011, .2, 10)) {
    z <- .msr_nb2_terms(y, eta, a, w, norm)
    ref <- if (a == 0) -sum(w * dpois(y, exp(eta), log = TRUE)) else
      -sum(w * dnbinom(y, mu = exp(eta), size = 1/a, log = TRUE))
    expect_equal(z$value, ref, tolerance = 1e-7)
    if (a > 1e-8) {
      h <- a * 1e-4
      numeric_score <- (.msr_nb2_terms(y, eta, a+h, w, norm)$value -
        .msr_nb2_terms(y, eta, a-h, w, norm)$value) / (2*h)
      expect_equal(z$alpha_score, numeric_score, tolerance = 2e-5)
    }
  }
  z <- .msr_nb2_terms(y, eta, 0, w, norm)
  expect_equal(z$alpha_score, sum(w * (y - (y-exp(eta))^2)) / 2, tolerance = 1e-10)
  # Positive dispersion permits extreme log means without exp(eta) overflow
  # contaminating the likelihood or its derivatives.
  z <- .msr_nb2_terms(c(0, 4), c(-1000, 1000), .2, c(1, 1),
                      .msr_nb2_normalizer(c(0, 4), c(1, 1)))
  expect_true(all(is.finite(c(z$value, z$score, z$alpha_score, z$curvature))))
})

test_that("joint NB2 matches profile fitting and observed uncertainty", {
  f <- profile_glm_fixture()
  j <- joint_nb2_fit(f)
  p <- suppressWarnings(profile_glm_fit(f))
  expect_identical(j$backend, "joint_nb2")
  expect_equal(j$optim_results$value, p$optim_results$value, tolerance = 1e-8)
  expect_equal(j$scale_est, p$scale_est, tolerance = 2e-5)
  expect_equal(j$dispersion_est$theta, p$dispersion_est$theta, tolerance = 2e-5)
  expect_equal(unname(vcov(j)), unname(vcov(p)), tolerance = 2e-4)
  expect_equal(predict(j), predict(p), tolerance = 2e-5)
  prepared <- prepare_prediction(j)
  expect_equal(predict(prepared), predict(j))
  expect_equal(predict(prepared, se.fit = TRUE, uncertainty = "conditional"),
               predict(j, se.fit = TRUE, uncertainty = "conditional"))
  if (j$uncertainty$available)
    expect_equal(predict(prepared, interval = "confidence"),
                 predict(j, interval = "confidence"))
  expect_equal(predict(j, se.fit = TRUE)$se.fit, predict(p, se.fit = TRUE)$se.fit,
               tolerance = 2e-4)
  expect_equal(unname(predict(j)), unname(predict(j$opt_mod, type = "response")), tolerance = 1e-7)
  expect_equal(nrow(j$profile_fit$start_results), 12L)
  expect_true(all(j$profile_fit$start_results$conditional_gap[
    j$profile_fit$start_results$eligible] <= 1e-5))
  expect_equal(-as.numeric(logLik(j$opt_mod)), j$optim_results$value, tolerance = 1e-10)
  expect_identical(eval(formals(multiScale_optim)$backend)[1L], "refit")
  expect_error(joint_nb2_fit(f, profile_control = list(theta_bounds = c(1e-4, 1e6))), "Inf")
  f$model <- glm(y ~ hab, data = f$data, family = poisson())
  expect_error(joint_nb2_fit(f), "requires an NB2")
})

test_that("exact Poisson fits preserve prediction and uncertainty contracts", {
  f <- profile_glm_fixture(n = 300, seed = 723, theta = Inf)
  # Deterministically underdispersed counts place the NB2 maximum at alpha zero.
  f$data$y <- round(f$mu)
  f$model <- glm(y ~ hab + z + group + offset(log(exposure)), data = f$data,
                 family = MASS::negative.binomial(2))
  j <- joint_nb2_fit(f)
  expect_identical(j$dispersion_est$exact_poisson, TRUE)
  expect_identical(j$dispersion_est$alpha, 0)
  prepared <- prepare_prediction(j)
  expect_equal(predict(prepared), predict(j))
  expect_equal(predict(prepared, interval = "confidence", uncertainty = "conditional"),
               predict(j, interval = "confidence", uncertainty = "conditional"))
  expect_error(predict(prepared, se.fit = TRUE), "Joint Wald uncertainty unavailable")
  expect_identical(j$dispersion_est$theta, Inf)
  expect_false(j$uncertainty$available)
  expect_true(j$uncertainty$boundary)
  expect_true(all(is.na(j$scale_est$SE)))
  expect_error(vcov(j), "exact Poisson")
  expect_error(predict(j, se.fit = TRUE), "exact Poisson")
  independent_data <- f$data
  independent_data$hab <- j$opt_mod$model$hab
  independent <- glm(y ~ hab + z + group + offset(log(exposure)),
    family = poisson(), data = independent_data, control = glm.control(epsilon = 1e-11))
  expect_equal(-as.numeric(logLik(independent)), j$optim_results$value, tolerance = 1e-9)
  # The fitted model is still the NB2 closure, with its dispersion parameter.
  expect_equal(attr(logLik(j$opt_mod), "df"), length(coef(independent)) + 1L)
  expect_equal(unname(predict(j)), unname(fitted(independent)), tolerance = 1e-7)
  columns <- names(j$profile_fit$contract$objective$evaluate(j$profile_fit$par)$beta)
  X <- model.matrix(independent)[, columns]
  independent_covariance <- solve(crossprod(X, X * fitted(independent)))
  expect_equal(unname(vcov(j, type = "conditional")), unname(independent_covariance),
               tolerance = 1e-9)
  nd <- f$data[c(7, 2, 91), ]; nd$hab <- 999
  expect_equal(predict(j, nd, f$inputs), predict(j)[rownames(nd)], tolerance = 1e-9)
  expect_true(all(is.finite(predict(j, nd, f$inputs, uncertainty = "conditional",
                                    interval = "confidence")$se.fit)))
  expect_error(predict(j, nd[, "z", drop = FALSE], f$inputs), "missing formula")
  nd$group <- "new"
  expect_error(predict(j, nd, f$inputs), "new level")
  restored <- unserialize(serialize(j, NULL))
  expect_equal(predict(restored), predict(j), tolerance = 1e-12)
  expect_equal(.msr_profile_at(j, 1, j$profile_fit$par[1]), j$optim_results$value,
               tolerance = 1e-7)
  expect_s3_class(summary(j), "summary_multiScaleR")
  expect_silent(capture.output(print(summary(j))))
  expect_silent(capture.output(print(j)))
  expect_equal(.msr_parameter_count(j$opt_mod), length(coef(j$opt_mod)) + 1L)
  a <- scale_inference(j, sigma_values = j$scale_est$Mean, n_starts = 2)
  expect_lte(a$reference$value, j$optim_results$value + 1e-6)
  expect_false(a$diagnostics$original_wald_available %in% TRUE)
  setup <- .msr_inference_setup(j, 2)
  ref <- .msr_inference_refit(setup)
  set.seed(741)
  boots <- .msr_scale_bootstrap(setup, ref, ref, 3)
  expect_true(all(boots$ok))
  expect_true(all(boots$theta > 0))
})

test_that("joint NB2 supports weighted two-kernel fits and rejects invalid starts", {
  f <- profile_glm_fixture(n = 180)
  f$inputs$binned$vsum$cover <- f$inputs$binned$vsum$hab[, 30:1]
  f$inputs$binned$csum$cover <- f$inputs$binned$csum$hab[, 30:1]
  engine <- .msr_prepare_kernel_engine(f$inputs$binned, covariates = c("cover", "hab"))
  kd <- .msr_eval_kernel_engine(engine, log(c(.4, .25)))
  f$data$hab <- kd$value[, "hab"]; f$data$cover <- kd$value[, "cover"]
  set.seed(951)
  f$data$y <- rnbinom(180, mu = exp(.5 + f$data$hab - .6*f$data$cover + .2*f$data$z), size = 2)
  f$model <- glm(y ~ cover + z + hab, data = f$data,
                 family = MASS::negative.binomial(2), weights = rep(2, 180))
  j <- joint_nb2_fit(f)
  p <- suppressWarnings(profile_glm_fit(f))
  expect_equal(j$optim_results$value, p$optim_results$value, tolerance = 1e-7)
  expect_equal(j$scale_est$Mean, p$scale_est$Mean, tolerance = 2e-4)
  expect_equal(predict(j), predict(p), tolerance = 2e-4)
  expect_error(joint_nb2_fit(f, par = c(.001, .25)), "Invalid starts")
  expect_error(joint_nb2_fit(f, profile_control = list(iter.max = 1)), "No analytic GLM start")
  # Expanded kernel interactions retain the native design contract.
  f$model <- glm(y ~ cover * hab, data = f$data, family = MASS::negative.binomial(2))
  expect_s3_class(.msr_profile_contract(f$model, f$inputs)$objective$response_design,
    "msr_response_design")
})

test_that("joint NB2 retains weak-scale safeguards and fits exponential-power kernels", {
  f <- profile_glm_fixture(n = 300, seed = 734007, effect = .05)
  j <- joint_nb2_fit(f)
  p <- suppressWarnings(profile_glm_fit(f))
  # A weak generating effect need not trigger local curvature diagnostics.
  # Preserve the existing decision; full profiles can still be bound limited.
  expect_identical(j$uncertainty$available, p$uncertainty$available)
  f <- profile_glm_fixture(n = 300, seed = 735004, effect = 0)
  j <- joint_nb2_fit(f)
  expect_false(j$uncertainty$available)
  expect_true(j$uncertainty$boundary || j$uncertainty$weak)
  expect_error(predict(j, interval = "confidence"), "Joint Wald uncertainty unavailable")
  expect_true(all(is.finite(predict(j))))
  f <- profile_glm_fixture(n = 180, kernel = "expow")
  j <- joint_nb2_fit(f)
  p <- suppressWarnings(profile_glm_fit(f))
  expect_equal(j$optim_results$value, p$optim_results$value, tolerance = 1e-7)
  expect_equal(predict(j), predict(p), tolerance = 2e-4)
  expect_true(is.finite(j$shape_est$Mean))
})

test_that("joint NB2 derivatives cover multiple kernels and shape parameters", {
  for (kernel in c("gaussian", "exp", "expow")) {
    f <- profile_glm_fixture(n = 100, kernel = kernel)
    b <- f$inputs$binned
    b$vsum$cover <- b$vsum$hab[, ncol(b$vsum$hab):1]
    b$csum$cover <- b$csum$hab[, ncol(b$csum$hab):1]
    e <- .msr_prepare_kernel_engine(b, kernel = kernel)
    obj <- .msr_joint_nb2_objective(e, f$data$y, X_fixed = model.matrix(~z, f$data),
      offset = log(f$data$exposure), weights = rep(2, 100))
    kp <- c(log(c(.2, .35)), if (kernel == "expow") log(c(1.3, 2.1)))
    par <- c(kp, .4, .2, .6, -.3, log1p(.4))
    numeric_gradient <- vapply(seq_along(par), function(i) {
      hi <- lo <- par; hi[i] <- hi[i]+1e-6; lo[i] <- lo[i]-1e-6
      (obj$joint$fn(hi)-obj$joint$fn(lo))/2e-6
    }, numeric(1))
    expect_equal(unname(obj$joint$gr(par)), numeric_gradient, tolerance = 2e-6)
    par[length(par)] <- -1e-16
    expect_true(is.finite(obj$joint$fn(par)))
    par[length(par)] <- -.01
    expect_error(obj$joint$fn(par), "Invalid NB2")
  }
})

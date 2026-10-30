weighted_glm_fixture <- function(kernel = "gaussian", surface_only = FALSE) {
  f <- profile_glm_fixture(n = 120, kernel = kernel)
  b <- f$inputs$binned
  f$inputs$scale_vars <- msr_vars(hab = kernel_var("hab"),
    rough = surface_var("elev", "sq", weighted = TRUE))
  f$inputs$d_list <- rep(list(b$bin_rep), nrow(f$data))
  f$inputs$raw_cov <- lapply(seq_len(nrow(f$data)), function(i) {
    v <- b$vsum$hab[i, ] / b$csum$hab[i, ]
    data.frame(hab = v, elev = v^2 + sin(i * b$bin_rep))
  })
  f$inputs$raw_cov[[1]]$elev[1] <- NA_real_
  f$inputs$kernel_dat <- data.frame(row.names = rownames(f$data))
  if (surface_only) {
    f$inputs$scale_vars <- f$inputs$scale_vars[2, , drop = FALSE]
    f$inputs$binned <- NULL
  }
  covs <- if (surface_only) "rough" else c("hab", "rough")
  e <- .msr_prepare_surface_engine(f$inputs, covs, seq_len(nrow(f$data)))
  p <- log(if (surface_only) .4 else c(.25, .4))
  kd <- .msr_eval_kernel_engine(e, p)
  for (nm in covs) f$data[[nm]] <- kd$value[, nm]
  f$eta <- .2 + .7 * f$data$rough + .2 * f$data$z +
    if (surface_only) 0 else .5 * f$data$hab
  f$parameters <- p
  f$formula <- reformulate(c(covs, "z", "offset(log(exposure))"), response = "y")
  f
}

test_that("analytic GLMs use weighted SQ values and profile derivatives", {
  for (kernel in c("gaussian", "exp")) for (only in c(FALSE, TRUE)) {
    f <- weighted_glm_fixture(kernel, only)
    for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
      set.seed(319)
      f$data$y <- switch(family, gaussian = f$eta + rnorm(120, sd = .2),
        poisson = rpois(120, exp(f$eta)), binomial = rbinom(120, 1, plogis(f$eta)),
        nbinom2 = rnbinom(120, mu = exp(f$eta), size = 2))
      fam <- switch(family, gaussian = gaussian(), poisson = poisson(),
        binomial = binomial(), nbinom2 = MASS::negative.binomial(2))
      model <- glm(f$formula, f$data, family = fam)
      a <- .msr_profile_contract(model, f$inputs)$objective
      s <- .msr_profile_contract(model, f$inputs, "stream")$objective
      p <- c(f$parameters, if (family == "nbinom2") log(2))
      numeric_score <- vapply(seq_along(p), function(j) {
        hi <- lo <- p; hi[j] <- hi[j] + 1e-5; lo[j] <- lo[j] - 1e-5
        (a$fn(hi) - a$fn(lo)) / 2e-5
      }, 0)
      expect_equal(a$gr(p), numeric_score, tolerance = 2e-5, ignore_attr = TRUE)
      expect_equal(s$fn(p), a$fn(p), tolerance = 1e-9)
      expect_equal(s$gr(p), a$gr(p), tolerance = 1e-8, ignore_attr = TRUE)
      fit <- a$evaluate(p); mu <- fitted(model)
      ref <- -sum(switch(family,
        gaussian = dnorm(f$data$y, mu, sqrt(mean((f$data$y-mu)^2)), log = TRUE),
        poisson = dpois(f$data$y, mu, log = TRUE),
        binomial = dbinom(f$data$y, 1, mu, log = TRUE),
        nbinom2 = dnbinom(f$data$y, mu = mu, size = 2, log = TRUE)))
      expect_equal(fit$value, ref, tolerance = 1e-7)
    }
  }
})

test_that("weighted SQ connects NB2 inference and mixed likelihood machinery", {
  f <- weighted_glm_fixture(surface_only = TRUE)
  set.seed(320); f$data$y <- rnbinom(120, mu = exp(f$eta), size = 3)
  model <- glm(f$formula, f$data, family = MASS::negative.binomial(3))
  for (backend in c("profile_glm", "joint_nb2")) {
    fit <- suppressWarnings(multiScale_optim(model, f$inputs, backend = backend,
      n_cores = 1, verbose = FALSE, profile_control = list(n_starts = 3)))
    expect_equal(predict(fit), predict(prepare_prediction(fit)))
    setup <- .msr_inference_setup(fit, 2L)
    expect_equal(setup$obj$fn(setup$starts[1, ]), fit$optim_results$value, tolerance = 1e-5)
    refit <- .msr_inference_refit(setup, fixed_scale = fit$profile_fit$par[1])
    expect_equal(refit$value, fit$optim_results$value, tolerance = 1e-5)
    profile <- scale_inference(fit, sigma_values = fit$scale_est$Mean,
      n_starts = 2, bootstrap = 0)
    expect_equal(profile$bootstrap, 0)
    expect_true(nrow(profile$grid) >= 3 && all(is.finite(profile$grid$nll)))
  }
  skip_if_not_installed("glmmTMB"); skip_if_not_installed("RTMB")
  f$data$site <- factor(rep(seq_len(20), each = 6))
  set.seed(322); group_effect <- rnorm(20, sd = .7)
  f$data$y <- rpois(120, exp(f$eta + group_effect[f$data$site]))
  m <- glmmTMB::glmmTMB(y ~ rough + z + (1|site) + offset(log(exposure)),
    data = f$data, family = poisson())
  fit <- multiScale_optim(m, f$inputs, backend = "joint_mixed", n_cores = 1,
    verbose = FALSE, profile_control = list(starts = matrix(c(.2,.4,.6), ncol = 1)))
  expect_equal(predict(fit), predict(fit$opt_mod, type = "response"))
  ids <- c(20, 1, 19)
  expect_equal(unname(predict(fit, f$data[ids, ], f$inputs)), unname(predict(fit)[ids]))
  expect_equal(predict(prepare_prediction(fit)), predict(fit))
  expect_error(predict(fit, se.fit = TRUE), "conditional")
  expect_equal(predict(fit, se.fit = TRUE, uncertainty = "conditional"),
    predict(fit$opt_mod, type = "response", se.fit = TRUE))
  units <- f$inputs; units$unit_conv <- units$unit_conv/2
  units$d_list <- lapply(units$d_list, function(d) d*2)
  expect_equal(predict(fit, f$data, units), predict(fit), tolerance = 1e-10)
  changed <- f$inputs; changed$scale_vars$source <- "different"
  expect_error(predict(fit, f$data, changed), "specifications")
  expect_error(predict(fit, f$data, within(f$inputs, raw_cov <- NULL)), "raw buffers")
})

test_that("weighted SQ fits connect to prediction, covariance and scale profiles", {
  for (only in c(FALSE, TRUE)) {
    f <- weighted_glm_fixture(surface_only = only)
    set.seed(321); f$data$y <- f$eta + rnorm(120, sd = .2)
    model <- glm(f$formula, f$data, family = gaussian())
    fit <- suppressWarnings(multiScale_optim(model, f$inputs, backend = "profile_glm",
      n_cores = 1, verbose = FALSE, profile_control = list(n_starts = 3)))
    expect_s3_class(fit$profile_fit$contract$objective$engine, "msr_surface_engine")
    expect_equal(predict(fit), fitted(fit$opt_mod), ignore_attr = TRUE)
    ids <- c(19, 2, 5)
    expect_equal(unname(predict(fit, f$data[ids, ], f$inputs)), unname(predict(fit)[ids]))
    expect_equal(predict(prepare_prediction(fit)), predict(fit))
    expect_equal(predict(unserialize(serialize(fit, NULL))), predict(fit))
    expect_equal(.msr_profile_at(fit, 1, fit$profile_fit$par[1]),
      fit$optim_results$value, tolerance = 1e-5)
    if (fit$uncertainty$available) {
      a <- predict(fit, se.fit = TRUE)
      b <- predict(prepare_prediction(fit), se.fit = TRUE)
      expect_equal(a, b)
      expect_true(all(is.finite(vcov(fit))))
    }
    bad <- f$inputs; bad$raw_cov <- NULL
    expect_error(.msr_profile_contract(model, bad), "raw buffers")
    bad <- f$inputs; bad$kernel <- "expow"
    expect_error(.msr_profile_contract(model, bad), "Surface derivatives")
    bad <- f$inputs; bad$scale_vars$weighted[bad$scale_vars$type == "surface"] <- FALSE
    expect_error(.msr_profile_contract(model, bad), "surface|spatial")
    bad <- f$inputs; bad$scale_vars$source[bad$scale_vars$type == "surface"] <- "changed"
    expect_error(predict(fit, f$data, bad), "specifications")
  }
})

test_that("weighted SQ mixed tapes agree with native Laplace likelihoods and scores", {
  skip_if_not_installed("RTMB"); skip_if_not_installed("glmmTMB")
  for (kernel in c("gaussian", "exp")) for (only in c(FALSE, TRUE)) {
    f <- weighted_glm_fixture(kernel, only)
    # A flat buffer has identically zero SQ and zero scale derivatives.
    f$inputs$raw_cov[[1]]$elev[] <- 5
    f$inputs$raw_cov[[2]]$elev[2] <- Inf
    covs <- if (only) "rough" else c("hab", "rough")
    e <- .msr_prepare_surface_engine(f$inputs, covs, 1:120)
    kd <- .msr_eval_kernel_engine(e, unname(f$parameters))
    for (nm in e$covariates) f$data[[nm]] <- kd$value[, nm]
    f$data$site <- factor(rep(1:20, each = 6))
    set.seed(324); re <- rnorm(20, sd = .7)
    eta <- .2 + .7*f$data$rough + .2*f$data$z + re[f$data$site] +
      log(f$data$exposure) + if (only) 0 else .5*f$data$hab
    for (family in c("poisson", "nbinom2")) {
      f$data$y <- if (family == "poisson") rpois(120, exp(eta)) else
        rnbinom(120, mu = exp(eta), size = 3)
      form <- update(f$formula, . ~ . + (1|site))
      native <- glmmTMB::glmmTMB(form, f$data,
        family = if (family == "poisson") poisson() else glmmTMB::nbinom2())
      design <- .msr_native_mixed_design(native, e, e$covariates,
        unname(f$parameters), rownames(f$data))
      obj <- do.call(.msr_mixed_rtmb, design)
      p <- obj$par
      expect_equal(as.numeric(obj$fn(p)), -as.numeric(logLik(native)), tolerance = 1e-5)
      score <- vapply(seq_along(p), function(j) {
        a <- b <- p; a[j] <- a[j]+1e-5; b[j] <- b[j]-1e-5
        (obj$fn(a)-obj$fn(b))/2e-5
      }, 0)
      expect_equal(unname(obj$gr(p)), score, tolerance = 2e-4, ignore_attr = TRUE)
    }
  }
})

test_that("surface-only conditional logistic models retain complete choice sets", {
  skip_if_not_installed("survival")
  coxph <- survival::coxph
  f <- weighted_glm_fixture(surface_only = TRUE)
  f$data$choice <- rep(1:40, each = 3)
  f$data$case_ <- rep(c(1,0,0), 40)
  form <- as.formula("case_ ~ rough + z + strata(choice)", env = asNamespace("survival"))
  model <- survival::clogit(form, f$data, model = TRUE, x = TRUE, y = TRUE)
  fit <- multiScale_optim(model, f$inputs, backend = "conditional_logit",
    n_cores = 1, verbose = FALSE, profile_control = list(n_starts = 3))
  expect_equal(predict(prepare_prediction(fit)), predict(fit))
  expect_equal(predict(fit, f$data, f$inputs), predict(fit))
  expect_equal(predict(fit, se.fit = TRUE, uncertainty = "conditional"),
    predict(prepare_prediction(fit), se.fit = TRUE, uncertainty = "conditional"))
  broken <- f$inputs; rownames(broken$kernel_dat)[1] <- "absent"
  expect_error(.msr_clogit_contract(model, broken), "row")
})

test_that("actual surface-only preparation reaches the analytic GLM backend", {
  skip_if_not_installed("terra"); skip_if_not_installed("sf")
  set.seed(325)
  r <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 1200,
    ymin = 0, ymax = 1200, crs = "EPSG:3857")
  terra::values(r) <- rnorm(144); names(r) <- "elev"
  xy <- expand.grid(x = seq(350, 850, length.out = 5), y = seq(350, 850, length.out = 4))
  rownames(xy) <- paste0("point", 1:20)
  pts <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 3857)
  inputs <- kernel_prep(pts, r, max_D = 300, sigma = 120, bin = FALSE,
    scale_vars = msr_vars(rough = surface_var("elev", "sq", weighted = TRUE)), verbose = FALSE)
  expect_null(inputs$binned)
  d <- inputs$kernel_dat; d$z <- rnorm(20); d$y <- d$rough + .2*d$z + rnorm(20, sd = .1)
  m <- glm(y ~ rough + z, data = d)
  o <- .msr_profile_contract(m, inputs)$objective
  expect_equal(o$fn(log(120/inputs$unit_conv)),
    -sum(dnorm(d$y, fitted(m), sqrt(mean(residuals(m)^2)), log = TRUE)), tolerance = 1e-8)
  fit <- suppressWarnings(multiScale_optim(m, inputs, backend = "profile_glm",
    n_cores = 1, verbose = FALSE, profile_control = list(n_starts = 3)))
  expect_equal(predict(fit, d, inputs), predict(fit))
  expect_equal(predict(prepare_prediction(fit)), predict(fit))
})

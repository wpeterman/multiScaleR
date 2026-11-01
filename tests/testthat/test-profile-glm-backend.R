testthat::skip_if_not_installed("MASS")

test_that("NB2 backend returns a compatible conditional model and joint information", {
  skip_if_not_installed("MASS")
  f <- profile_glm_fixture()
  x <- profile_glm_fit(f)
  expect_identical(x$backend, "profile_glm")
  expect_s3_class(x$opt_mod, "negbin")
  expect_error(update(x$opt_mod), "pre-optimization")
  expect_true(x$uncertainty$available)
  expect_true(all(x$profile_fit$start_results$eligible))
  ref <- MASS::glm.nb(y ~ hab + z + group + offset(log(exposure)),
                      data = transform(f$data, hab = x$opt_mod$model$hab))
  expect_equal(as.numeric(logLik(x$opt_mod)), as.numeric(logLik(ref)), tolerance = 1e-7)
  expect_equal(x$opt_mod$theta, ref$theta, tolerance = 1e-5)
  expect_equal(unname(coef(x$opt_mod)), unname(coef(ref)), tolerance = 1e-5)
  expect_equal(unname(predict(x)), unname(fitted(x$opt_mod)), tolerance = 1e-9)
  expect_equal(unname(predict(x, newdata = f$data, kernel_inputs = f$inputs)),
               unname(predict(x)), tolerance = 1e-9)
  V <- vcov(x)
  expect_identical(rownames(V), names(x$uncertainty$par))
  expect_true(all(eigen(V, symmetric = TRUE)$values > 0))
  expect_equal(unname(V[x$uncertainty$profile_index, x$uncertainty$profile_index]),
               unname(solve(x$profile_fit$hessian)), tolerance = 2e-4)
  s <- summary(x)
  expect_identical(unname(attr(s$opt_scale, "interval_method")), "log_wald")
  expect_true(s$opt_scale[1, "2.5%"] < x$scale_est$Mean)
  raw_engine <- x$profile_fit$contract$objective$engine
  raw_engine$standardize <- FALSE
  expect_equal(x$scl_params$sd, setNames(sd(
    .msr_eval_kernel_engine(raw_engine, x$profile_fit$par[1])$raw[, 1]), "hab"))
})

test_that("prediction uses training normalization, offsets, and joint derivatives", {
  f <- profile_glm_fixture()
  x <- profile_glm_fit(f)
  ids <- c(3L, 19L, 2L)
  nd <- f$data[ids, ]
  pred <- predict(x, nd, f$inputs, type = "link", interval = "confidence")
  expect_equal(unname(pred$fit), unname(predict(x, type = "link")[ids]), tolerance = 1e-9)
  one <- predict(x, nd[1, , drop = FALSE], f$inputs, se.fit = TRUE)
  expect_equal(one$fit, predict(x)[ids[1]], tolerance = 1e-9)
  conditional <- predict(x, nd, f$inputs, type = "link", uncertainty = "conditional", se.fit = TRUE)
  expect_true(all(pred$se.fit >= conditional$se.fit - 1e-10))
  response <- predict(x, nd, f$inputs, interval = "confidence")
  expect_equal(response$lower, exp(pred$lower))
  expect_equal(response$upper, exp(pred$upper))
  shifted <- nd
  shifted$exposure <- shifted$exposure * 2
  expect_equal(unname(predict(x, shifted, f$inputs)), 2 * unname(response$fit))
  converted <- f$inputs
  converted$unit_conv <- 1
  converted$binned$bin_rep <- converted$binned$bin_rep * 1000
  expect_equal(predict(x, nd, converted), response$fit, tolerance = 1e-9)
  expect_error(predict(x, nd), "matching smooth kernel_inputs")
  missing_fixed <- f$data
  missing_fixed$z <- NULL
  assign("z", f$data$z, envir = environment(formula(x$opt_mod)))
  expect_error(predict(x, missing_fixed, f$inputs), "missing formula variables: z")
  nd$group <- "unseen"
  expect_error(predict(x, nd, f$inputs), "new level")
  # Independent finite differences of the complete prediction map verify the
  # training centering/spread contribution to scale uncertainty.
  base <- x$uncertainty$par
  fn <- function(par) {
    xx <- x
    xx$uncertainty$par <- par
    xx$profile_fit$par[1] <- par[1]
    unname(predict(xx, f$data[ids, ], f$inputs, type = "link"))
  }
  J <- sapply(seq_along(base), function(j) {
    a <- b <- base
    a[j] <- a[j] + 1e-5
    b[j] <- b[j] - 1e-5
    (fn(a) - fn(b)) / 2e-5
  })
  expect_equal(unname(pred$se.fit), sqrt(rowSums((J %*% vcov(x)) * J)), tolerance = 1e-6)
})

test_that("retained rows and fixed transformed designs have explicit contracts", {
  f <- profile_glm_fixture()
  dat <- f$data
  dat$exposure[c(2, 9)] <- NA
  model <- glm(y ~ group + poly(z, 2) + hab + offset(log(exposure)),
               data = dat, subset = seq_len(nrow(dat)) %% 3 != 0,
               family = MASS::negative.binomial(2), na.action = na.exclude)
  f$model <- model
  x <- suppressWarnings(profile_glm_fit(f))
  rows <- as.integer(rownames(model$model))
  expect_identical(x$profile_fit$contract$rows, rows)
  expect_length(predict(x), length(rows))
  expect_equal(unname(predict(x, f$data[rows, ], f$inputs)),
               unname(predict(x)), tolerance = 1e-9)
  expect_equal(x$opt_mod$model$group, model$model$group)
  bad <- f
  bad$model <- glm(y ~ hab * z, data = f$data, family = poisson())
  expect_s3_class(.msr_profile_contract(bad$model, bad$inputs)$objective$response_design,
    "msr_response_design")
  bad$model <- glm(y ~ log(abs(hab) + 1), data = f$data, family = poisson())
  expect_error(profile_glm_fit(bad), "Unsupported dynamic function: abs")
  bad$model <- glm(y ~ hab + offset(hab), data = f$data, family = poisson())
  expect_error(profile_glm_fit(bad), "offsets")
  bad$model <- glm(y ~ hab, data = f$data, family = poisson(), offset = log(exposure))
  expect_error(profile_glm_fit(bad), "separate offset")
  bad$model <- glm(y ~ hab, data = f$data, family = poisson(), model = FALSE)
  expect_error(profile_glm_fit(bad), "model = TRUE")
  bad$model <- glm(y ~ hab, data = f$data, family = poisson(), weights = rep(.5, nrow(f$data)))
  expect_error(profile_glm_fit(bad), "frequency weights")
  collision <- transform(f$data, h = factor(rep(c("a", "ab"), length.out = nrow(f$data))))
  bad$model <- glm(y ~ h + hab, data = collision, family = poisson())
  expect_error(profile_glm_fit(bad), "Design column names must be distinct")
  bad <- f
  bad$inputs$binned$point_ids[1] <- "unmatched"
  expect_error(profile_glm_fit(bad), "match prepared point")
  expect_error(profile_glm_fit(f, profile_control = list(unknown = 1)), "Unknown")
  expect_error(multiScale_optim(f$model, f$inputs, backend = "refit",
                                profile_control = list(n_starts = 2)), "requires backend")
  expect_error(profile_glm_fit(f, start_strategy = "screen"), "serial full multi-start")
})

test_that("joint analytic scores cover all supported families", {
  f <- profile_glm_fixture(n = 90)
  e <- .msr_prepare_kernel_engine(f$inputs$binned)
  set.seed(11)
  for (family in c("gaussian", "poisson", "binomial", "nbinom2")) {
    y <- switch(family, gaussian = rnorm(90), binomial = rbinom(90, 5, .4), f$data$y)
    obj <- .msr_profile_glm_objective(e, y, family,
      trials = if (family == "binomial") rep(5, 90) else rep(1, 90))
    jp <- .msr_profile_joint(obj, c(log(.3), if (family == "nbinom2") log(2)))
    par <- jp$par + seq(-.05, .05, length.out = length(jp$par))
    numeric <- vapply(seq_along(par), function(j) {
      a <- b <- par
      a[j] <- a[j] + 1e-5
      b[j] <- b[j] - 1e-5
      (jp$fn(a) - jp$fn(b)) / 2e-5
    }, numeric(1))
    expect_equal(unname(jp$gr(par)), numeric, tolerance = 1e-5)
  }
})

test_that("Poisson-limit scores and predictions are stable and uncertainty is withheld", {
  f <- profile_glm_fixture(n = 200, theta = Inf)
  obj <- .msr_profile_contract(f$model, f$inputs)$objective
  for (theta in c(100, 1e4, 1e5, 1e6)) {
    par <- c(log(.25), log(theta))
    a <- b <- par
    a[2] <- a[2] + .001
    b[2] <- b[2] - .001
    numeric <- (obj$fn(a) - obj$fn(b)) / .002
    expect_lt(abs(obj$gr(par)[2] - numeric), 2e-5)
  }
  # Force the Poisson-boundary case with underdispersed counts.
  f$data$y <- round(f$mu)
  f$model <- glm(y ~ hab + z + group + offset(log(exposure)), data = f$data,
                 family = MASS::negative.binomial(2))
  expect_warning(x <- profile_glm_fit(f), "Joint Wald uncertainty unavailable")
  expect_true(x$uncertainty$poisson_limit)
  expect_true(x$uncertainty$boundary)
  expect_false(x$uncertainty$available)
  expect_true(is.na(x$opt_mod$SE.theta))
  expect_true(all(is.finite(predict(x))))
  expect_error(predict(x, se.fit = TRUE), "unavailable")
  expect_error(vcov(x), "unavailable")
  expect_true(all(is.finite(predict(x, se.fit = TRUE, uncertainty = "conditional")$se.fit)))
  expect_true(is.finite(kernel_dist(x)$Mean))
  expect_true(is.na(summary(x)$opt_scale$`2.5%`))
  fp <- f
  fp$model <- glm(y ~ hab + z + group + offset(log(exposure)), data = f$data,
                  family = poisson())
  xp <- profile_glm_fit(fp)
  expect_equal(unname(predict(x)), unname(predict(xp)), tolerance = 2e-3)
})

test_that("likelihood intervals reestimate nuisance parameters and expose search bounds", {
  f <- profile_glm_fixture(n = 120)
  x <- profile_glm_fit(f)
  pr <- profile_sigma_intervals(x)
  expect_true(pr$intervals$`2.5%` < x$scale_est$Mean)
  expect_true(pr$intervals$`97.5%` > x$scale_est$Mean)
  target <- x$optim_results$value + qchisq(.95, 1) / 2
  for (j in 1:2) if (!pr$bound_limited[1, j])
    expect_equal(.msr_profile_at(x, 1, log(pr$intervals[1, j] / 1000)),
                 target, tolerance = 1e-5)
  slice <- profile_sigma(x, sigma_values = c(100, 250, 500), metric = "LL", verbose = FALSE)
  expect_s3_class(slice, "sigma_profile")
  expect_true(all(is.finite(slice$profiles$LL)))
  # Narrow artificial bounds exercise reporting; this is not a refit.
  tight <- x
  tight$profile_fit$lower[1] <- x$profile_fit$par[1] - .001
  tight$profile_fit$upper[1] <- x$profile_fit$par[1] + .001
  tight$profile_fit$starts <- matrix(x$profile_fit$par, 1)
  bt <- profile_sigma_intervals(tight)
  expect_true(all(bt$bound_limited))
  expect_identical(unname(bt$method), "profile-bound-limited")
})

test_that("binary presence responses match GLM likelihood and probability prediction", {
  f <- profile_glm_fixture(n = 400, seed = 6221)
  set.seed(6222)
  f$data$presence <- rbinom(400, 1, plogis(-.4 + 1.2 * f$data$hab))
  f$model <- glm(presence ~ hab, family = binomial(), data = f$data)
  x <- profile_glm_fit(f)
  more <- profile_glm_fit(f, profile_control = list(n_starts = 24L))
  independent <- glm(presence ~ hab, family = binomial(), data = x$opt_mod$model,
                     control = glm.control(epsilon = 1e-11, maxit = 100L))
  expect_equal(unname(predict(x)), unname(fitted(independent)), tolerance = 1e-8)
  expect_equal(-x$optim_results$value, as.numeric(logLik(independent)), tolerance = 1e-7)
  expect_equal(x$optim_results$value, more$optim_results$value, tolerance = 1e-6)
  expect_true(x$uncertainty$available)
  prepared <- prepare_prediction(x)
  for (uncertainty in c("joint", "conditional")) {
    p <- predict(x, interval = "confidence", uncertainty = uncertainty)
    expect_equal(predict(prepared, interval = "confidence", uncertainty = uncertainty), p)
    expect_true(all(is.finite(p$se.fit)))
    expect_true(all(p$lower >= 0 & p$upper <= 1))
  }
})

test_that("grouped binomial predictions are probabilities and frequency weights preserve likelihood", {
  f <- profile_glm_fixture(n = 100)
  set.seed(91)
  f$data$success <- rbinom(100, 8, plogis(.1 + f$data$hab))
  f$model <- glm(cbind(success, 8 - success) ~ hab + z, family = binomial(), data = f$data)
  x <- suppressWarnings(profile_glm_fit(f))
  expect_equal(unname(predict(x)), unname(fitted(x$opt_mod)), tolerance = 1e-8)
  expect_true(all(predict(x) > 0 & predict(x) < 1))
  expect_equal(as.numeric(logLik(x$opt_mod)), -x$optim_results$value, tolerance = 1e-7)
  f$model <- glm(cbind(success, 8 - success) ~ hab + z - 1,
                 family = binomial(), data = f$data, weights = rep(2, 100))
  weighted <- suppressWarnings(profile_glm_fit(f))
  expect_equal(as.numeric(logLik(weighted$opt_mod)), -weighted$optim_results$value, tolerance = 1e-7)
  expect_identical(weighted$opt_mod$df.null, 100L)
  f$model <- glm(y ~ hab, family = MASS::negative.binomial(2), data = f$data,
                 weights = rep(2, 100))
  x <- suppressWarnings(profile_glm_fit(f))
  expect_equal(-x$optim_results$value,
               sum(2 * dnbinom(f$data$y, size = x$opt_mod$theta, mu = predict(x), log = TRUE)),
               tolerance = 1e-8)
})

test_that("kernel_prep inputs work and the refit default remains unchanged", {
  rs <- terra::subset(sim_rast(dim = 60, resolution = 10, user_seed = 7), "cont1")
  names(rs) <- "hab"
  terra::crs(rs) <- "EPSG:32617"
  bounds <- as.vector(terra::ext(rs))
  set.seed(33)
  pts <- terra::vect(cbind(runif(50, bounds[1] + 60, bounds[2] - 60),
                           runif(50, bounds[3] + 60, bounds[4] - 60)),
                     crs = terra::crs(rs))
  inputs <- kernel_prep(pts, rs, max_D = 90, bin = TRUE,
                        nbins = 60, store_cell_data = FALSE, verbose = FALSE)
  dat <- data.frame(inputs$kernel_dat)
  dat$y <- rpois(nrow(dat), exp(.5 + .8 * dat$hab))
  model <- glm(y ~ hab, data = dat, family = poisson())
  a <- suppressWarnings(multiScale_optim(model, inputs, n_cores = 1, verbose = FALSE))
  b <- suppressWarnings(multiScale_optim(model, inputs, n_cores = 1, verbose = FALSE, backend = "refit"))
  expect_equal(a$optim_results, b$optim_results)
  c <- suppressWarnings(multiScale_optim(model, inputs, n_cores = 1, verbose = FALSE,
                                         backend = "profile_glm"))
  expect_equal(unname(predict(c)), unname(predict(c, dat, inputs)), tolerance = 1e-8)
  expect_equal(as.numeric(logLik(c$opt_mod)), -c$optim_results$value, tolerance = 1e-8)
})

test_that("exponential-power new-location derivatives include shape and training scaling", {
  f <- profile_glm_fixture(n = 100, kernel = "expow")
  x <- suppressWarnings(profile_glm_fit(f, profile_control = list(n_starts = 4)))
  idx <- c(4L, 17L)
  nd <- f$data[idx, ]
  d <- .msr_profile_prediction_design(x, nd, f$inputs)
  for (j in 1:2) {
    a <- b <- x
    a$profile_fit$par[j] <- a$profile_fit$par[j] + 1e-5
    b$profile_fit$par[j] <- b$profile_fit$par[j] - 1e-5
    numeric <- (.msr_profile_prediction_design(a, nd, f$inputs)$X[, "hab"] -
                  .msr_profile_prediction_design(b, nd, f$inputs)$X[, "hab"]) / 2e-5
    analytic <- if (j == 1L) d$ds[, 1] else d$dh[, 1]
    expect_equal(unname(analytic), unname(numeric), tolerance = 1e-6)
  }
  expect_equal(unname(predict(x, nd, f$inputs)), unname(predict(x)[idx]), tolerance = 1e-8)
  # Isolate the distance-reporting contract with known, positive information.
  # Extra shape uncertainty must widen radius intervals without moving means.
  x$uncertainty$available <- TRUE
  x$uncertainty$vcov[1:2, 1:2] <- diag(c(.04, .001))
  narrow <- kernel_dist(x)
  x$uncertainty$vcov[2, 2] <- .2
  wide <- kernel_dist(x)
  expect_equal(wide$Mean, narrow$Mean)
  expect_true(wide$`2.5%` < narrow$`2.5%`)
  expect_true(wide$`97.5%` > narrow$`97.5%`)
})

test_that("Gaussian ML dispersion and kernel-only prediction survive serialization", {
  f <- profile_glm_fixture(n = 150)
  set.seed(70)
  f$data$y <- .4 + 1.2 * f$data$hab + rnorm(150, sd = .3)
  f$model <- glm(y ~ hab, data = f$data, family = gaussian())
  x <- profile_glm_fit(f)
  expect_equal(as.numeric(logLik(x$opt_mod)), -x$optim_results$value, tolerance = 1e-8)
  expect_equal(unname(predict(x)), unname(predict(x, kernel_inputs = f$inputs)), tolerance = 1e-9)
  restored <- unserialize(serialize(x, NULL))
  expect_equal(predict(restored, se.fit = TRUE), predict(x, se.fit = TRUE))
  prepared <- prepare_prediction(x, kernel_inputs = f$inputs)
  expect_equal(predict(prepared, se.fit = TRUE), predict(x, se.fit = TRUE))
  expect_equal(predict(unserialize(serialize(prepared, NULL)), interval = "confidence"),
               predict(x, interval = "confidence"))
  expect_identical(tail(rownames(vcov(x)), 1L), "log_residual_sd")
})

test_that("multiple kernel columns retain formula order through fitting and prediction", {
  f <- profile_glm_fixture(n = 150)
  other <- profile_glm_fixture(n = 150, seed = 912)
  f$inputs$binned$vsum$cover <- other$inputs$binned$vsum$hab
  f$inputs$binned$csum$cover <- other$inputs$binned$csum$hab
  engine <- .msr_prepare_kernel_engine(f$inputs$binned)
  values <- .msr_eval_kernel_engine(engine, log(c(.25, .4)))$value
  f$data$cover <- values[, "cover"]
  set.seed(71)
  f$data$y <- rnbinom(150, mu = exp(.5 + f$data$hab + .8 * f$data$cover + .2 * f$data$z),
                      size = 2)
  f$model <- glm(y ~ cover + z + hab, data = f$data, family = MASS::negative.binomial(2))
  reordered <- glm(y ~ hab + z + cover, data = f$data, family = MASS::negative.binomial(2))
  first <- .msr_profile_contract(f$model, f$inputs)$objective
  second <- .msr_profile_contract(reordered, f$inputs)$objective
  expect_equal(first$fn(log(c(.4, .25, 2))), second$fn(log(c(.25, .4, 2))), tolerance = 1e-9)
  x <- suppressWarnings(profile_glm_fit(f))
  expect_identical(rownames(x$scale_est), c("cover", "hab"))
  expect_identical(x$profile_fit$parameter_names, c("log_scale:cover", "log_scale:hab", "log_theta"))
  expect_equal(unname(predict(x)), unname(fitted(x$opt_mod)), tolerance = 1e-9)
  expect_equal(unname(predict(x, f$data[5:7, ], f$inputs)),
               unname(predict(x)[5:7]), tolerance = 1e-9)
})

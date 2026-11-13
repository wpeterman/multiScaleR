radius_grid_fixture <- function() {
  set.seed(148)
  n <- 120L
  distances <- c(.025, rep(c(.075, .125, .175, .225), each = 6L))
  values <- replicate(n, sample(1:3, length(distances), replace = TRUE), simplify = FALSE)
  ids <- sprintf("alternative_%03d", seq_len(n))
  names(values) <- ids
  raw <- lapply(values, function(x) cbind(cover = x))
  radius <- c(50, 100, 150, 200)
  # Independent composition calculation at each map-unit radius.
  entropy <- vapply(radius, function(r) vapply(values, function(x) {
    p <- prop.table(table(x[distances <= r / 1000]))
    -sum(p * log(p))
  }, numeric(1)), numeric(n))
  z <- as.numeric(scale(entropy[, 3L]))
  step <- rep(seq_len(n / 3L), each = 3L)
  used <- integer(n)
  for (idx in split(seq_len(n), step))
    used[sample(idx, 1L, prob = exp(.7 * z[idx]))] <- 1L
  data <- data.frame(used = used, step = step, diversity = as.numeric(scale(entropy[, 2L])),
    y = 2 * z + rnorm(n), row.names = ids)
  inputs <- structure(list(kernel_dat = data["diversity"], raw_cov = raw,
    d_list = stats::setNames(rep(list(distances), n), ids),
    scale_vars = msr_vars(diversity = landscape_var("cover", metric = "shdi")),
    kernel = "gaussian", sigma = .1, shape = NULL, min_D = 40, max_D = 240,
    unit_conv = 1000, cell_size = c(20, 20), resolution = 20, n_cols = 50L,
    n_covs = 1L, binned = NULL, cell_data_stored = TRUE,
    scl_params = list(mean = c(diversity = mean(entropy[, 2L])),
      sd = c(diversity = sd(entropy[, 2L])))), class = "multiScaleR_data")
  list(data = data, inputs = inputs, radius = radius, entropy = entropy)
}

test_that("a declared radius profile matches separate native likelihoods", {
  f <- radius_grid_fixture()
  d <- f$data
  initial <- glm(y ~ diversity, data = d)
  fit <- multiScale_optim(initial, f$inputs, radius_grid = rev(f$radius), verbose = FALSE)
  profile <- diagnostics(fit)$radius_profile
  ref <- vapply(2:4, function(i) {
    d$diversity <- as.numeric(scale(f$entropy[, i]))
    as.numeric(logLik(glm(y ~ diversity, data = d)))
  }, numeric(1))
  expect_equal(profile$radius, f$radius)
  expect_identical(profile$status, c("failed", rep("scored", 3L)))
  expect_true(nzchar(profile$message[1L]))
  expect_equal(profile$log_likelihood[-1L], ref, tolerance = 1e-10)
  expect_equal(fit$scale_est["diversity", "Mean"], f$radius[-1L][which.max(ref)])
  expect_equal(fit$optim_results$value, -max(ref), tolerance = 1e-10)
  expect_equal(as.numeric(predict(fit, f$data, f$inputs)), unname(fitted(fit$opt_mod)))
  expect_true(is.na(summary(fit)$opt_scale["diversity", "SE"]))
  expect_error(summary(fit, profile = TRUE), "radius_profile")
  expect_error(multiScale_optim(initial, f$inputs, radius_grid = c(50, 250)), "within min_D")
  expect_error(multiScale_optim(initial, f$inputs, backend = "profile_glm",
    radius_grid = f$radius), "requires backend")
  err <- tryCatch(multiScale_optim(initial, f$inputs, radius_grid = c(50, 60)),
    msr_radius_profile_error = identity)
  expect_s3_class(err, "msr_radius_profile_error")
  expect_identical(err$radius_profile$status, c("failed", "failed"))
  tie <- multiScale_optim(initial, f$inputs, radius_grid = c(105, 100, 105), verbose = FALSE)
  expect_equal(tie$scale_est["diversity", "Mean"], 100)
  bad_inputs <- f$inputs
  names(bad_inputs$d_list) <- rev(names(bad_inputs$d_list))
  expect_error(multiScale_optim(initial, bad_inputs, radius_grid = f$radius), "matching model row IDs")
})

test_that("choice CV profiles the hard radius using training sets only", {
  skip_if_not_installed("survival")
  f <- radius_grid_fixture()
  d <- f$data
  folds <- grouped_cv_folds(d$step, 3L, rownames(d))[1L]
  strata <- survival::strata
  coxph <- survival::coxph
  Surv <- survival::Surv
  form <- used ~ diversity + strata(step)
  run <- function(data, radius = f$radius) cross_validate_clogit(form, data, folds,
    f$inputs, backend = "refit", optim_args = list(radius_grid = radius),
    keep_fits = TRUE, keep_selection = TRUE)
  expect_error(cross_validate_clogit(form, d, folds, backend = "refit",
    optim_args = list(radius_grid = f$radius)), "requires kernel_inputs")
  cv <- run(d)
  expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = " "))
  fit <- cv$fits[[1L]]
  profile <- diagnostics(fit)$radius_profile
  analysis <- folds[[1L]]$analysis
  assessment <- folds[[1L]]$assessment
  ref <- vapply(2:4, function(i) {
    train <- d[analysis, , drop = FALSE]
    train$diversity <- as.numeric(scale(f$entropy[analysis, i]))
    as.numeric(logLik(survival::clogit(form, data = train, model = TRUE)))
  }, numeric(1))
  expect_equal(profile$log_likelihood[-1L], ref, tolerance = 1e-9)
  expect_equal(profile$radius[profile$selected], f$radius[-1L][which.max(ref)])
  best <- which(profile$selected)
  test <- d[assessment, , drop = FALSE]
  train_x <- f$entropy[analysis, best]
  test$diversity <- (f$entropy[assessment, best] - mean(train_x)) / sd(train_x)
  eta <- unname(coef(fit$opt_mod)) * test$diversity
  p <- exp(eta) / ave(exp(eta), test$step, FUN = sum)
  expect_equal(cv$predictions$prediction, unname(p), tolerance = 1e-9)
  expect_equal(cv$selection_records[[1L]]$data$K, 2L)
  expect_equal(cv$selection_records[[1L]]$data$nobs, length(analysis) / 3L)

  d[assessment, "used"] <- ave(d[assessment, "used"], d[assessment, "step"],
    FUN = function(x) x[c(2L, 3L, 1L)])
  changed <- run(d)
  expect_identical(diagnostics(changed$fits[[1L]])$radius_profile, profile)
  expect_equal(changed$predictions$prediction, cv$predictions$prediction)

  failed <- run(d, c(50, 60))
  expect_false(failed$complete)
  expect_identical(failed$errors[[1L]]$radius_profile$status, c("failed", "failed"))
})

test_that("Gaussian CV selects the radius and keeps training normalization", {
  f <- radius_grid_fixture()
  d <- f$data
  fold <- grouped_cv_folds(d$step, 3L, rownames(d))[1L]
  cv <- cross_validate(y ~ diversity, d, f$inputs, fold, family = "gaussian",
    backend = "refit", optim_args = list(radius_grid = f$radius), keep_fits = TRUE)
  expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = " "))
  fit <- cv$fits[[1L]]
  tab <- diagnostics(fit)$radius_profile
  train <- fold[[1L]]$analysis
  assessment <- fold[[1L]]$assessment
  best <- which(tab$selected)
  diversity <- f$entropy[train, best]
  native <- lm(y ~ diversity, data = data.frame(y = d[train, "y"], diversity = diversity))
  expected <- predict(native, data.frame(diversity = f$entropy[assessment, best]))
  expect_equal(cv$predictions$prediction, unname(expected), tolerance = 1e-10)
})

mixed_radius_grid_fixture <- function(kernel = "gaussian") {
  f <- radius_grid_fixture()
  d <- f$inputs$d_list[[1L]]
  set.seed(205)
  n <- 84L
  f$data <- f$data[seq_len(n), ]
  ids <- rownames(f$data)
  f$inputs$d_list <- f$inputs$d_list[ids]
  habitat <- matrix(rnorm(n * length(d)), nrow = n)
  elev <- matrix(rnorm(n * length(d)), nrow = n)
  f$inputs$raw_cov <- stats::setNames(lapply(seq_len(n), function(i)
    cbind(cover = f$inputs$raw_cov[[i]][, "cover"], habitat = habitat[i, ], elev = elev[i, ])), ids)
  f$inputs$kernel <- kernel
  # Preparation order intentionally differs from the formula and grid names.
  f$inputs$scale_vars <- msr_vars(rough = surface_var("elev", "sq"),
    diversity = landscape_var("cover", "shdi"), hab = kernel_var("habitat"))
  f$inputs$sigma <- c(.1, .15, .13)
  f$inputs$shape <- if (kernel == "expow") rep(2, 3L) else NULL
  f$inputs$n_covs <- 3L
  # Distinct distances occupy separate bins, preserving the raw-cell likelihood.
  f$inputs$binned <- .msr_build_kernel_bins(f$inputs$d_list, f$inputs$raw_cov,
    "hab", "habitat", nbins = 16L, point_ids = ids)
  entropy <- f$entropy[seq_len(n), , drop = FALSE]
  rough <- function(radius) apply(elev[, d <= radius / 1000, drop = FALSE], 1L, sd)
  mean_habitat <- function(sigma, shape = 2) {
    weights <- if (kernel == "expow") exp(-(d / sigma)^shape) else exp(-.5 * (d / sigma)^2)
    drop(habitat %*% (weights / sum(weights)))
  }
  prepared <- cbind(rough = rough(100), diversity = entropy[, 3L], hab = mean_habitat(.13))
  f$inputs$kernel_dat <- as.data.frame(scale(prepared))
  rownames(f$inputs$kernel_dat) <- ids
  f$inputs$scl_params <- list(mean = colMeans(prepared), sd = apply(prepared, 2L, sd))
  f$data <- cbind(f$data[setdiff(names(f$data), "diversity")], f$inputs$kernel_dat)
  signal <- .6 * scale(entropy[, 3L]) - .5 * scale(rough(200)) + 1.4 * scale(mean_habitat(.086, 2.7))
  f$data$y <- as.numeric(signal + rnorm(n, sd = .35))
  f$data$used <- 0L
  for (i in split(seq_len(n), f$data$step))
    f$data$used[sample(i, 1L, prob = exp(signal[i]))] <- 1L
  list(data = f$data, inputs = f$inputs, entropy = entropy, rough = rough,
    mean_habitat = mean_habitat, radius = list(rough = c(100, 200), diversity = c(100, 150)))
}

test_that("joint radii profile the Gaussian bandwidth against independent native fits", {
  f <- mixed_radius_grid_fixture()
  d <- f$data
  fit <- multiScale_optim(lm(y ~ hab + diversity + rough, data = d), f$inputs,
    radius_grid = f$radius, verbose = FALSE)
  tab <- diagnostics(fit)$radius_profile
  expect_equal(tab[["radius:diversity"]], rep(c(100, 150), each = 2L))
  expect_equal(tab[["radius:rough"]], rep(c(100, 200), 2L))
  expect_true(all(tab$status == "scored"), info = paste(tab$message, collapse = " "))
  reference <- lapply(seq_len(nrow(tab)), function(i) {
    train <- d
    train$diversity <- f$entropy[, match(tab[["radius:diversity"]][i], c(50, 100, 150, 200))]
    train$rough <- f$rough(tab[["radius:rough"]][i])
    optimize(function(log_sigma) {
      train$hab <- f$mean_habitat(exp(log_sigma))
      as.numeric(logLik(lm(y ~ hab + diversity + rough, data = train)))
    }, log(c(.04, .24)), maximum = TRUE, tol = 1e-10)
  })
  expect_equal(tab$log_likelihood, vapply(reference, function(x) x$objective, 0), tolerance = 1e-7)
  expect_equal(tab[["sigma:hab"]], 1000 * exp(vapply(reference, function(x) x$maximum, 0)),
    tolerance = 1e-4)
  expect_true(diff(range(tab[["sigma:hab"]])) > .5)
  expect_identical(which(tab$selected), which.max(tab$log_likelihood))
  expect_equal(as.numeric(predict(fit, d, f$inputs)), unname(fitted(fit$opt_mod)), tolerance = 1e-10)
  expect_true(is.finite(fit$scale_est["hab", "SE"]))
  expect_true(all(is.na(fit$scale_est[c("rough", "diversity"), "SE"])))
  expect_identical(diagnostics(fit)$refit_hessian$uncertainty,
    "outer_refit_wald_conditional_on_selected_hard_radii")
  expect_equal(.msr_selection_parameter_count(fit), 8)
  expect_error(multiScale_optim(lm(y ~ hab + diversity + rough, data = d),
    f$inputs, radius_grid = f$radius["rough"], verbose = FALSE), "every optimized hard-radius predictor")
})

test_that("choice CV trains joint radius, bandwidth and shape fits without assessment outcomes", {
  skip_if_not_installed("survival")
  f <- mixed_radius_grid_fixture("expow")
  d <- f$data
  folds <- grouped_cv_folds(d$step, 3L, rownames(d))[1L]
  strata <- survival::strata
  coxph <- survival::coxph
  Surv <- survival::Surv
  form <- used ~ hab + diversity + rough + strata(step)
  run <- function(data) cross_validate_clogit(form, data, folds, f$inputs,
    backend = "refit", optim_args = list(radius_grid = f$radius), keep_fits = TRUE, keep_selection = TRUE)
  cv <- run(d)
  expect_true(cv$complete, info = paste(unlist(cv$errors), collapse = " "))
  fit <- cv$fits[[1L]]
  tab <- diagnostics(fit)$radius_profile
  expect_true(all(tab$status == "scored"), info = paste(tab$message, collapse = " "))
  train <- match(folds[[1L]]$analysis, rownames(d))
  held <- match(folds[[1L]]$assessment, rownames(d))
  for (i in which(tab$status == "scored")) {
    train_d <- d[train, ]
    train_d$hab <- f$mean_habitat(tab[["sigma:hab"]][i] / 1000, tab[["shape:hab"]][i])[train]
    train_d$diversity <- f$entropy[train, match(tab[["radius:diversity"]][i], c(50, 100, 150, 200))]
    train_d$rough <- f$rough(tab[["radius:rough"]][i])[train]
    reference <- survival::clogit(form, data = train_d, model = TRUE)
    expect_equal(tab$log_likelihood[i], as.numeric(logLik(reference)), tolerance = 1e-8)
  }
  expect_true(all(is.finite(tab[["shape:hab"]])))
  expect_true(all(is.na(fit$shape_est[c("rough", "diversity"), "Mean"])))
  i <- which(tab$selected)
  nd <- d[held, ]
  raw <- cbind(hab = f$mean_habitat(tab[["sigma:hab"]][i] / 1000, tab[["shape:hab"]][i]),
    diversity = f$entropy[, match(tab[["radius:diversity"]][i], c(50, 100, 150, 200))],
    rough = f$rough(tab[["radius:rough"]][i]))
  nd[colnames(raw)] <- sweep(sweep(raw[held, , drop = FALSE], 2L,
    colMeans(raw[train, , drop = FALSE]), "-"), 2L, apply(raw[train, , drop = FALSE], 2L, sd), "/")
  eta <- drop(as.matrix(nd[names(coef(fit$opt_mod))]) %*% coef(fit$opt_mod))
  prob <- exp(eta) / ave(exp(eta), nd$step, FUN = sum)
  expect_equal(cv$predictions$prediction, unname(prob), tolerance = 1e-9)
  expect_equal(cv$selection_records[[1L]]$data$K, 7L)
  expect_equal(cv$selection_records[[1L]]$data$nobs, length(train) / 3L)
  d$used[held] <- ave(d$used[held], d$step[held], FUN = function(x) x[c(2L, 3L, 1L)])
  changed <- run(d)
  expect_identical(diagnostics(changed$fits[[1L]])$radius_profile, tab)
  expect_equal(changed$predictions$prediction, cv$predictions$prediction)
})

test_that("fixed kernels profile every estimated cutoff as a hard radius", {
  f <- mixed_radius_grid_fixture()
  d <- f$data
  inputs <- f$inputs
  inputs$kernel <- "fixed"
  inputs$binned <- NULL
  fit <- multiScale_optim(lm(y ~ hab + rough, data = d), inputs,
    radius_grid = list(rough = c(100, 200), hab = c(70, 150)), verbose = FALSE)
  tab <- diagnostics(fit)$radius_profile
  expect_true(all(tab$status == "scored"), info = paste(tab$message, collapse = " "))
  expect_false("sigma:hab" %in% names(tab))
  expect_true(all(is.na(fit$scale_est$SE)))
  expected <- vapply(seq_len(nrow(tab)), function(i) {
    d$rough <- f$rough(tab[["radius:rough"]][i])
    d$hab <- vapply(seq_len(nrow(d)), function(j)
      mean(inputs$raw_cov[[j]][inputs$d_list[[j]] < tab[["radius:hab"]][i] / 1000,
        "habitat"]), numeric(1))
    as.numeric(logLik(lm(y ~ hab + rough, data = d)))
  }, numeric(1))
  expect_equal(tab$log_likelihood, expected, tolerance = 1e-8)
})

test_that("named model data expressions preserve the native response audit", {
  f <- radius_grid_fixture()
  f$data$y <- drop(scale(f$data$y))
  initial <- lm(y ~ diversity, data = f$data)
  expect_false(is.null(attr(model.response(initial$model), "scaled:center")))
  fit <- multiScale_optim(initial, f$inputs, radius_grid = f$radius, verbose = FALSE)
  tab <- diagnostics(fit)$radius_profile
  expect_identical(tab$status, c("failed", rep("scored", 3L)))
  expect_identical(rownames(fit$opt_mod$model), rownames(initial$model))
  expect_identical(as.numeric(model.response(fit$opt_mod$model)),
    as.numeric(model.response(initial$model)))
})

test_that("a flat conditional bandwidth likelihood does not supply a Wald SE", {
  f <- mixed_radius_grid_fixture()
  inputs <- f$inputs
  for (i in seq_along(inputs$raw_cov))
    inputs$raw_cov[[i]][, "habitat"] <- f$data$hab[i]
  inputs$binned <- .msr_build_kernel_bins(inputs$d_list, inputs$raw_cov,
    "hab", "habitat", nbins = 16L)
  d <- f$data
  fit <- multiScale_optim(lm(y ~ hab + diversity + rough, data = d), inputs,
    radius_grid = f$radius, verbose = FALSE)
  diag <- diagnostics(fit)
  expect_true(all(diag$radius_profile$status == "scored"))
  expect_true(is.finite(fit$scale_est["hab", "Mean"]))
  expect_true(is.na(fit$scale_est["hab", "SE"]))
  expect_identical(diag$refit_hessian$status, "invalid_hessian")
})

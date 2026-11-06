testthat::skip_if_not_installed("MASS")

# Nine bounded scale fits within this test file. Raw cells and responses
# reuse deterministic CV fixtures; complete-set choice reuses the established seed.
cv_formula_spatial_inputs <- function(sq = FALSE) {
  i <- 1:72; d <- seq(.025, .99, length.out = 24); ids <- paste0("formulaSpatial", i)
  values <- outer(sin(i*.73), sin(7*d)) + outer(cos(i*.41), cos(13*d)) +
    outer(sin(i*.19), d) + .1*outer(sin(i*.37), sin(19*d))
  raw <- lapply(i, function(j) matrix(values[j, ], ncol = 1, dimnames = list(NULL, "terrain")))
  distances <- lapply(i, function(j) d*(1 - (j %% 5)*.005)); names(raw) <- names(distances) <- ids
  specs <- if (sq) msr_vars(hab = surface_var("terrain", "sq", weighted = TRUE)) else msr_vars(hab = kernel_var("terrain"))
  inputs <- structure(list(kernel_dat = data.frame(hab = rep(999, length(ids)), row.names = ids), raw_cov = raw,
    d_list = distances, scale_vars = specs, kernel = "gaussian", sigma = .3, shape = NULL,
    min_D = 35, max_D = 800, unit_conv = 800, n_covs = 1L, resolution = 20, n_cols = 40L,
    binned = if (sq) NULL else .msr_build_kernel_bins(distances, raw, "hab", "terrain", 32L, ids),
    scl_params = list(mean = c(hab = 999), sd = c(hab = 100)), cell_data_stored = TRUE), class = "multiScaleR_data")
  h <- vapply(i, function(j) .msr_eval_scale_vars(distances[[j]], raw[[j]], specs, .3, NULL,
    "gaussian", 800, 20, 40L)[1L], 0)
  h <- as.numeric(scale(h))
  data <- data.frame(hab = 999, z = sin(i*.31), group = factor(rep(c("a", "b"), 36)),
    exposure = 1 + (i %% 5)/10, row.names = ids)
  eta <- .9 + .55*h + .1*data$z + .1*(data$group == "b") + log(data$exposure)
  data$count <- round(exp(eta)*rep(c(.2, .5, 1, 2.5), 18))
  data$positive <- exp(.4 + .7*h + .15*data$z + log(data$exposure))*rep(c(.4, .8, 1.5, 2), 18)
  data$continuous <- .4 + .7*h + .15*data$z + .2*cos(i*.91)
  data$binary <- as.numeric(i %% 5 < 2)
  list(data = data, inputs = inputs, folds = grouped_cv_folds(rep(letters[1:12], each = 6), 3, ids))
}

# Independent raw-cell reduction, using the retained training bin reference and
# physical fitted scale. This oracle does not call the prediction evaluator.
cv_formula_spatial_data <- function(fixture, cv) {
  m <- cv$fits[[1]]; f <- fixture; sigma <- m$scale_est$Mean[1]/f$inputs$unit_conv
  reference <- cv$training[[1]]$bin_reference; ids <- rownames(f$data)
  raw <- vapply(ids, function(id) {
    distance <- f$inputs$d_list[[id]]; v <- f$inputs$raw_cov[[id]][, 1]
    if (!is.null(reference)) distance <- reference$bin_rep[.bincode(distance, reference$breaks, include.lowest = TRUE)]
    w <- exp(-distance^2/(2*sigma^2)); mean <- sum(w*v)/sum(w)
    if (isTRUE(f$inputs$scale_vars$metric == "sq")) sqrt(sum(w*(v - mean)^2)/sum(w)) else mean
  }, 0)
  analysis <- cv$folds[[1]]$analysis; center <- mean(raw[analysis]); spread <- sd(raw[analysis])
  expect_equal(unname(m$scl_params$mean["hab"]), center, tolerance = 1e-10)
  expect_equal(unname(m$scl_params$sd["hab"]), spread, tolerance = 1e-10)
  d <- f$data; d$hab <- (raw - center)/spread
  d$scaled_z <- (d$z - mean(d[analysis, "z"]))/sd(d[analysis, "z"])
  d$squared_z <- d$scaled_z^2; d$squared_hab <- d$hab^2
  d
}

test_that("spatial CV repairs final training metadata without changing native fits or gates", {
  skip_if_not_installed("MASS")
  planned <- 0L
  cases <- list(gaussian = list(family = "gaussian", backend = "profile_glm"),
    poisson = list(family = "poisson", backend = "profile_glm"),
    binomial = list(family = "binomial", backend = "profile_glm"),
    nb = list(family = "nbinom2", backend = "profile_glm"),
    gamma = list(family = "Gamma", backend = "refit", dispersion = .4),
    refit = list(family = "gaussian", backend = "refit"),
    sq = list(family = "gaussian", backend = "profile_glm"))
  for (label in names(cases)) {
    case <- cases[[label]]; f <- cv_formula_spatial_inputs(label == "sq")
    f$data$y <- f$data[[switch(case$family, gaussian = "continuous", poisson = "count",
      binomial = "binary", nbinom2 = "count", Gamma = "positive")]]
    held <- f$folds[[1]]$assessment; f$data[held, "z"] <- f$data[held, "z"] + 3
    form <- if (label == "refit") y ~ I(base::scale(hab)^2) + group + I(base::scale(z)^2) +
      offset(base::scale(z)) + offset(log(exposure)) else y ~ hab*group + I(base::scale(z)^2) +
      offset(base::scale(z)) + offset(log(exposure))
    args <- if (case$backend == "profile_glm") list(profile_control = list(starts = matrix(c(.08, .3, .65), ncol = 1))) else list()
    planned <- planned + 1L; stopifnot(planned <= 7L)
    cv <- do.call(cross_validate, c(list(formula = form, data = f$data, kernel_inputs = f$inputs,
      folds = f$folds[1], optim_args = args, keep_fits = TRUE, keep_selection = TRUE), case))
    expect_true(cv$complete, info = paste(label, paste(unlist(cv$errors), collapse = "; ")))
    if (!cv$complete) next
    d <- cv_formula_spatial_data(f, cv); fit <- cv$fits[[1]]; native <- fit$opt_mod
    manual <- if (label == "refit") y ~ squared_hab + group + squared_z + offset(scaled_z) + offset(log(exposure)) else
      y ~ hab*group + squared_z + offset(scaled_z) + offset(log(exposure))
    # Match the native conditional fit's retained family/theta and IRLS controls.
    # A second glm.nb optimization estimates another theta, a different oracle.
    fam <- native$family
    if (case$family == "nbinom2") fam$dispersion <- 1
    ref <- glm(manual, data = d[cv$folds[[1]]$analysis, ], family = fam,
      control = if (case$backend == "profile_glm") glm.control(epsilon = 1e-11, maxit = 100L) else glm.control())
    expected <- as.numeric(predict(ref, d[held, ], type = "response"))
    expect_equal(cv$predictions$prediction, expected, tolerance = 1e-7, info = label)
    expect_equal(unname(coef(native)), unname(coef(ref)), tolerance = 1e-7, info = label)
    expect_equal(unname(vcov(native)), unname(vcov(ref)), tolerance = 1e-7, info = label)
    expect_equal(as.numeric(logLik(native)), as.numeric(logLik(ref)), tolerance = 1e-7, info = label)
    nuisance <- as.integer(case$family == "nbinom2")
    expect_equal(attr(logLik(native), "df"), attr(logLik(ref), "df") + nuisance)
    expect_equal(cv$selection_records[[1]]$data$K, attr(logLik(ref), "df") + nuisance + 1L)
    expect_identical(cv$settings$parameter_uncertainty, "not_integrated")
    # Idempotent metadata repair cannot change native or ecological uncertainty.
    training_inputs <- .msr_cv_inputs(f$inputs, cv$folds[[1]]$analysis, "hab")
    again <- .msr_cv_prediction_metadata(fit, f$data[cv$folds[[1]]$analysis, ], training_inputs)
    expect_identical(coef(again$opt_mod), coef(native)); expect_identical(vcov(again$opt_mod), vcov(native))
    expect_identical(again$scale_est, fit$scale_est); expect_identical(again$diagnostics, fit$diagnostics)
    expect_identical(again$uncertainty, fit$uncertainty)
    if (case$backend == "profile_glm") {
      expect_gte(nrow(fit$profile_fit$start_results), 3L)
      expect_true(any(fit$profile_fit$start_results$eligible))
    }
  }
  expect_identical(planned, 7L)
})

test_that("complete-set choice retains learned offsets and expressions in ordinary and spatial CV", {
  skip_if_not_installed("survival")
  strata <- survival::strata; coxph <- survival::coxph
  # The existing profile fixture and seed are reused without any search/tuning.
  base <- profile_glm_fixture(); b <- base$inputs$binned; counts <- b$csum$hab
  ids <- paste0("formulaChoice", seq_len(nrow(base$data)))
  raw <- lapply(seq_along(ids), function(i) matrix(rep(b$vsum$hab[i, ]/counts[i, ], counts[i, ]), ncol = 1,
    dimnames = list(NULL, "hab")))
  distances <- lapply(seq_along(ids), function(i) rep(b$bin_rep, counts[i, ])); names(raw) <- names(distances) <- ids
  inputs <- structure(list(kernel_dat = data.frame(hab = base$data$hab, row.names = ids), raw_cov = raw,
    d_list = distances, scale_vars = msr_vars(hab = kernel_var("hab")), kernel = "gaussian", sigma = .25,
    shape = NULL, min_D = 50, max_D = 1000, unit_conv = 1000, cell_size = c(20, 20), resolution = 20,
    n_cols = 50L, n_covs = 1L, binned = .msr_build_kernel_bins(distances, raw, "hab", "hab", 32L, ids),
    cell_data_stored = TRUE, scl_params = list(mean = c(hab = 999), sd = c(hab = 100))), class = "multiScaleR_data")
  d <- base$data; rownames(d) <- ids; d$case_ <- rep(c(1, 0, 0), nrow(d)/3)
  d$choice <- rep(seq_len(nrow(d)/3), each = 3); folds <- grouped_cv_folds(d$choice, 3, ids)[1]
  d[folds[[1]]$assessment, "z"] <- d[folds[[1]]$assessment, "z"] + 3
  form <- case_ ~ hab + group + I(base::scale(z)^2) + offset(base::scale(z)) + strata(choice)
  for (route in c("ordinary_exact", "ordinary_efron", "conditional_logit", "refit")) {
    spatial <- route %in% c("conditional_logit", "refit")
    args <- if (route == "conditional_logit") list(profile_control = list(starts = matrix(c(.08, .25, .6), ncol = 1))) else list()
    cv <- cross_validate_clogit(form, d, folds, if (spatial) inputs else NULL,
      method = if (route == "ordinary_efron") "efron" else "exact",
      backend = if (spatial) route else "conditional_logit", optim_args = args,
      keep_fits = TRUE, keep_selection = TRUE)
    expect_true(cv$complete, info = paste(route, paste(unlist(cv$errors), collapse = "; ")))
    if (!cv$complete) next
    training <- d; fit <- cv$fits[[1]]
    if (spatial) training <- cv_formula_spatial_data(list(data = d, inputs = inputs), cv) else {
      training$scaled_z <- (d$z - mean(d[folds[[1]]$analysis, "z"]))/sd(d[folds[[1]]$analysis, "z"])
      training$squared_z <- training$scaled_z^2
    }
    native <- if (spatial) fit$opt_mod else fit
    ref <- survival::clogit(case_ ~ hab + group + squared_z + offset(scaled_z) + strata(choice),
      data = training[folds[[1]]$analysis, ], method = if (route == "ordinary_efron") "efron" else "exact", model = TRUE)
    assess <- training[folds[[1]]$assessment, ]; assess$choice <- training[folds[[1]]$analysis[1], "choice"]
    eta <- as.numeric(predict(ref, assess, type = "lp", reference = "zero"))
    expected <- numeric(length(eta))
    for (g in split(seq_along(eta), cv$predictions$choice)) expected[g] <- exp(eta[g] - max(eta[g]))/sum(exp(eta[g] - max(eta[g])))
    expect_equal(cv$predictions$prediction, expected, tolerance = 1e-7, info = route)
    expect_equal(unname(coef(native)), unname(coef(ref)), tolerance = 1e-7, info = route)
    expect_equal(unname(vcov(native)), unname(vcov(ref)), tolerance = 1e-7, info = route)
    expect_equal(as.numeric(logLik(native)), as.numeric(logLik(ref)), tolerance = 1e-7, info = route)
    expect_equal(cv$selection_records[[1]]$data$K, attr(logLik(ref), "df") + as.integer(spatial))
  }
})

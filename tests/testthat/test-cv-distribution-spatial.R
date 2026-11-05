testthat::skip_if_not_installed("MASS")

# Deterministic raw-cell plumbing fixture. No draws, seed search or recovery study.
cv_distribution_spatial_inputs <- function(kernel = "gaussian", sq = FALSE) {
  i <- 1:72; d <- seq(.025, .99, length.out = 24); ids <- paste0("distSpatial", i)
  values <- outer(sin(i*.73), sin(7*d)) + outer(cos(i*.41), cos(13*d)) +
    outer(sin(i*.19), d) + .1*outer(sin(i*.37), sin(19*d))
  raw <- lapply(i, function(j) matrix(values[j, ], ncol = 1, dimnames = list(NULL, "terrain")))
  distances <- lapply(i, function(j) d*(1 - (j %% 5)*.005))
  names(raw) <- names(distances) <- ids
  specs <- if (sq) msr_vars(hab = surface_var("terrain", "sq", weighted = TRUE)) else msr_vars(hab = kernel_var("terrain"))
  shape <- if (kernel == "expow") 1.8 else NULL
  inputs <- structure(list(kernel_dat = data.frame(hab = rep(999, 72), row.names = ids),
    raw_cov = raw, d_list = distances, scale_vars = specs, kernel = kernel, sigma = .3, shape = shape,
    min_D = 35, max_D = 800, unit_conv = 800, n_covs = 1L, resolution = 20, n_cols = 40L,
    binned = if (sq) NULL else .msr_build_kernel_bins(distances, raw, "hab", "terrain", 32L, ids),
    scl_params = list(mean = c(hab = 999), sd = c(hab = 100)), cell_data_stored = TRUE,
    raster_group = stats::setNames(rep(c("mapA", "mapB"), 36), ids)), class = "multiScaleR_data")
  h <- vapply(i, function(j) .msr_eval_scale_vars(distances[[j]], raw[[j]], specs, .3, shape,
    kernel, 800, 20, 40L)[1L], 0)
  h <- as.numeric(scale(h))
  data <- data.frame(hab = 999, z = sin(i*.31), group = factor(rep(c("a", "b"), 36)),
    exposure = 1 + (i %% 5)/10, row.names = ids)
  eta <- .9 + .55*h + .1*data$z + .1*(data$group == "b") + log(data$exposure)
  data$nb <- round(exp(eta)*rep(c(.2, .5, 1, 2.5), 18))
  data$positive <- exp(.4 + .7*h + .15*data$z + log(data$exposure))*rep(c(.4, .8, 1.5, 2), 18)
  list(data = data, inputs = inputs, mu = exp(eta), folds = grouped_cv_folds(rep(letters[1:12], each = 6), 3, ids))
}

cv_distribution_spatial_cache <- new.env(parent = emptyenv())
cv_distribution_spatial_fixture <- function() {
  if (exists("failure", cv_distribution_spatial_cache, inherits = FALSE)) stop(cv_distribution_spatial_cache$failure)
  if (!exists("value", cv_distribution_spatial_cache, inherits = FALSE)) {
    cv_distribution_spatial_cache$value <- tryCatch({
      g <- cv_distribution_spatial_inputs(); planned <- 0L; cases <- list()
      run <- function(name, f = g, family = "nbinom2", backend = "profile_glm", theta = NULL, dispersion = NULL) {
        planned <<- planned + 1L; stopifnot(planned <= 17L)
        form <- if (family == "nbinom2") nb ~ hab*group + I(z^2) + offset(log(exposure)) else
          positive ~ hab*group + I(z^2) + offset(log(exposure))
        args <- if (isTRUE(backend %in% c("profile_glm", "joint_nb2"))) list(profile_control = list(
          starts = if (f$inputs$kernel == "expow") cbind(c(.08, .3, .65), c(1.2, 1.8, 2.5)) else matrix(c(.08, .3, .65), ncol = 1))) else list()
        arguments <- list(formula = form, data = f$data, kernel_inputs = f$inputs, folds = f$folds[1L], family = family,
          optim_args = args, keep_fits = TRUE, keep_selection = TRUE, theta = theta, dispersion = dispersion)
        if (!is.null(backend)) arguments$backend <- backend
        cases[[name]] <<- list(cv = do.call(cross_validate, arguments), fixture = f)
      }
      run("nb_profile")
      for (kernel in c("exp", "expow")) run(paste0("nb_", kernel), cv_distribution_spatial_inputs(kernel))
      run("nb_sq", cv_distribution_spatial_inputs(sq = TRUE))
      run("nb_refit", backend = "refit")
      run("nb_known", backend = NULL, theta = 3)
      run("gamma", family = "Gamma", backend = "refit")
      run("gamma_known", family = "Gamma", backend = NULL, dispersion = .4)
      run("gamma_sq", cv_distribution_spatial_inputs(sq = TRUE), "Gamma", "refit")
      run("gamma_expow", cv_distribution_spatial_inputs("expow"), "Gamma", "refit")
      run("gamma_exp", cv_distribution_spatial_inputs("exp"), "Gamma", "refit")
      run("nb_fixed", cv_distribution_spatial_inputs("fixed"), backend = "refit")
      run("gamma_fixed", cv_distribution_spatial_inputs("fixed"), "Gamma", "refit")
      run("nb_joint", backend = "joint_nb2")
      limit <- g; limit$data$nb <- round(g$mu)
      run("nb_poisson_limit", limit, backend = "joint_nb2")
      changed <- g; held <- g$folds[[1]]$assessment
      changed$data[held, "nb"] <- changed$data[held, "nb"]*3 + 2
      changed$data[held, "positive"] <- changed$data[held, "positive"]*4
      changed$data[held, "z"] <- changed$data[held, "z"] + 3
      changed$inputs$kernel_dat$hab <- -1e9
      changed$inputs$scl_params <- list(mean = c(hab = -1e9), sd = c(hab = .001))
      changed$inputs$binned$bin_rep[] <- 1e-5
      for (id in held) {
        changed$inputs$raw_cov[[id]] <- changed$inputs$raw_cov[[id]]*7
        changed$inputs$d_list[[id]] <- changed$inputs$d_list[[id]]*.6
      }
      run("nb_leak", changed)
      run("gamma_leak", changed, "Gamma", "refit")
      list(cases = cases, planned = planned)
    }, error = function(e) { cv_distribution_spatial_cache$failure <- e; stop(e) })
  }
  cv_distribution_spatial_cache$value
}

test_that("spatial NB/Gamma folds preserve frozen mean/SQ predictions, native density and ecological penalties", {
  skip_if_not_installed("MASS")
  z <- cv_distribution_spatial_fixture()
  expect_identical(z$planned, 17L)
  for (label in names(z$cases)) {
    case <- z$cases[[label]]; x <- case$cv; f <- case$fixture
    expect_true(x$complete, info = paste(label, paste(unlist(x$errors), collapse = "; ")))
    if (!x$complete) next
    m <- x$fits[[1]]; p <- x$predictions; d <- x$distribution[[1]]; r <- x$selection_records[[1]]$data
    if (label %in% c("nb_known", "gamma_known")) expect_identical(x$settings$backend, "refit")
    expect_equal(r$K, attr(logLik(m$opt_mod), "df") + nrow(m$scale_est) + if (f$inputs$kernel == "expow") nrow(m$shape_est) else 0)
    expect_equal(r$native_K, m$opt_mod$rank + d$estimated)
    expect_identical(x$raster_groups[[1]], lapply(x$folds[[1]], function(ids) f$inputs$raster_group[ids]))
    # Independent mean/SQ kernel reduction with the training-only representative
    # distances and optimized transform; raw assessment data learn nothing.
    sigma <- m$scale_est$Mean[1]/f$inputs$unit_conv
    shape <- if (f$inputs$kernel == "expow") m$shape_est$Mean[1] else NULL
    reference <- x$training[[1]]$bin_reference
    raw <- vapply(p$row_id, function(id) {
      distance <- f$inputs$d_list[[id]]; v <- f$inputs$raw_cov[[id]][, 1]
      if (!is.null(reference)) distance <- reference$bin_rep[.bincode(distance, reference$breaks, include.lowest = TRUE)]
      w <- switch(f$inputs$kernel, gaussian = exp(-distance^2/(2*sigma^2)),
        exp = exp(-distance/sigma), fixed = as.numeric(distance < sigma), expow = exp(-(distance/sigma)^shape))
      mean <- sum(w*v)/sum(w)
      if (grepl("sq", label)) sqrt(sum(w*(v - mean)^2)/sum(w)) else mean
    }, 0)
    nd <- f$data[p$row_id, ]; nd$hab <- (raw - m$scl_params$mean[["hab"]])/m$scl_params$sd[["hab"]]
    expect_equal(p$prediction, as.numeric(predict(m$opt_mod, nd, type = "response")), tolerance = 1e-8, info = label)
    expected <- if (x$settings$family == "nbinom2") -dnbinom(p$observed, size = d$theta, mu = p$prediction, log = TRUE) else
      -dgamma(p$observed, shape = 1/d$dispersion, scale = p$prediction*d$dispersion, log = TRUE)
    expect_equal(p$neg_log_score, expected, tolerance = 1e-10)
    expect_identical(x$settings$parameter_uncertainty, "not_integrated")
    compact <- x; compact$fits <- NULL
    for (criterion in c("AIC", "AICc", "BIC")) {
      a <- select_cv(list(candidate = x), criterion); b <- select_cv(list(candidate = compact), criterion)
      expect_true(a$complete); expect_true(b$complete)
      expect_identical(a$selection, b$selection); expect_identical(a$predictions, b$predictions)
    }
    if (f$inputs$kernel == "fixed") {
      expect_true(all(is.na(m$scale_est$SE)))
      expect_true(all(is.finite(vcov(m$opt_mod))))
      expect_identical(m$diagnostics$refit_hessian$status, "no_active_smooth_parameters")
      expect_false(m$diagnostics$refit_hessian$triggered)
      expect_identical(unname(m$diagnostics$refit_hessian$active), FALSE)
    }
    if (.msr_analytic_backend(m)) {
      expect_gte(nrow(m$profile_fit$start_results), 3L)
      expect_true(any(m$profile_fit$start_results$eligible))
      if (isTRUE(m$uncertainty$available)) expect_true(all(is.finite(m$scale_est$SE) & m$scale_est$SE > 0)) else
        expect_true(all(is.na(m$scale_est$SE)))
    }
  }
})

test_that("spatial NB/Gamma training nuisance, scales, transforms and starts resist assessment perturbation", {
  z <- cv_distribution_spatial_fixture()
  for (pair in list(c("nb_profile", "nb_leak"), c("gamma", "gamma_leak"))) {
    a <- z$cases[[pair[1]]]$cv; b <- z$cases[[pair[2]]]$cv
    expect_true(a$complete); expect_true(b$complete)
    expect_identical(a$distribution, b$distribution)
    untimed <- function(training) lapply(training, function(fold) {
      if (!is.null(fold$start_results)) fold$start_results$seconds <- NULL
      fold
    })
    expect_identical(untimed(a$training), untimed(b$training))
    if (a$complete && b$complete) {
      expect_identical(coef(a$fits[[1]]$opt_mod), coef(b$fits[[1]]$opt_mod))
      expect_identical(a$fits[[1]]$scl_params, b$fits[[1]]$scl_params)
    }
  }
})

test_that("unqualified Gamma and fixed-theta analytic scale routes reject before fitting", {
  f <- cv_distribution_spatial_inputs()
  expect_error(cross_validate(positive ~ hab, f$data, f$inputs, f$folds, "Gamma", "profile_glm"), "require backend")
  expect_error(cross_validate(nb ~ hab, f$data, f$inputs, f$folds, "nbinom2", "profile_glm", theta = 3), "require backend")
  expect_error(cross_validate(nb ~ hab, f$data, f$inputs, f$folds, "nbinom2", "joint_nb2", theta = 3), "require backend")
  expect_error(cross_validate(nb ~ hab, f$data, f$inputs, f$folds, "poisson", "joint_nb2"), "requires family")
})

test_that("joint-NB2 CV preserves the exact Poisson boundary and estimated-nuisance penalty", {
  x <- cv_distribution_spatial_fixture()$cases$nb_poisson_limit$cv
  expect_true(x$complete, info = paste(unlist(x$errors), collapse = "; "))
  if (!x$complete) return(invisible(NULL))
  m <- x$fits[[1]]; d <- x$distribution[[1]]
  expect_true(m$dispersion_est$exact_poisson)
  expect_true(d$poisson_limit); expect_identical(d$theta, Inf)
  expect_true(d$estimated)
  expect_identical(x$training[[1]]$initial_distribution$theta, 1)
  expect_false(x$training[[1]]$initial_distribution$estimated)
  expect_equal(x$predictions$neg_log_score, -dpois(x$predictions$observed, x$predictions$prediction, log = TRUE))
  expect_equal(x$selection_records[[1]]$data$native_K, m$opt_mod$rank + 1)
  compact <- x; compact$fits <- NULL
  expect_identical(select_cv(list(candidate = x), "AIC")$selection,
    select_cv(list(candidate = compact), "AIC")$selection)
  bad <- compact; bad$settings$backend <- "profile_glm"
  expect_error(cv_tab(list(bad = bad), "neg_log_score"), "nuisance")
})

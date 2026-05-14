test_that("multiScale_optim returns a valid optimized object", {
  fix <- make_core_fixture()
  opt <- fix$opt

  expect_s3_class(opt, "multiScaleR")
  expect_equal(rownames(opt$scale_est), "cont1")
  expect_equal(colnames(opt$scale_est), c("Mean", "SE"))
  expect_true(inherits(opt$opt_mod, "glm"))
  expect_true(is.list(opt$scl_params))
  expect_true(all(c("mean", "sd") %in% names(opt$scl_params)))
  expect_true(opt$scale_est[1, "Mean"] > 0)
  expect_true(is.list(opt$diagnostics))
  expect_true(is.list(opt$next_run))
  expect_equal(unname(opt$next_run$start_sigma),
               unname(opt$scale_est[, "Mean"]))
  expect_equal(opt$next_run$max_D, opt$max_D)
  expect_equal(unname(opt$next_run$start_par),
               unname(opt$next_run$start_sigma / opt$next_run$max_D))
})

test_that("build_opt_context stores complete-case indices from model data", {
  fake_mod <- structure(list(), class = "glm")
  dat <- data.frame(y = c(1, 2, 3), x = c(0.1, NA, 0.3), z = c(1, 2, 3))

  ctx <- with_mocked_bindings(
    build_opt_context(
      fitted_mod = fake_mod,
      cov_df = list(data.frame(x = 1, other = 2))
    ),
    find_predictors = function(x) list(c("x")),
    get_data = function(x, ...) dat,
    .package = "multiScaleR"
  )

  expect_equal(ctx$complete_idx, c(1, 3))
  expect_equal(nrow(ctx$data_template), 2)
  expect_equal(ctx$data_template$y, c(1, 3))
})

test_that("build_opt_context preserves original row indices from model frames", {
  fake_mod <- structure(list(), class = "glm")
  dat <- data.frame(
    y = c(1, 3, 4, 6),
    x = c(0.1, 0.3, 0.4, 0.6),
    z = c(1, 3, 4, 6),
    row.names = c("1", "3", "4", "6")
  )

  ctx <- with_mocked_bindings(
    build_opt_context(
      fitted_mod = fake_mod,
      cov_df = replicate(6, data.frame(x = 1, other = 2), simplify = FALSE)
    ),
    find_predictors = function(x) list(c("x")),
    get_data = function(x, ...) dat,
    .package = "multiScaleR"
  )

  expect_equal(ctx$data_idx, 1:4)
  expect_equal(ctx$complete_idx, c(1, 3, 4, 6))
  expect_equal(nrow(ctx$data_template), 4)
  expect_equal(ctx$data_template$y, c(1, 3, 4, 6))
})

test_that("kernel_scale_fn aligns scaled covariates to model complete cases", {
  fake_mod <- structure(list(), class = "glm")
  captured <- new.env(parent = emptyenv())
  captured$seen <- numeric(0)
  captured$refit_data <- NULL

  opt_context <- list(
    fitted_mod = fake_mod,
    mod_class = "glm",
    covs = "x",
    n_covs = 1,
    data_template = data.frame(y = c(2, 4), x = c(0, 0),
                               row.names = c("2", "4")),
    complete_idx = c(2, 4)
  )

  out <- with_mocked_bindings(
    kernel_scale_fn(
      par = 0.2,
      d_list = as.list(1:4),
      cov_df = list(
        data.frame(x = 1),
        data.frame(x = 2),
        data.frame(x = 3),
        data.frame(x = 4)
      ),
      kernel = "gaussian",
      fitted_mod = fake_mod,
      opt_context = opt_context
    ),
    scale_type = function(d, kernel, sigma, shape, r_stack.df) {
      captured$seen <- c(captured$seen, r_stack.df$x[[1]])
      r_stack.df$x[[1]]
    },
    .refit_model = function(model, data, opt_context) {
      captured$refit_data <- data
      structure(list(), class = "mock_fit")
    },
    .neg_loglik_model = function(model, mod_class) 0,
    .package = "multiScaleR"
  )

  expected_x <- as.numeric(scale(matrix(c(2, 4), ncol = 1)))
  expect_equal(captured$seen, c(2, 4))
  expect_equal(captured$refit_data$y, c(2, 4))
  expect_equal(captured$refit_data$x, expected_x)
  expect_equal(out, 0)
})

test_that("multiScale_optim builds and forwards cached optimization context", {
  fix <- make_core_fixture()
  captured <- new.env(parent = emptyenv())
  refit_fn <- function(model, data, context) model

  out <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = fix$kernel_inputs,
      par = 0.2,
      verbose = FALSE,
      refit_fn = refit_fn
    ),
    optim = function(..., opt_context) {
      captured$optim_context <- opt_context
      list(par = 0.2, hessian = matrix(1, 1, 1))
    },
    kernel_scale_fn = function(..., opt_context, mod_return = NULL) {
      captured$final_context <- opt_context
      if (isTRUE(mod_return)) {
        list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
      } else {
        0
      }
    },
    kernel_dist = function(...) data.frame(Mean = 20, low = 10, high = 30),
    .package = "multiScaleR"
  )

  expect_s3_class(out, "multiScaleR")
  expect_equal(captured$optim_context$mod_class, "glm")
  expect_equal(captured$optim_context$covs, "cont1")
  expect_equal(captured$optim_context$n_covs, 1)
  expect_identical(captured$optim_context$refit_fn, refit_fn)
  expect_true(is.data.frame(captured$optim_context$data_template))
  expect_identical(captured$final_context$covs, captured$optim_context$covs)
  expect_identical(captured$final_context$cov_idx, captured$optim_context$cov_idx)
  expect_false(isTRUE(out$diagnostics$max_distance$triggered))
  expect_equal(unname(out$diagnostics$max_distance$effective_distance[[1]]), 20)
  expect_equal(out$diagnostics$max_distance$suggested_max_D, 40)
  expect_equal(out$next_run$max_D, fix$kernel_inputs$max_D)
  expect_equal(unname(out$next_run$start_sigma), 0.2 * fix$kernel_inputs$unit_conv)
  expect_equal(unname(out$next_run$start_par), 0.2)
})

test_that("multiScale_optim can prescreen sigma starts with marginal scans only", {
  fix <- make_core_fixture()
  captured <- new.env(parent = emptyenv())
  sigma_grid <- exp(seq(
    log(fix$kernel_inputs$min_D / fix$kernel_inputs$unit_conv),
    log(fix$kernel_inputs$max_D / fix$kernel_inputs$unit_conv),
    length.out = 5
  ))
  target <- sigma_grid[4]

  out <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = fix$kernel_inputs,
      start_strategy = "screen",
      screen_n_sigma = 5,
      screen_n_jitter = 0,
      verbose = FALSE
    ),
    optim = function(par, fn, hessian, ...) {
      captured$full_par <- par
      list(
        par = par,
        hessian = matrix(1, 1, 1),
        value = fn(par, ...)
      )
    },
    kernel_scale_fn = function(par, mod_return = NULL, ...) {
      if (isTRUE(mod_return)) {
        list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
      } else {
        (par[1] - target)^2
      }
    },
    kernel_dist = function(...) data.frame(Mean = 20, low = 10, high = 30),
    .package = "multiScaleR"
  )

  expect_s3_class(out, "multiScaleR")
  expect_equal(captured$full_par[1], target, tolerance = 1e-10)
})

test_that("multiScale_optim screening keeps exactly one full optimization", {
  fix <- make_core_fixture()
  captured <- new.env(parent = emptyenv())
  captured$screen_calls <- 0L
  captured$full_calls <- 0L
  target <- 0.3

  out <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = fix$kernel_inputs,
      start_strategy = "screen",
      screen_n_sigma = 4,
      screen_n_jitter = 4,
      screen_maxit = 3,
      verbose = FALSE
    ),
    optim = function(par, fn, hessian, ...) {
      if (isTRUE(hessian)) {
        captured$full_calls <- captured$full_calls + 1L
        captured$full_par <- par
        return(list(
          par = par,
          hessian = matrix(1, 1, 1),
          value = fn(par, ...)
        ))
      }

      captured$screen_calls <- captured$screen_calls + 1L
      list(
        par = c(target),
        value = fn(c(target), ...),
        convergence = 0L
      )
    },
    kernel_scale_fn = function(par, mod_return = NULL, ...) {
      if (isTRUE(mod_return)) {
        list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
      } else {
        (par[1] - target)^2
      }
    },
    kernel_dist = function(...) data.frame(Mean = 20, low = 10, high = 30),
    .package = "multiScaleR"
  )

  expect_s3_class(out, "multiScaleR")
  expect_gt(captured$screen_calls, 0L)
  expect_equal(captured$full_calls, 1L)
  expect_equal(captured$full_par[1], target, tolerance = 1e-10)
})

test_that("kernel_scale_fn can use a custom refit function", {
  fake_mod <- structure(list(), class = "custom_model")
  captured <- new.env(parent = emptyenv())
  captured$data <- NULL
  captured$model_class <- NULL

  opt_context <- list(
    fitted_mod = fake_mod,
    mod_class = "other",
    covs = "x",
    n_covs = 1,
    data_template = data.frame(y = c(1, 2), x = c(0, 0)),
    complete_idx = c(1, 2),
    refit_fn = function(model, data, context) {
      captured$data <- data
      structure(list(data = data), class = "custom_fit")
    }
  )

  out <- with_mocked_bindings(
    kernel_scale_fn(
      par = 0.2,
      d_list = list(1, 1),
      cov_df = list(data.frame(x = 2), data.frame(x = 4)),
      kernel = "gaussian",
      fitted_mod = fake_mod,
      opt_context = opt_context
    ),
    scale_type = function(d, kernel, sigma, shape, r_stack.df) r_stack.df$x[[1]],
    .neg_loglik_model = function(model, mod_class) {
      captured$model_class <- class(model)
      1.23
    },
    .package = "multiScaleR"
  )

  expect_equal(out, 1.23)
  expect_equal(captured$model_class, "custom_fit")
  expect_equal(names(captured$data), c("y", "x"))
  expect_false(all(captured$data$x == 0))
})

test_that("kernel_scale_fn can refit survival models through stats generics", {
  testthat::skip_if_not_installed("survival")

  fix <- make_core_fixture()
  n <- nrow(fix$kernel_inputs$kernel_dat)
  dat <- data.frame(
    time = seq_len(n),
    status = rep(c(1, 1, 0), length.out = n),
    cont1 = fix$kernel_inputs$kernel_dat$cont1
  )

  mod <- survival::coxph(survival::Surv(time, status) ~ cont1, data = dat)
  opt_context <- build_opt_context(fitted_mod = mod, cov_df = fix$kernel_inputs$raw_cov)

  neg_ll <- kernel_scale_fn(
    par = 40 / fix$kernel_inputs$unit_conv,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = "gaussian",
    fitted_mod = mod,
    opt_context = opt_context
  )

  fit <- kernel_scale_fn(
    par = 40 / fix$kernel_inputs$unit_conv,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = "gaussian",
    fitted_mod = mod,
    opt_context = opt_context,
    mod_return = TRUE
  )

  expect_true(is.finite(neg_ll))
  expect_s3_class(fit$mod, "coxph")
})

test_that("nested clogit wrappers can be optimized through their analysis model", {
  testthat::skip_if_not_installed("survival")

  fix <- make_core_fixture()
  n <- nrow(fix$kernel_inputs$kernel_dat)
  dat <- data.frame(
    case_ = rep(c(TRUE, FALSE, FALSE), length.out = n),
    cont1 = fix$kernel_inputs$kernel_dat$cont1,
    site = seq_len(n) / 10,
    stratum = factor(rep(seq_len(ceiling(n / 3)), each = 3, length.out = n))
  )
  coxph <- survival::coxph
  strata <- survival::strata

  inner <- suppressWarnings(
    survival::clogit(
      case_ ~ cont1 + site + strata(stratum),
      data = dat,
      model = TRUE
    )
  )
  wrapped <- structure(
    list(model = inner, sl_ = NULL, ta_ = NULL, more = NULL),
    class = c("fit_clogit", "list")
  )

  opt_context <- build_opt_context(
    fitted_mod = wrapped,
    cov_df = fix$kernel_inputs$raw_cov
  )

  neg_ll <- suppressWarnings(
    kernel_scale_fn(
      par = 40 / fix$kernel_inputs$unit_conv,
      d_list = fix$kernel_inputs$d_list,
      cov_df = fix$kernel_inputs$raw_cov,
      kernel = "gaussian",
      fitted_mod = wrapped,
      opt_context = opt_context
    )
  )

  fit <- suppressWarnings(
    kernel_scale_fn(
      par = 40 / fix$kernel_inputs$unit_conv,
      d_list = fix$kernel_inputs$d_list,
      cov_df = fix$kernel_inputs$raw_cov,
      kernel = "gaussian",
      fitted_mod = wrapped,
      opt_context = opt_context,
      mod_return = TRUE
    )
  )

  expect_equal(opt_context$covs, "cont1")
  expect_true(all(c("case_", "cont1", "site", "stratum") %in%
                    names(opt_context$data_template)))
  expect_true(is.finite(neg_ll))
  expect_s3_class(fit$mod, "fit_clogit")
  expect_s3_class(fit$mod$model, "coxph")
})

test_that("single-covariate models work when kernel inputs contain extra raster layers", {
  set.seed(555)

  pts <- terra::vect(cbind(c(5, 7, 9, 11, 13), c(13, 11, 9, 7, 5)))
  rast_stack <- terra::rast(list(
    r1 = terra::rast(matrix(rnorm(20^2), nrow = 20)),
    r2 = terra::rast(matrix(rnorm(20^2), nrow = 20))
  ))

  kernel_inputs <- kernel_prep(
    pts = pts,
    raster_stack = rast_stack,
    max_D = 5,
    kernel = "gaussian",
    verbose = FALSE
  )

  dat <- data.frame(y = rnorm(5), kernel_inputs$kernel_dat)
  mod1 <- glm(y ~ r1, data = dat)

  expect_s3_class(
    suppressWarnings(
      suppressMessages(
        multiScale_optim(
          fitted_mod = mod1,
          kernel_inputs = kernel_inputs,
          par = 0.5,
          n_cores = NULL,
          verbose = FALSE
        )
      )
    ),
    "multiScaleR"
  )
})

test_that("summary and distance methods return structured outputs", {
  fix <- make_core_fixture()

  sum_opt <- summary(fix$opt, prob = 0.95)
  dist_opt <- kernel_dist(fix$opt, prob = 0.95)

  expect_s3_class(sum_opt, "summary_multiScaleR")
  expect_equal(colnames(sum_opt$opt_scale), c("Mean", "SE", "2.5%", "97.5%"))
  expect_equal(colnames(dist_opt), c("Mean", "2.5%", "97.5%"))
  expect_true(dist_opt[1, "Mean"] > 0)
  expect_identical(sum_opt$opt_distance, sum_opt$opt_dist)
  expect_identical(sum_opt$diagnostics, fix$opt$diagnostics)
})

test_that("diagnostics accessor returns structured warning metadata", {
  fix <- make_core_fixture()

  expect_identical(diagnostics(fix$opt), fix$opt$diagnostics)

  legacy_obj <- structure(list(warn_message = c(0, 1)), class = "multiScaleR")
  expect_identical(
    diagnostics(legacy_obj),
    list(
      max_distance = NULL,
      sigma_precision = NULL,
      shape_precision = NULL
    )
  )
})

test_that("print methods emit readable summaries", {
  fix <- make_core_fixture()
  warn_obj <- fix$opt
  warn_obj$warn_message <- c(0, 1)
  warn_obj$diagnostics$max_distance <- list(
    triggered = TRUE,
    variables = "cont1",
    suggested_max_D = 123.45
  )

  expect_output(print(fix$opt), "Optimized Scale of Effect")
  expect_output(print(summary(fix$opt)), "Fitted Model Summary")
  expect_output(print(warn_obj), "123.45")
})

test_that("multiScale_optim validates kernel inputs and parameter lengths", {
  fix <- make_core_fixture()

  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = list(), verbose = FALSE),
    "kernel_inputs must be a list"
  )

  expect_error(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = fix$kernel_inputs,
      par = c(0.1, 0.2),
      verbose = FALSE
    ),
    "length of par"
  )

  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = fix$kernel_inputs, join_by = 1, verbose = FALSE),
    "join_by must be a data frame"
  )

  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = fix$kernel_inputs, par = "a", verbose = FALSE),
    "par must be numeric"
  )

  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = fix$kernel_inputs, n_cores = 0, verbose = FALSE),
    "n_cores must be a positive integer"
  )
  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = fix$kernel_inputs, PSOCK = NA, verbose = FALSE),
    "PSOCK"
  )
  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = fix$kernel_inputs, verbose = NA),
    "verbose"
  )

  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = fix$kernel_inputs, refit_fn = 1, verbose = FALSE),
    "`refit_fn` must be a function"
  )

  bad_ki <- fix$kernel_inputs
  colnames(bad_ki$raw_cov[[1]]) <- "not_in_model"
  expect_error(
    multiScale_optim(fitted_mod = fix$fitted_mod, kernel_inputs = bad_ki, verbose = FALSE),
    "do not match the variables used in your fitted model"
  )
})

test_that("multiScale_optim covers expow and warning branches via mocks", {
  fix <- make_core_fixture()
  expow_ki <- fix$kernel_inputs
  expow_ki$kernel <- "expow"
  expow_ki$shape <- 2

  out <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = expow_ki,
      par = NULL,
      verbose = FALSE
    ),
    optim = function(...) {
      list(par = c(0.2, 2), hessian = diag(c(0.04, 0.04)))
    },
    kernel_scale_fn = function(..., mod_return = NULL) {
      if (isTRUE(mod_return)) {
        list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
      } else {
        0
      }
    },
    kernel_dist = function(...) data.frame(Mean = 200, low = 100, high = 300),
    .package = "multiScaleR"
  )

  expect_s3_class(out, "multiScaleR")
  expect_false(is.null(out$shape_est))
  expect_true(all(c(1, 2, 3) %in% out$warn_message))
  expect_true(isTRUE(out$diagnostics$max_distance$triggered))
  expect_true(isTRUE(out$diagnostics$sigma_precision$triggered))
  expect_true(isTRUE(out$diagnostics$shape_precision$triggered))
  expect_equal(out$next_run$max_D, 400)
  expect_equal(unname(out$next_run$start_sigma), 0.2 * fix$kernel_inputs$unit_conv)
  expect_equal(unname(out$next_run$start_shape), 2)
  expect_equal(unname(out$next_run$start_par), c(50 / 400, 2))
  expect_equal(out$next_run$flags,
               c("max_distance", "sigma_precision", "shape_precision"))
})

test_that("multiScale_optim covers singular hessian and failure branches via mocks", {
  fix <- make_core_fixture()

  singular_out <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = fix$kernel_inputs,
      par = 0.2,
      verbose = FALSE
    ),
    optim = function(...) {
      list(par = 0.2, hessian = matrix(0, 1, 1))
    },
    kernel_scale_fn = function(..., mod_return = NULL) {
      if (isTRUE(mod_return)) {
        list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
      } else {
        0
      }
    },
    kernel_dist = function(...) data.frame(Mean = 20, low = 10, high = 30),
    .package = "multiScaleR"
  )

  expect_s3_class(singular_out, "multiScaleR")
  expect_type(singular_out$scale_est$SE, "double")
  expect_true(is.infinite(singular_out$scale_est$SE[[1]]))
  singular_summary <- summary(singular_out)
  expect_s3_class(singular_summary, "summary_multiScaleR")
  expect_type(singular_summary$opt_scale$SE, "double")

  expow_ki <- fix$kernel_inputs
  expow_ki$kernel <- "expow"
  expow_ki$shape <- 2

  singular_expow <- with_mocked_bindings(
      multiScale_optim(
        fitted_mod = fix$fitted_mod,
        kernel_inputs = expow_ki,
        par = c(0.2, 2),
        verbose = FALSE
      ),
      optim = function(...) {
        list(par = c(0.2, 2), hessian = matrix(0, 2, 2))
      },
      kernel_scale_fn = function(..., mod_return = NULL) {
        if (isTRUE(mod_return)) {
          list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
        } else {
          0
        }
      },
      kernel_dist = function(...) data.frame(Mean = 20, low = 10, high = 30),
      .package = "multiScaleR"
  )

  expect_s3_class(singular_expow, "multiScaleR")
  expect_type(singular_expow$scale_est$SE, "double")
  expect_type(singular_expow$shape_est$SE, "double")
  expect_true(is.infinite(singular_expow$scale_est$SE[[1]]))
  expect_true(is.infinite(singular_expow$shape_est$SE[[1]]))
  expow_summary <- summary(singular_expow)
  expect_s3_class(expow_summary, "summary_multiScaleR")
  expect_type(expow_summary$opt_scale$SE, "double")
  expect_type(expow_summary$opt_shape$SE, "double")

  expect_error(
    with_mocked_bindings(
      multiScale_optim(
        fitted_mod = fix$fitted_mod,
        kernel_inputs = fix$kernel_inputs,
        par = 0.2,
        verbose = TRUE
      ),
      optim = function(...) {
        structure(
          "boom",
          class = "try-error",
          condition = simpleError("mock standard failure")
        )
      },
      .package = "multiScaleR"
    ),
    "mock standard failure"
  )
})

test_that("multiScale_optim covers parallel and unmarked branches via mocks", {
  fix <- make_core_fixture()

  par_out <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = fix$fitted_mod,
      kernel_inputs = fix$kernel_inputs,
      par = 0.2,
      n_cores = 2,
      verbose = FALSE
    ),
    makeCluster = function(...) "cl",
    setDefaultCluster = function(cl = NULL) invisible(NULL),
    cluster_prep = function(model, cl) "stats",
    optimParallel = function(...) list(par = 0.2, hessian = matrix(1, 1, 1)),
    stopCluster = function(cl) invisible(NULL),
    kernel_scale_fn = function(..., mod_return = NULL) {
      if (isTRUE(mod_return)) {
        list(mod = fix$fitted_mod, scl_params = fix$kernel_inputs$scl_params)
      } else {
        0
      }
    },
    kernel_dist = function(...) data.frame(Mean = 20, low = 10, high = 30),
    .package = "multiScaleR"
  )

  expect_s3_class(par_out, "multiScaleR")

  expect_error(
    with_mocked_bindings(
      multiScale_optim(
        fitted_mod = fix$fitted_mod,
        kernel_inputs = fix$kernel_inputs,
        par = 0.2,
        n_cores = 2,
        verbose = TRUE
      ),
      makeCluster = function(...) "cl",
      setDefaultCluster = function(cl = NULL) invisible(NULL),
      cluster_prep = function(model, cl) "stats",
      optimParallel = function(...) {
        structure(
          "boom",
          class = "try-error",
          condition = simpleError("mock parallel failure")
        )
      },
      stopCluster = function(cl) invisible(NULL),
      .package = "multiScaleR"
    ),
    "mock parallel failure"
  )

  unmark <- make_unmarked_fixture()
  mocked_unmark <- with_mocked_bindings(
    multiScale_optim(
      fitted_mod = unmark$mod,
      kernel_inputs = unmark$kernel_inputs,
      par = c(0.2, 0.2),
      verbose = FALSE
    ),
    optim = function(...) list(par = c(0.2, 0.2), hessian = diag(2)),
    kernel_scale_fn = function(..., mod_return = NULL) {
      if (isTRUE(mod_return)) {
        list(mod = unmark$mod, scl_params = unmark$kernel_inputs$scl_params)
      } else {
        0
      }
    },
    kernel_dist = function(...) data.frame(Mean = c(20, 25), low = c(10, 15), high = c(30, 35)),
    .package = "multiScaleR"
  )

  expect_s3_class(mocked_unmark, "multiScaleR")
  expect_equal(rownames(mocked_unmark$scale_est), c("bin1", "cont1"))
})

test_that("PSOCK workers load the same multiScaleR namespace as the main session", {
  skip_on_cran()
  skip_if_not_installed("pkgload")

  fix <- make_core_fixture()
  cl <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  cluster_prep(fix$fitted_mod, cl)

  master_path <- normalizePath(getNamespaceInfo("multiScaleR", "path"),
                               winslash = "/", mustWork = FALSE)
  master_system_path <- normalizePath(system.file(package = "multiScaleR"),
                                      winslash = "/", mustWork = FALSE)
  master_load_all <- nzchar(master_system_path) &&
    !identical(master_path, master_system_path)
  master_version <- as.character(utils::packageVersion("multiScaleR"))

  worker_info <- parallel::clusterEvalQ(
    cl,
    list(
      version = as.character(utils::packageVersion("multiScaleR")),
      path = normalizePath(getNamespaceInfo("multiScaleR", "path"),
                           winslash = "/", mustWork = FALSE)
    )
  )[[1]]

  expect_equal(worker_info$version, master_version)
  if (isTRUE(master_load_all)) {
    expect_equal(worker_info$path, master_path)
  }
})

test_that("cluster_prep exports library paths for PSOCK workers", {
  fix <- make_core_fixture()
  captured <- new.env(parent = emptyenv())

  with_mocked_bindings(
    cluster_prep(fix$fitted_mod, cl = "cl"),
    clusterExport = function(cl, varlist, envir) {
      captured$varlist <- varlist
      captured$lib_paths <- get("lib_paths", envir = envir)
      captured$r_libs_user <- get("r_libs_user", envir = envir)
      captured$r_libs <- get("r_libs", envir = envir)
      invisible(NULL)
    },
    clusterEvalQ = function(cl, expr) {
      list(list(
        version = as.character(utils::packageVersion("multiScaleR")),
        path = normalizePath(getNamespaceInfo("multiScaleR", "path"),
                             winslash = "/", mustWork = FALSE)
      ))
    },
    .package = "multiScaleR"
  )

  expect_true(all(c("lib_paths", "r_libs_user", "r_libs") %in% captured$varlist))
  expect_true(length(captured$lib_paths) >= 1)
  expect_equal(captured$r_libs_user, Sys.getenv("R_LIBS_USER", unset = ""))
  expect_equal(captured$r_libs, Sys.getenv("R_LIBS", unset = ""))
})

test_that("PSOCK optimization works for unqualified MASS model calls", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("pkgload")

  fix <- make_core_fixture()
  mass_attached <- "package:MASS" %in% search()
  if (!isTRUE(mass_attached)) {
    library(MASS)
    on.exit(detach("package:MASS", unload = FALSE), add = TRUE)
  }

  nb_mod <- glm.nb(y ~ cont1 + site, data = fix$df)

  serial <- suppressWarnings(
    suppressMessages(
      multiScale_optim(
        fitted_mod = nb_mod,
        kernel_inputs = fix$kernel_inputs,
        par = 40 / fix$kernel_inputs$unit_conv,
        n_cores = NULL,
        verbose = FALSE
      )
    )
  )

  parallel <- suppressWarnings(
    suppressMessages(
      multiScale_optim(
        fitted_mod = nb_mod,
        kernel_inputs = fix$kernel_inputs,
        par = 40 / fix$kernel_inputs$unit_conv,
        n_cores = 2,
        PSOCK = TRUE,
        verbose = FALSE
      )
    )
  )

  expect_equal(parallel$scale_est$Mean, serial$scale_est$Mean,
               tolerance = 1e-5)
  expect_true(is.finite(parallel$scale_est$SE[[1]]))
})

test_that("PSOCK optimization works for nested clogit wrappers", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("pkgload")

  fix <- make_core_fixture()
  n <- nrow(fix$kernel_inputs$kernel_dat)
  dat <- data.frame(
    case_ = rep(c(TRUE, FALSE, FALSE), length.out = n),
    cont1 = fix$kernel_inputs$kernel_dat$cont1,
    site = seq_len(n) / 10,
    stratum = factor(rep(seq_len(ceiling(n / 3)), each = 3, length.out = n))
  )
  coxph <- survival::coxph
  strata <- survival::strata

  inner <- suppressWarnings(
    survival::clogit(
      case_ ~ cont1 + site + strata(stratum),
      data = dat,
      model = TRUE
    )
  )
  wrapped <- structure(
    list(model = inner, sl_ = NULL, ta_ = NULL, more = NULL),
    class = c("fit_clogit", "list")
  )

  serial <- suppressWarnings(
    suppressMessages(
      multiScale_optim(
        fitted_mod = wrapped,
        kernel_inputs = fix$kernel_inputs,
        par = 40 / fix$kernel_inputs$unit_conv,
        n_cores = NULL,
        verbose = FALSE
      )
    )
  )

  parallel <- suppressWarnings(
    suppressMessages(
      multiScale_optim(
        fitted_mod = wrapped,
        kernel_inputs = fix$kernel_inputs,
        par = 40 / fix$kernel_inputs$unit_conv,
        n_cores = 2,
        PSOCK = TRUE,
        verbose = FALSE
      )
    )
  )

  expect_equal(parallel$scale_est$Mean, serial$scale_est$Mean,
               tolerance = 1e-5)
  expect_s3_class(parallel$opt_mod, "fit_clogit")
  expect_s3_class(parallel$opt_mod$model, "coxph")
})

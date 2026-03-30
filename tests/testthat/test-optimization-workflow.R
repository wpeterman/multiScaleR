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
})

test_that("summary and distance methods return structured outputs", {
  fix <- make_core_fixture()

  sum_opt <- summary(fix$opt, prob = 0.95)
  dist_opt <- kernel_dist(fix$opt, prob = 0.95)

  expect_s3_class(sum_opt, "summary_multiScaleR")
  expect_equal(colnames(sum_opt$opt_scale), c("Mean", "SE", "2.5%", "97.5%"))
  expect_equal(colnames(dist_opt), c("Mean", "2.5%", "97.5%"))
  expect_true(dist_opt[1, "Mean"] > 0)
})

test_that("print methods emit readable summaries", {
  fix <- make_core_fixture()

  expect_output(print(fix$opt), "Optimized Scale of Effect")
  expect_output(print(summary(fix$opt)), "Fitted Model Summary")
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
})

test_that("multiScale_optim covers singular hessian and failure branches via mocks", {
  fix <- make_core_fixture()

  expect_error(
    with_mocked_bindings(
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
    ),
    "non-numeric argument"
  )

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
    "Standard optimization failed"
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
    "Parallel optimization failed"
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

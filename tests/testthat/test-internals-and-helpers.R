test_that("ci_func_r covers finite and non-finite branches", {
  x <- matrix(c(5, 1, 10, Inf), ncol = 2, byrow = TRUE)
  out <- multiScaleR:::ci_func_r(x, df = 20, min_D = 2, names = c("a", "b"))

  expect_equal(rownames(out), c("a", "b"))
  expect_true(all(out["a", c("2.5%", "97.5%")] >= 2))
  expect_true(all(is.nan(out["b", c("2.5%", "97.5%")])))
})

test_that("scale_type_cpp covers kernels and validates arguments", {
  d <- c(0, 1, 2)
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)

  expect_length(multiScaleR:::scale_type_cpp(d, kernel = "gaussian", sigma_ = c(1, 2), r_stack_df = mat), 2)
  expect_length(multiScaleR:::scale_type_cpp(d, kernel = "exp", sigma_ = c(1, 2), r_stack_df = mat), 2)
  expect_length(multiScaleR:::scale_type_cpp(d, kernel = "fixed", sigma_ = c(1, 2), r_stack_df = mat), 2)
  expect_length(multiScaleR:::scale_type_cpp(d, kernel = "expow", sigma_ = c(1, 2), shape_ = c(2, 3), r_stack_df = mat), 2)
  expect_length(multiScaleR:::scale_type_cpp(d, kernel = "gaussian", sigma_ = 2, r_stack_df = mat[, 1, drop = FALSE], output = 1L), 3)

  expect_warning(
    multiScaleR:::scale_type_cpp(d, kernel = "fixed", sigma_ = 0, r_stack_df = mat[, 1, drop = FALSE]),
    "Sum of weights"
  )

  expect_error(multiScaleR:::scale_type_cpp(d, kernel = "gaussian"), "sigma")
  expect_error(multiScaleR:::scale_type_cpp(d, kernel = "gaussian", sigma_ = 1), "r_stack_df")
  expect_error(multiScaleR:::scale_type_cpp(d, kernel = "expow", sigma_ = c(1, 2), r_stack_df = mat), "shape must be provided")
  expect_error(multiScaleR:::scale_type_cpp(d, kernel = "expow", sigma_ = c(1, 2), shape_ = 1, r_stack_df = matrix(1, nrow = 3, ncol = 3)), "same length")
  expect_error(multiScaleR:::scale_type_cpp(d, kernel = "bogus", sigma_ = c(1, 2), r_stack_df = mat), "Invalid kernel")
})

test_that("extract_namespace identifies namespaced calls and plain calls", {
  data(quine, package = "MASS")
  nb_mod <- MASS::glm.nb(Days ~ Sex/(Age + Eth * Lrn), data = quine)

  expect_equal(multiScaleR:::extract_namespace(nb_mod), "MASS")

  glm_mod <- glm(mpg ~ wt, data = mtcars)
  expect_null(multiScaleR:::extract_namespace(glm_mod))
})

test_that("cluster_prep covers direct, fallback, and error branches", {
  data(quine, package = "MASS")
  nb_mod <- MASS::glm.nb(Days ~ Sex/(Age + Eth * Lrn), data = quine)
  cl1 <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl1), add = TRUE)

  pkg1 <- multiScaleR:::cluster_prep(nb_mod, cl1)
  expect_equal(pkg1, "MASS")

  glm_mod <- glm(mpg ~ wt, data = mtcars)
  cl2 <- parallel::makeCluster(1)
  on.exit(parallel::stopCluster(cl2), add = TRUE)

  pkg2 <- multiScaleR:::cluster_prep(glm_mod, cl2)
  expect_equal(pkg2, "stats")

  expect_error(
    multiScaleR:::cluster_prep(structure(list(), class = "mystery_model"), cl2),
    "Could not determine the package"
  )
})

test_that("safe_predict and namespace handle supported and unsupported models", {
  glm_mod <- glm(mpg ~ wt, data = mtcars)
  pred <- multiScaleR:::safe_predict(glm_mod, data.frame(wt = c(2, 3)))
  expect_true(is.list(pred))

  expect_equal(multiScaleR:::namespace(glm_mod), "stats")
  expect_error(multiScaleR:::safe_predict(structure(list(), class = "mystery_model"), data.frame(x = 1)), "Could not determine")
  expect_error(multiScaleR:::namespace(structure(list(), class = "mystery_model")), "Could not determine")

  unmark <- make_unmarked_fixture()
  expect_equal(multiScaleR:::namespace(unmark$mod), "unmarked")
})

test_that("extract_model_data covers fallback paths", {
  assign("model.frame.fakeDataSlot", function(formula, ...) NULL, envir = .GlobalEnv)
  assign("model.frame.fakeCallData", function(formula, ...) NULL, envir = .GlobalEnv)
  assign("model.frame.HLfit", function(formula, ...) NULL, envir = .GlobalEnv)
  assign("model.frame.fakeWarn", function(formula, ...) NULL, envir = .GlobalEnv)
  assign("model.frame.fakeErr", function(formula, ...) stop("boom"), envir = .GlobalEnv)
  on.exit(
    rm(
      list = c("model.frame.fakeDataSlot", "model.frame.fakeCallData",
               "model.frame.HLfit", "model.frame.fakeWarn", "model.frame.fakeErr"),
      envir = .GlobalEnv
    ),
    add = TRUE
  )

  expect_s3_class(multiScaleR:::extract_model_data(glm(mpg ~ wt, data = mtcars)), "data.frame")

  fake1 <- structure(list(data = data.frame(x = 1:3)), class = "fakeDataSlot")
  expect_equal(nrow(multiScaleR:::extract_model_data(fake1)), 3)

  stored_df <- data.frame(y = 1:3)
  fake2 <- structure(list(call = list(data = as.name("stored_df"))), class = "fakeCallData")
  expect_equal(nrow(multiScaleR:::extract_model_data(fake2)), 3)

  fake3 <- structure(list(fr = data.frame(z = 3:1)), class = "HLfit")
  expect_equal(nrow(multiScaleR:::extract_model_data(fake3)), 3)

  fake3b <- structure(list(data = data.frame(w = 1:2)), class = "HLfit")
  expect_equal(nrow(multiScaleR:::extract_model_data(fake3b)), 2)

  fake4 <- structure(list(), class = "fakeWarn")
  expect_warning(expect_null(multiScaleR:::extract_model_data(fake4)), "Could not extract data")

  fake5 <- structure(list(), class = "fakeErr")
  expect_warning(expect_null(multiScaleR:::extract_model_data(fake5)), "Failed to extract data")
})

test_that("kernel_scale_fn covers glm, other, unmarked, fallback, and negative sigma paths", {
  fix <- make_core_fixture()

  expect_equal(
    multiScaleR:::kernel_scale_fn(
      par = -0.1,
      d_list = fix$kernel_inputs$d_list,
      cov_df = fix$kernel_inputs$raw_cov,
      kernel = fix$kernel_inputs$kernel,
      fitted_mod = fix$fitted_mod
    ),
    1e6^10
  )

  glm_obj <- multiScaleR:::kernel_scale_fn(
    par = 40 / fix$kernel_inputs$unit_conv,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = fix$fitted_mod,
    mod_return = TRUE
  )
  expect_true(is.list(glm_obj))
  expect_true(inherits(glm_obj$mod, "glm"))

  zi <- make_zeroinfl_fixture()
  other_obj <- multiScaleR:::kernel_scale_fn(
    par = 40 / fix$kernel_inputs$unit_conv,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = zi$mod
  )
  expect_true(is.numeric(other_obj))

  unmark <- make_unmarked_fixture()
  unmarked_obj <- multiScaleR:::kernel_scale_fn(
    par = c(40 / unmark$kernel_inputs$unit_conv, 60 / unmark$kernel_inputs$unit_conv),
    d_list = unmark$kernel_inputs$d_list,
    cov_df = unmark$kernel_inputs$raw_cov,
    kernel = unmark$kernel_inputs$kernel,
    fitted_mod = unmark$mod,
    join_by = data.frame(site = seq_len(nrow(unmark$kernel_inputs$kernel_dat)))
  )
  expect_true(is.numeric(unmarked_obj))

  gls_mod <- nlme::gls(y ~ cont1 + site, data = fix$df)
  gls_obj <- multiScaleR:::kernel_scale_fn(
    par = 40 / fix$kernel_inputs$unit_conv,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = gls_mod,
    mod_return = TRUE
  )
  expect_true(inherits(gls_obj$mod, "gls"))

  unmarked_no_join <- multiScaleR:::kernel_scale_fn(
    par = c(40 / unmark$kernel_inputs$unit_conv, 60 / unmark$kernel_inputs$unit_conv),
    d_list = unmark$kernel_inputs$d_list,
    cov_df = unmark$kernel_inputs$raw_cov,
    kernel = unmark$kernel_inputs$kernel,
    fitted_mod = unmark$mod
  )
  expect_true(is.numeric(unmarked_no_join))

  expect_equal(
    with_mocked_bindings(
      multiScaleR:::kernel_scale_fn(
        par = 40 / fix$kernel_inputs$unit_conv,
        d_list = fix$kernel_inputs$d_list,
        cov_df = fix$kernel_inputs$raw_cov,
        kernel = fix$kernel_inputs$kernel,
        fitted_mod = fix$fitted_mod
      ),
      logLik = function(x) stop("boom"),
      get_loglikelihood = function(x) 7,
      .package = "multiScaleR"
    ),
    -7
  )

  fake_other <- structure(list(), class = "fake_other_model")
  expect_error(
    with_mocked_bindings(
      multiScaleR:::kernel_scale_fn(
        par = 0.1,
        d_list = list(c(0, 1)),
        cov_df = list(Matrix::Matrix(matrix(c(1, 2), ncol = 1), sparse = TRUE)),
        kernel = "gaussian",
        fitted_mod = fake_other
      ),
      find_predictors = function(x) list("x"),
      get_data = function(x, ...) NULL,
      extract_model_data = function(model) NULL,
      .package = "multiScaleR"
    ),
    "Data from original model not saved"
  )
})

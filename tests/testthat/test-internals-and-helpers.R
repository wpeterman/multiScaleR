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
  expect_warning(
    out4 <- multiScaleR:::extract_model_data(fake4),
    "Could not extract data from the model object"
  )
  expect_null(out4)

  fake5 <- structure(list(), class = "fakeErr")
  expect_warning(
    out5 <- multiScaleR:::extract_model_data(fake5),
    "Failed to extract model data needed for refitting"
  )
  expect_null(out5)
})

test_that("build_opt_context caches reusable model metadata", {
  fix <- make_core_fixture()

  glm_context <- multiScaleR:::build_opt_context(
    fitted_mod = fix$fitted_mod,
    cov_df = fix$kernel_inputs$raw_cov
  )

  expect_equal(glm_context$mod_class, "glm")
  expect_equal(glm_context$covs, "cont1")
  expect_equal(glm_context$n_covs, 1)
  expect_true(is.data.frame(glm_context$data_template))
  expect_equal(glm_context$cov_idx, match("cont1", colnames(glm_context$data_template)))

  unmark <- make_unmarked_fixture()
  join_by <- data.frame(site = seq_len(nrow(unmark$kernel_inputs$kernel_dat)))
  um_context <- multiScaleR:::build_opt_context(
    fitted_mod = unmark$mod,
    cov_df = unmark$kernel_inputs$raw_cov,
    join_by = join_by
  )

  expect_equal(um_context$mod_class, "unmarked")
  expect_equal(um_context$covs, c("bin1", "cont1"))
  expect_equal(um_context$n_covs, 2)
  expect_equal(um_context$join_cols, "site")
  expect_true(all(c("site", "nuisance") %in% colnames(um_context$umf_template@siteCovs)))
  expect_false(any(um_context$covs %in% colnames(um_context$umf_template@siteCovs)))
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
  other_obj <- with_mocked_bindings(
    multiScaleR:::kernel_scale_fn(
      par = 40 / fix$kernel_inputs$unit_conv,
      d_list = fix$kernel_inputs$d_list,
      cov_df = fix$kernel_inputs$raw_cov,
      kernel = fix$kernel_inputs$kernel,
      fitted_mod = zi$mod
    ),
    find_predictors = function(x) list("cont1"),
    get_data = function(x, effects = "all") fix$df,
    .package = "multiScaleR"
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
      .neg_loglik_model = function(model, mod_class) -7,
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
        cov_df = list({
          m <- Matrix::Matrix(matrix(c(1, 2), ncol = 1), sparse = TRUE)
          colnames(m) <- "x"
          m
        }),
        kernel = "gaussian",
        fitted_mod = fake_other
      ),
      find_predictors = function(x) list("x"),
      get_data = function(x, ...) NULL,
      extract_model_data = function(model) NULL,
      .package = "multiScaleR"
    ),
    "Could not recover the original model data"
  )
})

test_that("build_opt_context and kernel_scale_fn report user-facing setup failures clearly", {
  unmark <- make_unmarked_fixture()

  expect_error(
    multiScaleR:::build_opt_context(
      fitted_mod = unmark$mod,
      cov_df = unmark$kernel_inputs$raw_cov,
      join_by = data.frame(not_site = seq_len(nrow(unmark$kernel_inputs$kernel_dat)))
    ),
    "join_by"
  )

  fix <- make_core_fixture()
  expect_error(
    with_mocked_bindings(
      multiScaleR:::kernel_scale_fn(
        par = 40 / fix$kernel_inputs$unit_conv,
        d_list = fix$kernel_inputs$d_list,
        cov_df = fix$kernel_inputs$raw_cov,
        kernel = fix$kernel_inputs$kernel,
        fitted_mod = fix$fitted_mod,
        mod_return = TRUE
      ),
      .refit_model = function(...) stop("mock refit failure"),
      .package = "multiScaleR"
    ),
    "Failed to refit the model with the optimized covariates"
  )
})

test_that("kernel_scale_fn returns consistent results with cached optimization context", {
  fix <- make_core_fixture()
  par <- 40 / fix$kernel_inputs$unit_conv
  glm_context <- multiScaleR:::build_opt_context(
    fitted_mod = fix$fitted_mod,
    cov_df = fix$kernel_inputs$raw_cov
  )

  glm_obj <- multiScaleR:::kernel_scale_fn(
    par = par,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = fix$fitted_mod
  )
  glm_obj_cached <- multiScaleR:::kernel_scale_fn(
    par = par,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = fix$fitted_mod,
    opt_context = glm_context
  )

  expect_equal(glm_obj_cached, glm_obj, tolerance = 1e-8)

  glm_mod <- multiScaleR:::kernel_scale_fn(
    par = par,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = fix$fitted_mod,
    mod_return = TRUE
  )
  glm_mod_cached <- multiScaleR:::kernel_scale_fn(
    par = par,
    d_list = fix$kernel_inputs$d_list,
    cov_df = fix$kernel_inputs$raw_cov,
    kernel = fix$kernel_inputs$kernel,
    fitted_mod = fix$fitted_mod,
    mod_return = TRUE,
    opt_context = glm_context
  )

  expect_equal(unname(stats::coef(glm_mod_cached$mod)),
               unname(stats::coef(glm_mod$mod)),
               tolerance = 1e-8)
  expect_equal(glm_mod_cached$scl_params$mean, glm_mod$scl_params$mean)
  expect_equal(glm_mod_cached$scl_params$sd, glm_mod$scl_params$sd)

  unmark <- make_unmarked_fixture()
  join_by <- data.frame(site = seq_len(nrow(unmark$kernel_inputs$kernel_dat)))
  um_context <- multiScaleR:::build_opt_context(
    fitted_mod = unmark$mod,
    cov_df = unmark$kernel_inputs$raw_cov,
    join_by = join_by
  )
  um_par <- c(40 / unmark$kernel_inputs$unit_conv, 60 / unmark$kernel_inputs$unit_conv)

  um_obj <- multiScaleR:::kernel_scale_fn(
    par = um_par,
    d_list = unmark$kernel_inputs$d_list,
    cov_df = unmark$kernel_inputs$raw_cov,
    kernel = unmark$kernel_inputs$kernel,
    fitted_mod = unmark$mod,
    join_by = join_by
  )
  um_obj_cached <- multiScaleR:::kernel_scale_fn(
    par = um_par,
    d_list = unmark$kernel_inputs$d_list,
    cov_df = unmark$kernel_inputs$raw_cov,
    kernel = unmark$kernel_inputs$kernel,
    fitted_mod = unmark$mod,
    join_by = join_by,
    opt_context = um_context
  )

  expect_equal(um_obj_cached, um_obj, tolerance = 1e-8)
})

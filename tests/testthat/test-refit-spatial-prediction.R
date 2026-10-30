test_that("native refit prediction supports weighted surfaces and all kernels", {
  for (kernel in c("gaussian", "exp", "expow")) {
    for (metric in c("sa", "sq", "ssk", "sku")) {
      f <- refit_spatial_fixture(kernel, metric)
      rows <- c(19, 3, 28, 9)
      nd <- f$data[rows, , drop = FALSE]
      expected <- predict(f$object$opt_mod, nd, type = "response")
      nd[c("mean", "rough", "fixed")] <- -999
      actual <- predict(f$object, nd, f$inputs, include_fixed = TRUE)
      expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-11)
      expect_identical(attr(actual, "prediction_contract")$uncertainty$kernel_parameters,
                       "fixed_at_estimate")
      expect_equal(as.numeric(predict(f$object, type = "link")),
        as.numeric(predict(f$object$opt_mod, type = "link")), tolerance = 1e-12)
    }
  }
})

test_that("hard-radius surface features retain native refit predictions", {
  for (metric in c("sa", "sq", "ssk", "sku", "sdq", "sdr")) {
    f <- refit_spatial_fixture(metric = metric, weighted = FALSE)
    rows <- c(2, 25, 14)
    expected <- predict(f$object$opt_mod, f$data[rows, ], type = "response")
    expect_equal(as.numeric(predict(f$object, f$data[rows, ], f$inputs,
      include_fixed = TRUE)), as.numeric(expected), tolerance = 1e-10)
  }
})

test_that("refit prediction preserves units, scenarios, and conditional uncertainty", {
  f <- refit_spatial_fixture()
  rows <- c(9, 3, 9, 18)
  nd <- f$data[rows, , drop = FALSE]
  keys <- data.frame(site_id = rownames(f$data))
  nd$site_id <- keys$site_id[rows]
  rownames(nd) <- paste0("scenario", seq_along(rows))
  expected <- predict(f$object$opt_mod, nd, type = "response", se.fit = TRUE)
  inputs <- f$inputs
  inputs$unit_conv <- inputs$unit_conv * 2
  inputs$d_list <- lapply(inputs$d_list, function(d) d / 2)
  if (!is.null(inputs$binned)) inputs$binned$bin_rep <- inputs$binned$bin_rep / 2
  inputs$scale_vars <- inputs$scale_vars[c(4, 2, 1, 3), ]
  actual <- predict(f$object, nd, inputs, join_by = keys, include_fixed = TRUE,
    se.fit = TRUE, uncertainty = "conditional")
  expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-11)
  expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-11)
  expect_equal(unname(actual$fit[1]), unname(actual$fit[3]))
  expect_error(predict(f$object, nd, inputs, se.fit = TRUE), "conditional")
  expect_error(predict(f$object, nd, inputs, interval = "confidence"), "conditional")
  expect_error(predict(f$object, nd), "together")
  nd$exposure <- NULL
  expect_error(predict(f$object, nd, inputs, join_by = keys), "explicitly.*exposure")
  bad <- inputs
  bad$scale_vars$source[bad$scale_vars$covariate == "rough"] <- "hab"
  expect_error(predict(f$object, f$data, bad), "definitions")
  bad <- inputs
  bad$max_D <- bad$max_D + 1
  expect_error(predict(f$object, f$data, bad), "buffer extent")
  bad_obj <- f$object
  bad_obj$optim_results$par[1] <- bad_obj$optim_results$par[1] + 1
  expect_error(predict(bad_obj, f$data, inputs), "disagree")
})

test_that("ordinary linear and GAM models keep native prediction semantics", {
  f <- refit_spatial_fixture()
  for (kind in c("lm", "gam")) {
    if (kind == "gam") skip_if_not_installed("mgcv")
    obj <- f$object
    obj$opt_mod <- if (kind == "lm") lm(y ~ mean + rough + fixed + group, f$data) else
      mgcv::gam(y ~ mean + rough + fixed + group, family = poisson(), data = f$data)
    rows <- c(3, 24, 11)
    nd <- f$data[rows, ]
    native_type <- if (kind == "lm") "response" else "link"
    expected <- predict(obj$opt_mod, nd, type = native_type, se.fit = TRUE)
    actual <- predict(obj, nd, f$inputs, type = "link", se.fit = TRUE,
                      uncertainty = "conditional", include_fixed = TRUE)
    expect_equal(as.numeric(actual$fit), as.numeric(expected$fit), tolerance = 1e-11)
    expect_equal(as.numeric(actual$se.fit), as.numeric(expected$se.fit), tolerance = 1e-11)
  }
})

test_that("a real refit optimizer supports point prediction without a second fit", {
  f <- make_core_fixture()
  rows <- c(8, 2, 13)
  nd <- .model_data(f$opt$opt_mod)[rows, , drop = FALSE]
  value <- predict(f$opt, nd, f$kernel_inputs)
  expect_equal(as.numeric(value), as.numeric(predict(f$opt$opt_mod, nd,
               type = "response")), tolerance = 1e-10)
})

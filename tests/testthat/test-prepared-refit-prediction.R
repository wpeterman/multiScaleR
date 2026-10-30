test_that("native refit snapshots freeze mixed-metric outputs without native fits", {
  f <- refit_spatial_fixture("expow")
  nd <- f$data[c(3, 25, 11), , drop = FALSE]
  snap <- prepare_prediction(f$object, nd, f$inputs, include_fixed = TRUE)
  expect_s3_class(snap, "multiScaleR_prediction")
  expect_identical(snap$native_kind, "refit")
  expect_named(snap$results, c("response", "link"))
  expect_null(snap$opt_mod)
  expect_null(snap$kernel_inputs)
  expect_null(snap$design)
  expect_identical(snap$provenance$row_names, rownames(nd))
  expect_identical(snap$provenance$spatial$source_specs, f$inputs$scale_vars)
  for (type in c("response", "link")) {
    expect_equal(predict(snap, type = type),
      predict(f$object, nd, f$inputs, type = type, include_fixed = TRUE))
    expect_equal(predict(snap, type = type, se.fit = TRUE, uncertainty = "conditional"),
      predict(f$object, nd, f$inputs, type = type, include_fixed = TRUE,
        se.fit = TRUE, uncertainty = "conditional"))
  }
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(snap, path)
  restored <- readRDS(path)
  expect_equal(predict(restored), predict(snap))
  f$object$opt_mod <- NULL
  f$inputs$raw_cov <- NULL
  expect_equal(predict(restored), predict(snap))
  expect_error(predict(snap, se.fit = TRUE), "conditional")
  expect_error(predict(snap, newdata = nd), "fixed inputs")
})

test_that("unavailable native SEs are retained explicitly without losing point output", {
  skip_if_not_installed("nlme")
  f <- refit_spatial_fixture()
  f$object$opt_mod <- nlme::gls(y ~ mean + rough + fixed, data = f$data,
                                method = "ML")
  nd <- f$data[c(7, 26), ]
  snap <- prepare_prediction(f$object, nd, f$inputs, component = "response")
  expect_equal(as.numeric(predict(snap)), as.numeric(predict(f$object$opt_mod, nd)))
  expect_match(snap$results$response$se$error, "did not supply se.fit")
  expect_false(snap$prediction_contract$response$uncertainty$conditional_se_available)
  expect_error(predict(snap, se.fit = TRUE, uncertainty = "conditional"),
                "did not supply se.fit")
  expect_error(prepare_prediction(f$object, nd, f$inputs, native_args = list(newdata = nd)),
                "supported native")
  expect_error(prepare_prediction(f$object, nd, f$inputs, max_bytes = 1), "max_bytes")
})

test_that("survival refit snapshots preserve named native targets", {
  skip_if_not_installed("survival")
  f <- refit_spatial_fixture()
  dat <- f$data
  dat$time <- seq_len(nrow(dat)) + as.numeric(dat$mean > 0) / 2
  dat$status <- rep(c(1, 1, 0), length.out = nrow(dat))
  f$object$opt_mod <- survival::coxph(survival::Surv(time, status) ~
    mean + rough + fixed, data = dat, model = TRUE)
  nd <- dat[c(15, 29, 2), ]
  snap <- prepare_prediction(f$object, nd, f$inputs, include_fixed = TRUE,
                             native_args = list(reference = "zero"))
  expect_named(snap$results, c("lp", "risk"))
  expect_equal(as.numeric(predict(snap, type = "lp")),
    as.numeric(predict(f$object$opt_mod, nd, type = "lp", reference = "zero")),
    tolerance = 1e-10)
  expect_equal(as.numeric(predict(snap, type = "risk")),
    as.numeric(predict(f$object$opt_mod, nd, type = "risk", reference = "zero")),
    tolerance = 1e-10)
  expect_error(predict(f$object, nd, f$inputs, type = "response"),
                "no generic response probability")
  expect_identical(snap$prediction_contract$lp$quantity, "native_component")
})

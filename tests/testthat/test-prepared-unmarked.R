test_that("N-mixture snapshots retain native tables and component controls", {
  skip_if_not_installed("unmarked")
  f <- make_unmarked_fixture()
  object <- f$obj
  for (component in c("state", "det")) {
    prepared <- prepare_prediction(object, component = component)
    native <- unmarked::predict(object$opt_mod, type = component)
    expect_equal(predict(prepared), native)
    expect_equal(predict(prepared, type = component), native)
    expect_equal(predict(unserialize(serialize(prepared, NULL))), native)
    expect_null(prepared$opt_mod)
    expect_null(prepared$kernel_inputs)
    expect_error(predict(prepared, se.fit = TRUE), "native prediction tables")
    expect_error(predict(prepared, level = .9), "native prediction tables")
    expect_error(predict(prepared, backTransform = FALSE), "fixed inputs")
  }
  nd <- object$opt_mod@data@siteCovs
  nd <- nd[nrow(nd):1, , drop = FALSE]
  prepared <- prepare_prediction(object, nd, component = "state",
    native_args = list(backTransform = FALSE, level = .9, appendData = TRUE))
  expected <- unmarked::predict(object$opt_mod, type = "state", newdata = nd,
    backTransform = FALSE, level = .9, appendData = TRUE)
  expect_equal(predict(prepared), expected)
  nd$cont1 <- 999
  object$opt_mod@data@siteCovs$cont1 <- 999
  expect_equal(predict(prepared), expected)
  expect_error(prepare_prediction(object), "component explicitly")
  expect_error(prepare_prediction(object, kernel_inputs = f$kernel_inputs,
                                 component = "state"), "newdata must")
  expect_error(prepare_prediction(object, component = "state", max_bytes = 1), "max_bytes")
  expect_error(prepare_prediction(object, component = "state",
                                 native_args = list(unknown = TRUE)), "supported")
})

test_that("occupancy and dynamic occupancy snapshots retain observation and year contracts", {
  skip_if_not_installed("unmarked")
  e <- new.env()
  utils::data("crossbill", package = "unmarked", envir = e)
  dat <- e$crossbill
  y <- as.matrix(dat[, grep("^det", names(dat))])
  date <- as.matrix(dat[, grep("^date", names(dat))]) / 100
  sc <- data.frame(forest = as.numeric(scale(dat$forest)))
  annual <- unmarked::unmarkedMultFrame(y, siteCovs = sc,
    obsCovs = list(date = date), numPrimary = ncol(y) / 3)
  dynamic <- unmarked::colext(~forest, ~1, ~1, ~date, data = annual)
  single <- unmarked::unmarkedFrameOccu(y[, 1:3], siteCovs = sc,
    obsCovs = list(date = date[, 1:3]))
  expect_warning(occupancy <- unmarked::occu(~date ~forest, data = single),
                 "dropped due to missing values")
  for (model in list(occupancy, dynamic)) {
    object <- structure(list(opt_mod = model), class = "multiScaleR")
    components <- if (methods::is(model, "unmarkedFitColExt"))
      c("psi", "col", "ext", "det") else c("state", "det")
    for (component in components) {
      expected <- unmarked::predict(model, type = component, na.rm = FALSE)
      prepared <- prepare_prediction(object, component = component,
                                     native_args = list(na.rm = FALSE))
      expect_equal(predict(prepared), expected)
      expect_equal(predict(unserialize(serialize(prepared, NULL))), expected)
    }
  }
})

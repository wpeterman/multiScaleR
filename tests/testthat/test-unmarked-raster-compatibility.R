test_that("unmarked S4 raster preparation discovers formulas without generic dollar access", {
  skip_if_not_installed("unmarked")
  f <- refit_unmarked_spatial_fixture()
  formulas <- .unmarked_model_formulas(f$object$opt_mod)
  expect_true(all(vapply(formulas, inherits, logical(1), "formula")))
  expect_setequal(.unmarked_model_predictors(f$object$opt_mod), c("bin1", "cont1"))
  legacy <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, verbose = FALSE)
  expect_s4_class(legacy, "SpatRaster")
  inputs <- kernel_scale.raster(f$rasters, multiScaleR = f$object,
    scale_center = TRUE, match_buffer = TRUE, verbose = FALSE)
  expect_s4_class(inputs, "SpatRaster")
  expect_error(predict(f$object, inputs), "component explicitly")
  data <- as.data.frame(inputs, na.rm = FALSE)
  expected <- unmarked::predict(f$object$opt_mod, newdata = data, type = "state")
  actual <- predict(f$object, inputs, type = "state")
  expect_equal(terra::values(actual)[, 1L], expected$Predicted, tolerance = 1e-10)
  uncertain <- predict(f$object, inputs, type = "state", se.fit = TRUE,
                        uncertainty = "conditional")
  expect_equal(unname(terra::values(uncertain)), as.matrix(unname(expected)),
    ignore_attr = TRUE, tolerance = 1e-10)
})

test_that("double unmarked formulas retain state and detection offset inputs", {
  offsets <- .msr_raster_offset_terms(~offset(log(effort)) ~
    habitat + offset(log(area * duration)))
  expect_length(offsets, 2L)
  expect_identical(lapply(offsets, all.vars), list("effort", c("area", "duration")))
})

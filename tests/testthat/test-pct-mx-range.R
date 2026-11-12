# `pct_mx` moves each clamping limit by a fraction of the observed covariate
# range: lower = min - pct_mx * (max - min), upper = max + pct_mx * (max - min).
# The earlier rule multiplied each limit by (1 + pct_mx), which is wrong for a
# range that spans zero (the standardized covariates used for clamping always do).

test_that("pct_mx clamping limits follow the observed range for a range spanning zero", {
  fix <- make_core_fixture()

  # Standardized covariate at the survey points, the scale clamping works on.
  obs <- range(fix$opt$opt_mod$model$cont1)
  expect_lt(obs[1], 0)
  expect_gt(obs[2], 0)
  width <- diff(obs)

  # Large blocks far below and far above the surveyed values, so the smoothed,
  # standardized raster reaches both clamping limits for every `pct_mx`.
  r <- fix$rs
  vals <- rep(c(-1e4, 1e4), each = terra::ncell(r) / 2)
  terra::values(r) <- vals

  limits <- function(pct_mx) {
    out <- kernel_scale.raster(raster_stack = r,
                               multiScaleR = fix$opt,
                               scale_center = TRUE,
                               clamp = TRUE,
                               pct_mx = pct_mx,
                               verbose = FALSE)
    range(terra::values(out[["cont1"]], mat = FALSE), na.rm = TRUE)
  }

  # Zero keeps the exact observed range.
  expect_equal(limits(0), obs)
  # Positive values move both limits outward by the same amount.
  expect_equal(limits(0.1), c(obs[1] - 0.1 * width, obs[2] + 0.1 * width))
  # Negative values move both limits inward by the same amount.
  expect_equal(limits(-0.1), c(obs[1] + 0.1 * width, obs[2] - 0.1 * width))
})

test_that("pct_mx applies the same rule to prepared (kernel_prep) objects", {
  fix <- make_core_fixture()
  obs <- range(fix$kernel_inputs$kernel_dat$cont1)
  width <- diff(obs)

  r <- fix$rs
  terra::values(r) <- rep(c(-1e4, 1e4), each = terra::ncell(r) / 2)
  out <- kernel_scale.raster(raster_stack = r,
                             multiScaleR = fix$kernel_inputs,
                             sigma = 40,
                             scale_center = TRUE,
                             clamp = TRUE,
                             pct_mx = 0.1,
                             verbose = FALSE)
  expect_equal(range(terra::values(out[["cont1"]], mat = FALSE), na.rm = TRUE),
               c(obs[1] - 0.1 * width, obs[2] + 0.1 * width))
})

test_that("pct_mx that would cross the clamping limits is rejected", {
  fix <- make_core_fixture()

  for (bad in c(-0.5, -0.75)) {
    expect_error(
      kernel_scale.raster(raster_stack = fix$rs,
                          multiScaleR = fix$opt,
                          scale_center = TRUE,
                          clamp = TRUE,
                          pct_mx = bad,
                          verbose = FALSE),
      "`pct_mx` must be greater than -0.5"
    )
  }
  expect_error(
    multiScaleR:::scale_center_raster(fix$rs, fix$opt, clamp = TRUE, pct_mx = -0.5),
    "`pct_mx` must be greater than -0.5"
  )
  expect_error(
    kernel_scale.raster(raster_stack = fix$rs, sigma = 40, pct_mx = NA_real_,
                        verbose = FALSE),
    "pct_mx"
  )
  # Just inside the bound is accepted.
  expect_no_error(
    kernel_scale.raster(raster_stack = fix$rs,
                        multiScaleR = fix$opt,
                        scale_center = TRUE,
                        clamp = TRUE,
                        pct_mx = -0.49,
                        verbose = FALSE)
  )
})

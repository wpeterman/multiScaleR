# `base` does not apply to the modified Simpson metrics, and relative patch
# richness needs `classes_max`.

fix_cover_layer <- function() {
  land <- terra::rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
  cover <- terra::classify(land$land3,
    rcl = matrix(c(-Inf, .45, 1, .45, .55, 2, .55, Inf, 3), ncol = 3, byrow = TRUE))
  names(cover) <- "cover"
  cover
}

fix_cover_mean <- function(spec) {
  data("surv_pts", package = "multiScaleR", envir = environment())
  prep <- kernel_prep(pts = surv_pts, raster_stack = fix_cover_layer(), max_D = 800,
    scale_vars = msr_vars(m = spec), verbose = FALSE)
  unlist(prep$scl_params$mean, use.names = FALSE)
}

test_that("landscape_var warns that base does not apply to msidi and msiei", {
  for (metric in c("msidi", "msiei")) {
    expect_warning(
      spec <- landscape_var("cover", metric = metric, radius = 500, base = 2),
      "always uses the natural log"
    )
    # The stored specification records the base actually used.
    expect_equal(spec$base, exp(1))
    expect_silent(landscape_var("cover", metric = metric, radius = 500))
    expect_silent(landscape_var("cover", metric = metric, radius = 500, base = exp(1)))
  }
  # Metrics that use the base are unaffected.
  expect_silent(spec <- landscape_var("cover", metric = "shdi", radius = 500, base = 2))
  expect_equal(spec$base, 2)
})

test_that("msidi values are unchanged and shdi still responds to base", {
  skip_on_cran()
  default <- fix_cover_mean(landscape_var("cover", metric = "msidi", radius = 500))
  with_base <- fix_cover_mean(suppressWarnings(
    landscape_var("cover", metric = "msidi", radius = 500, base = 2)))
  expect_identical(default, with_base)

  shdi_e <- fix_cover_mean(landscape_var("cover", metric = "shdi", radius = 500))
  shdi_2 <- fix_cover_mean(landscape_var("cover", metric = "shdi", radius = 500, base = 2))
  expect_equal(shdi_2, shdi_e / log(2))
})

test_that("rpr requires classes_max when the covariate is defined", {
  expect_error(landscape_var("cover", metric = "rpr", radius = 500), "`classes_max` is required")
  expect_error(landscape_var("cover", metric = "rpr"), "`classes_max` is required")
  spec <- landscape_var("cover", metric = "rpr", radius = 500, classes_max = 5)
  expect_equal(spec$classes_max, 5)
  # With a maximum supplied the values are defined. Every 500 m buffer on this
  # layer holds all three classes, so the existing zero-variance guard is the
  # correct outcome, not the former "all NA" failure.
  expect_error(fix_cover_mean(spec), "zero variance")
  # Other metrics still do not need it.
  expect_silent(landscape_var("cover", metric = "pr", radius = 500))
})

test_that("landscape_var warns when base has no effect on the metric", {
  # Metrics without a logarithm ignore base.
  for (metric in c("sidi", "pr", "ed", "ai", "pladj")) {
    expect_warning(
      spec <- landscape_var("cover", metric = metric, radius = 500, base = 2),
      "does not apply to"
    )
    expect_equal(spec$base, exp(1))
    expect_silent(landscape_var("cover", metric = metric, radius = 500))
  }
  expect_warning(landscape_var("cover", metric = "pland", radius = 500, base = 10, class = 1),
                 "does not apply to `pland`")
  # Entropy metrics use base, and log-ratio metrics accept it silently.
  for (metric in c("shdi", "ent", "shei", "iji", "contag", "relmutinf")) {
    expect_silent(spec <- landscape_var("cover", metric = metric, radius = 500, base = 2))
    expect_equal(spec$base, 2)
  }
})

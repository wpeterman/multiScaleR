test_that("simulation helpers return expected structures", {
  fix <- make_simulation_fixture()

  expect_true(inherits(fix$rs, "SpatRaster"))
  expect_equal(names(fix$rs), c("bin1", "bin2", "cont1", "cont2"))

  rs_vals <- terra::values(fix$rs, mat = FALSE)
  expect_gte(min(rs_vals, na.rm = TRUE), 0)
  expect_lte(max(rs_vals, na.rm = TRUE), 1)

  expect_equal(names(fix$sim), c("obs", "df", "pts"))
  expect_equal(dim(fix$sim$df), c(15, 3))
  expect_false(anyNA(fix$sim$df))

  expect_equal(names(fix$sim_umf), c("y", "df", "pts"))
  expect_equal(dim(fix$sim_umf$y), c(15, 4))
  expect_false(anyNA(fix$sim_umf$y))
  expect_false(anyNA(fix$sim_umf$df))
})

test_that("aic_tab and bic_tab compare optimized and standard models", {
  fix <- make_core_fixture()

  aic_res <- aic_tab(list(fix$opt, fix$plain_mod))
  bic_res <- bic_tab(list(fix$opt, fix$plain_mod))

  expect_s3_class(aic_res, "aictab")
  expect_s3_class(bic_res, "bictab")
  expect_equal(nrow(aic_res), 2)
  expect_equal(nrow(bic_res), 2)
})

test_that("aic_tab and bic_tab count one optimized parameter per spatial term", {
  fix <- make_core_fixture()

  opt_multi <- fix$opt
  opt_multi$scale_est <- data.frame(
    Mean = c(40, 80),
    SE = c(1, 1),
    row.names = c("cont1", "site")
  )

  expected_k <- nrow(insight::get_parameters(opt_multi$opt_mod)) +
    nrow(opt_multi$scale_est)

  aic_res <- aic_tab(
    list(opt_multi, fix$plain_mod),
    AICc = FALSE,
    mod_names = c("multi", "plain")
  )
  bic_res <- bic_tab(
    list(opt_multi, fix$plain_mod),
    mod_names = c("multi", "plain")
  )

  expect_equal(aic_res$K[aic_res$Modnames == "multi"], expected_k)
  expect_equal(bic_res$K[bic_res$Modnames == "multi"], expected_k)
})

test_that("aic_tab rejects models fitted to different sample sizes", {
  fix <- make_core_fixture()
  smaller_mod <- glm(y ~ site, family = poisson(), data = fix$df[-1, ])

  expect_error(
    aic_tab(list(fix$opt, smaller_mod)),
    "different number of sample locations"
  )
})

test_that("aic_tab and bic_tab reject models fit to different observation sets", {
  fix <- make_core_fixture()
  mod_a <- glm(y ~ site, family = poisson(), data = fix$df[-1, ])
  mod_b <- glm(y ~ site, family = poisson(), data = fix$df[-nrow(fix$df), ])

  expect_error(
    aic_tab(list(mod_a, mod_b)),
    "different observation sets"
  )
  expect_error(
    bic_tab(list(mod_a, mod_b)),
    "different observation sets"
  )
})

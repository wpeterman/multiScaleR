test_that("unmarked site keys preserve detection-history row order", {
  fix <- make_unmarked_fixture()
  n <- nrow(fix$site_covs)
  order <- c(seq(n, 1), seq_len(n), seq(n, 1))
  site_covs <- fix$site_covs[order, , drop = FALSE]
  y <- fix$umf@y[order, , drop = FALSE]
  umf <- unmarked::unmarkedFramePCount(y = y, siteCovs = site_covs)
  mod <- unmarked::pcount(~1 ~ bin1 + cont1, data = umf, K = 50)
  join_by <- data.frame(site = seq_len(n))
  ctx <- multiScaleR:::build_opt_context(mod, fix$kernel_inputs$raw_cov,
                                          join_by = join_by)
  expect_equal(ctx$unmarked_site_idx, order)
  expect_error(multiScaleR:::build_opt_context(mod,
               fix$kernel_inputs$raw_cov), "Supply `join_by`")
  expect_error(multiScale_optim(mod, fix$kernel_inputs,
               par = c(40, 60) / fix$kernel_inputs$unit_conv,
               verbose = FALSE), "Supply `join_by`")

  par <- c(40, 60) / fix$kernel_inputs$unit_conv
  raw <- multiScaleR:::.msr_kernel_cov_w(par, fix$kernel_inputs$d_list,
                                         fix$kernel_inputs$raw_cov,
                                         fix$kernel_inputs$kernel, ctx)
  result <- multiScaleR:::kernel_scale_fn(par, fix$kernel_inputs$d_list,
                fix$kernel_inputs$raw_cov, fix$kernel_inputs$kernel,
                mod, join_by = join_by, opt_context = ctx, mod_return = TRUE)
  expect_equal(result$mod@data@y, y)
  expect_equal(result$mod@data@siteCovs$site, site_covs$site)
  expect_equal(result$mod@data@siteCovs$nuisance, site_covs$nuisance)
  expect_equal(as.matrix(result$mod@data@siteCovs[c("bin1", "cont1")]),
               unname(scale(raw)[order, , drop = FALSE]),
               ignore_attr = TRUE)
})

test_that("unmarked join keys reject ambiguous and incomplete matches", {
  fix <- make_unmarked_fixture()
  site_covs <- fix$site_covs
  n <- nrow(site_covs)
  map <- multiScaleR:::.unmarked_site_map
  keys <- data.frame(site = seq_len(n))
  expect_equal(map(site_covs, NULL, n), seq_len(n))
  expect_error(map(site_covs, NULL, n - 1L), "Supply `join_by`")
  reordered <- site_covs[c(2, 1, seq.int(3, n)), , drop = FALSE]
  expect_error(map(reordered, NULL, n, prepared_ids = as.character(seq_len(n))),
               "row IDs do not match")
  expect_error(map(site_covs, data.frame(), n), "one or more")
  expect_error(map(site_covs, keys[-1, , drop = FALSE], n), "exactly one row")
  bad <- keys
  bad$site[2] <- bad$site[1]
  expect_error(map(site_covs, bad, n), "unique key")
  bad <- keys
  bad$site[1] <- NA_integer_
  expect_error(map(site_covs, bad, n), "missing values")
  bad <- site_covs
  bad$site[1] <- NA_integer_
  expect_error(map(bad, keys, n), "missing values")
  bad <- site_covs
  bad$site[1] <- n + 1L
  expect_error(map(bad, keys, n), "must match exactly one")
  expect_error(map(site_covs, data.frame(other = seq_len(n)), n),
               "also be present")
  multi <- data.frame(site = seq_len(n), group = rep(c("a", "b"), length.out = n))
  site_covs$group <- multi$group
  expect_equal(map(site_covs, multi, n), seq_len(n))
  expect_error(multiScaleR:::build_opt_context(fix$mod,
               fix$kernel_inputs$raw_cov,
               join_by = data.frame(bin1 = site_covs$bin1)),
               "key columns cannot be optimized")
})

test_that("site maps handle year-major, site-major, and missing site-years", {
  prepared <- data.frame(site = c("A", "B", "C"))
  year_major <- data.frame(site = c("A", "B", "C", "A", "B", "C"))
  site_major <- data.frame(site = c("A", "A", "B", "B", "C", "C"))
  incomplete <- data.frame(site = c("B", "A", "C", "B"))
  map <- multiScaleR:::.unmarked_site_map
  expect_equal(map(year_major, prepared, 3L), c(1L, 2L, 3L, 1L, 2L, 3L))
  expect_equal(map(site_major, prepared, 3L), c(1L, 1L, 2L, 2L, 3L, 3L))
  expect_equal(map(incomplete, prepared, 3L), c(2L, 1L, 3L, 2L))
})

test_that("unmarked open-population models use component formulas and static site covariates", {
  set.seed(2)
  y <- matrix(stats::rpois(60, 3), nrow = 10)
  site_covs <- data.frame(site = seq_len(10), x = as.numeric(scale(seq_len(10))))
  umf <- unmarked::unmarkedFramePCO(y = y, siteCovs = site_covs,
                                    numPrimary = 3)
  mod <- unmarked::pcountOpen(~x, ~1, ~1, ~1, data = umf, K = 20,
                              se = FALSE)
  expect_equal(multiScaleR:::.unmarked_model_predictors(mod), "x")
  nested <- mod
  nested@formlist <- list(state = list(lambda = ~x), detection = list(~1))
  expect_equal(multiScaleR:::.unmarked_model_predictors(nested), "x")
  cov_df <- lapply(seq_len(10), function(i) {
    out <- matrix(c(i, i + 1), ncol = 1)
    colnames(out) <- "x"
    out
  })
  ctx <- multiScaleR:::build_opt_context(mod, cov_df)
  expect_equal(ctx$covs, "x")
  result <- multiScaleR:::kernel_scale_fn(1, rep(list(c(0, 1)), 10),
                cov_df, "gaussian", mod, opt_context = ctx, mod_return = TRUE)
  expect_s4_class(result$mod, "unmarkedFitPCO")
  expect_equal(result$mod@data@y, y)
  expect_equal(result$mod@data@siteCovs$site, site_covs$site)
})

test_that("unmarked dynamic occupancy models use component formulas", {
  set.seed(42)
  y <- matrix(stats::rbinom(60, 1, 0.5), nrow = 10)
  site_covs <- data.frame(site = seq_len(10),
                          x = as.numeric(scale(seq_len(10))))
  umf <- unmarked::unmarkedMultFrame(y = y, siteCovs = site_covs,
                                     numPrimary = 3)
  mod <- unmarked::colext(~x, ~1, ~1, ~1, data = umf, se = FALSE)
  expect_equal(multiScaleR:::.unmarked_model_predictors(mod), "x")
  cov_df <- lapply(seq_len(10), function(i) {
    out <- matrix(c(i, i + 1), ncol = 1)
    colnames(out) <- "x"
    out
  })
  ctx <- multiScaleR:::build_opt_context(mod, cov_df)
  result <- multiScaleR:::kernel_scale_fn(1, rep(list(c(0, 1)), 10),
                cov_df, "gaussian", mod, opt_context = ctx, mod_return = TRUE)
  expect_s4_class(result$mod, "unmarkedFitColExt")
  expect_equal(result$mod@data@y, y)
})

test_that("unmarked refits cannot silently change the modeled site set", {
  fix <- make_unmarked_fixture()
  par <- c(40, 60) / fix$kernel_inputs$unit_conv
  changed_sites <- function(model, data, context) {
    model@sitesRemoved <- 1L
    model
  }
  ctx <- multiScaleR:::build_opt_context(fix$mod,
                fix$kernel_inputs$raw_cov, refit_fn = changed_sites)
  expect_equal(multiScaleR:::kernel_scale_fn(par,
                 fix$kernel_inputs$d_list, fix$kernel_inputs$raw_cov,
                 fix$kernel_inputs$kernel, fix$mod, opt_context = ctx),
               1e6^10)
  expect_error(multiScaleR:::kernel_scale_fn(par,
                 fix$kernel_inputs$d_list, fix$kernel_inputs$raw_cov,
                 fix$kernel_inputs$kernel, fix$mod, opt_context = ctx,
                 mod_return = TRUE), "changed which sites")
})

test_that("site covariates must contain the optimized unmarked predictor", {
  set.seed(11)
  y <- matrix(stats::rpois(30, 2), nrow = 10)
  umf <- unmarked::unmarkedFramePCount(
    y = y, siteCovs = data.frame(site = seq_len(10)),
    obsCovs = list(z = matrix(seq_len(30), nrow = 10)))
  mod <- unmarked::pcount(~z ~1, data = umf, K = 20, se = FALSE)
  cov_df <- rep(list(matrix(1:2, ncol = 1,
                            dimnames = list(NULL, "z"))), 10)
  expect_error(multiScaleR:::build_opt_context(mod, cov_df),
               "must be static `siteCovs`")
})

test_that("multiScale_optim completes an open-population unmarked workflow", {
  fix <- make_unmarked_fixture()
  set.seed(2)
  y <- matrix(stats::rpois(90, 3), nrow = 15)
  umf <- unmarked::unmarkedFramePCO(y = y, siteCovs = fix$site_covs,
                                    numPrimary = 3)
  mod <- unmarked::pcountOpen(~cont1, ~1, ~1, ~1, data = umf,
                              K = 20, se = FALSE)
  expect_warning(
    result <- multiScale_optim(mod, fix$kernel_inputs,
                 par = 40 / fix$kernel_inputs$unit_conv, verbose = FALSE),
    NA
  )
  expect_s4_class(result$opt_mod, "unmarkedFitPCO")
  expect_equal(result$opt_mod@data@y, y)
  expect_equal(multiScaleR:::.msr_model_nobs(result$opt_mod), 15L)
  expect_equal(diagnostics(result)$sample_size$prepared_sites, 15L)
  expect_no_error(summary(result))
})

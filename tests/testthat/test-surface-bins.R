# Linearly binned power sums for kernel-weighted surface metrics (sq, ssk,
# sku). See `.msr_build_surface_bins()` in R/kernel_bins.R.

# Exact weighted metric from every cell, the reference for the binned path.
exact_surface_metric <- function(d, v, sigma, metric, kernel = "gaussian", shape = NULL) {
  .surface_weighted_by_buffer(d, matrix(v, ncol = 1, dimnames = list(NULL, "elev")),
                              kernel, sigma, shape, metric, validate = FALSE)[[1]]
}

test_that("binned power sums reproduce exact weighted moments on the node grid", {
  # Every distance sits on a node, so linear binning is exact and only the
  # algebra from centered power sums to central moments is tested.
  set.seed(11)
  nodes <- 0:20
  d_list <- list(a = sample(0:20, 400, TRUE), b = sample(0:20, 300, TRUE))
  values <- list(matrix(rnorm(400, 1e4, 3), ncol = 1, dimnames = list(NULL, "elev")),
                 matrix(rexp(300), ncol = 1, dimnames = list(NULL, "elev")))
  values[[2]][c(5, 9)] <- NA_real_
  specs <- data.frame(covariate = c("rough", "skew", "kurt"), source = "elev",
                      metric = c("sq", "ssk", "sku"))
  sb <- .msr_build_surface_bins(d_list, values, specs, nodes = nodes)
  expect_identical(sb$nbins, 20L)
  expect_identical(sb$point_ids, c("a", "b"))
  expect_length(sb$sums$elev, 5L)
  for (kernel in c("gaussian", "exp", "expow")) for (sigma in c(1.5, 4, 12)) {
    shape <- if (kernel == "expow") 1.6 else NULL
    got <- .msr_binned_surface_metrics(sb, specs$covariate, specs$covariate,
      rep(sigma, 3), rep(1.6, 3), kernel, 1:2)
    for (j in seq_len(3)) for (i in 1:2) {
      want <- exact_surface_metric(d_list[[i]], values[[i]][, 1], sigma,
                                   specs$metric[j], kernel, shape)
      expect_equal(unname(got[i, j]), unname(want), tolerance = 1e-10,
                   info = paste(kernel, sigma, specs$metric[j], i))
    }
  }
  # Zero spread gives SQ 0 and missing skewness and kurtosis, as the exact path.
  flat <- .msr_build_surface_bins(list(z = c(0, 1, 2)),
    list(matrix(5, 3, 1, dimnames = list(NULL, "elev"))), specs, nodes = 0:3)
  expect_equal(unname(.msr_binned_surface_metrics(flat, specs$covariate,
    specs$covariate, rep(1, 3), NULL, "gaussian", 1L)[1, ]), c(0, NA, NA))
  # Rows are matched by point ID; unknown IDs fall back to the exact path.
  expect_identical(.msr_surface_bin_rows(sb, c("b", "a"), 2L), 2:1)
  expect_null(.msr_surface_bin_rows(sb, c("b", "zz"), 2L))
  expect_identical(.msr_surface_bin_rows(sb, NULL, 2L), 1:2)
  expect_null(.msr_surface_bin_rows(sb, NULL, 1L))
  # Subsetting keeps only the sources still in use.
  sub <- .msr_subset_surface_bins(sb, "rough")
  expect_identical(sub$covariates, "rough")
  expect_null(.msr_subset_surface_bins(sb, "other"))
})

surface_bin_fixture <- function(metric = "sq", bin = TRUE, kernel = "gaussian") {
  r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 400,
                   ymin = 0, ymax = 400, crs = "EPSG:3857")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- 50 + xy[, 1] / 9 + sin(xy[, 2] / 13) * (1 + xy[, 1] / 60) +
    cos(xy[, 1] / 7 + xy[, 2] / 17)
  names(r) <- "elev"
  points <- expand.grid(x = seq(110, 290, 30), y = seq(110, 290, 30))
  pts <- terra::vect(as.matrix(points), crs = terra::crs(r))
  kernel_prep(pts, r, max_D = 100, kernel = kernel, bin = bin, verbose = FALSE,
              scale_vars = msr_vars(rough = surface_var("elev", metric, weighted = TRUE)))
}

test_that("kernel_prep() bins weighted sq, ssk and sku within kernel-mean accuracy", {
  skip_on_cran()
  for (metric in c("sq", "ssk", "sku")) {
    inputs <- surface_bin_fixture(metric)
    sb <- inputs$surface_binned
    expect_false(is.null(sb))
    expect_identical(sb$covariates, "rough")
    expect_identical(sb$point_ids, names(inputs$d_list))
    n <- length(inputs$d_list)
    # Kernel scales from about 2 to 8 cells, in the internal units of d_list.
    for (sigma_m in c(20, 40, 80)) {
      sigma <- sigma_m / inputs$unit_conv
      got <- .msr_binned_surface_metrics(sb, "rough", "rough", sigma, NULL,
                                         "gaussian", seq_len(n))[, 1]
      want <- vapply(seq_len(n), function(i) exact_surface_metric(inputs$d_list[[i]],
        inputs$raw_cov[[i]][, "elev"], sigma, metric), 0)
      # Scaled by the across-point spread, as the standardized covariate is.
      expect_lt(max(abs(got - want)) / stats::sd(want), 0.01)
    }
  }
  expect_null(surface_bin_fixture("sq", bin = FALSE)$surface_binned)
  expect_null(surface_bin_fixture("sa")$surface_binned)
  expect_null(surface_bin_fixture("sq", kernel = "fixed")$surface_binned)
  expect_output(print(surface_bin_fixture("sq")), "Binned weighted surface metrics: rough")
})

test_that("refit optimization uses binned surface moments and agrees with exact fits", {
  skip_on_cran()
  fits <- lapply(c(TRUE, FALSE), function(bin) {
    inputs <- surface_bin_fixture("sq", bin = bin)
    dat <- inputs$kernel_dat
    dat$y <- c(4, 1, 0, 3, 6, 2, 1)[(seq_len(nrow(dat)) %% 7L) + 1L] +
      as.integer(dat$rough > 0) * 3L
    model <- stats::glm(y ~ rough, data = dat, family = stats::poisson())
    fit <- suppressWarnings(multiScale_optim(model, inputs, par = 30 / inputs$unit_conv,
                                             n_cores = 1L, verbose = FALSE))
    list(inputs = inputs, fit = fit)
  })
  binned <- fits[[1]]$fit
  exact <- fits[[2]]$fit
  expect_false(is.null(binned$opt_context$surface_binned))
  expect_null(exact$opt_context$surface_binned)
  expect_equal(binned$scale_est$Mean, exact$scale_est$Mean, tolerance = 0.02)
  expect_equal(as.numeric(logLik(binned$opt_mod)), as.numeric(logLik(exact$opt_mod)),
               tolerance = 1e-3)

  # Evaluating a reordered subset of points gives the same rows as the full set.
  inputs <- fits[[1]]$inputs
  ctx <- binned$opt_context
  par <- binned$optim_results$par
  full <- .msr_kernel_cov_w(par, inputs$d_list, inputs$raw_cov, inputs$kernel, ctx)
  ids <- rev(names(inputs$d_list))[1:10]
  part <- .msr_kernel_cov_w(par, inputs$d_list[ids], inputs$raw_cov[ids],
                            inputs$kernel, ctx)
  expect_equal(unname(part[, "rough"]), unname(full[ids, "rough"]))
  # Cell lists whose IDs the bins do not know use the exact path instead.
  renamed <- stats::setNames(inputs$d_list[ids], paste0("new", seq_along(ids)))
  cells <- stats::setNames(inputs$raw_cov[ids], names(renamed))
  exact_part <- .msr_kernel_cov_w(par, renamed, cells, inputs$kernel, ctx)
  want <- vapply(ids, function(i) exact_surface_metric(inputs$d_list[[i]],
    inputs$raw_cov[[i]][, "elev"], par[[1]], "sq"), 0)
  expect_equal(unname(exact_part[, "rough"]), unname(want), tolerance = 1e-12)

  # Cross-validation folds rebuild the sums on the prepared node grid.
  fold <- .msr_cv_inputs(inputs, ids, "rough")
  expect_identical(fold$surface_binned$nodes, inputs$surface_binned$nodes)
  expect_identical(fold$surface_binned$point_ids, ids)
  rows <- match(ids, inputs$surface_binned$point_ids)
  expect_equal(fold$surface_binned$sums$elev[[3]],
               inputs$surface_binned$sums$elev[[3]][rows, , drop = FALSE])

  # Prediction inputs must be binned the same way as the fit.
  expect_error(predict(binned, binned$opt_mod$data, fits[[2]]$inputs),
               "binned versus exact")
})

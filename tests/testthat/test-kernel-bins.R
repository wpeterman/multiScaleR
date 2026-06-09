test_that("binned kernel-weighted mean matches the exact per-cell result", {
  set.seed(1)
  d <- runif(2000, 0, 1)
  v <- rnorm(2000, 5, 2)
  sparse <- methods::as(matrix(v, ncol = 1, dimnames = list(NULL, "hab")),
                        "sparseMatrix")

  binned <- .msr_build_kernel_bins(
    d_list = list(`1` = d),
    value_list = list(sparse),
    covariates = "hab",
    sources = "hab",
    nbins = 256L,
    point_ids = "1"
  )

  for (sig in c(0.05, 0.1, 0.25, 0.5)) {
    exact <- scale_type(d, kernel = "gaussian", sigma = sig, shape = NULL,
                        r_stack.df = sparse)[[1]]
    approx <- .msr_binned_kernel_means(
      binned = binned, covs = "hab", param_covariates = "hab",
      sigma = sig, shape = NULL, kernel = "gaussian", eval_idx = NULL
    )[1, 1]
    expect_equal(unname(approx), unname(exact), tolerance = 1e-3,
                 info = paste("sigma =", sig))
  }
})

test_that("binned optimization recovers the same sigma as the exact path", {
  rs <- terra::subset(sim_rast(dim = 60, resolution = 10, user_seed = 7), "cont1")
  names(rs) <- "cov"
  terra::crs(rs) <- "EPSG:32617"
  e <- as.vector(terra::ext(rs))
  lo <- e[1] + 0.12 * (e[2] - e[1]); hi <- e[2] - 0.12 * (e[2] - e[1])
  set.seed(3)
  pts <- terra::vect(cbind(stats::runif(40, lo, hi), stats::runif(40, lo, hi)),
                     type = "points")

  ki_exact <- kernel_prep(pts, rs, max_D = 150, kernel = "gaussian",
                          verbose = FALSE, bin = FALSE)
  ki_bin <- kernel_prep(pts, rs, max_D = 150, kernel = "gaussian",
                        verbose = FALSE, bin = TRUE, nbins = 256)

  expect_null(ki_exact$binned)
  expect_false(is.null(ki_bin$binned))
  expect_true(ki_bin$cell_data_stored)

  y <- stats::rpois(40, exp(0.4 + 0.6 * as.numeric(scale(ki_exact$kernel_dat$cov))))
  fit <- function(ki) stats::glm(y ~ cov, family = stats::poisson(),
                                 data = data.frame(y = y, ki$kernel_dat))

  oe <- suppressWarnings(suppressMessages(
    multiScale_optim(fit(ki_exact), ki_exact, verbose = FALSE)))
  ob <- suppressWarnings(suppressMessages(
    multiScale_optim(fit(ki_bin), ki_bin, verbose = FALSE)))

  expect_equal(ob$scale_est[1, 1], oe$scale_est[1, 1],
               tolerance = 0.02 * oe$scale_est[1, 1])
  expect_equal(as.numeric(stats::logLik(ob$opt_mod)),
               as.numeric(stats::logLik(oe$opt_mod)), tolerance = 1e-2)
})

test_that("lean mode drops cell data, still optimizes, and matches binned", {
  rs <- terra::subset(sim_rast(dim = 60, resolution = 10, user_seed = 7), "cont1")
  names(rs) <- "cov"
  terra::crs(rs) <- "EPSG:32617"
  e <- as.vector(terra::ext(rs))
  lo <- e[1] + 0.12 * (e[2] - e[1]); hi <- e[2] - 0.12 * (e[2] - e[1])
  set.seed(5)
  pts <- terra::vect(cbind(stats::runif(40, lo, hi), stats::runif(40, lo, hi)),
                     type = "points")

  ki_bin <- kernel_prep(pts, rs, max_D = 150, kernel = "gaussian",
                        verbose = FALSE, bin = TRUE, nbins = 256)
  ki_lean <- kernel_prep(pts, rs, max_D = 150, kernel = "gaussian",
                         verbose = FALSE, bin = TRUE, nbins = 256,
                         store_cell_data = FALSE)

  expect_null(ki_lean$d_list)
  expect_null(ki_lean$raw_cov)
  expect_false(is.null(ki_lean$binned))
  expect_false(ki_lean$cell_data_stored)
  expect_lt(as.numeric(object.size(ki_lean)),
            as.numeric(object.size(ki_bin)))

  y <- stats::rpois(40, exp(0.4 + 0.6 * as.numeric(scale(ki_bin$kernel_dat$cov))))
  fit <- function(ki) stats::glm(y ~ cov, family = stats::poisson(),
                                 data = data.frame(y = y, ki$kernel_dat))

  ob <- suppressWarnings(suppressMessages(
    multiScale_optim(fit(ki_bin), ki_bin, verbose = FALSE)))
  ol <- suppressWarnings(suppressMessages(
    multiScale_optim(fit(ki_lean), ki_lean, verbose = FALSE)))

  expect_equal(ol$scale_est[1, 1], ob$scale_est[1, 1], tolerance = 1e-6)
})

test_that("store_cell_data = FALSE is ignored when landscape metrics are present", {
  rs <- sim_rast(dim = 60, resolution = 10, user_seed = 7)
  land <- terra::classify(terra::subset(rs, "cont1"),
                          rcl = matrix(c(-Inf, 0.5, 1, 0.5, Inf, 2), ncol = 3, byrow = TRUE))
  names(land) <- "lc"
  terra::crs(land) <- "EPSG:32617"
  e <- as.vector(terra::ext(land))
  lo <- e[1] + 0.15 * (e[2] - e[1]); hi <- e[2] - 0.15 * (e[2] - e[1])
  set.seed(9)
  pts <- terra::vect(cbind(stats::runif(25, lo, hi), stats::runif(25, lo, hi)),
                     type = "points")

  vars <- msr_vars(lc_ed = landscape_var("lc", metric = "ed", radius = 80))
  ki <- kernel_prep(pts, land, max_D = 200, kernel = "gaussian",
                    scale_vars = vars, verbose = FALSE,
                    store_cell_data = FALSE)

  # Landscape metrics require cell-level data, so it must be retained.
  expect_false(is.null(ki$d_list))
  expect_true(ki$cell_data_stored)
})

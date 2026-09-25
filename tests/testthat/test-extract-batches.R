test_that("kernel_prep() extracts points in batches without mixing them up", {
  set.seed(11)
  r <- terra::rast(matrix(stats::rnorm(60 * 60), 60, 60),
                   extent = terra::ext(0, 600, 0, 600), crs = "EPSG:32617")
  names(r) <- "cov"
  # 60 points span three extraction batches of 25
  xy <- cbind(stats::runif(60, 150, 450), stats::runif(60, 150, 450))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                      coords = c("x", "y"), crs = "EPSG:32617")
  ki <- kernel_prep(pts, r, max_D = 100, kernel = "gaussian", verbose = FALSE)

  expect_length(ki$d_list, 60)
  expect_length(ki$raw_cov, 60)
  for (i in c(1, 25, 26, 51, 60)) {
    one <- exactextractr::exact_extract(r, sf::st_buffer(pts[i, ], 100),
                                        include_xy = TRUE, progress = FALSE)[[1]]
    d <- sqrt((one$x - xy[i, 1])^2 + (one$y - xy[i, 2])^2) / 100
    expect_equal(unname(ki$d_list[[i]]), d)
    expect_equal(unname(as.numeric(ki$raw_cov[[i]][, "cov"])), one$value)
  }
})

test_that("batch boundaries preserve multiple layers and lean binning", {
  set.seed(12)
  a <- terra::rast(matrix(stats::rnorm(40 * 40), 40, 40),
                   extent = terra::ext(0, 400, 0, 400), crs = "EPSG:32617")
  b <- a
  terra::values(b) <- stats::rnorm(terra::ncell(b))
  rasters <- c(a, b)
  names(rasters) <- c("a", "b")
  xy <- cbind(stats::runif(27, 100, 300), stats::runif(27, 100, 300))
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]),
                      coords = c("x", "y"), crs = "EPSG:32617")

  full <- kernel_prep(pts, rasters, max_D = 60, bin = TRUE,
                      store_cell_data = TRUE, verbose = FALSE)
  lean <- kernel_prep(pts, rasters, max_D = 60, bin = TRUE,
                      store_cell_data = FALSE, verbose = FALSE)
  expect_equal(lean$kernel_dat, full$kernel_dat)
  expect_equal(lean$binned, full$binned)
  expect_null(lean$d_list)
  expect_null(lean$raw_cov)

  for (i in c(25, 26, 27)) {
    extracted <- exactextractr::exact_extract(
      rasters, sf::st_buffer(pts[i, ], 60), include_xy = TRUE,
      progress = FALSE
    )[[1]]
    expect_equal(as.numeric(full$raw_cov[[i]][, "a"]), extracted$a)
    expect_equal(as.numeric(full$raw_cov[[i]][, "b"]), extracted$b)
    expected_d <- sqrt((extracted$x - xy[i, 1])^2 +
                       (extracted$y - xy[i, 2])^2) / 60
    expect_equal(unname(full$d_list[[i]]), expected_d)
  }
})

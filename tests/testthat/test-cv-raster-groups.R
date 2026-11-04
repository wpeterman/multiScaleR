test_that("native grouped-map preparation subsets buffers and row labels together", {
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 200, ymin = 0, ymax = 200, crs = "EPSG:32617")
  v <- seq_len(terra::ncell(r))
  terra::values(r) <- sin(v / 7) + cos(v / 13); names(r) <- "terrain"
  other <- r; terra::values(other) <- cos(v / 9) + sin(v / 17) + .1
  xy <- expand.grid(x = c(60, 100, 140), y = c(50, 80, 110, 140))
  points <- sf::st_as_sf(xy, coords = c("x", "y"), crs = 32617)
  ids <- paste0("native", seq_len(nrow(points))); rownames(points) <- ids
  labels <- rep(c("A", "B"), length.out = length(ids))
  cases <- list(gaussian = FALSE, exp = FALSE, expow = FALSE, weighted_sq = TRUE)
  for (name in names(cases)) {
    sq <- cases[[name]]; kernel <- if (sq) "gaussian" else name
    vars <- if (sq) msr_vars(hab = surface_var("terrain", metric = "sq", weighted = TRUE)) else
      msr_vars(hab = kernel_var("terrain"))
    inputs <- kernel_prep_by_group(points, list(A = r, B = other), labels,
      max_D = 35, kernel = kernel, sigma = 12, shape = if (kernel == "expow") 1.8 else NULL,
      scale_vars = vars, nbins = 16L, verbose = FALSE)
    expect_equal(inputs$cell_size, terra::res(r))
    original <- serialize(inputs, NULL)
    analysis <- ids[c(8, 1:7)]; assessment <- ids[9:12]
    train <- .msr_cv_inputs(inputs, analysis, "hab")
    held <- .msr_cv_inputs(inputs, assessment, "hab",
      if (is.null(train$binned)) NULL else train$binned$cv_reference)
    expect_identical(inputs$raster_group, stats::setNames(labels, ids))
    expect_identical(train$raster_group, inputs$raster_group[analysis])
    expect_identical(held$raster_group, inputs$raster_group[assessment])
    expect_identical(train$raw_cov, inputs$raw_cov[analysis])
    expect_identical(held$d_list, inputs$d_list[assessment])
    expect_identical(train$scale_vars, inputs$scale_vars)
    expect_equal(.msr_cv_initial_values(train), .msr_cv_initial_values(inputs)[analysis, , drop = FALSE], tolerance = 0)
    expect_equal(.msr_cv_initial_values(held), .msr_cv_initial_values(inputs)[assessment, , drop = FALSE], tolerance = 0)
    if (!is.null(train$binned)) expect_identical(train$binned$cv_reference, held$binned$cv_reference)
    expect_identical(serialize(inputs, NULL), original)
  }
})


test_that("native rectangular cells use a fixed geometric CV domain without clipping", {
  r <- terra::rast(nrows = 10, ncols = 20, xmin = 0, xmax = 200, ymin = 0, ymax = 200, crs = "EPSG:32617")
  terra::values(r) <- sin(seq_len(terra::ncell(r)) / 7); names(r) <- "hab"
  points <- sf::st_as_sf(expand.grid(x = c(60, 100, 140), y = c(60, 100, 140)), coords = c("x", "y"), crs = 32617)
  rownames(points) <- paste0("rect", seq_len(nrow(points)))
  inputs <- kernel_prep(points, r, max_D = 35, sigma = 12, nbins = 16L, verbose = FALSE)
  ids <- rownames(inputs$kernel_dat)
  expect_equal(inputs$cell_size, c(10, 20))
  expect_true(any(unlist(inputs$d_list) > 1))
  train <- .msr_cv_inputs(inputs, ids[1:6], "hab")
  reference <- train$binned$cv_reference
  expect_equal(reference$domain, (35 + sqrt(10^2 + 20^2) / 2) / 35)
  expect_identical(reference$cell_size, c(10, 20))
  expect_identical(train$d_list, inputs$d_list[ids[1:6]])
  expect_equal(as.numeric(rowSums(train$binned$csum$hab)), as.numeric(lengths(train$d_list)))
  held <- .msr_cv_inputs(inputs, ids[7:9], "hab", reference)
  expect_identical(held$binned$cv_reference, reference)
  expect_equal(as.numeric(rowSums(held$binned$csum$hab)), as.numeric(lengths(held$d_list)))
  bad <- inputs; bad$d_list[[1]][1] <- reference$domain * (1 + 1e-4)
  expect_error(.msr_cv_inputs(bad, ids[1:6], "hab"), "extraction envelope")
  bad <- inputs; bad$cell_size <- c(10, 0)
  expect_error(.msr_cv_inputs(bad, ids[1:6], "hab"), "finite positive")
  bad <- inputs; bad$cell_size <- rev(inputs$cell_size)
  expect_error(.msr_cv_inputs(bad, ids[7:9], "hab", reference), "frozen training bin")
  legacy <- inputs; legacy$cell_size <- NULL
  expect_error(.msr_cv_inputs(legacy, ids[1:6], "hab"), "rebuild legacy")
})

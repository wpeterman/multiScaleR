test_that("grouped folds keep whole sampling units and preserve the RNG", {
  ids <- paste0("id", 1:12); g <- rep(letters[1:6], each = 2)
  before <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  f <- grouped_cv_folds(g, 3, ids)
  after <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  expect_identical(before, after)
  expect_s3_class(f, "multiScaleR_folds")
  expect_setequal(unlist(lapply(f, `[[`, "assessment")), ids)
  expect_identical(lengths(lapply(f, `[[`, "assessment")), c(fold1 = 4L, fold2 = 4L, fold3 = 4L))
  for (s in f) {
    expect_length(intersect(s$analysis, s$assessment), 0L)
    expect_length(intersect(g[match(s$analysis, ids)], g[match(s$assessment, ids)]), 0L)
  }
  reverse <- grouped_cv_folds(rev(g), 3, rev(ids))
  expect_equal(lapply(f, function(s) sort(s$assessment)), lapply(reverse, function(s) sort(s$assessment)), ignore_attr = TRUE)
  expect_output(print(f), "grouped")
  expect_error(grouped_cv_folds(c("a", NA), 2), "nonmissing")
  expect_error(grouped_cv_folds(c("a", ""), 2), "empty")
  expect_error(grouped_cv_folds(rep("a", 4), 2), "at least k")
  expect_error(grouped_cv_folds(g, 1), "at least two")
  expect_error(grouped_cv_folds(g, 2, rep("x", 12)), "unique")
  expect_error(grouped_cv_folds(stats::setNames(g, rev(ids)), 2, ids), "aligned")
})

test_that("projected spatial blocks and cross-block groups stay indivisible", {
  ids <- paste0("s", 1:12)
  xy <- cbind(x = rep(seq(0, 500, 100), each = 2), y = rep(c(0, 10), 6))
  f <- spatial_cv_folds(xy, 50, 3, ids, crs = 32617)
  info <- attr(f, "fold_info")
  expect_equal(info$block_size, 50)
  expect_identical(info$units, "metre")
  expect_true(all(lengths(split(info$assignment$fold, info$assignment$block)) == 2L))
  expect_true(all(vapply(split(info$assignment$fold, info$assignment$block), function(z) length(unique(z)) == 1L, logical(1))))
  reverse <- spatial_cv_folds(xy[12:1, ], 50, 3, rev(ids), crs = 32617)
  expect_equal(lapply(f, function(s) sort(s$assessment)), lapply(reverse, function(s) sort(s$assessment)), ignore_attr = TRUE)
  pts <- sf::st_as_sf(data.frame(x = xy[, 1], y = xy[, 2], row.names = ids), coords = c("x", "y"), crs = 32617)
  rownames(pts) <- ids # sf construction resets data-frame row names.
  expect_equal(spatial_cv_folds(pts, 50, 3), f)
  expect_equal(lapply(spatial_cv_folds(terra::vect(pts), 50, 3, ids), `[[`, "assessment"), lapply(f, `[[`, "assessment"))
  g <- as.character(seq_len(12)); g[c(1, 3)] <- "a"; g[c(4, 5)] <- "b"
  joined <- spatial_cv_folds(xy, 50, 3, ids, group = g, crs = 32617)
  a <- attr(joined, "fold_info")$assignment
  expect_length(unique(a$unit[1:6]), 1L)
  expect_length(unique(a$fold[1:6]), 1L)
  expect_error(spatial_cv_folds(xy, 50, 5, ids, group = g, crs = 32617), "at least k")
  expect_error(spatial_cv_folds(xy, 50, 3), "explicit projected")
  expect_error(spatial_cv_folds(xy, 50, 3, crs = 4326), "projected CRS")
  expect_error(spatial_cv_folds(pts, 50, 3, crs = 32618), "must match")
  expect_error(spatial_cv_folds(pts, 50, 3, rev(ids)), "row names")
  expect_error(spatial_cv_folds(xy, 0, 3, crs = 32617), "positive finite")
  expect_error(spatial_cv_folds(xy, 50, 3, buffer = -1, crs = 32617), "nonnegative")
  edge <- spatial_cv_folds(cbind(c(-10, -.1, 0, 9.9, 10, 19.9), 0), 10, 3, crs = 32617)
  expect_identical(attr(edge, "fold_info")$assignment$block, rep(1:3, each = 2))
})

test_that("distance buffers exclude whole units and agree with a tiny distance oracle", {
  xy <- cbind(x = seq(0, 230, 10), y = 0); ids <- paste0("p", seq_len(nrow(xy)))
  f <- spatial_cv_folds(xy, 20, 4, ids, buffer = 10, crs = 32617)
  info <- attr(f, "fold_info"); a <- info$assignment
  for (nm in names(f)) {
    s <- f[[nm]]; train <- match(s$analysis, ids); assess <- match(s$assessment, ids)
    # Dense distances are limited to this 24-point software oracle.
    d <- abs(outer(xy[, 1], xy[assess, 1], "-"))
    near <- a$unit[apply(d, 1, min) <= 10]
    expected <- ids[!ids %in% s$assessment & a$unit %in% near]
    expect_setequal(info$excluded[[nm]], expected)
    expect_true(all(d[train, , drop = FALSE] > 10))
    expect_length(intersect(s$analysis, s$assessment), 0L)
    for (unit in unique(a$unit)) expect_true(all(ids[a$unit == unit] %in% s$analysis) ||
      all(!ids[a$unit == unit] %in% s$analysis))
  }
  expect_error(spatial_cv_folds(xy, 20, 3, ids, buffer = 1000, crs = 32617), "every training")
  expect_true(all(table(unlist(lapply(f, `[[`, "assessment"))) == 1L))
})

# Use the established spatial fixture at fixed, declared scales. No optimizer
# or uncertainty restriction is changed to test the prediction-data adapter.
unmarked_spatial_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- make_unmarked_fixture()
    x <- f$kernel_inputs
    ctx <- multiScaleR:::build_opt_context(f$mod, x$raw_cov,
      scale_vars = x$scale_vars, unit_conv = x$unit_conv,
      resolution = x$resolution, n_cols = x$n_cols, binned = x$binned)
    result <- multiScaleR:::kernel_scale_fn(c(40, 60) / x$unit_conv,
      x$d_list, x$raw_cov, x$kernel, f$mod, opt_context = ctx, mod_return = TRUE)
    object <- structure(list(opt_mod = result$mod, scl_params = result$scl_params,
      kernel_inputs = x, min_D = x$min_D, max_D = x$max_D,
      scale_est = data.frame(Mean = c(40, 60), row.names = c("bin1", "cont1"))),
      class = "multiScaleR")
    cache <<- list(f = f, object = object)
    cache
  }
})

test_that("keyed spatial prediction preserves native rows and training normalization", {
  f <- unmarked_spatial_fixture()
  object <- f$object
  x <- f$f$kernel_inputs
  n <- nrow(x$kernel_dat)
  order <- c(n, 2, 2, 1)
  nd <- object$opt_mod@data@siteCovs[order, , drop = FALSE]
  nd$occasion <- factor(c("b", "a", "b", "a"), levels = c("a", "b"))
  nd$effort <- c(1, 2, 4, 3)
  nd$bin1 <- nd$cont1 <- -999
  keys <- data.frame(site = seq_len(n))
  result <- prepare_unmarked_data(object, nd, x, keys)
  expected <- object$opt_mod@data@siteCovs[order, c("bin1", "cont1")]
  expect_equal(result[c("bin1", "cont1")], expected)
  expect_identical(result[setdiff(names(nd), c("bin1", "cont1"))],
                   nd[setdiff(names(nd), c("bin1", "cont1"))])
  expect_identical(row.names(result), row.names(nd))
  for (component in c("state", "det")) {
    native <- unmarked::predict(object$opt_mod, component, newdata = result,
                                appendData = TRUE, level = .9)
    snapshot <- prepare_prediction(object, nd, x, component = component,
      native_args = list(appendData = TRUE, level = .9), join_by = keys)
    expect_equal(predict(snapshot), native)
    expect_equal(predict(unserialize(serialize(snapshot, NULL))), native)
    expect_null(snapshot$kernel_inputs)
    expect_identical(snapshot$provenance$spatial$site_keys$site, nd$site)
    expect_identical(snapshot$provenance$spatial$row_names, rownames(nd))
  }
  converted <- x
  converted$unit_conv <- x$unit_conv * 1000
  converted$d_list <- lapply(x$d_list, function(d) d / 1000)
  converted$binned$bin_rep <- x$binned$bin_rep / 1000
  expect_equal(prepare_unmarked_data(object, nd, converted, keys), result)
  lean <- x
  lean$d_list <- lean$raw_cov <- NULL
  expect_equal(prepare_unmarked_data(object, nd, lean, keys), result)
  # Later source/model mutation cannot change cached numerical outputs.
  x$binned$vsum[[1]][] <- 99
  object$scl_params$mean[] <- 99
  expect_equal(predict(snapshot), native)
})

test_that("new spatial sites are evaluated once and joined to explicit scenarios", {
  f <- unmarked_spatial_fixture()
  object <- f$object
  selected <- c(2L, 4L, 7L)
  x <- kernel_prep(f$f$sim$pts[selected, ], f$f$rs_two,
    max_D = object$max_D, kernel = object$kernel_inputs$kernel, verbose = FALSE)
  keys <- data.frame(site = c("B", "D", "G"))
  nd <- data.frame(site = c("G", "B", "D", "B"), year = c(2, 1, 2, 2))
  out <- prepare_unmarked_data(object, nd, x, keys)
  raw <- multiScaleR:::.msr_binned_kernel_means(x$binned,
    c("bin1", "cont1"), c("bin1", "cont1"), c(40, 60) / x$unit_conv,
    NULL, x$kernel)
  manual <- sweep(sweep(raw, 2, object$scl_params$mean),
                  2, object$scl_params$sd, "/")
  expect_equal(as.matrix(out[c("bin1", "cont1")]),
    unname(manual[c(3, 1, 2, 1), ]), ignore_attr = TRUE)
  # A subset of sites must not be standardized as a new sample.
  expect_gt(max(abs(out$cont1 - as.numeric(scale(out$cont1)))), .01)
  expect_equal(unmarked::predict(object$opt_mod, "state", newdata = out),
    predict(prepare_prediction(object, nd, x, component = "state", join_by = keys)))
})

test_that("spatial input and key mismatches fail before native prediction", {
  f <- unmarked_spatial_fixture()
  object <- f$object
  x <- f$f$kernel_inputs
  nd <- object$opt_mod@data@siteCovs
  keys <- data.frame(site = nd$site)
  expect_error(prepare_unmarked_data(object, nd, x), "explicit join_by")
  bad <- keys
  bad$site[2] <- bad$site[1]
  expect_error(prepare_unmarked_data(object, nd, x, bad), "unique key")
  bad <- nd
  bad$site[1] <- -1
  expect_error(prepare_unmarked_data(object, bad, x, keys), "match exactly one")
  expect_error(prepare_unmarked_data(object, nd, x, keys[-1, , drop = FALSE]), "exactly one row")
  expect_error(prepare_unmarked_data(object, nd, x, nd["cont1"]), "cannot be spatial")
  bad <- x
  bad$scale_vars$source[1] <- "cont1"
  expect_error(prepare_unmarked_data(object, nd, bad, keys), "source and metric")
  bad <- x
  bad$kernel <- "exp"
  expect_error(prepare_unmarked_data(object, nd, bad, keys), "fitted kernel")
  bad <- x
  bad$max_D <- x$max_D * 2
  expect_error(prepare_unmarked_data(object, nd, bad, keys), "buffer extent")
  bad <- x
  bad$resolution <- x$resolution * 2
  expect_error(prepare_unmarked_data(object, nd, bad, keys), "raster resolution")
  bad <- x
  bad$binned <- NULL
  expect_error(prepare_unmarked_data(object, nd, bad, keys), "binned versus exact")
  bad <- object
  bad$scl_params$sd[1] <- 0
  expect_error(prepare_unmarked_data(bad, nd, x, keys), "standard deviations")
  expect_error(prepare_prediction(object, nd, component = "state", join_by = keys),
               "requires explicit kernel_inputs")
  bad <- nd
  bad$site[1] <- NA_integer_
  expect_error(prepare_unmarked_data(object, bad, x, keys), "missing values")
  # Missing support stays missing rather than changing native row membership.
  bad <- x
  bad$binned$csum$cont1[1, ] <- 0
  bad$binned$vsum$cont1[1, ] <- 0
  out <- prepare_unmarked_data(object, nd, bad, keys)
  expect_true(is.na(out$cont1[1]))
  expect_equal(predict(prepare_prediction(object, nd, bad, component = "state",
    native_args = list(na.rm = FALSE), join_by = keys)),
    unmarked::predict(object$opt_mod, "state", newdata = out, na.rm = FALSE))
})

test_that("weighted SQ spatial prediction retains its fitted metric definition", {
  f <- unmarked_spatial_fixture()
  specs <- msr_vars(bin1 = kernel_var("bin1"),
                    cont1 = surface_var("cont1", "sq", weighted = TRUE))
  x <- kernel_prep(f$f$sim$pts, f$f$rs_two, max_D = f$object$max_D,
    kernel = "gaussian", scale_vars = specs, verbose = FALSE)
  object <- f$object
  object$kernel_inputs <- x
  # Fixed known scales test the same metric evaluator used by the refit backend.
  context <- list(covs = c("bin1", "cont1"), n_covs = 2L,
    scale_vars = x$scale_vars, unit_conv = x$unit_conv,
    resolution = x$resolution, n_cols = x$n_cols, binned = x$binned)
  raw <- multiScaleR:::.msr_kernel_cov_w(c(40, 60) / x$unit_conv,
    x$d_list, x$raw_cov, x$kernel, context)
  object$scl_params <- list(mean = colMeans(raw), sd = apply(raw, 2, sd))
  nd <- data.frame(site = seq_len(nrow(raw)))
  expected <- scale(raw)
  out <- prepare_unmarked_data(object, nd, x, nd)
  expect_equal(as.matrix(out[c("bin1", "cont1")]), unname(expected), ignore_attr = TRUE)
  bad <- x
  bad$scale_vars$weighted[2] <- FALSE
  expect_error(prepare_unmarked_data(object, nd, bad, nd), "source and metric")
  bad <- x
  bad$d_list <- bad$raw_cov <- NULL
  expect_error(prepare_unmarked_data(object, nd, bad, nd), "Cell-level data")
})

test_that("exponential and exponential-power data use fitted scales and shapes", {
  f <- unmarked_spatial_fixture()
  for (kernel in c("exp", "expow")) {
    x <- f$f$kernel_inputs
    x$kernel <- kernel
    par <- c(40, 60) / x$unit_conv
    if (kernel == "expow") par <- c(par, 1.5, 3)
    ctx <- multiScaleR:::build_opt_context(f$f$mod, x$raw_cov,
      scale_vars = x$scale_vars, unit_conv = x$unit_conv,
      resolution = x$resolution, n_cols = x$n_cols, binned = x$binned)
    refit <- multiScaleR:::kernel_scale_fn(par, x$d_list, x$raw_cov, kernel,
      f$f$mod, opt_context = ctx, mod_return = TRUE)
    object <- f$object
    object$opt_mod <- refit$mod
    object$scl_params <- refit$scl_params
    object$kernel_inputs <- x
    if (kernel == "expow") object$shape_est <- data.frame(
      Mean = c(1.5, 3), row.names = c("bin1", "cont1"))
    nd <- refit$mod@data@siteCovs
    keys <- nd["site"]
    result <- prepare_unmarked_data(object, nd, x, keys)
    expect_equal(result[c("bin1", "cont1")], nd[c("bin1", "cont1")])
    if (kernel == "expow") {
      object$shape_est$Mean[1] <- NA_real_
      expect_error(prepare_unmarked_data(object, nd, x, keys), "kernel shapes")
    }
  }
})

test_that("native unmarked raster projection preserves cells, masks and output columns", {
  f <- unmarked_spatial_fixture()
  object <- f$object
  # Already transformed native covariates isolate the terra/S4 prediction bridge.
  # The guide uses kernel_scale.raster() to prepare these layers at fitted scales.
  grid <- terra::rast(nrows = 3, ncols = 5, nlyrs = 2,
    xmin = 0, xmax = 50, ymin = 0, ymax = 30, crs = "EPSG:32617")
  names(grid) <- c("bin1", "cont1")
  dat <- object$opt_mod@data@siteCovs[c("bin1", "cont1")]
  dat$cont1[c(2, 8)] <- NA_real_
  terra::values(grid) <- as.matrix(dat)
  # Explicit S4 dispatch works without attaching unmarked to the search path.
  native_unmarked <- function(model, data, type, ...) {
    unmarked::predict(model, type = type, newdata = data, ...)
  }
  path <- tempfile(fileext = ".tif")
  on.exit(unlink(path), add = TRUE)
  projected <- terra::predict(grid, object$opt_mod, fun = native_unmarked,
    type = "state", na.rm = TRUE, filename = path, overwrite = TRUE,
    wopt = list(datatype = "FLT8S"))
  good <- stats::complete.cases(dat)
  expected <- matrix(NA_real_, nrow(dat), 4,
    dimnames = list(NULL, c("Predicted", "SE", "lower", "upper")))
  expected[good, ] <- as.matrix(unmarked::predict(object$opt_mod, "state",
    newdata = dat[good, , drop = FALSE]))
  expect_equal(terra::values(projected), expected, tolerance = 1e-12)
  expect_identical(names(projected), colnames(expected))
  expect_true(terra::compareGeom(grid, projected, stopOnError = FALSE))
  expect_equal(terra::values(terra::rast(path)), expected, tolerance = 1e-12)
})

test_that("stored fit parameters avoid reporting round trips without accepting stale scales", {
  f <- unmarked_spatial_fixture()
  x <- f$f$kernel_inputs
  candidates <- seq(.03, .4, length.out = 101)
  internal <- head(candidates[candidates != (candidates * x$unit_conv) / x$unit_conv], 2)
  expect_length(internal, 2L)
  ctx <- multiScaleR:::build_opt_context(f$f$mod, x$raw_cov,
    scale_vars = x$scale_vars, unit_conv = x$unit_conv,
    resolution = x$resolution, n_cols = x$n_cols, binned = x$binned)
  refit <- multiScaleR:::kernel_scale_fn(internal, x$d_list, x$raw_cov, x$kernel,
    f$f$mod, opt_context = ctx, mod_return = TRUE)
  object <- f$object
  object$opt_mod <- refit$mod
  object$scl_params <- refit$scl_params
  object$optim_results <- list(par = internal)
  object$scale_est$Mean <- internal * x$unit_conv
  nd <- refit$mod@data@siteCovs
  out <- prepare_unmarked_data(object, nd, x, nd["site"])
  expect_identical(out[c("bin1", "cont1")], nd[c("bin1", "cont1")])
  object$scale_est$Mean[1] <- object$scale_est$Mean[1] + 1
  expect_error(prepare_unmarked_data(object, nd, x, nd["site"]),
               "Stored fitted parameters disagree")
})

unmarked_fixed_prediction_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- make_unmarked_fixture()
    specs <- msr_vars(bin1 = kernel_var("bin1"),
      rough = surface_var("cont1", "sq", radius = 90),
      edge = landscape_var("bin1", "ed", radius = 90))
    x <- kernel_prep(f$sim$pts, f$rs_two, max_D = 250,
      scale_vars = specs, kernel = "gaussian", verbose = FALSE)
    site_covs <- data.frame(site = seq_len(nrow(x$kernel_dat)), x$kernel_dat)
    umf <- unmarked::unmarkedFramePCount(y = f$sim$y, siteCovs = site_covs)
    mod <- unmarked::pcount(~1 ~bin1 + rough + edge, data = umf, K = 50)
    context <- multiScaleR:::build_opt_context(mod, x$raw_cov,
      scale_vars = x$scale_vars, unit_conv = x$unit_conv,
      resolution = x$resolution, n_cols = x$n_cols, binned = x$binned)
    result <- multiScaleR:::kernel_scale_fn(40 / x$unit_conv,
      x$d_list, x$raw_cov, x$kernel, mod, opt_context = context,
      mod_return = TRUE)
    object <- structure(list(opt_mod = result$mod,
      kernel_inputs = x, scl_params = result$scl_params,
      scale_est = data.frame(Mean = 40, row.names = "bin1"),
      min_D = x$min_D, max_D = x$max_D), class = "multiScaleR")
    cache <<- list(object = object, inputs = x, context = context)
    cache
  }
})

test_that("fixed-radius preparation is opt-in and matches native inputs", {
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  n <- nrow(f$inputs$kernel_dat)
  keys <- data.frame(site = seq_len(n))
  order <- c(n, 2, 2, 1)
  nd <- object$opt_mod@data@siteCovs[order, , drop = FALSE]
  expected <- nd
  nd$rough <- nd$edge <- -999
  manual <- prepare_unmarked_data(object, nd, f$inputs, keys)
  expect_identical(manual$rough, nd$rough)
  expect_identical(manual$edge, nd$edge)
  result <- prepare_unmarked_data(object, nd, f$inputs, keys, include_fixed = TRUE)
  expect_equal(result, expected, tolerance = 0)
  expect_identical(rownames(result), rownames(nd))
  for (type in c("state", "det")) {
    snapshot <- prepare_prediction(object, nd, f$inputs,
      component = type, join_by = keys, include_fixed = TRUE)
    native <- unmarked::predict(object$opt_mod, type, newdata = expected)
    expect_identical(predict(snapshot), native)
    expect_identical(predict(unserialize(serialize(snapshot, NULL))), native)
    expect_true(snapshot$provenance$spatial$include_fixed)
    expect_identical(snapshot$provenance$spatial$fixed_covariates, c("rough", "edge"))
  }
  converted <- f$inputs
  converted$unit_conv <- converted$unit_conv * 1000
  converted$d_list <- lapply(converted$d_list, function(d) d / 1000)
  converted$binned$bin_rep <- converted$binned$bin_rep / 1000
  expect_equal(prepare_unmarked_data(object, nd, converted, keys, include_fixed = TRUE),
               result)
})

test_that("fixed-radius preparation refuses unverified normalization and definitions", {
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  x <- f$inputs
  nd <- object$opt_mod@data@siteCovs
  keys <- nd["site"]
  bad <- x
  bad$scale_vars$radius[2] <- 100
  expect_error(prepare_unmarked_data(object, nd, bad, keys, include_fixed = TRUE),
               "Fixed-radius source and metric")
  bad <- object
  bad$opt_mod@data@siteCovs$rough <- bad$opt_mod@data@siteCovs$rough * 2
  expect_error(prepare_unmarked_data(bad, nd, x, keys, include_fixed = TRUE),
               "does not match")
  bad <- object
  bad$kernel_inputs$scl_params$mean["rough"] <-
    bad$kernel_inputs$scl_params$mean["rough"] + 1
  expect_error(prepare_unmarked_data(bad, nd, x, keys, include_fixed = TRUE),
               "does not match")
  bad <- object
  bad$kernel_inputs$scl_params$sd["edge"] <- 0
  expect_error(prepare_unmarked_data(bad, nd, x, keys, include_fixed = TRUE),
               "positive standard deviations")
  bad <- object
  bad$kernel_inputs$raw_cov <- bad$kernel_inputs$d_list <- NULL
  expect_error(prepare_unmarked_data(bad, nd, x, keys, include_fixed = TRUE),
               "Cell-level data")
  expect_error(prepare_prediction(object, nd, component = "state", include_fixed = TRUE),
               "requires explicit kernel_inputs")
  expect_error(prepare_unmarked_data(object, nd, x, keys, include_fixed = NA),
               "include_fixed")
  expect_error(prepare_prediction(object, include_fixed = NA), "include_fixed")
  expect_error(prepare_prediction(structure(list(), class = "multiScaleR"),
    include_fixed = TRUE), "only to unmarked")
})

test_that("fixed native training rows use their original explicit site keys", {
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  nd <- object$opt_mod@data@siteCovs
  n <- nrow(nd)
  keys <- data.frame(site = seq_len(n))
  order <- c(seq.int(n, 1L), seq_len(n))
  repeated_frame <- unmarked::unmarkedFramePCount(
    y = object$opt_mod@data@y[order, , drop = FALSE],
    siteCovs = nd[order, , drop = FALSE])
  object$opt_mod <- unmarked::pcount(~1 ~bin1 + rough + edge,
    data = repeated_frame, K = 50)
  object$join_by <- keys
  scenarios <- data.frame(site = c(3, 1, n, 3))
  result <- prepare_unmarked_data(object, scenarios, f$inputs, keys,
                                   include_fixed = TRUE)
  expect_equal(result$rough, nd$rough[c(3, 1, n, 3)], tolerance = 0)
  expect_equal(result$edge, nd$edge[c(3, 1, n, 3)], tolerance = 0)
  object$join_by <- NULL
  expect_error(prepare_unmarked_data(object, scenarios, f$inputs, keys,
    include_fixed = TRUE), "Supply `join_by`")
})

test_that("native fixed-metric raster projection retains predictors and verified scaling", {
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  rs <- make_unmarked_fixture()$rs_two
  before <- object$scl_params
  native_data <- object$opt_mod@data@siteCovs
  expect_identical(multiScaleR:::.msr_scale_vars_in_model(f$inputs$scale_vars,
    object)$covariate, c("bin1", "rough", "edge"))
  expect_error(kernel_scale.raster(rs, multiScaleR = object,
    scale_center = TRUE, verbose = FALSE),
    "Rasters scaled in model absent|Scaling parameters are missing")
  raw <- kernel_scale.raster(rs, multiScaleR = object, verbose = FALSE)
  expect_identical(names(raw), c("bin1", "rough", "edge"))
  scaled <- kernel_scale.raster(rs, multiScaleR = object,
    scale_center = TRUE, include_fixed = TRUE, verbose = FALSE)
  center <- c(object$scl_params$mean, f$inputs$scl_params$mean[c("rough", "edge")])
  spread <- c(object$scl_params$sd, f$inputs$scl_params$sd[c("rough", "edge")])
  expected <- sweep(sweep(terra::values(raw), 2L, center), 2L, spread, "/")
  expect_equal(terra::values(scaled), expected, tolerance = 1e-12)
  expect_identical(object$scl_params, before)
  expect_identical(object$opt_mod@data@siteCovs, native_data)
  native_unmarked <- function(model, data, ...) {
    unmarked::predict(model, newdata = data, ...)
  }
  mapped <- terra::predict(scaled, object$opt_mod, fun = native_unmarked,
                           type = "state", na.rm = TRUE)
  dat <- as.data.frame(terra::values(scaled))
  good <- stats::complete.cases(dat)
  expected_native <- unmarked::predict(object$opt_mod, "state",
    newdata = dat[good, , drop = FALSE])
  expect_equal(terra::values(mapped)[good, , drop = FALSE],
    as.matrix(expected_native), tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(names(mapped), names(expected_native))
  expect_error(kernel_scale.raster(rs, multiScaleR = object, include_fixed = TRUE,
    verbose = FALSE), "scale_center = TRUE")
  bad <- object$kernel_inputs$scale_vars
  bad$radius[2] <- 100
  expect_error(kernel_scale.raster(rs, multiScaleR = object, scale_vars = bad,
    scale_center = TRUE, include_fixed = TRUE, verbose = FALSE),
    "Fixed-radius source and metric")
  bad <- object
  bad$opt_mod@data@siteCovs$edge <- bad$opt_mod@data@siteCovs$edge + 1
  expect_error(kernel_scale.raster(rs, multiScaleR = bad,
    scale_center = TRUE, include_fixed = TRUE, verbose = FALSE), "does not match")
})

test_that("native raster clamping ignores unrelated character site identifiers", {
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  object$opt_mod@data@siteCovs$site <- paste0("site", object$opt_mod@data@siteCovs$site)
  object$opt_mod@data@siteCovs$group <- factor(rep(c("A", "B"), length.out =
    nrow(object$opt_mod@data@siteCovs)))
  rs <- make_unmarked_fixture()$rs_two
  mapped <- kernel_scale.raster(rs, multiScaleR = object,
    scale_center = TRUE, clamp = TRUE, include_fixed = TRUE, verbose = FALSE)
  values <- terra::values(mapped)
  for (nm in names(mapped)) {
    bounds <- range(object$opt_mod@data@siteCovs[[nm]])
    expect_gte(min(values[, nm], na.rm = TRUE), bounds[1])
    expect_lte(max(values[, nm], na.rm = TRUE), bounds[2])
  }
})

test_that("verified fixed raster scaling supports retained native GLM frames", {
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  dat <- object$opt_mod@data@siteCovs
  dat$y <- make_unmarked_fixture()$sim$y[, 1L]
  object$opt_mod <- stats::glm(y ~bin1 + rough + edge,
    data = dat, family = poisson(), model = TRUE)
  expect_true(all(is.finite(stats::coef(object$opt_mod))))
  rs <- make_unmarked_fixture()$rs_two
  layers <- kernel_scale.raster(rs, multiScaleR = object,
    scale_center = TRUE, include_fixed = TRUE, verbose = FALSE)
  mapped <- terra::predict(layers, object$opt_mod, type = "response", na.rm = TRUE)
  nd <- as.data.frame(terra::values(layers))
  good <- stats::complete.cases(nd)
  expect_equal(terra::values(mapped)[good, 1L],
    as.numeric(stats::predict(object$opt_mod, newdata = nd[good, , drop = FALSE],
                              type = "response")), tolerance = 1e-12)
  # A retained native subset maps by original row IDs, not new positions.
  object$opt_mod <- stats::glm(y ~bin1 + rough + edge,
    data = dat[c(15, 3, 7, 2, 4, 5, 9, 10, 11, 12), ],
    family = poisson(), model = TRUE)
  expect_equal(multiScaleR:::.msr_fixed_training_contract(object)$center,
    f$inputs$scl_params$mean[c("rough", "edge")], tolerance = 0)
  bad <- object
  rownames(bad$opt_mod$model)[1] <- "unknown"
  expect_error(multiScaleR:::.msr_fixed_training_contract(bad), "row IDs")
  bad <- object
  bad$opt_mod$model <- NULL
  expect_error(multiScaleR:::.msr_fixed_training_contract(bad), "retained model frame")
  bad <- object
  bad$opt_mod$model$rough <- bad$opt_mod$model$rough * 2
  expect_error(multiScaleR:::.msr_fixed_training_contract(bad), "does not match")
})

test_that("conditional-logistic fixed rasters retain native relative-score meaning", {
  skip_if_not_installed("survival")
  f <- unmarked_fixed_prediction_fixture()
  object <- f$object
  dat <- object$opt_mod@data@siteCovs
  dat$case_ <- rep(c(TRUE, FALSE, FALSE), length.out = nrow(dat))
  dat$choice <- rep(seq_len(nrow(dat) / 3L), each = 3L)
  coxph <- survival::coxph
  form <- stats::as.formula("case_ ~bin1 + rough + strata(choice)",
                            env = asNamespace("survival"))
  object$opt_mod <- survival::clogit(form, data = dat, model = TRUE)
  expect_true(all(is.finite(stats::coef(object$opt_mod))))
  expect_gt(min(eigen(stats::vcov(object$opt_mod), symmetric = TRUE,
                      only.values = TRUE)$values), 0)
  rs <- make_unmarked_fixture()$rs_two
  expect_warning(layers <- kernel_scale.raster(rs, multiScaleR = object,
    scale_center = TRUE, include_fixed = TRUE, verbose = FALSE), "categorical")
  expect_identical(names(layers), c("bin1", "rough"))
  # This placeholder supplies native formula structure. reference='zero' avoids
  # arbitrary choice-set centering; the score is not a choice probability.
  mapped <- terra::predict(layers, object$opt_mod, type = "risk", reference = "zero",
    const = data.frame(choice = dat$choice[1L]), na.rm = TRUE)
  nd <- as.data.frame(terra::values(layers))
  nd$choice <- dat$choice[1L]
  good <- stats::complete.cases(nd)
  expected <- stats::predict(object$opt_mod, nd[good, , drop = FALSE],
                             type = "risk", reference = "zero")
  expect_equal(terra::values(mapped)[good, 1L], as.numeric(expected), tolerance = 1e-12)
})

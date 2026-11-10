# No new outcomes. Reuse the declared native mixed fixture and deterministic
# node coordinates. Cache one accepted native field, never a seed-selected fit.
gdr_native_field_fixture <- local({
  value <- NULL; failure <- NULL
  function() {
    if (!is.null(value)) return(value)
    if (!is.null(failure)) stop(paste("Shared native field setup failed:", failure), call. = FALSE)
    f <- gdr_tmb_refit_fixture("mixed")
    block <- .msr_gdr_field_map(f$model@TMB$env$data$Z_lambda, 160L)
    nodes <- cbind(430000 + 40000 * (0:7), 4400000)
    xy <- nodes[block, , drop = FALSE]
    rownames(xy) <- rownames(f$model@data@siteCovs)
    before <- gdr_guard_state(f$model)
    before_threads <- TMB::openmp(DLL = "unmarked_TMBExports")
    fit <- tryCatch(testthat::with_mocked_bindings(
      fit_gdr_response_field(f$model, xy, 32617, 50000, f$starts,
        gdr_test_audit_path("native-field-fit-"), f$settings),
      .msr_gdr_gaussian_context = function(...) stop("Native fitting must not create RTMB callbacks.")),
      error = function(e) { failure <<- conditionMessage(e); stop(e) })
    stopifnot(identical(gdr_guard_state(f$model), before))
    stopifnot(identical(TMB::openmp(DLL = "unmarked_TMBExports"), before_threads))
    value <<- list(model = fit, template = f, coords = xy, nodes = nodes, block = block)
    value
  }
})

# A small actual source raster has a linear feature inside each response node's
# neighborhood. Physical kernel focal points vary within nodes, as coarse
# response groups and ecological predictor locations need not be identical.
gdr_native_field_projection_fixture <- local({
  value <- NULL
  function() {
    if (!is.null(value)) return(value)
    f <- gdr_native_field_fixture(); frame <- f$model@data
    raster <- terra::rast(nrows = 21L, ncols = 320L,
      xmin = 410000, xmax = 730000, ymin = 4389500, ymax = 4410500,
      crs = "EPSG:32617")
    x <- terra::xyFromCell(raster, seq_len(terra::ncell(raster)))[, 1L]
    center <- 430000 + 40000 * pmin(7, pmax(0, round((x - 430000)/40000)))
    terra::values(raster) <- (x - center) * 1.5/9500; names(raster) <- "hab"
    xy <- f$coords
    xy[, 1L] <- xy[, 1L] + frame@siteCovs$hab * 9500/1.5
    points <- terra::vect(xy, crs = "EPSG:32617")
    inputs <- kernel_prep(points, raster, max_D = 6000, kernel = "gaussian",
      scale_vars = msr_vars(hab = kernel_var("hab")), bin = FALSE, verbose = FALSE)
    callback <- make_gdr_refitter(f$model, f$template$starts,
      gdr_test_audit_path("native-field-kernel-"), tmb_settings = f$template$settings)
    context <- build_opt_context(f$model, inputs$raw_cov, refit_fn = callback,
      scale_vars = inputs$scale_vars, unit_conv = inputs$unit_conv,
      resolution = inputs$resolution, n_cols = inputs$n_cols,
      binned = inputs$binned, d_list = inputs$d_list)
    context$covariate_transform <- "raw"
    par <- 2000/inputs$unit_conv
    actual <- kernel_scale_fn(par, inputs$d_list, inputs$raw_cov, "gaussian",
      f$model, opt_context = context, mod_return = TRUE)
    object <- structure(list(opt_mod = actual$mod, backend = "refit", kernel_inputs = inputs,
      scale_est = data.frame(Mean = 2000, SE = NA_real_, row.names = "hab"),
      optim_results = list(par = par), min_D = inputs$min_D, max_D = inputs$max_D,
      scl_params = .msr_merge_scl_params(actual$scl_params, inputs$scl_params)), class = "multiScaleR")
    value <<- list(field = f, object = object, inputs = inputs, raster = raster,
      points = points, result = actual, context = context)
    value
  }
})

# Bounded adapter fixture at declared physical scales. This exercises native
# model behavior, not scale-estimation accuracy or statistical calibration.
refit_spatial_fixture <- function(kernel = "gaussian", metric = "sq",
                                  weighted = TRUE, bin = FALSE) {
  r <- terra::rast(nrows = 24, ncols = 24, xmin = 0, xmax = 240,
                   ymin = 0, ymax = 240, crs = "EPSG:3857")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- sin(xy[, 1] / 17) + cos(xy[, 2] / 23) + xy[, 1] / 150
  names(r) <- "hab"
  elevation <- r
  terra::values(elevation) <- 30 + xy[, 1] / 8 + sin(xy[, 2] / 11) *
    (1 + xy[, 1] / 30) + cos(xy[, 1] / 13 + xy[, 2] / 19)
  names(elevation) <- "elev"
  points <- expand.grid(x = seq(45, 195, 30), y = seq(45, 195, 30))
  pts <- terra::vect(as.matrix(points), crs = terra::crs(r))
  specs <- msr_vars(mean = kernel_var("hab"),
    rough = surface_var("elev", metric, weighted = weighted),
    fixed = surface_var("elev", "sa", radius = 25), unused = kernel_var("hab"))
  inputs <- kernel_prep(pts, c(r, elevation), max_D = 40, kernel = kernel,
    scale_vars = specs, bin = bin, verbose = FALSE)
  covs <- c("mean", "rough")
  par <- c(19, 29) / inputs$unit_conv
  if (kernel == "expow") par <- c(par, 1.3, 1.7)
  context <- list(covs = covs, n_covs = 2L, scale_vars = inputs$scale_vars,
    unit_conv = inputs$unit_conv, resolution = inputs$resolution,
    n_cols = inputs$n_cols, binned = inputs$binned)
  raw <- .msr_kernel_cov_w(par, inputs$d_list, inputs$raw_cov, kernel, context)
  scaled <- scale(raw)
  dat <- inputs$kernel_dat
  dat[covs] <- as.data.frame(scaled)
  n <- nrow(dat)
  dat$group <- factor(rep(c("a", "b", "c"), length.out = n))
  contrasts(dat$group) <- stats::contr.sum(3)
  dat$site <- factor(rep(seq_len(6), each = 6))
  dat$z <- sin(seq_len(n) / 3)
  dat$exposure <- 1 + seq_len(n) / (2 * n)
  dat$y <- (seq_len(n) %% 5L) + as.integer(dat$rough > 0)
  model <- stats::glm(y ~ mean + rough + fixed + group + z + offset(log(exposure)),
                      data = dat, family = poisson())
  # Native glm prediction restores the saved contrast matrix. Its model.frame
  # warns if newdata redundantly carries an explicit contrast attribute.
  attr(dat$group, "contrasts") <- NULL
  obj <- structure(list(opt_mod = model, kernel_inputs = inputs,
    scale_est = data.frame(Mean = par[1:2] * inputs$unit_conv,
                          SE = NA_real_, row.names = covs),
    shape_est = if (kernel == "expow") data.frame(Mean = par[3:4],
                         SE = NA_real_, row.names = covs) else NULL,
    optim_results = list(par = par), min_D = inputs$min_D, max_D = inputs$max_D,
    scl_params = .msr_merge_scl_params(list(mean = attr(scaled, "scaled:center"),
      sd = attr(scaled, "scaled:scale")), inputs$scl_params)), class = "multiScaleR")
  list(object = obj, data = dat, inputs = inputs, rasters = c(r, elevation),
       points = pts, raw = raw)
}

# One retained native refit at declared scales, shared by point/raster tests.
refit_unmarked_spatial_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    f <- make_unmarked_fixture()
    # The original synthetic fixture has no CRS. Declare projected toy units
    # and prepare its buffers afresh for checked raster/point interoperability.
    rasters <- f$rs_two
    terra::crs(rasters) <- "EPSG:3857"
    points <- terra::vect(f$sim$pts)
    terra::crs(points) <- terra::crs(rasters)
    x <- kernel_prep(points, rasters, max_D = 250, verbose = FALSE)
    context <- build_opt_context(f$mod, x$raw_cov, scale_vars = x$scale_vars,
      unit_conv = x$unit_conv, resolution = x$resolution, n_cols = x$n_cols,
      binned = x$binned, d_list = x$d_list)
    par <- rep(40 / x$unit_conv, context$n_covs)
    final <- kernel_scale_fn(par, x$d_list, x$raw_cov, x$kernel, f$mod,
                              opt_context = context, mod_return = TRUE)
    object <- structure(list(opt_mod = final$mod, kernel_inputs = x,
      scale_est = data.frame(Mean = par * x$unit_conv, SE = NA_real_, row.names = context$covs),
      optim_results = list(par = par), min_D = x$min_D, max_D = x$max_D,
      scl_params = .msr_merge_scl_params(final$scl_params, x$scl_params)), class = "multiScaleR")
    cache <<- list(object = object, data = object$opt_mod@data@siteCovs,
                   inputs = x, rasters = rasters, points = points)
    cache
  }
})

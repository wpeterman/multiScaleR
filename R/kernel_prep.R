#' @title Prepare Kernel Inputs for Scale Optimization
#'
#' @description Extracts the raster cells around each survey point, records
#' the distance from the point to each cell, and computes a first set of
#' kernel-weighted covariates. Fit your starting model to those covariates,
#' then pass the model and the returned object to [multiScale_optim()]. This
#' is the first step of every `multiScaleR` workflow, shown in
#' `vignette("v01-quickstart", package = "multiScaleR")`.
#'
#' @param pts Survey point locations as a `SpatVector` (terra) or `sf` object,
#'   in a projected coordinate system whose units match `max_D`. Row names
#'   are kept as point IDs.
#' @param raster_stack A `SpatRaster` (terra) with one or more layers, on the
#'   same projection as `pts`. Without `scale_vars`, each layer becomes a
#'   covariate with the layer's name, so the names must match the model
#'   formula. Layer names cannot be `"x"`, which is reserved for internal
#'   coordinate data. The raster must cover a buffer of `max_D` around every
#'   point.
#' @param max_D Positive numeric scalar. Radius, in map units, of the circle
#'   around each point within which raster cells are extracted. It is also
#'   the upper bound of the sigma search, and kernels are cut off at this
#'   distance. Set it above the largest plausible scale of effect, at least
#'   about 3.3 times the largest plausible sigma for a Gaussian kernel. If an
#'   estimate comes back near `max_D`, increase it and run `kernel_prep()`
#'   again. See the section on choosing `max_D` in [multiScaleR-kernels].
#' @param kernel Character. Distance-weighting function used for every
#'   covariate: `"gaussian"` (default), `"exp"` (negative exponential),
#'   `"expow"` (exponential power, which also estimates a shape), or
#'   `"fixed"` (equal weight within a radius). The weight functions are given
#'   in [multiScaleR-kernels]. To compare kernels, run `kernel_prep()` once
#'   per kernel.
#' @param scale_vars Optional covariate specifications created with
#'   [msr_vars()]. Omit it to make each raster layer one kernel-weighted mean.
#'   Supply it to define derived covariates, such as a kernel mean and a
#'   landscape or surface metric from the same layer.
#' @param sigma Optional numeric vector of starting sigma values, in map
#'   units, one per covariate whose scale is estimated (one per raster layer
#'   when `scale_vars` is omitted). Default: `NULL`, which starts every sigma
#'   at `max_D / 2`. The default is recommended. These values set the
#'   covariates in `kernel_dat` and are not the starting values of
#'   [multiScale_optim()], which chooses its own.
#' @param shape Optional numeric vector of starting shape values for
#'   `kernel = "expow"`, with the same length as `sigma`. Default: `NULL`,
#'   which starts every shape at 2 (the Gaussian form). Ignored for other
#'   kernels.
#' @param projected Logical. Whether `pts` and `raster_stack` use projected
#'   (planar) coordinates. Only projected coordinates are supported, so
#'   `FALSE` is an error. Default: `TRUE`.
#' @param progress Logical. Show progress bars during extraction and weighting.
#'   Default: `FALSE`.
#' @param verbose Logical. Print status messages. Default: `TRUE`.
#' @param bin Logical. If `TRUE` (default), summarize each kernel mean's cell
#'   values into distance bins. [multiScale_optim()] then evaluates kernel
#'   means from the bins instead of every cell, which is much faster and no
#'   longer depends on `max_D`, at the cost of a negligible approximation.
#'   The analytic backends and [multiScale_bayes()] require it for kernel
#'   means. With the default `refit` backend it also bins kernel-weighted
#'   surface metrics (`surface_var(weighted = TRUE)` with `"sq"`, `"ssk"`, or
#'   `"sku"`): each cell's centered value powers are split between the two
#'   nearest distance nodes, which keeps their approximation as small as that
#'   of a binned kernel mean. Weighted `"sa"`, `kernel = "fixed"`, landscape
#'   metrics, and hard-radius surface metrics are never binned, and the
#'   analytic backends and [multiScale_bayes()] evaluate weighted surface
#'   metrics from every cell. `FALSE` forces exact cell-by-cell evaluation.
#' @param nbins Integer, at least 2. Number of equal-width distance bins when
#'   `bin = TRUE`. More bins increase accuracy and memory. The default
#'   (`256`) reproduces exact optimization results to several significant
#'   figures. Ignored when `bin = FALSE`.
#' @param store_cell_data Logical. If `TRUE` (default), the per-cell distances
#'   (`d_list`) and values (`raw_cov`) are kept. If `FALSE`, they are dropped
#'   after binning, often shrinking the object by one to two orders of
#'   magnitude. This lean mode needs `bin = TRUE` and only kernel means. With
#'   landscape or surface metrics, or `bin = FALSE`, the cell data are kept
#'   and a message says so. The bins suffice for optimization,
#'   [profile_sigma()], and the final refit.
#'
#' @return An object of class `multiScaleR_data`, a list. Its full structure
#' is described in [multiScaleR_data-class]. The parts you use directly:
#' * `kernel_dat`: data frame with one row per point and one column per
#'   covariate. Values are centered and scaled (mean 0, standard deviation 1)
#'   at the starting sigma. Combine it with your response data to fit the
#'   starting model.
#' * `scl_params`: list with the `mean` and `sd` used for that scaling, in
#'   the covariates' own units.
#' * `kernel`, `max_D`, `min_D`, and `scale_vars`: the settings passed on to
#'   [multiScale_optim()] and [kernel_scale.raster()].
#'
#' Printing the object reports the number of points and covariates, the
#' kernel, the number of cells extracted around the first point ("Number of
#' elements"), whether binning is on, the minimum distance (`min_D`, about one
#' cell width and the lower bound of the sigma search), the maximum distance
#' (`max_D`), and the unit conversion (distances are stored divided by
#' `max_D`).
#'
#' @details
#' # Prepare the points and rasters
#'
#' Points and rasters must share a defined, projected CRS. Every distance,
#' including `max_D` and any `sigma` you supply, is in the linear unit of
#' that projection, usually meters. Landscape and surface metrics also need
#' square cells. Row names of `pts` are kept throughout, so model data can be
#' joined back to the points.
#'
#' Each point's buffer is the circle of radius `max_D`, and every cell it
#' intersects is extracted, including partly covered boundary cells, which
#' keep their full kernel weight. Missing cells drop out of each weighted
#' mean, and the remaining weights renormalize. Extraction runs in batches of
#' 25 points to limit the peak memory of temporary tables. Point order,
#' returned values, and distances are not affected.
#'
#' # Manage memory
#'
#' A larger `max_D` extracts more cells per point. Check the size of the
#' object and the memory needed for parallel optimization with
#' [estimate_multiscale_ram()]. When every covariate is a kernel mean,
#' `store_cell_data = FALSE` keeps only the bins.
#'
#' # Related preparation functions
#'
#' To assign each observation to a different annual raster while estimating
#' one pooled landscape relationship, use [kernel_prep_by_group()]. To
#' prepare inputs for `unmarked` models, see [prepare_unmarked_data()].
#'
#' @seealso [multiScale_optim()], [msr_vars()], [multiScaleR-kernels],
#'   [multiScaleR_data-class], [estimate_multiscale_ram()].
#'   `vignette("v02-user-guide", package = "multiScaleR")` explains each
#'   argument on a worked example.
#'
#' @examples
#' ## Survey points and a binary habitat raster bundled with the package
#' data("pts")
#' data("count_data")
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#'
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = hab,
#'                              max_D = 250,          # map units (meters)
#'                              kernel = "gaussian",
#'                              verbose = FALSE)
#' # 75 points, one covariate, cells within 250 m of each point
#' kernel_inputs
#'
#' # Standardized kernel-weighted habitat at the starting sigma (125 m)
#' head(kernel_inputs$kernel_dat)
#'
#' # Fit the starting model that multiScale_optim() will optimize
#' df <- data.frame(count_data["y"], kernel_inputs$kernel_dat)
#' mod <- glm(y ~ hab, family = poisson, data = df)
#'
#' ## Lean mode keeps only the distance bins
#' lean_inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                            store_cell_data = FALSE, verbose = FALSE)
#' # Object sizes in megabytes
#' c(full = as.numeric(object.size(kernel_inputs)),
#'   lean = as.numeric(object.size(lean_inputs))) / 1e6
#' @rdname kernel_prep
#' @export
#' @importFrom exactextractr exact_extract
#' @importFrom terra nlyr vect
#' @importFrom sf st_as_sf st_buffer st_coordinates st_crs st_crs<-
#' @importFrom dplyr bind_rows
#' @importFrom fields rdist rdist.earth
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix as.matrix
#' @importFrom methods as

kernel_prep <- function(pts,
                        raster_stack,
                        max_D,
                        kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                        scale_vars = NULL,
                        sigma = NULL,
                        shape = NULL,
                        projected = TRUE,
                        progress = FALSE,
                        verbose = TRUE,
                        bin = TRUE,
                        nbins = 256L,
                        store_cell_data = TRUE){
  unit_conv <- max_D

  kernel <- match.arg(kernel)
  validate_scalar_numeric(max_D, "max_D", positive = TRUE)
  validate_scalar_logical(projected, "projected")
  validate_scalar_logical(progress, "progress")
  validate_scalar_logical(verbose, "verbose")
  validate_scalar_logical(bin, "bin")
  validate_scalar_logical(store_cell_data, "store_cell_data")
  if (isTRUE(bin)) {
    validate_scalar_numeric(nbins, "nbins", integerish = TRUE, lower = 2)
  }

  if(!inherits(raster_stack, "SpatRaster")){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }
  validate_reserved_variable_names(names(raster_stack), "raster layer name")

  scale_vars <- .msr_validate_scale_vars(scale_vars = scale_vars,
                                         raster_stack = raster_stack,
                                         kernel = kernel)
  validate_reserved_variable_names(scale_vars$covariate, "derived covariate name")
  opt_scale_vars <- .msr_optimized_scale_vars(scale_vars)
  n_optimized <- nrow(opt_scale_vars)

  if (any(scale_vars$type %in% c("landscape", "surface"))) {
    raster_res <- terra::res(raster_stack)
    if (!isTRUE(all.equal(raster_res[[1]], raster_res[[2]]))) {
      stop("Landscape and surface metric variables require square raster cells.",
           call. = FALSE)
    }
  }

  if(is.null(sigma)){
    sigma <- rep(max_D/2, n_optimized) ## Need to set to minimum distance...no scale
    if(verbose && n_optimized > 0){
      cat(paste("\nNo sigma values provided...",
                "Creating necessary elements to optimize sigma\n", sep = '\n'))
    }
  }

  if(is.null(shape) & kernel == 'expow' & n_optimized > 0){
    shape <- rep(2, n_optimized) ## Need to set shape
    if(verbose){
      cat(paste("\nNo shape values provided...",
                "Creating necessary elements to optimize shape\n", sep = '\n'))
    }
  }


  if(length(sigma) != n_optimized){
    stop("Number of sigma values must equal the number of optimized covariates!!!")
  }
  if(n_optimized > 0){
    validate_numeric_vector(sigma,
                            "sigma",
                            length_ = n_optimized,
                            positive = TRUE)
  }

  if(kernel == 'expow' && n_optimized > 0){
    validate_numeric_vector(shape,
                            "shape",
                            length_ = n_optimized,
                            positive = TRUE)
  }

  if(isFALSE(projected)){
    stop("Currently, spatial points must be projected")
    r_pts <- st_as_sf(pts, coords = c(1,2))
    st_crs(r_pts) <- 4326

    buff_poly <- st_buffer(r_pts,
                           dist = max_D)
    spat_poly <- vect(buff_poly)

    if(verbose){
      cat(paste0("\nExtracting values from raster layers...\n"))
    }
    r_ext <- exact_extract(raster_stack,
                           buff_poly,
                           include_cell = T,
                           include_xy = T)


    # Convert to per-point value matrices (dense or sparse, whichever is smaller)
    cell_list <- lapply(r_ext, df_to_values)

    names(r_ext) <- 1:length(r_ext)

    r_ext_ <- bind_rows(r_ext, .id = "id")

    ## Progress bar
    if(isTRUE(progress)){
      cat(paste0("\nCalculating distances...\n"))
      pb = txtProgressBar(min = 0,
                          max = dim(r_pts)[1],
                          initial = 0,
                          char = "*",
                          style = 3)
    }

    D <- vector('list', dim(r_pts)[1])
    for (i in 1:dim(r_pts)[1]) {
      if(isTRUE(progress)){
        setTxtProgressBar(pb,i)
      }

      D[[i]] <- rdist.earth(st_coordinates(r_pts[i,]),
                            # r_ext[r_ext$ID == 1,c("x","y")],
                            r_ext[[i]][,c("x","y")],
                            miles = F)[1,] * 1000
    }

    if(isTRUE(progress)){
      close(pb)
    }
    min_D <- floor(rdist.earth(st_coordinates(r_ext[[1]][1:2,c("x","y")]),
                               miles = F)[1,2] * 1000)

  } else { ## Projected points
    if (!any(grep(paste(c("SpatVector",
                          "sf"), collapse = "|"), class(pts))))
      stop("`pts` must be terra or sf class object")

    # if (class(pts)[1] == "sf") {
    if(any(grep(paste(c("sfc", "sfc_MULTIPOINT", "sf"), collapse = "|"),
                class(pts)))){
      buff_poly <- st_buffer(pts,
                             dist = max_D)
    }

    if(class(pts)[1] == "SpatVector"){
      pts <- st_as_sf(pts)
      buff_poly <- st_buffer(pts,
                             dist = max_D)
    }

    if(nrow(pts) == 0){
      stop("`pts` must contain at least one point.")
    }

    if(verbose){
      cat(paste0("\nExtracting values from raster layers...\n"))
    }

    # Extract, convert, and compute distances in batches of points so the full
    # extraction table (values, x, y, coverage fraction, and cell for every cell
    # in every buffer) never sits in memory at once. Each point's buffer is
    # extracted independently, so results are identical to a single call.
    n_pts_all <- dim(pts)[1]
    batch_size <- 25L
    batches <- split(seq_len(n_pts_all),
                     ceiling(seq_len(n_pts_all) / batch_size))
    need_cells <- .msr_needs_cells(scale_vars)
    single_layer <- nlyr(raster_stack) == 1
    cell_list <- vector("list", n_pts_all)
    D <- vector("list", n_pts_all)
    # The lower bound of the sigma search is one cell width, as in
    # kernel_prep_by_group(). It used to be the distance between the first two
    # cells extracted for the first point. Those are neighbors in one row
    # unless the top row of that buffer holds a single cell, when the second
    # cell lies on the next row and min_D became several cells wide.
    min_D <- floor(terra::res(raster_stack)[1])

    if(isTRUE(progress)){
      pb = txtProgressBar(min = 0,
                          max = n_pts_all,
                          initial = 0,
                          char = "*",
                          style = 3)
    }

    for (b in batches) {
      r_ext <- exact_extract(raster_stack,
                             buff_poly[b, ],
                             include_xy = T,
                             include_cell = need_cells,
                             progress = FALSE)
      for (j in seq_along(b)) {
        i <- b[j]
        df <- r_ext[[j]]
        if (single_layer) {
          c_names <- colnames(df)
          c_names[1] <- names(raster_stack)
          colnames(df) <- c_names
        }
        xy <- df[, c("x", "y")]
        cell_list[[i]] <- df_to_values(df)
        D[[i]] <- rdist(st_coordinates(pts[i,]), xy)[1,] / unit_conv
        r_ext[j] <- list(NULL)
      }
      rm(r_ext, df, xy)
      if(isTRUE(progress)){
        setTxtProgressBar(pb, max(b))
      }
    }
    if(isTRUE(progress)){
      close(pb)
    }
  } ## End ifelse for projected points

  n_pts <- dim(pts)[1]
  point_ids <- row.names(pts)
  if (is.null(point_ids) || length(point_ids) != n_pts ||
      anyNA(point_ids) || any(!nzchar(point_ids)) || anyDuplicated(point_ids)) {
    point_ids <- as.character(seq_len(n_pts))
  }

  cov.w <- matrix(NA_real_,
                  nrow = n_pts,
                  ncol = nrow(scale_vars))
  colnames(cov.w) <- scale_vars$covariate
  rownames(cov.w) <- point_ids
  names(D) <- point_ids
  names(cell_list) <- point_ids
  sigma <- sigma / unit_conv

  if(isTRUE(progress)){
    cat(paste0("\nCalculating weights...\n"))


    pb = txtProgressBar(min = 0,
                        max = dim(pts)[1],
                        initial = 0,
                        char = "*",
                        style = 3)
  }
  # system.time(
  for(i in seq_len(n_pts)){
    if(isTRUE(progress)){
      setTxtProgressBar(pb,i)
    }

    cov.w[i,] <- .msr_eval_scale_vars(
      d = D[[i]],
      cov_df = cell_list[[i]],
      scale_vars = scale_vars,
      sigma = sigma,
      shape = shape,
      kernel = kernel,
      unit_conv = unit_conv,
      resolution = terra::res(raster_stack)[[1]],
      n_cols = terra::ncol(raster_stack)
    )

  }
  if(isTRUE(progress)){
    close(pb)
  }

  validate_covariates_before_scale(cov.w,
                                   context = "`kernel_prep()` covariates",
                                   scale_vars = scale_vars)
  scl_df <- scale(cov.w)
  kernel_dat <- as.data.frame(scl_df)

  # Precompute distance-binned summaries for kernel-type covariates so the
  # optimizer can evaluate kernel-weighted means in O(nbins) per point instead
  # of O(n_cells). See R/kernel_bins.R.
  binned <- NULL
  kernel_specs <- scale_vars[scale_vars$type == "kernel", , drop = FALSE]
  if (isTRUE(bin) && nrow(kernel_specs) > 0) {
    binned <- .msr_build_kernel_bins(
      d_list = D,
      value_list = cell_list,
      covariates = kernel_specs$covariate,
      sources = kernel_specs$source,
      nbins = nbins,
      point_ids = point_ids
    )
  }

  # Kernel-weighted surface moments (sq, ssk, sku) are linearly binned power
  # sums, so the optimizer evaluates them in O(nbins) as well. The cell data
  # stay stored for these covariates. See R/kernel_bins.R.
  surface_binned <- NULL
  if (isTRUE(bin)) {
    surface_binned <- .msr_build_surface_bins(
      d_list = D,
      value_list = cell_list,
      specs = .msr_surface_bin_specs(scale_vars, kernel),
      nbins = nbins,
      point_ids = point_ids
    )
  }

  # Cell-level data (`d_list`/`raw_cov`) can be dropped only when every
  # covariate is kernel-type (landscape metrics need the full per-cell spatial
  # detail) and binned summaries are available to drive the optimizer.
  has_landscape <- any(scale_vars$type == "landscape")
  drop_cells <- isFALSE(store_cell_data) && !is.null(binned) &&
    !has_landscape && all(scale_vars$type == "kernel")
  if (isFALSE(store_cell_data) && !drop_cells && verbose) {
    cat(paste0("\n`store_cell_data = FALSE` was ignored: cell-level data are ",
               "still required for this configuration.\n"))
  }

  out <- list(kernel_dat = kernel_dat,
              d_list = if (drop_cells) NULL else D,
              raw_cov = if (drop_cells) NULL else cell_list,
              kernel = kernel,
              shape = shape,
              min_D = min_D,
              max_D = max_D,
              n_covs = n_optimized,
              unit_conv = unit_conv,
              sigma = sigma,
              scale_vars = scale_vars,
              resolution = terra::res(raster_stack)[[1]],
              cell_size = terra::res(raster_stack),
              n_cols = terra::ncol(raster_stack),
              binned = binned,
              surface_binned = surface_binned,
              cell_data_stored = !drop_cells,
              scl_params = list(mean = attr(scl_df, "scaled:center"),
                                sd = attr(scl_df, "scaled:scale")))

  class(out) <- 'multiScaleR_data'
  return(out)
}

# Convert an extracted data frame to a compact per-point matrix of raster
# values. Dense storage is faster to weight (`scale_type()` evaluates it in
# vectorized R instead of via per-cell sparse random access) and avoids the
# repeated sparse-to-dense conversion in the landscape-metric C++ path. Sparse
# storage is smaller only when most cells are zero (e.g. low-prevalence binary
# rasters). We keep whichever representation is smaller per point, so memory is
# never worse than the previous sparse-only storage while continuous and
# categorical layers (the common cases) get the faster dense path. Raster cell
# ids are preserved as an attribute for the landscape edge/adjacency metrics;
# `scale_type()` and the landscape path dispatch on the stored type.
df_to_values <- function(df) {
  cells <- if ("cell" %in% names(df)) df$cell else NULL
  value_cols <- setdiff(names(df), c("x", "y", "cell", "coverage_fraction"))
  dense <- as.matrix(df[, value_cols, drop = FALSE])
  storage.mode(dense) <- "double"

  out <- dense
  sparse <- tryCatch(as(dense, "sparseMatrix"), error = function(e) NULL)
  if (!is.null(sparse) &&
      as.numeric(utils::object.size(sparse)) <
        as.numeric(utils::object.size(dense))) {
    out <- sparse
  }
  if (!is.null(cells)) {
    attr(out, "cell") <- cells
  }
  out
}

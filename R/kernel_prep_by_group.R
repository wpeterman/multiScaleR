#' Prepare Kernel Inputs from Group-Specific Rasters
#'
#' Assign each observation to one raster map, then prepare one pooled data
#' object for \code{\link{multiScale_optim}}. Use this when land cover changes
#' across years but the model estimates shared landscape effects and scales.
#'
#' @param pts Projected point locations as an \code{sf} or \code{SpatVector}
#'   object. Point order matches \code{group} and the returned covariate rows.
#' @param raster_stacks Named list of \code{SpatRaster} maps. Each map has the
#'   same extent, grid, CRS, layer names, and class coding. Its name identifies
#'   a year or other observation group.
#' @param group Character or factor vector, one value per point, identifying
#'   the corresponding name in \code{raster_stacks}. Every point uses only its
#'   assigned map, including cells affected by changes outside the point's
#'   immediate location.
#' @param max_D Maximum extraction radius in the projected map units.
#' @param kernel Kernel weighting function; see \code{\link{kernel_prep}}.
#' @param scale_vars Optional \code{\link{msr_vars}} specification. By default,
#'   each raster layer becomes a kernel-weighted mean covariate.
#' @param sigma Optional starting scales in map units. Defaults to
#'   \code{max_D / 2} for each optimized covariate.
#' @param shape Optional starting shapes for the exponential power kernel.
#'   Defaults to 2 for each optimized covariate.
#' @param bin Logical. Build pooled distance bins for kernel covariates.
#'   Default: \code{TRUE}. All observations use one set of bin definitions.
#' @param nbins Number of distance bins when \code{bin = TRUE}; default 256.
#' @param store_cell_data Logical. Retain per-cell values and distances. Default:
#'   \code{TRUE}. Set \code{FALSE} for a lean kernel-only binned object; metrics
#'   requiring cell detail retain it regardless.
#' @param progress Logical. Display extraction progress. Default \code{FALSE}.
#' @param verbose Logical. Print preparation messages. Default \code{TRUE}.
#'
#' @return A \code{multiScaleR_data} list containing:
#' \describe{
#'   \item{\code{kernel_dat}}{Initial covariate data frame, with one row per
#'     point. Columns are centered and divided by their pooled standard
#'     deviations. Use this in the starting model.}
#'   \item{\code{d_list}, \code{raw_cov}}{Named per-point distance vectors
#'     (divided by \code{unit_conv}) and raster-value matrices. Both are
#'     \code{NULL} in lean mode.}
#'   \item{\code{binned}}{Pooled distance-bin summaries for kernel covariates,
#'     or \code{NULL} when \code{bin = FALSE}. Each covariate has finite-value
#'     sums and counts by point and distance bin.}
#'   \item{\code{raster_group}}{Named character vector recording the map
#'     assigned to each point, in input order.}
#'   \item{\code{scl_params}}{Pooled means and standard deviations used to
#'     create \code{kernel_dat}; the fitted model replaces these with scaling
#'     parameters at the optimized scales.}
#'   \item{\code{kernel}, \code{sigma}, \code{shape}}{Kernel name, initial
#'     scales divided by \code{unit_conv}, and optional initial shape values.}
#'   \item{\code{min_D}, \code{max_D}, \code{unit_conv}}{Lower scale bound,
#'     maximum extraction radius, and internal distance scaling factor. These
#'     are in map units except the internally scaled \code{sigma}.}
#'   \item{\code{n_covs}, \code{scale_vars}}{Number of optimized covariates and
#'     the validated covariate definitions.}
#'   \item{\code{resolution}, \code{n_cols}}{Raster cell size and column count,
#'     used to calculate landscape configuration metrics.}
#'   \item{\code{cell_data_stored}}{Whether per-cell values and distances are
#'     retained.}
#' }
#' Pass the full object to \code{\link{multiScale_optim}}.
#'
#' @details Raster maps must be aligned because landscape configuration metrics
#'   use cell adjacency and the resulting model assumes common map units.
#'   Group-specific maps may contain different missing cells. For kernel
#'   covariates, pooled bins track both finite-cell value sums and counts.
#'   The fitted model has common slopes and scales unless its formula explicitly
#'   includes group interactions. Choose the response likelihood and any group
#'   effect according to the sampling design. A pooled model does not establish
#'   that effects are stable across years.
#'
#' @examples
#' library(terra)
#' r1 <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 200,
#'            ymin = 0, ymax = 200, crs = "EPSG:26915")
#' values(r1) <- seq_len(ncell(r1)) / ncell(r1)
#' names(r1) <- "habitat"
#' r2 <- r1
#' r2[seq_len(30)] <- 0
#' pts <- vect(cbind(c(50, 70, 90, 110, 130, 150),
#'                   c(60, 90, 120, 140, 100, 70)),
#'             type = "points", crs = crs(r1))
#' prep <- kernel_prep_by_group(pts, list(y1 = r1, y2 = r2),
#'                              group = rep(c("y1", "y2"), 3),
#'                              max_D = 25, verbose = FALSE)
#' head(prep$kernel_dat)
#'
#' @seealso \code{\link{kernel_prep}}, \code{\link{multiScale_optim}}
#' @export
kernel_prep_by_group <- function(pts, raster_stacks, group, max_D,
                                 kernel = c("gaussian", "exp", "expow", "fixed"),
                                 scale_vars = NULL, sigma = NULL, shape = NULL,
                                 bin = TRUE, nbins = 256L,
                                 store_cell_data = TRUE,
                                 progress = FALSE, verbose = TRUE) {
  kernel <- match.arg(kernel)
  validate_scalar_numeric(max_D, "max_D", positive = TRUE)
  validate_scalar_logical(bin, "bin")
  validate_scalar_logical(store_cell_data, "store_cell_data")
  validate_scalar_logical(progress, "progress")
  validate_scalar_logical(verbose, "verbose")
  if (isTRUE(bin)) {
    validate_scalar_numeric(nbins, "nbins", integerish = TRUE, lower = 2)
  }
  if (!inherits(pts, c("sf", "SpatVector"))) {
    stop("`pts` must be an sf or terra SpatVector of projected points.", call. = FALSE)
  }
  pts <- sf::st_as_sf(pts)
  if (nrow(pts) == 0 || !all(sf::st_geometry_type(pts) == "POINT") ||
      is.na(sf::st_crs(pts)) || isTRUE(sf::st_is_longlat(pts))) {
    stop("`pts` must contain projected POINT geometries with a defined CRS.",
         call. = FALSE)
  }
  if (!is.list(raster_stacks) || length(raster_stacks) == 0 ||
      is.null(names(raster_stacks)) ||
      anyNA(names(raster_stacks)) || any(!nzchar(names(raster_stacks))) ||
      anyDuplicated(names(raster_stacks)) ||
      !all(vapply(raster_stacks, inherits, logical(1), "SpatRaster"))) {
    stop("`raster_stacks` must be a named list of SpatRaster maps with unique names.",
         call. = FALSE)
  }
  if (length(group) != nrow(pts) || anyNA(group) ||
      !all(as.character(group) %in% names(raster_stacks))) {
    stop("`group` must assign every point to a name in `raster_stacks`.",
         call. = FALSE)
  }
  group <- as.character(group)
  ref <- raster_stacks[[1]]
  if (isTRUE(terra::is.lonlat(ref)) ||
      !isTRUE(sf::st_crs(pts) == sf::st_crs(terra::crs(ref)))) {
    stop("Points and rasters must share the same projected CRS.", call. = FALSE)
  }
  for (r in raster_stacks) {
    if (!isTRUE(terra::compareGeom(ref, r, stopOnError = FALSE)) ||
        !identical(names(ref), names(r))) {
      stop("All group rasters must have identical grid geometry and layer names.",
           call. = FALSE)
    }
  }
  validate_reserved_variable_names(names(ref), "raster layer name")
  scale_vars <- .msr_validate_scale_vars(scale_vars, ref, kernel)
  validate_reserved_variable_names(scale_vars$covariate,
                                   "derived covariate name")
  if (any(scale_vars$type %in% c("landscape", "surface")) &&
      !isTRUE(all.equal(terra::res(ref)[1], terra::res(ref)[2]))) {
    stop("Landscape and surface metrics require square raster cells.",
         call. = FALSE)
  }
  n_optimized <- nrow(.msr_optimized_scale_vars(scale_vars))
  if (is.null(sigma)) sigma <- rep(max_D / 2, n_optimized)
  if (n_optimized > 0) {
    validate_numeric_vector(sigma, "sigma", length_ = n_optimized,
                            positive = TRUE)
  } else if (length(sigma)) {
    stop("`sigma` must have one value per optimized covariate.", call. = FALSE)
  }
  if (kernel == "expow" && n_optimized > 0) {
    if (is.null(shape)) shape <- rep(2, n_optimized)
    validate_numeric_vector(shape, "shape", length_ = n_optimized,
                            positive = TRUE)
  }
  point_ids <- row.names(pts)
  if (is.null(point_ids) || anyNA(point_ids) || anyDuplicated(point_ids) ||
      any(!nzchar(point_ids))) point_ids <- as.character(seq_len(nrow(pts)))
  distances <- values <- vector("list", nrow(pts))
  names(distances) <- names(values) <- point_ids
  for (g in unique(group)) {
    idx <- which(group == g)
    if (verbose) message("Extracting raster group: ", g)
    extracted <- exactextractr::exact_extract(
      raster_stacks[[g]], sf::st_buffer(pts[idx, ], dist = max_D),
      include_xy = TRUE, include_cell = .msr_needs_cells(scale_vars),
      progress = progress
    )
    for (j in seq_along(idx)) {
      i <- idx[[j]]
      cells <- extracted[[j]]
      if (nrow(cells) < 2L) {
        stop(sprintf("Point %s has fewer than two raster cells in its buffer.",
                     point_ids[[i]]), call. = FALSE)
      }
      values[[i]] <- df_to_values(cells)
      if (terra::nlyr(ref) == 1L) colnames(values[[i]])[1] <- names(ref)
      distances[[i]] <- fields::rdist(
        sf::st_coordinates(pts[i, ]), cells[, c("x", "y")]
      )[1, ] / max_D
    }
  }
  cov_raw <- matrix(NA_real_, nrow = nrow(pts), ncol = nrow(scale_vars),
                    dimnames = list(point_ids, scale_vars$covariate))
  for (i in seq_len(nrow(pts))) {
    cov_raw[i, ] <- .msr_eval_scale_vars(
      d = distances[[i]], cov_df = values[[i]],
      scale_vars = scale_vars, sigma = sigma / max_D, shape = shape,
      kernel = kernel, unit_conv = max_D,
      resolution = terra::res(ref)[1], n_cols = terra::ncol(ref)
    )
  }
  validate_covariates_before_scale(cov_raw,
                                   context = "pooled group covariates",
                                   scale_vars = scale_vars)
  scaled <- scale(cov_raw)
  kernel_specs <- scale_vars[scale_vars$type == "kernel", , drop = FALSE]
  binned <- if (isTRUE(bin) && nrow(kernel_specs)) {
    .msr_build_kernel_bins(distances, values, kernel_specs$covariate,
                           kernel_specs$source, nbins, point_ids)
  } else NULL
  drop_cells <- isFALSE(store_cell_data) && !is.null(binned) &&
    all(scale_vars$type == "kernel")
  if (isFALSE(store_cell_data) && !drop_cells && verbose) {
    message("Cell data retained because the requested covariates need it.")
  }
  out <- list(kernel_dat = as.data.frame(scaled),
              d_list = if (drop_cells) NULL else distances,
              raw_cov = if (drop_cells) NULL else values,
              kernel = kernel, shape = shape,
              min_D = floor(terra::res(ref)[1]), max_D = max_D,
              n_covs = n_optimized, unit_conv = max_D,
              sigma = sigma / max_D, scale_vars = scale_vars,
              resolution = terra::res(ref)[1], n_cols = terra::ncol(ref),
              binned = binned, cell_data_stored = !drop_cells,
              scl_params = list(mean = attr(scaled, "scaled:center"),
                                sd = attr(scaled, "scaled:scale")),
              raster_group = stats::setNames(group, point_ids))
  class(out) <- "multiScaleR_data"
  out
}

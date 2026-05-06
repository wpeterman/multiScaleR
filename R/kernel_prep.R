#' @title Kernel Scale Preparation
#'
#' @description Prepares the data inputs required for multiscale kernel
#' optimization. Extracts raster values within a buffer around each point,
#' computes pairwise distances, and calculates initial kernel-weighted covariate
#' values. The result is passed directly to \code{\link{multiScale_optim}}.
#'
#' @param pts Spatial point locations as a \code{SpatVector} (terra) or \code{sf}
#'   object. Points must be in a projected coordinate system that shares units
#'   with \code{max_D}.
#' @param raster_stack One or more raster layers of class \code{SpatRaster}
#'   (terra). Layer names must match the variable names in the model formula
#'   used with \code{\link{multiScale_optim}}, unless \code{scale_vars} is
#'   provided to define derived covariates.
#' @param max_D Positive numeric. The maximum radius (in the same units as the
#'   projection of \code{pts} and \code{raster_stack}) around each point to
#'   sample from the raster. Should be at least 2–3 times the largest expected
#'   sigma value. Raster layers must extend far enough beyond the points to
#'   fully cover buffers of this radius.
#' @param kernel Kernel function to be used for weighting raster values by
#'   distance. One of:
#'   \describe{
#'     \item{\code{"gaussian"}}{(Default) Gaussian (normal) decay: weights fall
#'       off as a normal distribution with standard deviation \code{sigma}.}
#'     \item{\code{"exp"}}{Negative exponential decay: weights fall off as
#'       \code{exp(-d / sigma)}.}
#'     \item{\code{"fixed"}}{Fixed-radius (step) buffer: all cells within
#'       \code{sigma} receive equal weight; cells beyond receive zero.}
#'     \item{\code{"expow"}}{Exponential power kernel parameterized by both
#'       \code{sigma} (scale) and \code{shape} (shape). Requires \code{shape}
#'       to be specified.}
#'   }
#' @param scale_vars Optional variable specifications created with
#'   \code{\link{msr_vars}}. When omitted, each raster layer becomes one
#'   kernel-weighted covariate with the same name, preserving historical
#'   behavior. Provide \code{scale_vars} to define derived model covariates —
#'   for example, combining a kernel-weighted mean with landscape composition
#'   or edge metrics from the same source layer. See \code{\link{kernel_var}}
#'   and \code{\link{landscape_var}}.
#' @param sigma Optional numeric vector of initial sigma values for
#'   optimization. Length must equal the number of covariates with
#'   \code{optimize = TRUE} in \code{scale_vars} (or the number of raster
#'   layers when \code{scale_vars} is not provided). Values must be in the same
#'   units as the projection. Default: \code{NULL} — initial values are
#'   generated automatically as \code{max_D / 2}, which is recommended.
#' @param shape Optional numeric vector of initial shape values when
#'   \code{kernel = "expow"}. Length requirements same as \code{sigma}. Default:
#'   \code{NULL} — starting values of 2 are generated automatically.
#' @param projected Logical. Whether \code{pts} and \code{raster_stack} are in a
#'   projected (planar) coordinate system. Currently only projected coordinates
#'   are supported. Default: \code{TRUE}.
#' @param progress Logical. Print progress bars to the console during extraction
#'   and distance calculations. Default: \code{FALSE}.
#' @param verbose Logical. Print status messages during preparation. Default:
#'   \code{TRUE}.
#'
#' @return A list of class \code{"multiScaleR_data"} containing:
#' \describe{
#'   \item{\code{kernel_dat}}{Data frame of scaled (mean-centered, unit-variance)
#'     initial kernel-weighted covariate values, one row per point and one column
#'     per covariate. Row names match \code{pts} row names or sequential integers.
#'     Use this directly to fit the initial model passed to
#'     \code{\link{multiScale_optim}}.}
#'   \item{\code{d_list}}{Named list (one element per point) of numeric distance
#'     vectors from the point to every raster cell within the buffer.}
#'   \item{\code{raw_cov}}{Named list (one element per point) of sparse matrices
#'     containing raw raster cell values within the buffer, aligned with
#'     \code{d_list}.}
#'   \item{\code{kernel}}{Character string identifying the kernel used.}
#'   \item{\code{sigma}}{Numeric vector of initial sigma values on the internal
#'     (scaled) parameter space.}
#'   \item{\code{shape}}{Numeric vector of initial shape values, or \code{NULL}.}
#'   \item{\code{min_D}}{Numeric. Approximate minimum inter-cell distance,
#'     used as the lower bound for sigma during optimization.}
#'   \item{\code{max_D}}{Numeric. The \code{max_D} value supplied by the user.}
#'   \item{\code{n_covs}}{Integer. Number of covariates with scale being
#'     optimized.}
#'   \item{\code{unit_conv}}{Numeric. Internal distance scaling factor (equals
#'     \code{max_D}).}
#'   \item{\code{scale_vars}}{Data frame of class \code{"multiScaleR_vars"}
#'     describing all covariate specifications.}
#'   \item{\code{resolution}}{Numeric. Raster cell resolution.}
#'   \item{\code{n_cols}}{Integer. Number of raster columns (used for adjacency
#'     landscape metrics).}
#'   \item{\code{scl_params}}{Named list with elements \code{mean} and \code{sd}
#'     — the centering and scaling parameters applied to \code{kernel_dat}.
#'     Stored for use by \code{\link{kernel_scale.raster}} when
#'     \code{scale_center = TRUE}.}
#' }
#'
#' @details
#' Point locations and raster layers must share a defined CRS and both must be
#' projected. All units (including \code{max_D} and any user-supplied
#' \code{sigma}) must be in the same linear unit as the projection (typically
#' metres or feet).
#'
#' The \code{max_D} buffer should be large enough to encompass the plausible
#' range of the true scale of effect. A common rule of thumb is to set
#' \code{max_D} to at least 2–3 times the largest expected sigma. After running
#' \code{\link{multiScale_optim}}, the \code{max_distance} diagnostic will warn
#' if the estimated scale approaches \code{max_D}.
#'
#' Initial \code{sigma} values do not need to be precise — the optimizer will
#' refine them. Provide explicit starting values only if the default
#' (\code{max_D / 2}) leads to convergence problems.
#'
#' Row names from \code{pts} are preserved throughout the returned object so
#' that downstream model data frames can be joined back to the original point
#' order.
#' @examples
#' library(terra)
#' pts <- vect(cbind(c(3,5,7),
#'                   c(7,5,3)))
#'
#' mat_list <- list(r1 = rast(matrix(rnorm(100),
#'                                   nrow = 10)),
#'                  r2 = rast(matrix(rnorm(100),
#'                                   nrow = 10)))
#' rast_stack <- rast(mat_list)
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = rast_stack,
#'                              max_D = 2,
#'                              kernel = 'gaussian',
#'                              sigma = NULL)
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
                        verbose = TRUE){
  unit_conv <- max_D

  kernel <- match.arg(kernel)
  validate_scalar_numeric(max_D, "max_D", positive = TRUE)
  validate_scalar_logical(projected, "projected")
  validate_scalar_logical(progress, "progress")
  validate_scalar_logical(verbose, "verbose")

  if(!inherits(raster_stack, "SpatRaster")){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }

  scale_vars <- .msr_validate_scale_vars(scale_vars = scale_vars,
                                         raster_stack = raster_stack,
                                         kernel = kernel)
  opt_scale_vars <- .msr_optimized_scale_vars(scale_vars)
  n_optimized <- nrow(opt_scale_vars)

  if (any(scale_vars$type == "landscape")) {
    raster_res <- terra::res(raster_stack)
    if (!isTRUE(all.equal(raster_res[[1]], raster_res[[2]]))) {
      stop("Landscape metric variables require square raster cells.",
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


    # Convert to list of sparse matrices
    sparse_list <- lapply(r_ext, df_to_sparse)

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
    r_ext <- exact_extract(raster_stack,
                           buff_poly,
                           # full_colnames = T,
                           # force_df = T,
                           include_xy = T,
                           include_cell = .msr_needs_cells(scale_vars),
                           progress = progress)

    # Convert to list of sparse matrices
    sparse_list <- lapply(r_ext, df_to_sparse)



    if(nlyr(raster_stack) == 1){
      re_name <- function(x){
        c_names <- colnames(x)
        c_names[1] <- names(raster_stack)
        colnames(x) <- c_names
        return(x)
      }

      r_ext <- lapply(r_ext, re_name)
      sparse_list <- lapply(sparse_list, re_name)
    }

    ## Progress bar
    D <- vector('list', dim(pts)[1])

    if(isTRUE(progress)){

      cat(paste0("\nCalculating distances...\n"))

      pb = txtProgressBar(min = 0,
                          max = dim(pts)[1],
                          initial = 0,
                          char = "*",
                          style = 3)

    }

    for (i in 1:dim(pts)[1]) {
      if(isTRUE(progress)){
        setTxtProgressBar(pb,i)
      }

      D[[i]] <- rdist(st_coordinates(pts[i,]),
                      # r_ext[r_ext$id == i, c("x","y")])[1,]
                      r_ext[[i]][, c("x","y")])[1,] / unit_conv

    }
    if(isTRUE(progress)){
      close(pb)
    }
    min_D <- floor(rdist(r_ext[[1]][1:2,c("x","y")])[1,2])
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
  names(sparse_list) <- point_ids
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
      cov_df = sparse_list[[i]],
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

  scl_df <- scale(cov.w)
  kernel_dat <- as.data.frame(scl_df)

  out <- list(kernel_dat = kernel_dat,
              d_list = D,
              raw_cov = sparse_list,
              kernel = kernel,
              shape = shape,
              min_D = min_D,
              max_D = max_D,
              n_covs = n_optimized,
              unit_conv = unit_conv,
              sigma = sigma,
              scale_vars = scale_vars,
              resolution = terra::res(raster_stack)[[1]],
              n_cols = terra::ncol(raster_stack),
              scl_params = list(mean = attr(scl_df, "scaled:center"),
                                sd = attr(scl_df, "scaled:scale")))

  class(out) <- 'multiScaleR_data'
  return(out)
}

# Function to convert a data frame to a sparse matrix
df_to_sparse <- function(df) {
  cells <- if ("cell" %in% names(df)) df$cell else NULL
  value_cols <- setdiff(names(df), c("x", "y", "cell", "coverage_fraction"))
  out <- as(as.matrix(df[, value_cols, drop = FALSE]), "sparseMatrix")
  if (!is.null(cells)) {
    attr(out, "cell") <- cells
  }
  out
}

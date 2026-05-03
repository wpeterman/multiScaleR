#' @title Kernel Scale Preparation
#' @description Function to prepare data inputs for kernel scale analysis
#' @param pts Point locations provided as `SpatVector` or `sf` objects
#' @param raster_stack Raster layer(s) of class `SpatRaster`
#' @param max_D The maximum distance to consider during the scale optimization
#' @param kernel Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow'; Default: 'gaussian')
#' @param scale_vars Optional variable specifications created with `msr_vars()`.
#' If omitted, each raster layer is treated as one optimized kernel-weighted
#' covariate, preserving the historical behavior. Use `kernel_var()` and
#' `landscape_var()` inside `msr_vars()` to explicitly define derived model
#' covariates from source raster layers.
#' @param sigma Initial values for optimizing the scale parameter. Default: NULL, initial values will be automatically generated. This is recommended.
#' @param shape Initial values for optimizing the shape parameter if using exponential power kernel. Default: NULL, starting values will be automatically generated. This is recommended.
#' @param projected Logical. Are `pts` and `raster_stack` projected. Function currently requires that both are projected. Default: TRUE
#' @param progress Should progress bars be printed to console. Default: FALSE
#' @param verbose Logical. Print preparation information to the console. Default: TRUE
#' @return A list of class `multiscaleR` with necessary elements to conduct scale optimization using the `multiScale_optim` function
#' @details Spatial point locations and raster layers should have a defined projection and be the same CRS. If providing starting values for `sigma` or `shape`, it must be a vector of length equal to the number of raster layers for which scale is being assessed and should be provided in the unit of the used projection. When specifying `max_D`, ensure that your raster layers adequately extend beyond the points provided so that the surrounding landscape can be meaningfully sampled during scale optimization. Row names from `pts` are preserved in the returned kernel data, distance list, and raw covariate list so downstream model data can be aligned to the original point order.
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

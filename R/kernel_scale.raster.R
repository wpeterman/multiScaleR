#' @title Create scaled rasters
#' @description Function to create scaled rasters
#' @param raster_stack Stack of combined `SpatRaster` layers
#' @param sigma Vector of parameters listed in order to scale each raster
#' @param scale_opt If scale optimization with `multiScale_optim` has been completed, provide the `multiscaleR` object here. Default: NULL
#' @param shape Vector of parameters listed in order to scale each raster if using 'expow' kernel. Default: NULL
#' @param kernel Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow'; Default: 'gaussian')
#' @param pct_wt The percentage of the weighted density to include when applying the kernel smoothing function, Default: 0.95
#' @return `SpatRaster` object containing scaled rasters
#' @details NA
#' @examples
#' ## Not Run
#' r1 <- rast(matrix(rnorm(25^2),
#'                   nrow = 25))
#'
#' r1_s <- kernel_scale.raster(r1,
#'                             sigma = 4,
#'                             kernel = 'gaussian')
#' plot(c(r1, r1_s))
#'
#' @usage
#' kernel_scale.raster(raster_stack,
#'                     sigma,
#'                     scale_opt = NULL,
#'                     shape = NULL,
#'                     kernel = 'gaussian',
#'                     pct_wt = 0.95)
#'
#' @rdname kernel_scale.raster
#' @export
#' @importFrom Hmisc wtd.Ecdf
#' @importFrom terra crs rast subset cellFromRowCol crds focal res xyFromCell nlyr

kernel_scale.raster <- function(raster_stack,
                                sigma,
                                scale_opt = NULL,
                                shape = NULL,
                                kernel = 'gaussian',
                                pct_wt = 0.95){

  if(class(raster_stack) != 'SpatRaster'){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }

  if(!is.null(scale_opt) & class(scale_opt) == 'multiScaleR'){
    covs <- rownames(scale_opt$scale_est)
    sigma <- scale_opt$scale_est[,1]
    shape <- scale_opt$shape_est[,1]
    kernel <- scale_opt$kernel_inputs$kernel

    if(isFALSE(covs %in% names(raster_stack))){
      stop('optimized covariate is not present in the provided SpatRaster!')
    } else {
      raster_stack <- subset(raster_stack, covs)
    }
  } else {
    covs <- names(raster_stack)
  }

  # if(kernel != 'gaussian'){
  #   stop('Only gaussian kernel is currently supported with this function!')
  # }

  if(length(sigma) != nlyr(raster_stack)){
    warning("Number of sigma values must equal the number of raster layers!!!  \n  All raster layers will be smoothed using the same sigma value")
    sigma <- rep(sigma[1], nlyr(raster_stack))
  }


  smooth_list <- wt_list <-  vector('list', length(sigma))

  for(i in 1:length(sigma)){
    lyr <- covs[i]
    d <- seq(1, 1e6,
             length.out = 1e6)
    wt <- scale_type(d = d,
                     kernel = kernel,
                     sigma = sigma[i],
                     shape = shape[i],
                     output = 'wts')

    mx <- Hmisc::wtd.Ecdf(d, weights = wt)
    mx <- round(mx$x[which(mx$ecdf > pct_wt)[1]], digits = -1)

    r_res <- res(raster_stack)[1]
    focal_d <- ceiling(mx / r_res) * 2

    if((focal_d %% 2) == 0) {
      focal_d <- focal_d + 1
    }

    mat <- matrix(rnorm(focal_d^2), focal_d, focal_d)
    r_wt <- rast(mat)
    terra::crs(r_wt) <- crs(raster_stack)
    cntr_crd <- cellFromRowCol(r_wt, focal_d/2, focal_d/2)
    cntr_crd <- xyFromCell(r_wt, ceiling(length(mat)/2))
    cell_crds <- crds(r_wt)
    r_wt[] <- rdist(cntr_crd, cell_crds)[1,] * r_res
    r_wt[] <- scale_type(d = as.vector(r_wt),
                         kernel = kernel,
                         sigma = sigma[i],
                         shape = shape[i],
                         output = 'wts')


    wt_mat <- as.matrix(r_wt, wide = T)

    cat(paste0("\nSmoothing SpatRaster ",i, " of ", length(sigma), ": ",lyr,"\n"))

    smooth_list[[i]] <- focal(raster_stack[[lyr]],
                              wt_mat,
                              mean,
                              expand = T)
  }
  smooth_stack <- rast(smooth_list)
  names(smooth_stack) <- names(raster_stack)

  # class(smooth_stack) <- c('multiScaleR_SpatRaster')
  return(smooth_stack)
}

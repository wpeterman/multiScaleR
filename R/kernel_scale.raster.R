#' @title Create scaled rasters
#' @description Function to create scaled rasters
#' @param raster_stack Stack of combined `SpatRaster` layers
#' @param sigma Vector of parameters listed in order to scale each raster
#' @param multiScaleR If scale optimization with `multimultiScaleRim` has been completed, provide the `multiscaleR` object here. You can also pass an object of class "multiScaleR_data" created using `kernel_prep`.  Default: NULL
#' @param shape Vector of parameters listed in order to scale each raster if using 'expow' kernel. Default: NULL
#' @param kernel Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow'; Default: 'gaussian')
#' @param pct_wt The percentage of the weighted density to include when applying the kernel smoothing function, Default: 0.975
#' @param fft Logical. If TRUE (Default), a fast Fourier transformation will be used to smooth the raster surface. See details.
#' @param scale_center Logical. If `TRUE`, raster values are scaled and centered accordingly to the data used to fit the model. Necessary when predicting model results across the landscape.
#' @param clamp Logical. If `TRUE`, scaled values are clamped to the covariate range in the model data.
#' @param pct_mx Numeric. If `clamp` is `TRUE`, this value specifies the amount (percentage; positive or negative) by which to expand/contract the min/max range when clamping. Can range from -0.99–0.99 (Default = 0).
#' @param na.rm Logical. If TRUE (Default), NA values are removed from the weighted mean calculation.
#' @param ... Not used
#' @return `SpatRaster` object containing scaled rasters
#' @details The fast Fourier transformation is substantially faster when scaling large raster surfaces with large kernel areas. There will be some edge effects on the outer boundaries.
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
#'
#' @rdname kernel_scale.raster
#' @export
#' @importFrom terra crs rast subset cellFromRowCol crds focal res xyFromCell nlyr clamp nrow ncol values as.matrix

kernel_scale.raster <- function(raster_stack,
                                sigma = NULL,
                                multiScaleR = NULL,
                                shape = NULL,
                                kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                                pct_wt = 0.975,
                                fft = TRUE,
                                scale_center = FALSE,
                                clamp = FALSE,
                                pct_mx = 0,
                                na.rm = TRUE,
                                ...){

  args <- list(...)

  if ("scale_opt" %in% names(args)) {
    warning("Argument 'scale_opt' is deprecated. Use 'multiScaleR' instead.", call. = FALSE)
    if (is.null(multiScaleR) & !is.null(args$scale_opt)) {
      multiScaleR <- args$scale_opt
    }
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")){
    covs <- rownames(multiScaleR$scale_est)
    sigma <- multiScaleR$scale_est[,1]
    shape <- multiScaleR$shape_est[,1]
    kernel <- multiScaleR$kernel_inputs$kernel

    if(!any(covs %in% names(raster_stack))){
      stop('optimized covariate(s) not present in the provided SpatRaster!')
    } else {
      var <- intersect(covs, names(raster_stack))
      c <- which(covs %in% var)
      raster_stack <- subset(raster_stack, var)
      sigma <- sigma[c]
      shape <- shape[c]
    }
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")){
    covs <- colnames(multiScaleR$kernel_dat)
    var <- intersect(covs, names(raster_stack))
    c <- which(var %in% covs)

    sigma <- multiScaleR$sigma[c] * multiScaleR$unit_conv
    shape <- multiScaleR$shape[c]
    kernel <- multiScaleR$kernel

    if(isFALSE(covs %in% names(raster_stack))){
      stop('optimized covariate is not present in the provided SpatRaster!')
    } else {
      raster_stack <- subset(raster_stack, covs)
    }
  } else {
    covs <- var <- names(raster_stack)
  }

  kernel <- match.arg(kernel)

  if(is.null(sigma)){
    stop("sigma values must be specified\n")
  }

  if(!inherits(raster_stack, "SpatRaster")){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }

  if(!is.logical(na.rm)){
    stop("`na.rm` must be TRUE or FALSE")
  }

  if(!is.logical(fft)){
    stop("`fft` must be TRUE or FALSE")
  }

  if(length(sigma) != nlyr(raster_stack)){
    warning("Number of sigma values must equal the number of raster layers!!!  \n  All raster layers will be smoothed using the same sigma value")
    sigma <- rep(sigma[1], nlyr(raster_stack))
  }


  smooth_list <- wt_list <-  vector('list', length(sigma))
  out <- rast(raster_stack[[1]])

  for(i in 1:length(sigma)){
    # lyr <- covs[i]
    lyr <- var[i]

    mx <- kernel_dist(kernel = kernel,
                      sigma = sigma[i],
                      shape = shape[i],
                      prob = pct_wt)

    r_res <- res(raster_stack)[1]
    focal_d <- ceiling(mx / r_res) * 2

    if((focal_d %% 2) == 0) {
      focal_d <- focal_d + 1
    }

    mat <- matrix(rnorm(focal_d^2), focal_d, focal_d)
    r_wt <- rast(mat)
    crs(r_wt) <- crs(raster_stack)
    cntr_crd <- cellFromRowCol(r_wt, focal_d/2, focal_d/2)
    cntr_crd <- xyFromCell(r_wt, ceiling(length(mat)/2))
    cell_crds <- crds(r_wt)
    r_wt[] <- rdist(cntr_crd, cell_crds)[1,] * r_res
    r_wt[] <- scale_type_r(d = as.vector(r_wt),
                           kernel = kernel,
                           sigma = sigma[i],
                           shape = shape[i],
                           output = 'wts')


    wt_mat <- as.matrix(r_wt, wide = T)

    cat(paste0("\nSmoothing spatRaster ",i, " of ", length(sigma), ": ",lyr," at sigma = ",floor(sigma[i]),"\n"))

    if(isTRUE(fft)){
      mat <- as.matrix(raster_stack[[lyr]], wide = T)

      fft_mat <- fft_convolution(mat,
                                 wt_mat,
                                 fun = "mean",
                                 na.rm = na.rm)
      values(out) <- as.vector(fft_mat)

      smooth_list[[i]] <- out
    } else {
      smooth_list[[i]] <- focal(raster_stack[[lyr]],
                                w = wt_mat,
                                fun = "mean",
                                na.rm = na.rm,
                                expand = F)
    }

  }
  smooth_stack <- rast(smooth_list)
  names(smooth_stack) <- names(raster_stack)

  if(isTRUE(scale_center) &&
     (inherits(multiScaleR, "multiScaleR_data") || inherits(multiScaleR, "multiScaleR"))){
    smooth_stack <- scale_center_raster(r = smooth_stack,
                                        multiScaleR = multiScaleR,
                                        clamp = clamp,
                                        pct_mx = pct_mx)
  }

  return(smooth_stack)
}

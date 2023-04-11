# Raster simulation -------------------------------------------------------
#' Function to simulate raster surfaces
#' @description
#' FUnction to create four spatRaster surfaces, simulated using the \code{\link[nlm_gaussianfield]{NLMR}} function from the `NLMR` package.
#'  description
#' @param dim Dimension (number of cells) on a side a square raster (Default = 100)
#' @param resolution Resolution of raster cells (Default = 10)
#' @param autocorr_range1 Optional, range in cells for autocorrelation of first binary raster. If NULL (default), autocorrelation range will be 0.03 of specified dimension.
#' @param autocorr_range2 Optional, range in cells for autocorrelation of second second raster. If NULL (default), autocorrelation range will be 0.15 of specified dimension.
#' @param mag_var Magnitude of variation over the entire landscape (Default = 10)
#' @param nug Magnitude of variation in the scale of autocorr_range, smaller values lead to more homogeneous landscapes. (Default = 3)
#' @param plot Logical. If TRUE, the spatRaster stack will be plotted following the simulation
#' @param user_seed Optional seed to replicate simulated surfaces
#' @param ... Additional arguments. Not currently used
#' @usage
#' sim_rast(dim = 100,
#'          resolution = 10,
#'          autocorr_range1 = NULL,
#'          autocorr_range2 = NULL,
#'          mag_var = 10,
#'          nug = 3,
#'          plot = FALSE,
#'          user_seed = NULL,
#'          ...)
#'
#' @return
#' Four spatRaster surfaces. Two 1/0 binary surfaces and two continuous surfaces.
#' @export
#'
#' @seealso
#' \code{\link[nlm_gaussianfield]{NLMR}}
#' @examples
#'
#' sim1 <- sim_rast()
#'
#' sim2 <- sim_rast(dim = 150,
#'                  resolution = 25)
#'
#'
#' @details
#' Requires `NLMR` and `fields` packages to be installed. This is a simple wrapped to create four different raster surfaces. Surfaces differ in the range of autocorrelation. Binary surfaces are created by thresholding continuous values of the Gaussian random surface.
#'
#' @rdname sim_rast
#' @importFrom terra as.int rast plot

sim_rast <- function(dim = 100,
                     resolution = 10,
                     autocorr_range1 = NULL,
                     autocorr_range2 = NULL,
                     mag_var = 10,
                     nug = 3,
                     plot = FALSE,
                     user_seed = NULL,
                     ...){

  if(is.null(autocorr_range1)){
    autocorr_range1 <- floor(0.03 * dim)
  }

  if(is.null(autocorr_range2)){
    autocorr_range2 <- floor(0.15 * dim)
  }

  ## Make binary surface
  if(!is.null(user_seed)){
    user_seed1 <- user_seed
    user_seed2 <- user_seed * 9
    user_seed3 <- user_seed * 99
    user_seed4 <- user_seed * 55
  } else {
    user_seed1 <- user_seed2 <- user_seed3 <- user_seed4 <- user_seed
  }

  bin1 <- NLMR::nlm_gaussianfield(ncol = dim,
                                  nrow = dim,
                                  resolution = resolution,
                                  autocorr_range = autocorr_range1,
                                  mag_var = mag_var,
                                  nug = nug,
                                  user_seed = user_seed1)
  bin1 <- (bin1 >= 0.55)

  bin2 <- NLMR::nlm_gaussianfield(ncol = dim,
                                  nrow = dim,
                                  resolution = resolution,
                                  autocorr_range = autocorr_range2,
                                  mag_var = mag_var,
                                  nug = nug,
                                  user_seed = user_seed2)
  bin2 <- (bin2 < 0.4)

  ## Make continuous surface
  cont1 <- NLMR::nlm_gaussianfield(ncol = dim,
                                   nrow = dim,
                                   resolution = resolution,
                                   autocorr_range = floor(autocorr_range1*0.75),
                                   mag_var = mag_var,
                                   nug = nug,
                                   user_seed = user_seed3)

  cont2 <- NLMR::nlm_gaussianfield(ncol = dim,
                                   nrow = dim,
                                   resolution = resolution,
                                   autocorr_range = floor(autocorr_range2*1.25),
                                   mag_var = mag_var,
                                   nug = nug,
                                   user_seed = user_seed4)

  # r_stack <- stack(list(bin1 = bin1,
  #                       bin2 = bin2,
  #                       cont1 = cont1,
  #                       cont2 = cont2))

  r_stack <- rast(list(bin1 = as.int(rast(bin1)),
                       bin2 = as.int(rast(bin2)),
                       cont1 = rast(cont1),
                       cont2 = rast(cont2)))

  if(isTRUE(plot)){
    plot(r_stack)
  }

  return(r_stack)
}



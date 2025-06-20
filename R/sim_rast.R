# Raster simulation -------------------------------------------------------
#' Function to simulate raster surfaces
#' @description
#' Function to create four spatRaster surfaces
#' @param dim Dimension (number of cells) on a side a square raster (Default = 100)
#' @param resolution Resolution of raster cells (Default = 10)
#' @param autocorr_range1 Optional, Numeric. Spatial correlation range in map cells. Controls the decay of the exponential covariance. If NULL (default), autocorrelation range will be 5\% of specified dimension.
#' @param autocorr_range2 Optional, Numeric. Spatial correlation range in map cells. Controls the decay of the exponential covariance. If NULL (default), autocorrelation range will be 25\% of specified dimension.
#' @param sill Numeric. Variance (partial sill) of the random field (default = 10).
#' @param plot Logical. If TRUE, the spatRaster stack will be plotted following the simulation
#' @param user_seed Optional seed to replicate simulated surfaces
#' @param ... Additional arguments. Not currently used
#' @return
#' Four spatRaster surfaces. Two 1/0 binary surfaces and two continuous surfaces.
#' @export
#' @examples
#'
#' sim1 <- sim_rast()
#'
#' sim2 <- sim_rast(dim = 150,
#'                  resolution = 25)
#'
#'
#' @details
#' This is a simple wrapper to create four different raster surfaces. Surfaces differ in the range of autocorrelation. Binary surfaces are created by thresholding continuous values of the Gaussian random surface.
#'
#' @rdname sim_rast
#' @importFrom terra as.int rast plot minmax ext<- values<- crs<-
#' @importFrom stats rpois rnbinom plogis rbinom dist predict

sim_rast <- function(dim = 100,
                     resolution = 10,
                     autocorr_range1 = NULL,
                     autocorr_range2 = NULL,
                     sill = 10,
                     plot = FALSE,
                     user_seed = NULL,
                     ...){

  if(is.null(autocorr_range1)){
    autocorr_range1 <- floor(0.05 * dim)
  }

  if(is.null(autocorr_range2)){
    autocorr_range2 <- floor(0.25 * dim)
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

  r <- simulate_fft_grf(dim = dim,
                        resolution = resolution,
                        range = autocorr_range1,
                        sill = sill,
                        user_seed = user_seed1)

  r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])
  bin1 <- as.numeric((r >= 0.55))

  r <- simulate_fft_grf(dim = dim,
                        resolution = resolution,
                        range = autocorr_range2,
                        sill = sill,
                        user_seed = user_seed2)
  r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])
  bin2 <- as.numeric(r < 0.4)

  ## Make continuous surface
  r <- simulate_fft_grf(dim = dim,
                        resolution = resolution,
                        range = floor(autocorr_range1*0.75),
                        sill = sill,
                        user_seed = user_seed3)
  r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])
  cont1 <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])


  r <- simulate_fft_grf(dim = dim,
                        resolution = resolution,
                        range = floor(autocorr_range2*1.25),
                        sill = sill,
                        user_seed = user_seed4)
  r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])
  cont2 <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])

  r_stack <- rast(list(bin1 = (bin1),
                       bin2 = (bin2),
                       cont1 = (cont1),
                       cont2 = (cont2)))

  if(isTRUE(plot)){
    plot(r_stack)
  }

  return(r_stack)
}


simulate_fft_grf <- function(dim = 100, resolution = 10, range = 15, sill = 10, user_seed = NULL) {
  if (!is.null(user_seed)) set.seed(user_seed)

  n <- dim
  # Create grid of distances from center
  x <- seq(-n / 2, n / 2 - 1)
  y <- seq(-n / 2, n / 2 - 1)
  d2 <- outer(x^2, y^2, "+") # squared Euclidean distances

  # Exponential covariance function
  cov_mat <- sill * exp(-sqrt(d2) / range)

  # FFT of covariance (circulant embedding assumption)
  eigvals <- Re(fft(cov_mat))

  # Clamp eigenvalues to non-negative to avoid NaNs in sqrt
  eigvals <- pmax(eigvals, 0)

  # White noise in Fourier space
  Z <- matrix(rnorm(n * n), nrow = n)
  Z_fft <- fft(Z)

  # Element-wise multiply and invert FFT
  sim_complex <- Z_fft * sqrt(eigvals)
  sim <- Re(fft(sim_complex, inverse = TRUE)) / (n * n)

  # Rescale to [0, 1] and convert to SpatRaster
  sim <- (sim - min(sim)) / (max(sim) - min(sim))
  rast(sim, extent = ext(0, dim * resolution, 0, dim * resolution))
}

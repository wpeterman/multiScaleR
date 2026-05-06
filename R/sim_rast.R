# Raster simulation -------------------------------------------------------
#' Simulate spatially autocorrelated raster surfaces
#'
#' @description
#' Creates a \code{SpatRaster} stack of four simulated landscape surfaces for
#' use with \code{\link{sim_dat}} and \code{\link{sim_dat_unmarked}}. The stack
#' contains two binary (0/1) and two continuous (0–1) surfaces with differing
#' spatial autocorrelation ranges, allowing simulation of multiscale ecological
#' processes.
#'
#' @param dim Positive integer. Number of cells on each side of the square
#'   raster. The output will be a \code{dim x dim} grid. Default: \code{100}.
#' @param resolution Positive numeric. Cell size (edge length) in map units.
#'   The raster extent will be \code{dim * resolution} in each direction.
#'   Default: \code{10}.
#' @param autocorr_range1 Optional positive numeric. Spatial autocorrelation
#'   range (in map cells) for surfaces 1 and 3 (\code{bin1}, \code{cont1}).
#'   Controls the decay rate of the exponential covariance: larger values
#'   produce smoother, more broadly correlated patterns. If \code{NULL}
#'   (default), set to 5\% of \code{dim} (fine-scale autocorrelation).
#' @param autocorr_range2 Optional positive numeric. Autocorrelation range for
#'   surfaces 2 and 4 (\code{bin2}, \code{cont2}). If \code{NULL} (default),
#'   set to 25\% of \code{dim} (broad-scale autocorrelation).
#' @param sill Positive numeric. Partial sill (variance) of the Gaussian
#'   random field used for all four surfaces. Default: \code{10}.
#' @param plot Logical. If \code{TRUE}, the raster stack is plotted using
#'   \code{terra::plot} after simulation. Default: \code{FALSE}.
#' @param user_seed Optional integer seed for reproducibility. Different
#'   transformations of \code{user_seed} are applied to each of the four
#'   surfaces so they are independent but fully reproducible. Default:
#'   \code{NULL}.
#' @param ... Additional arguments. Not currently used.
#'
#' @return A \code{SpatRaster} with four named layers:
#' \describe{
#'   \item{\code{bin1}}{Binary (0/1) surface with fine-scale autocorrelation
#'     (\code{autocorr_range1}). Values of 1 where the underlying continuous
#'     field exceeds its 55th percentile.}
#'   \item{\code{bin2}}{Binary (0/1) surface with broad-scale autocorrelation
#'     (\code{autocorr_range2}). Values of 1 where the underlying continuous
#'     field falls below its 40th percentile.}
#'   \item{\code{cont1}}{Continuous surface rescaled to [0, 1] with fine-scale
#'     autocorrelation (75\% of \code{autocorr_range1}).}
#'   \item{\code{cont2}}{Continuous surface rescaled to [0, 1] with broad-scale
#'     autocorrelation (125\% of \code{autocorr_range2}).}
#' }
#' All layers span the same spatial extent:
#' \code{[0, dim * resolution] x [0, dim * resolution]}.
#'
#' @details
#' Each surface is generated as a Gaussian random field using a fast Fourier
#' transform (FFT) circulant embedding approach with an exponential covariance
#' function. The underlying continuous fields are normalized to [0, 1] before
#' thresholding (binary surfaces) or returning directly (continuous surfaces).
#'
#' The two autocorrelation ranges allow simulation of covariates that operate
#' at different spatial scales — a common scenario in landscape ecology where
#' some resources are patchily distributed at fine scales and others vary
#' broadly across the study area.
#'
#' When \code{user_seed} is provided, independent but reproducible seeds are
#' derived for each surface as multiples of \code{user_seed}.
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

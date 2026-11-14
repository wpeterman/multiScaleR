# Raster simulation -------------------------------------------------------
#' Simulate Spatially Autocorrelated Raster Surfaces
#'
#' @description
#' Creates a `SpatRaster` with four simulated landscape layers, two binary
#' (0 or 1) and two continuous (0 to 1), for testing and teaching. Each pair
#' has one layer with a fine and one with a broad autocorrelation range. Use the layers with
#' [sim_dat()] or [sim_dat_unmarked()] to simulate responses at known scales
#' of effect, then check that [multiScale_optim()] recovers them.
#'
#' @param dim Positive integer. Number of cells along each side of the square
#'   raster, so the output is a `dim` x `dim` grid. With the default ranges,
#'   use at least 40, because smaller values round the fine-scale range down
#'   to zero cells and `cont1` (and below 20, `bin1`) returns `NaN`. See
#'   Details. Default: `100`.
#' @param resolution Positive numeric. Cell size (edge length) in map units.
#'   The raster spans 0 to `dim * resolution` in each direction. Default:
#'   `10`.
#' @param autocorr_range1 Optional positive numeric. Autocorrelation range, in
#'   cells, of the fields behind `bin1` and `cont1`. Larger values give
#'   smoother, more broadly correlated patterns. Default `NULL` uses
#'   `floor(0.05 * dim)` (fine-scale autocorrelation).
#' @param autocorr_range2 Optional positive numeric. Autocorrelation range, in
#'   cells, of the fields behind `bin2` and `cont2`. Default `NULL` uses
#'   `floor(0.25 * dim)` (broad-scale autocorrelation).
#' @param sill Positive numeric. Variance (partial sill) of the exponential
#'   covariance used for all four fields. Each field is rescaled to 0 to 1
#'   afterward, so the sill does not change the returned values. Default:
#'   `10`.
#' @param plot Logical. If `TRUE`, plot the four layers with `terra::plot()`.
#'   Default: `FALSE`.
#' @param user_seed Optional integer seed for reproducibility. Each layer uses
#'   a different multiple of the seed (1, 9, 99, and 55 times `user_seed`), so
#'   the layers differ but are fully reproducible. The seed also sets the
#'   global random number stream. Default: `NULL` (not reproducible).
#' @param ... Not used.
#'
#' @return A `SpatRaster` with `dim` rows and columns, cell size `resolution`,
#' extent 0 to `dim * resolution` in x and y, no coordinate reference system,
#' and four layers:
#' * `bin1`: binary, fine-scale (`autocorr_range1`). 1 where the field,
#'   rescaled to 0 to 1, is at least 0.55, and 0 elsewhere.
#' * `bin2`: binary, broad-scale (`autocorr_range2`). 1 where the rescaled
#'   field is below 0.4, and 0 elsewhere.
#' * `cont1`: continuous, rescaled to 0 to 1, with a range of
#'   `floor(0.75 * autocorr_range1)` cells.
#' * `cont2`: continuous, rescaled to 0 to 1, with a range of
#'   `floor(1.25 * autocorr_range2)` cells.
#'
#' The thresholds are on the rescaled values, not percentiles, so the share
#' of cells equal to 1 in `bin1` and `bin2` varies from one simulation to the
#' next.
#'
#' @details
#' Each layer comes from a Gaussian random field generated with a fast
#' Fourier transform (FFT) method from the exponential covariance function
#' `sill * exp(-d / range)`, with `d` in cells. Each field is rescaled to run
#' from 0 to 1 before it is thresholded (binary layers) or returned
#' (continuous layers). A range rounded down to zero cells leaves the
#' covariance undefined and fills the layer with `NaN`. This happens for
#' `cont1` whenever `autocorr_range1` is below 2, and for `bin1` whenever
#' `autocorr_range1` is 0. With the default ranges, that means `dim` below 40
#' and below 20, respectively.
#'
#' Two ranges let you simulate covariates that act at different spatial
#' scales, as when one resource is patchy at fine scales and another varies
#' broadly across the study area. `vignette("v02-user-guide", package =
#' "multiScaleR")` uses these layers to simulate zero-inflated, N-mixture,
#' and occupancy data.
#'
#' @seealso [sim_dat()], [sim_dat_unmarked()], [kernel_prep()]
#'
#' @examples
#' ## Four 50 x 50 layers with 20 m cells (1,000 m on a side)
#' r_sim <- sim_rast(dim = 50, resolution = 20, user_seed = 555)
#' r_sim
#' # Binary layers in the top row, continuous below. Layers ending in 1 are
#' # fine-scale and layers ending in 2 broad-scale
#' terra::plot(r_sim)
#' # Share of cells equal to 1 in each binary layer
#' terra::global(r_sim[[c("bin1", "bin2")]], "mean")
#'
#' ## Broader patterns: raise the fine-scale range from 2 to 10 cells
#' r_smooth <- sim_rast(dim = 50, resolution = 20, autocorr_range1 = 10,
#'                      user_seed = 555)
#' terra::plot(r_smooth[[c("bin1", "cont1")]])
#'
#' @export
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

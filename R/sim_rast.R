# Raster simulation -------------------------------------------------------
#' Function to simulate raster surfaces
#' @description
#' Function to create four spatRaster surfaces, simulated using the \code{\link[gstat]{gstat}} package.
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
#' @importFrom gstat gstat vgm
#' @examples
#'
#' sim1 <- sim_rast()
#'
#' sim2 <- sim_rast(dim = 150,
#'                  resolution = 25)
#'
#'
#' @details
#' Requires `NLMR` package to be installed (https://github.com/ropensci/NLMR). This is a simple wrapper to create four different raster surfaces. Surfaces differ in the range of autocorrelation. Binary surfaces are created by thresholding continuous values of the Gaussian random surface.
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

  grid <- expand.grid(
    x = seq(1, dim, length.out = dim),  # x-coordinates
    y = seq(1, dim, length.out = dim)   # y-coordinates
  )

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

  r <- terra::rast(nrows = dim, ncols = dim,
                   xmin = 0, xmax = dim*resolution,
                   ymin = 0, ymax = dim*resolution,
                   resolution = resolution)
  # terra::values(r) <- as.vector(terra::rast(sim_data[, c("x", "y", "sim1")]))
  # r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])

  # plot(r, main = "Gaussian Random Field")

  set.seed(user_seed1)
  g <- gstat(
    formula = z ~ 1,          # Simulate a mean-zero field
    locations = ~x + y,       # Coordinates
    dummy = TRUE,             # Generate a dummy variable
    beta = 0,                 # Mean value (0 for centered Gaussian)
    model = vgm(              # Variogram model
      psill = mag_var,              # Partial sill (variance)
      model = "Exp",          # Exponential covariance
      range = autocorr_range1,  # Range of spatial correlation
      nugget = nug            # Nugget effect (added uncorrelated noise)
    ),
    nmax = 30                 # Number of nearest neighbors to consider
  )

  # Simulate the Gaussian field
  sim_data <- predict(g, newdata = grid, nsim = 1)
  terra::values(r) <- as.vector(terra::rast(sim_data[, c("x", "y", "sim1")]))
  r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])
  bin1 <- as.numeric((r >= 0.55))

  set.seed(user_seed2)
  g <- gstat(
    formula = z ~ 1,          # Simulate a mean-zero field
    locations = ~x + y,       # Coordinates
    dummy = TRUE,             # Generate a dummy variable
    beta = 0,                 # Mean value (0 for centered Gaussian)
    model = vgm(              # Variogram model
      psill = mag_var,              # Partial sill (variance)
      model = "Exp",          # Exponential covariance
      range = autocorr_range2,  # Range of spatial correlation
      nugget = nug            # Nugget effect (added uncorrelated noise)
    ),
    nmax = 30                 # Number of nearest neighbors to consider
  )
  # Simulate the Gaussian field
  sim_data <- predict(g, newdata = grid, nsim = 1)
  terra::values(r) <- as.vector(terra::rast(sim_data[, c("x", "y", "sim1")]))
  r <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])
  bin2 <- as.numeric(r < 0.4)

  ## Make continuous surface
  set.seed(user_seed3)
  g <- gstat(
    formula = z ~ 1,          # Simulate a mean-zero field
    locations = ~x + y,       # Coordinates
    dummy = TRUE,             # Generate a dummy variable
    beta = 0,                 # Mean value (0 for centered Gaussian)
    model = vgm(              # Variogram model
      psill = mag_var,              # Partial sill (variance)
      model = "Exp",          # Exponential covariance
      range = floor(autocorr_range1*0.75),  # Range of spatial correlation
      nugget = nug            # Nugget effect (added uncorrelated noise)
    ),
    nmax = 30                 # Number of nearest neighbors to consider
  )
  # Simulate the Gaussian field
  sim_data <- predict(g, newdata = grid, nsim = 1)
  terra::values(r) <- as.vector(terra::rast(sim_data[, c("x", "y", "sim1")]))
  cont1 <- (r - minmax(r)[1]) / (minmax(r)[2] - minmax(r)[1])


  set.seed(user_seed3)
  g <- gstat(
    formula = z ~ 1,          # Simulate a mean-zero field
    locations = ~x + y,       # Coordinates
    dummy = TRUE,             # Generate a dummy variable
    beta = 0,                 # Mean value (0 for centered Gaussian)
    model = vgm(              # Variogram model
      psill = mag_var,              # Partial sill (variance)
      model = "Exp",          # Exponential covariance
      range = floor(autocorr_range2*1.25),  # Range of spatial correlation
      nugget = nug            # Nugget effect (added uncorrelated noise)
    ),
    nmax = 30                 # Number of nearest neighbors to consider
  )
  # Simulate the Gaussian field
  sim_data <- predict(g, newdata = grid, nsim = 1)
  terra::values(r) <- as.vector(terra::rast(sim_data[, c("x", "y", "sim1")]))
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



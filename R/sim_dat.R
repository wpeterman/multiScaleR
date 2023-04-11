# Data sim function -------------------------------------------------------
#' @title Simulate data for optimizing scales of effect
#' @description
#'  Function to simulate data with known scales of effect from spatial spatRaster variables
#'
#' @param alpha Intercept term for GLM (Default = 1)
#' @param beta Slope term(s) for GLM. Should be vector equal in length to number of spatRaster surfaces provided
#' @param kernel Type of kernel transformation. Valid options are 'gaussian', 'exp' (negative exponential), 'expow' (exponential power), and 'fixed' fixed width buffer. (Default = 'gaussian')
#' @param type Type of response data to simulate. Valid options are 'count' for Poisson distributed count; 'count_nb' for negative binomial counts; 'occ' for binomial response; and 'gaussian' for normally distributed response.
#' 'count' for normally distributed response (Default = 'count')
#' @param StDev If specifying 'count_nb' or 'gaus' for type, this is the dispersion term for those respective processes (Default = 0.5)
#' @param n_points Number of spatial sample points (Default = 50).
#' @param min_D Minimum distance between points. Function will attempt to create the number of sample points specified while honoring this minimum distance.
#' @param raster_stack A spatRaster object
#' @param sigma The scale term dictating the rate of decay with distance
#' @param shape If using an exponential power function, the shape parameter must also be specified. Values between 1-50 are generally valid
#' @param max_D The maximum distance surrounding spatial points to consider. This typically needs to be >= 2.5x greater than sigma
#' @param user_seed Optional seed to reproduce simulation
#' @param ... Additional arguments. Not currently used
#'
#' @return
#' Returns a list containing:
#' \tabular{ll}{
#' \tab * obs --> The simulated response variable \cr
#' \tab * df --> A data frame with the simulated response (obs) as well as the true kernel weighted mean values for each raster surface included \cr
#' \tab * pts --> An `sf` object with the simulated spatial point locations \cr
#'   }
#' @export
#'
#' @details
#' This function distributes sample points across the landscape on a hexagonal grid, then subsamples to the specified number. The weighted values of each landscape are determined according to the simulation parameters, then the specified response is generated.
#'
#' @examples
#' rs <- sim_rast()
#' rs <- terra::subset(rs, c(1,3))
#' s_dat <- sim_dat(alpha = 0.5,
#'                  beta = c(0.75,-0.75),
#'                  kernel = 'gaussian',
#'                  sigma = c(75, 150),
#'                  type = 'count',
#'                  raster_stack = rs,
#'                  max_D = 400)
#'
#' plot(s_dat$df$y ~ s_dat$df$bin1)
#' plot(s_dat$df$y ~ s_dat$df$cont1)
#'
#' @importFrom terra as.polygons ext nlyr rast
#' @importFrom sf  st_as_sf st_make_grid

sim_dat <- function(alpha = 1,
                    beta = NULL,
                    kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                    type = c('count', 'count_nb', 'occ', 'gaussian'),
                    StDev = 0.5,
                    n_points = 50,
                    min_D = NULL,
                    raster_stack = NULL,
                    sigma = NULL,
                    shape = NULL,
                    max_D = NULL,
                    user_seed = NULL)
{

  kernel <- match.arg(kernel)
  type <- match.arg(type)

  if(!is.null(user_seed)){
    set.seed(user_seed)
  }

  if(is.null(min_D)){
    min_D <- 1.55 * max(sigma)
  }

  if(isFALSE((length(beta) <= nlyr(raster_stack)) &
             (length(beta) <= length(sigma)))){
    stop("The number of beta coefficients, sigma values and/or raster layers must be equal!!!")
  }

  if(kernel == 'expow' & is.null(shape)){
    stop("Shape parameter(s) must be specified when using the `expow` kernel!")
  }

  s_ext <- as.vector(ext(raster_stack[[1]]))
  min_x <- min_y <- floor(s_ext[1] + (s_ext[2] * 0.15))
  max_x <- max_y <- floor(s_ext[2] - (s_ext[2] * 0.15))
  r <- rast()
  ext(r) <- c(min_x, max_x, min_y, max_y)
  poly <- terra::as.polygons(ext(c(min_x, max_x, min_y, max_y)))
  poly_sf <- st_as_sf(poly)

  pts <- 0
  while(length(pts) < n_points){
    min_D <- min_D * 0.97
    pts <- sf::st_make_grid(poly_sf,
                            cellsize = min_D,
                            # n = n_points,
                            what = 'centers')
  }

  # pts <- spatSample(raster_stack[[1]],
  #                   size = n_points,
  #                   ext = ext(r),
  #                   as.points = TRUE)

  # names(pts) <- 'sample_points'
  pts <- sf::st_as_sf(pts)

  pts <- pts[sample(dim(pts)[1], n_points),]

  kernel_out <- kernel_prep(pts = pts,
                            sigma = sigma,
                            shape = shape,
                            kernel = kernel,
                            max_D = max_D,
                            raster_stack = raster_stack,
                            projected = T)

  if(type =='count'){

    obs <- exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    # alpha + kernel_out$kernel_dat$bin1*beta[1] + kernel_out$kernel_dat$cont2*beta[2]

    obs <- rpois(n_points,
                 obs)
  } else if(type =='count_nb'){
    obs <- exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    obs <- rnbinom(n_points,
                   mu = obs,
                   size = StDev)

  } else if(type == 'occ') {

    # obs <- exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*'))) /
    #   (1 + exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*'))))

    obs <- plogis(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    obs <- rbinom(n_points,
                  1,
                  obs)
  } else {

    obs <- (alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    obs <- rnorm(n_points,
                 obs,
                 sd = StDev)
  }

  pts$obs <- obs

  # df <- data.frame(y = obs,
  #                  kernel_out$kernel_dat)
  ## Scaled
  df <- data.frame(y = obs,
                   scale(kernel_out$kernel_dat))
  names(df) <- c('y', names(raster_stack))

  return(list(obs = obs,
              df = df,
              pts = pts))
}

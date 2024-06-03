# Data sim function -------------------------------------------------------
#' @title Simulate data for optimizing scales of effect with `unmarked`
#' @description
#'  Function to simulate data with known scales of effect from spatial spatRaster variables for analysis with the R package `unmarked`
#'
#' @param alpha Intercept term for GLM (Default = 1)
#' @param beta Slope term(s) for GLM. Should be vector equal in length to number of spatRaster surfaces provided
#' @param kernel Type of kernel transformation. Valid options are 'gaussian', 'exp' (negative exponential), 'expow' (exponential power), and 'fixed' fixed width buffer. (Default = 'gaussian')
#' @param type Type of response data to simulate in `unmarked`. Valid options are 'count' for Poisson distributed count; 'count_nb' for negative binomial counts; and 'occ' for binomial response.(Default = 'count')
#' @param StDev If specifying 'count_nb' or 'gaus' for type, this is the dispersion term for those respective processes (Default = 0.5)
#' @param n_points Number of spatial sample points (Default = 50).
#' @param n_surv Number of surveys to simulate in `unmarked` (Default = 3).
#' @param det The probability of detection. (Default = 0.5)
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
#' \tab * y --> The simulated observation matrix for use in an unmarkedFrame \cr
#' \tab * df --> A data frame with the simulated response (obs) as well as the true kernel weighted mean values for each raster surface included \cr
#' \tab * pts --> An `sf` object with the simulated spatial point locations \cr
#'   }
#' @export
#'
#' @details
#' This function distributes sample points across the landscape on a hexagonal grid, then subsamples to the specified number. The weighted values of each landscape are determined according to the simulation parameters, then the specified response is generated.
#'
#' @examples
#' rs <- sim_rast(user_seed = 123)
#' rs <- terra::subset(rs, c(1,3))
#' s_dat <- sim_dat_unmarked(alpha = 1,
#'                           beta = c(0.75,-0.75),
#'                           kernel = 'gaussian',
#'                           sigma = c(75, 150),
#'                           n_points = 75,
#'                           n_surv = 5,
#'                           det = 0.5,
#'                           type = 'count',
#'                           raster_stack = rs,
#'                           max_D = 550,
#'                           user_seed = 123)
#' plot(s_dat$df$y ~ s_dat$df$bin1)
#' plot(s_dat$df$y ~ s_dat$df$cont1)
#' ## unmarked analysis
#' library(unmarked)
#' kernel_inputs <- kernel_prep(pts = s_dat$pts,
#'                              raster_stack = rs,
#'                              max_D = 550,
#'                              kernel = 'gaus')
#'
#' umf <- unmarkedFramePCount(y = s_dat$y,
#'                            siteCovs = kernel_inputs$kernel_dat)
#'
#' ## Base unmarked model
#' mod0 <- pcount(~1 ~bin1 + cont1,
#'                data = umf,
#'                K = 100)
#'
#' ## `multiscale_optim`
#' opt1 <- multiScale_optim(fitted_mod = mod0,
#'                          kernel_inputs = kernel_inputs)
#'
#' summary(opt1)
#'
#' @importFrom terra as.polygons ext nlyr rast
#' @importFrom sf  st_as_sf st_make_grid

sim_dat_unmarked <- function(alpha = 1,
                             beta = NULL,
                             kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                             type = c('count', 'count_nb', 'occ'),
                             StDev = 0.5,
                             n_points = 50,
                             n_surv = 3,
                             det = 0.5,
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

  max_D <- kernel_dist(kernel = kernel,
                       prob = 0.99,
                       sigma = max(sigma)) * 1.1

  kernel_out <- kernel_prep(pts = pts,
                            sigma = sigma,
                            shape = shape,
                            kernel = kernel,
                            max_D = max_D,
                            raster_stack = raster_stack,
                            projected = T)

  if(type =='count'){

    obs <- exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    ## Multiply by detection?
    ## C[,i] <- rbinom(n = M, size = N, prob = p[,i])

    obs_ <- rpois(n_points,obs)

    obs_mat <- matrix(NA, n_points, n_surv)
    for(i in 1:n_surv){
      obs_mat[,i] <- rbinom(n_points, obs_, prob = det)
    }



  } else if(type =='count_nb'){
    obs <- exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    obs_ <- rnbinom(n_points,
                    mu = obs,
                    size = StDev)

    obs_mat <- matrix(NA, n_points, n_surv)
    for(i in 1:n_surv){
      obs_mat[,i] <- rbinom(n_points, obs_, prob = det)
    }

  } else if(type == 'occ') {

    # obs <- exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*'))) /
    #   (1 + exp(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*'))))

    obs <- plogis(alpha + rowSums(sweep(kernel_out$kernel_dat, MARGIN = 2, beta, '*')))

    obs_ <- rbinom(n_points, 1, obs)

    obs_mat <- matrix(NA, n_points, n_surv)
    for(i in 1:n_surv){
      obs_mat[,i] <- rbinom(n_points, 1, obs_*det)
    }
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
  df <- data.frame(y = obs_,
                   scale(kernel_out$kernel_dat))
  names(df) <- c('y', names(raster_stack))


  # package unmarked data ------------------------------------------------------
  # umf <- unmarkedFrame(y = obs_mat,
  #                      siteCovs = kernel_out$kernel_dat)

  return(list(y = obs_mat,
              df = df,
              pts = pts))
}

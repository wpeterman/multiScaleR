# Data sim function -------------------------------------------------------
#' @title Simulate data for optimizing scales of effect
#'
#' @description Generates simulated response data with known kernel-weighted
#' landscape covariates at a controlled scale of effect. Useful for testing and
#' demonstrating \code{\link{multiScale_optim}} with data where the true
#' parameters are known.
#'
#' @param alpha Numeric scalar. Intercept term for the linear predictor.
#'   Default: \code{1}.
#' @param beta Numeric vector of slope coefficients, one per raster layer in
#'   \code{raster_stack}. Length must equal \code{nlyr(raster_stack)}.
#' @param kernel Character. Kernel function used to weight raster values by
#'   distance when generating the true covariate values. One of
#'   \code{"gaussian"} (default), \code{"exp"} (negative exponential),
#'   \code{"expow"} (exponential power), or \code{"fixed"} (fixed-radius
#'   buffer).
#' @param type Character. Distribution of the simulated response variable. One
#'   of:
#'   \describe{
#'     \item{\code{"count"}}{(Default) Poisson-distributed counts.}
#'     \item{\code{"count_nb"}}{Negative binomial counts with dispersion
#'       \code{StDev}.}
#'     \item{\code{"occ"}}{Bernoulli (0/1) occupancy data via a logistic
#'       link.}
#'     \item{\code{"gaussian"}}{Normally distributed continuous response with
#'       standard deviation \code{StDev}.}
#'   }
#' @param StDev Positive numeric. Dispersion parameter for \code{"count_nb"}
#'   (the \code{size} argument of \code{rnbinom}) or the standard deviation for
#'   \code{"gaussian"} responses. Default: \code{0.5}.
#' @param n_points Positive integer or a \code{SpatVector}/\code{sf} point
#'   object. When an integer, that many points are placed on a hexagonal grid
#'   covering the raster extent (trimmed to 85\% of the x/y range) and then
#'   randomly subsampled. When a spatial object is supplied, those points are
#'   used directly. Default: \code{50}.
#' @param min_D Positive numeric. Minimum inter-point spacing on the hexagonal
#'   grid. If \code{NULL} (default), set automatically to
#'   \code{1.55 * max(sigma)}.
#' @param raster_stack A \code{SpatRaster} object. Layer names become covariate
#'   names in the returned data frame.
#' @param sigma Positive numeric vector. True kernel scale parameters, one per
#'   raster layer. These are the values that \code{\link{multiScale_optim}}
#'   should recover. Must be in the same units as the raster projection.
#' @param shape Positive numeric vector. Shape parameters for the exponential
#'   power kernel (\code{kernel = "expow"}), one per raster layer. Values
#'   between 1 and 50 are typical. Required when \code{kernel = "expow"}.
#' @param max_D Positive numeric. Maximum buffer radius for
#'   \code{\link{kernel_prep}} during simulation. If \code{NULL} (default), set
#'   automatically to 110\% of the distance enclosing 99\% of the kernel weight
#'   at \code{max(sigma)}.
#' @param user_seed Optional integer seed for reproducibility. Default:
#'   \code{NULL}.
#' @param ... Additional arguments. Not currently used.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{obs}}{Numeric vector of length \code{n_points} containing the
#'     simulated response values.}
#'   \item{\code{df}}{Data frame with \code{n_points} rows and one column per
#'     raster layer plus a \code{y} column. The raster columns contain the
#'     true kernel-weighted mean values, scaled to zero mean and unit variance.
#'     This data frame can be used to fit a model for
#'     \code{\link{multiScale_optim}}.}
#'   \item{\code{pts}}{An \code{sf} POINT object with \code{n_points} rows.
#'     An \code{obs} column is appended containing the simulated response
#'     values. Pass this to \code{\link{kernel_prep}} as the \code{pts}
#'     argument.}
#' }
#'
#' @details
#' Points are distributed on a hexagonal grid across the interior of the raster
#' extent (85\% of the range in each direction to avoid edge effects), then
#' randomly subsampled to \code{n_points}. The \code{min_D} spacing controls
#' the grid resolution; if the grid produces fewer points than \code{n_points},
#' \code{min_D} is reduced iteratively by 3\% until enough points are generated.
#'
#' Kernel-weighted covariate values are computed using
#' \code{\link{kernel_prep}} at the specified \code{sigma} (and \code{shape})
#' values. These represent the "true" covariate values that the optimization
#' should recover.
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
#' @export
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
                    user_seed = NULL,
                    ...)
{

  kernel <- match.arg(kernel)
  type <- match.arg(type)
  validate_scalar_numeric(alpha, "alpha")
  validate_scalar_numeric(StDev, "StDev", positive = TRUE)

  if(is.null(raster_stack) || !inherits(raster_stack, "SpatRaster")){
    stop("`raster_stack` must be provided as a `SpatRaster` object.")
  }

  validate_numeric_vector(beta,
                          "beta",
                          length_ = nlyr(raster_stack))
  validate_numeric_vector(sigma,
                          "sigma",
                          length_ = nlyr(raster_stack),
                          positive = TRUE)

  if(!is.null(user_seed)){
    validate_scalar_numeric(user_seed, "user_seed", integerish = TRUE)
    set.seed(user_seed)
  }

  if(is.null(min_D)){
    min_D <- 1.55 * max(sigma)
  } else {
    validate_scalar_numeric(min_D, "min_D", positive = TRUE)
  }

  if(kernel == 'expow' & is.null(shape)){
    stop("Shape parameter(s) must be specified when using the `expow` kernel!")
  }
  if(kernel == 'expow'){
    validate_numeric_vector(shape,
                            "shape",
                            length_ = nlyr(raster_stack),
                            positive = TRUE)
  }

  if(is.numeric(n_points)){
    validate_scalar_numeric(n_points, "n_points", integerish = TRUE, positive = TRUE)
    s_ext <- as.vector(ext(raster_stack[[1]]))
    min_x <- min_y <- floor(s_ext[1] + (s_ext[2] * 0.15))
    max_x <- max_y <- floor(s_ext[2] - (s_ext[2] * 0.15))
    r <- rast()
    ext(r) <- c(min_x, max_x, min_y, max_y)
    poly <- as.polygons(ext(c(min_x, max_x, min_y, max_y)))
    poly_sf <- st_as_sf(poly)

    pts <- 0
    while(length(pts) < n_points){
      min_D <- min_D * 0.97
      pts <- st_make_grid(poly_sf,
                          cellsize = min_D,
                          # n = n_points,
                          what = 'centers')
    }

    pts <- st_as_sf(pts)
    pts <- pts[sample(dim(pts)[1], n_points),]
  } else {
    pts <- n_points
    if(!is_spatial(pts)){
      stop("A spatVector or sf point object is required if `n_points` is not numeric")
    }

    n_points <- length(pts)


  }

  if(is.null(max_D)){
    max_D <- kernel_dist(kernel = kernel,
                         prob = 0.99,
                         sigma = max(sigma)) * 1.1
  } else {
    validate_scalar_numeric(max_D, "max_D", positive = TRUE)
  }

  kernel_out <- kernel_prep(pts = pts,
                            sigma = sigma,
                            shape = shape,
                            kernel = kernel,
                            max_D = max_D,
                            raster_stack = raster_stack,
                            projected = T,
                            verbose = F)

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

is_spatial <- function(x) {
  # Check for spatVector (terra package)
  spatvec_check <- inherits(x, "SpatVector")

  # Check for sf object (sf package)
  sf_check <- inherits(x, "sf") || inherits(x, "sfc")

  # Return TRUE if either check passes
  spatvec_check || sf_check
}

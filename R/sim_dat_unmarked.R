# Data sim function -------------------------------------------------------
#' @title Simulate data for optimizing scales of effect with `unmarked`
#'
#' @description Generates simulated replicated detection/non-detection or count
#' data formatted for use with the \code{unmarked} package, with known
#' kernel-weighted landscape covariates at a controlled scale of effect.
#'
#' @param alpha Numeric scalar. Intercept for the abundance/occupancy linear
#'   predictor. Default: \code{1}.
#' @param beta Numeric vector of slope coefficients, one per raster layer in
#'   \code{raster_stack}. Length must equal \code{nlyr(raster_stack)}.
#' @param kernel Character. Kernel function used to weight raster values by
#'   distance. One of \code{"gaussian"} (default), \code{"exp"} (negative
#'   exponential), \code{"expow"} (exponential power), or \code{"fixed"}
#'   (fixed-radius buffer).
#' @param type Character. Response type to simulate. One of:
#'   \describe{
#'     \item{\code{"count"}}{(Default) Poisson-distributed abundance with
#'       binomial thinning to simulate detection probability \code{det}.}
#'     \item{\code{"count_nb"}}{Negative binomial abundance with dispersion
#'       \code{StDev}, then binomial thinning.}
#'     \item{\code{"occ"}}{Bernoulli occupancy data; occupancy is generated
#'       from a logistic model and detections follow
#'       \code{Bernoulli(occ * det)}.}
#'   }
#' @param StDev Positive numeric. Dispersion (size) parameter for
#'   \code{"count_nb"}. Default: \code{0.5}.
#' @param n_points Positive integer. Number of spatial sample sites. Points are
#'   placed on a hexagonal grid and randomly subsampled. Default: \code{50}.
#' @param n_surv Positive integer. Number of repeated surveys per site, forming
#'   the columns of the returned observation matrix. Default: \code{3}.
#' @param det Numeric between 0 and 1. Per-survey detection probability applied
#'   via binomial thinning of the true abundance/occupancy. Default: \code{0.5}.
#' @param min_D Positive numeric. Minimum inter-point spacing on the hexagonal
#'   grid. If \code{NULL} (default), set automatically to
#'   \code{1.55 * max(sigma)}.
#' @param raster_stack A \code{SpatRaster} object. Layer names become covariate
#'   names in the returned data frame.
#' @param sigma Positive numeric vector. True kernel scale parameters, one per
#'   raster layer. Must be in the same units as the raster projection.
#' @param shape Positive numeric vector. Shape parameters for
#'   \code{kernel = "expow"}, one per raster layer. Required when using the
#'   exponential power kernel.
#' @param max_D Positive numeric. Maximum buffer radius for
#'   \code{\link{kernel_prep}} during simulation. If \code{NULL} (default),
#'   set automatically to 110\% of the distance enclosing 99\% of the kernel
#'   weight at \code{max(sigma)}.
#' @param user_seed Optional integer seed for reproducibility. Default:
#'   \code{NULL}.
#' @param ... Additional arguments. Not currently used.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{y}}{Integer matrix of dimensions \code{n_points x n_surv}
#'     containing the simulated replicated observations. Pass as the \code{y}
#'     argument when constructing an \code{unmarkedFrame}.}
#'   \item{\code{df}}{Data frame with \code{n_points} rows. The \code{y} column
#'     contains the simulated true abundance/occupancy (before detection
#'     thinning). Remaining columns are the true kernel-weighted covariate
#'     values, scaled to zero mean and unit variance.}
#'   \item{\code{pts}}{An \code{sf} POINT object with \code{n_points} rows.
#'     An \code{obs} column is appended containing the true
#'     abundance/occupancy values. Pass to \code{\link{kernel_prep}} as
#'     \code{pts}.}
#' }
#'
#' @details
#' Sites are distributed on a hexagonal grid across the interior of the raster
#' extent (85\% of the x/y range) and then randomly subsampled. The true
#' abundance or occupancy at each site is a function of the kernel-weighted
#' landscape covariates at the specified scale. Repeated surveys introduce
#' imperfect detection controlled by \code{det}.
#'
#' Use the returned \code{y} matrix and \code{df} (as \code{siteCovs}) to
#' construct an \code{unmarkedFrame}, fit a model, and then pass both the
#' fitted model and fresh \code{\link{kernel_prep}} output to
#' \code{\link{multiScale_optim}}.
#'
#' @examples
#' \donttest{
#' rs <- sim_rast(dim = 50, user_seed = 123)
#' rs <- terra::subset(rs, c(1, 3))
#' s_dat <- sim_dat_unmarked(alpha = 1,
#'                           beta = c(0.75, -0.75),
#'                           kernel = 'gaussian',
#'                           sigma = c(40, 80),
#'                           n_points = 20,
#'                           n_surv = 3,
#'                           det = 0.5,
#'                           type = 'count',
#'                           raster_stack = rs,
#'                           max_D = 220,
#'                           user_seed = 123)
#' plot(s_dat$df$y ~ s_dat$df$bin1)
#' plot(s_dat$df$y ~ s_dat$df$cont1)
#' ## unmarked analysis
#' kernel_inputs <- kernel_prep(pts = s_dat$pts,
#'                              raster_stack = rs,
#'                              max_D = 220,
#'                              kernel = 'gaus',
#'                              verbose = FALSE)
#'
#' umf <- unmarked::unmarkedFramePCount(y = s_dat$y,
#'                                      siteCovs = kernel_inputs$kernel_dat)
#'
#' ## Base unmarked model
#' mod0 <- unmarked::pcount(~1 ~ bin1 + cont1,
#'                          data = umf,
#'                          K = 30)
#'
#' ## The unmarked optimization is slower than the simulation and base model
#' ## fit, so it is skipped in automated example checks.
#' if (interactive()) {
#'   opt1 <- multiScale_optim(fitted_mod = mod0,
#'                            kernel_inputs = kernel_inputs)
#'
#'   summary(opt1)
#' }
#'}
#' @export
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
                             user_seed = NULL,
                             ...)
{

  kernel <- match.arg(kernel)
  type <- match.arg(type)
  validate_scalar_numeric(alpha, "alpha")
  validate_scalar_numeric(StDev, "StDev", positive = TRUE)
  validate_scalar_numeric(n_points, "n_points", integerish = TRUE, positive = TRUE)
  validate_scalar_numeric(n_surv, "n_surv", integerish = TRUE, positive = TRUE)
  validate_scalar_numeric(det, "det", lower = 0, upper = 1)

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
                            what = 'centers')
  }

  pts <- st_as_sf(pts)

  pts <- pts[sample(dim(pts)[1], n_points),]

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

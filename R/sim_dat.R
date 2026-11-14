# Data sim function -------------------------------------------------------
#' @title Simulate Responses at Known Scales of Effect
#'
#' @description Simulates survey points and a response that depends on
#' kernel-weighted landscape covariates at scales of effect you choose. Fit
#' the simulated data with [kernel_prep()] and [multiScale_optim()] to check
#' whether the true scales and coefficients are recovered, or to explore how
#' sample size and kernel choice affect estimation.
#'
#' @param alpha Numeric scalar. Intercept of the linear predictor, on the link
#'   scale of `type` (log for counts, logit for occurrence, identity for
#'   `"gaussian"`). Default: `1`.
#' @param beta Numeric vector of slopes, one per layer of `raster_stack` and
#'   in the same order, on the link scale. Each slope is the effect of a
#'   one standard deviation change in that standardized covariate. Required.
#' @param kernel Character. Kernel used to weight raster values by distance
#'   when computing the true covariates: `"gaussian"` (default), `"exp"`
#'   (negative exponential), `"expow"` (exponential power), or `"fixed"`
#'   (fixed-radius buffer). See [multiScaleR-kernels].
#' @param type Character. Distribution of the simulated response:
#'   * `"count"` (default): Poisson counts with mean `exp(eta)`.
#'   * `"count_nb"`: negative binomial counts with mean `exp(eta)` and size
#'     (dispersion) `StDev`.
#'   * `"occ"`: binary (0 or 1) occurrence with probability `plogis(eta)`.
#'   * `"gaussian"`: normal values with mean `eta` and standard deviation
#'     `StDev`.
#' @param StDev Positive numeric scalar. The `size` (dispersion) of the
#'   negative binomial for `"count_nb"`, where smaller values give more
#'   overdispersion, or the residual standard deviation for `"gaussian"`.
#'   Unused for `"count"` and `"occ"`. Default: `0.5`.
#' @param n_points Either a positive integer number of points to simulate, or
#'   a `SpatVector` of points to use directly. An integer places points on a
#'   regular square grid over the interior of the raster and samples
#'   `n_points` of them at random (see Details). Default: `50`.
#' @param min_D Positive numeric scalar. Starting spacing, in map units, of
#'   the point grid when `n_points` is an integer. Default `NULL` uses
#'   `1.55 * max(sigma)`. The spacing shrinks until the grid holds at least
#'   `n_points` points.
#' @param raster_stack A `SpatRaster` whose layers are the landscape
#'   variables, such as the output of [sim_rast()]. Layer names become the
#'   covariate names. Required.
#' @param sigma Positive numeric vector of true scale parameters, in map
#'   units, one per layer of `raster_stack`. These are the values that
#'   [multiScale_optim()] should recover. Required.
#' @param shape Positive numeric vector of true shape parameters, one per
#'   layer. Required when `kernel = "expow"` and ignored otherwise.
#' @param max_D Positive numeric scalar. Radius, in map units, of the buffer
#'   used to compute the true covariates with [kernel_prep()]. Default `NULL`
#'   uses 1.1 times the line-based (1D) distance holding 99% of the kernel
#'   weight at `max(sigma)`. That is shorter than the planar radius holding
#'   99%, so supply a larger value if you need less truncation. Supply it
#'   when `kernel = "expow"`, because the default cannot be computed for
#'   that kernel.
#' @param user_seed Optional integer seed. When supplied, `set.seed()` is
#'   called with it, which makes the simulation reproducible and also resets
#'   the global random number stream. Default: `NULL`.
#' @param ... Not used.
#'
#' @return A named list with three elements:
#' * `obs`: numeric vector with the simulated response, one value per point.
#' * `df`: data frame with the response in column `y` and one column per
#'   raster layer. The layer columns hold the true kernel-weighted means,
#'   standardized to mean 0 and standard deviation 1. A model fitted to this
#'   data frame estimates the coefficients at the true scales.
#' * `pts`: the survey points with an `obs` column holding the response.
#'   With an integer `n_points`, an `sf` POINT object with no coordinate
#'   reference system. Otherwise the supplied points. Pass it to
#'   [kernel_prep()] as `pts`.
#'
#' @details
#' **Model.** For each point, `sim_dat()` computes the kernel-weighted mean of
#' each layer at its true `sigma` (and `shape`) with [kernel_prep()], and
#' standardizes it across points. The linear predictor is
#' `eta = alpha + sum(beta * x)`, where `x` holds those standardized
#' covariates, and the response is drawn from the distribution set by `type`.
#'
#' **Point placement.** With an integer `n_points`, points start as the cell
#' centers of a regular square grid over a window inset from the raster
#' edges, which reduces edge effects. The inset is 15% of the raster's
#' maximum x coordinate on each side, and the same limits are used for x and
#' y. This assumes a square raster whose extent starts at 0, as produced by
#' [sim_rast()]. Rasters with other extents can fail or place points poorly,
#' so supply your own points through `n_points` for them. The grid spacing
#' starts at 97% of `min_D` and shrinks by 3% at a time until the grid holds
#' at least `n_points` points, and then `n_points` of them are sampled at
#' random.
#'
#' **Recovering the scales.** Prepare `pts` and the same `raster_stack` with
#' [kernel_prep()], using a `max_D` comfortably above the largest true
#' scale. Fit a starting model to the standardized covariates, run
#' [multiScale_optim()], and compare the estimated sigma values with `sigma`.
#' The example below shows these steps.
#'
#' @seealso [sim_rast()], [sim_dat_unmarked()], [kernel_prep()],
#'   [multiScale_optim()]. `vignette("v02-user-guide", package =
#'   "multiScaleR")` simulates zero-inflated counts with known scales.
#'
#' @examples
#' ## Four simulated layers. Keep one binary and one continuous layer
#' rs <- sim_rast(user_seed = 555)
#' rs <- terra::subset(rs, c("bin1", "cont1"))
#'
#' ## Counts with true scales of 75 m (bin1) and 150 m (cont1)
#' s_dat <- sim_dat(alpha = 0.5,
#'                  beta = c(0.75, -0.75),   # log-scale effect per SD
#'                  kernel = "gaussian",
#'                  sigma = c(75, 150),
#'                  type = "count",
#'                  n_points = 75,
#'                  raster_stack = rs,
#'                  max_D = 400,
#'                  user_seed = 999)
#' # Response and the true standardized covariates
#' head(s_dat$df)
#' # Counts rise with bin1 and fall with cont1
#' plot(y ~ bin1, data = s_dat$df)
#' plot(y ~ cont1, data = s_dat$df)
#'
#' \donttest{
#' ## Check whether multiScale_optim() recovers the true scales
#' # max_D well above the largest true scale
#' kernel_inputs <- kernel_prep(pts = s_dat$pts, raster_stack = rs,
#'                              max_D = 750, verbose = FALSE)
#' mod <- glm(y ~ bin1 + cont1, family = poisson,
#'            data = data.frame(y = s_dat$obs, kernel_inputs$kernel_dat))
#' fit <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs,
#'                         verbose = FALSE)
#' # Compare the sigma estimates with the true 75 and 150 m
#' summary(fit)
#' }
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
                         dimension = "1d",
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

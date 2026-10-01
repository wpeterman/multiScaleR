#' @title Scale Distance
#'
#' @description Converts a kernel scale parameter (\code{sigma}) into a distance
#' in map units: the distance within which a specified proportion of the kernel
#' weight falls. Distances can be measured along a line through the point
#' (one-dimensional) or as the radius of a circle around it (two-dimensional).
#' Can be used with a fitted \code{multiScaleR} object to report optimized scale
#' distances, or with manually supplied kernel parameters to explore kernel
#' behavior.
#'
#' @param model A fitted \code{multiScaleR} object from \code{\link{multiScale_optim}}.
#'   When provided, scale distances and 95\% confidence intervals are returned for
#'   all optimized covariates. When omitted, \code{sigma}, \code{kernel}, and
#'   (for \code{expow}) \code{beta} must be supplied via \code{...}.
#' @param prob Numeric between 0 and 1 (exclusive). Cumulative proportion
#'   used to define the distance. Default: \code{0.9}. The proportion refers
#'   to weight along a line for \code{dimension = "1d"} and weight on a plane
#'   for \code{dimension = "2d"}.
#' @param dimension Character. How the distance is measured. \code{"1d"}
#'   (default) returns the distance at which the kernel's cumulative density,
#'   read along a line through the point, reaches \code{prob}; this is the value
#'   \code{kernel_dist()} has always returned. \code{"2d"} returns the radius of
#'   the circle around the point that holds \code{prob} of the total weight the
#'   ideal kernel gives to locations on an unbounded plane. Because raster cells
#'   lie across a plane, the \code{"2d"} radius answers "within what
#'   distance does the ideal kernel assign 90\% of its total weight?". For
#'   distance-weighted kernels it is larger than the \code{"1d"} distance;
#'   hard-radius covariates report the same radius for both. See Details.
#' @param ... Additional parameters used when \code{model} is not supplied:
#'   \describe{
#'     \item{\code{sigma}}{Numeric (positive). The kernel scale parameter in the
#'       same units as the projection of \code{pts} and \code{raster_stack} passed
#'       to \code{\link{kernel_prep}}. For Gaussian kernels this is the standard
#'       deviation; for negative exponential kernels this is the decay rate.}
#'     \item{\code{kernel}}{Character. The kernel function to use. One of
#'       \code{"gaussian"}, \code{"exp"} (negative exponential),
#'       \code{"fixed"} (fixed-radius buffer), or \code{"expow"} (exponential
#'       power). Required when \code{model} is not provided.}
#'     \item{\code{beta}}{Numeric (positive). Shape parameter for the exponential
#'       power kernel. Required when \code{kernel = "expow"} and \code{model} is
#'       not provided. Ignored for all other kernels.}
#'   }
#'
#' @return
#' When \code{model} is provided: a data frame with one row per optimized
#' covariate and three columns, all in map units: \code{Mean} (distance at the
#' estimated sigma), \code{2.5\%} (distance at the lower 95\% limit of sigma),
#' and \code{97.5\%} (distance at the upper 95\% limit of sigma). Values are
#' rounded to two decimal places. For hard-radius covariates (landscape metrics
#' and unweighted surface metrics) the row reports the radius itself, whichever
#' \code{dimension} is requested.
#'
#' When kernel parameters are supplied directly via \code{...}: a single numeric
#' value giving the distance, in map units, at which the cumulative kernel
#' weight under the requested dimensional definition reaches \code{prob}.
#'
#' @details
#' \strong{One or two dimensions.} A kernel assigns each raster cell a weight
#' that declines with the cell's distance from the point. Two distances can
#' summarize how far that weight reaches:
#' \itemize{
#'   \item \code{dimension = "1d"} treats the kernel as a curve along a single
#'     line (a transect) through the point and returns the distance on either
#'     side that encloses \code{prob} of the area under that curve. For a
#'     Gaussian kernel the 90\% distance is about 1.65 sigma.
#'   \item \code{dimension = "2d"} accounts for locations lying on a plane, so
#'     the number of locations at distance \emph{d} grows with the circumference
#'     2 pi \emph{d}. It returns the radius of the circle that holds \code{prob}
#'     of the ideal kernel weight. For a Gaussian kernel the 90\% radius is about
#'     2.15 sigma.
#' }
#' With sigma = 300 m, a Gaussian kernel gives a 90\% distance of about 493 m in
#' one dimension and 644 m in two. On an unbounded plane, locations within
#' 493 m carry only about 74\% of the total weight. When you describe the area
#' associated with kernel weighting, or compare with values defined on the plane
#' (for example, a radius holding 90\% of the kernel weight), use
#' \code{"2d"}. Always state which definition you report.
#'
#' \strong{Kernel formulas.} In one dimension:
#' \itemize{
#'   \item \strong{Gaussian}: the inverse normal CDF, \code{qnorm((1 + prob) / 2) * sigma}.
#'   \item \strong{Negative exponential}: \code{-sigma * log(1 - prob)}.
#'   \item \strong{Fixed buffer}: \code{sigma * prob} (the fraction of
#'     the buffer radius).
#'   \item \strong{Exponential power}: integrates the density numerically; both
#'     \code{sigma} and \code{beta} (shape) must be specified.
#' }
#' In two dimensions:
#' \itemize{
#'   \item \strong{Gaussian}: \code{sigma * sqrt(-2 * log(1 - prob))}.
#'   \item \strong{Negative exponential}: \code{sigma * qgamma(prob, shape = 2)},
#'     about 3.89 sigma at 90\%.
#'   \item \strong{Fixed buffer}: \code{sigma * sqrt(prob)}, the radius of the
#'     circle holding \code{prob} of the buffer's area.
#'   \item \strong{Exponential power}:
#'     \code{sigma * qgamma(prob, shape = 2 / beta)^(1 / beta)}.
#' }
#' Both definitions describe the ideal continuous kernel. They do not account
#' for truncation at \code{max_D}, raster boundaries, missing cells, or cell
#' geometry. The realized share of weight among extracted raster cells can
#' therefore differ. If the two-dimensional radius approaches or exceeds
#' \code{max_D}, consider increasing \code{max_D} in
#' \code{\link{kernel_prep}} and refitting the model.
#'
#' Confidence intervals for the fitted \code{model} case are derived from the
#' Hessian-based standard errors (or profile-likelihood intervals when
#' \code{\link{profile_sigma}} has been run and stored on the object). Each
#' limit of sigma is converted to a distance with the same formula.
#'
#' \code{summary()} of a fitted model reports both distances; see
#' \code{\link{summary.multiScaleR}}.
#' @seealso \code{\link{summary.multiScaleR}}, \code{\link{plot.multiScaleR}},
#'   \code{\link{plot_kernel}}
#' @examples
#' \donttest{
#' ## Using package data
#' data('pts')
#' data('count_data')
#' hab <- terra::rast(system.file('extdata',
#'                    'hab.tif', package = 'multiScaleR'))
#'
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = hab,
#'                              max_D = 250,
#'                              kernel = 'gaussian')
#'
#' mod <- glm(y ~ hab,
#'            family = poisson,
#'            data = count_data)
#'
#' ## Optimize scale
#' opt <- multiScale_optim(fitted_mod = mod,
#'                         kernel_inputs = kernel_inputs)
#'
#' ## Uses of `kernel_dist`
#' kernel_dist(model = opt)                     # along a line (1D)
#' kernel_dist(model = opt, dimension = "2d")   # radius of a circle (2D)
#' kernel_dist(model = opt, prob = 0.95)
#' kernel_dist(sigma = 500, kernel = 'gaussian', prob = 0.95)
#' kernel_dist(sigma = 100, prob = 0.975, kernel = "exp")
#' kernel_dist(sigma = 100, prob = 0.95, kernel = "expow", beta = 1.5)
#' kernel_dist(sigma = 100, kernel = "fixed")
#' }
#'
#' ## The two definitions for a Gaussian kernel with sigma = 300
#' kernel_dist(sigma = 300, kernel = "gaussian")                     # about 493
#' kernel_dist(sigma = 300, kernel = "gaussian", dimension = "2d")   # about 644
#'
#' @rdname kernel_dist
#' @export
#' @importFrom insight get_df

kernel_dist <- function(model,
                        prob = 0.9,
                        dimension = c("1d", "2d"),
                        ...){
  param_list <- list(...)
  dimension <- match.arg(dimension)
  validate_scalar_numeric(prob,
                          "prob",
                          lower = 0,
                          upper = 1,
                          inclusive_lower = FALSE,
                          inclusive_upper = FALSE)

  if(!missing("model")){
    if(!inherits(model, "multiScaleR")){
      stop("Provide a fitted `multiScaleR` model object")
    }
  }

  if(!missing("model")){
    if(length(param_list) >= 1){
      warning("Calculating fitted scale relationship; Ignoring specified `sigma` and/or `shape` parameters")
    }

    if(!missing("model")){
      # ci_ <- summary(model)$opt_scale

      opt_mod <- .analysis_model(model$opt_mod)

      if(any(class(opt_mod) == 'gls')){
        df <- opt_mod$dims$N - opt_mod$dims$p
        names <- all.vars(formula(opt_mod)[-2])

      } else if(any(grepl("^unmarked", class(opt_mod)))){
        df <- dim(opt_mod@data@y)[1]
        names <- .unmarked_model_predictors(opt_mod)

      } else {
        df <- get_df(opt_mod, type = "residual")
        names <- all.vars(formula(opt_mod)[-2])
      }

      if(!is.null(model$profile_scale_est)){
        ci_ <- model$profile_scale_est
      } else {
        ci_ <- ci_func(model$scale_est,
                       df = df,
                       min_D = model$min_D,
                       names = row.names(model$scale_est))
      }

      # browser()

      dist_list <- vector('list', nrow(ci_))
      hard_radius <- .msr_hard_radius_covariates(model$kernel_inputs$scale_vars)
      for(i in 1:nrow(ci_)){
        is_landscape <- rownames(ci_)[[i]] %in% hard_radius
        if (isTRUE(is_landscape) && is.finite(ci_[i, 1])) {
            scale_mn <- ci_[i, 1]
            scale_l <- ci_[i, 3]
            scale_u <- ci_[i, 4]
        } else if (is.finite(ci_[i, 2])) {
            shape_i <- if (!is.null(model$shape_est)) model$shape_est[i,1] else NULL
          # wt_mn <- scale_type_r(d = d,
          #                       kernel = model$kernel_inputs$kernel,
          #                       sigma = ci_[i, 1],
          #                       shape = model$shape_est[i,1],
          #                       output = 'wts')
          #
          #
          # wt_l <- scale_type_r(d = d,
          #                      kernel = model$kernel_inputs$kernel,
          #                      sigma = ci_[i,3],
          #                      shape = model$shape_est[i,1],
          #                      output = 'wts')
          #
          # wt_u <- scale_type_r(d = d,
          #                      kernel = model$kernel_inputs$kernel,
          #                      sigma = ci_[i,4],
          #                      shape = model$shape_est[i,1],
          #                      output = 'wts')
          #
          # scale_mn <- wtd.Ecdf(d, weights = wt_mn)
          # scale_mn <- round(scale_mn$x[which(scale_mn$ecdf > prob)[1]], digits = 2)
          #
          # scale_l <- wtd.Ecdf(d, weights = wt_l)
          # scale_l <- round(scale_l$x[which(scale_l$ecdf > prob)[1]], digits = 2)
          #
          # scale_u <- wtd.Ecdf(d, weights = wt_u)
          # scale_u <- round(scale_u$x[which(scale_u$ecdf > prob)[1]], digits = 2)

          scale_mn <- k_dist(sigma = ci_[i, 1],
                             prob = prob,
                             kernel = model$kernel_inputs$kernel,
                             beta = shape_i,
                             dimension = dimension)
          scale_l <- k_dist(sigma = ci_[i, 3],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = shape_i,
                            dimension = dimension)
          scale_u <- k_dist(sigma = ci_[i, 4],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = shape_i,
                            dimension = dimension)
        } else {
          scale_mn <- NaN
          scale_l <- NaN
          scale_u <- NaN
        }


        dist_list[[i]] <- data.frame(mn = round(scale_mn, digits = 2),
                                     l = round(scale_l, digits = 2),
                                     u = round(scale_u, digits = 2))
      }
      dist_out <- do.call(rbind, dist_list)
      rownames(dist_out) <- rownames(ci_)
      colnames(dist_out) <- colnames(ci_)[c(1,3,4)]
    }
    return(dist_out)
  } else if(length(param_list) >= 1){
    sig_ <- param_list$sigma
    shp_ <- param_list$beta
    kern <- param_list$kernel

    if(is.null(sig_)){
      stop('\nA value for `sigma` must be provided!\n')
    }
    validate_scalar_numeric(sig_, "sigma", positive = TRUE)
    if(is.null(kern)){
      stop('\nYou must specify `kernel` function; See Details\n')
    }
    kern <- match.arg(kern, c("gaussian", "exp", "expow", "fixed"))
    if(kern == 'expow' & is.null(shp_)){
      stop('\nBoth a `sigma` and `shape` parameter must be specified when using the `expow` kernel; See Details\n')
    }

    if (kern == "expow") {
      validate_scalar_numeric(shp_, "beta", positive = TRUE)
    }

    # d <- seq(1, round(sig_*1000,0), length.out = round(sig_*1000,0))
    # wt <- scale_type_r(d = d,
    #                    kernel = kern,
    #                    sigma = sig_,
    #                    shape = shp_,
    #                    output = 'wts')
    #
    # mx <- wtd.Ecdf(d, weights = wt)
    # mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)
    #
    # d <- seq(1, mx, length.out = 100)
    # wt <- scale_type_r(d = d,
    #                    kernel = kern,
    #                    sigma = sig_,
    #                    shape = shp_,
    #                    output = 'wts')
    #
    # scale_d <- round(d[which(wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], 2)

    scale_d <- k_dist(sigma = sig_,
                      prob = prob,
                      kernel = kern,
                      beta = shp_,
                      dimension = dimension)

    return(round(scale_d, digits = 2))
  } else {
    stop("Parameters not correctly specified to calculate distance. See Details and try again.")
  }
}


#' @title Distance at Cumulative Kernel Proportion
#' @description Compute the distance at which a given cumulative proportion of
#'   kernel weight is reached for several kernel types, measured either along a
#'   line through the point (one-dimensional) or as the radius of a circle
#'   around it (two-dimensional).
#' @param sigma Numeric. Scale parameter. For Gaussian and exponential, this is standard deviation or decay rate. For expow, this is the kernel bandwidth.
#' @param prob Numeric. Desired cumulative proportion (e.g., 0.95).
#' @param kernel Character. One of "gaussian", "exp", "expow", or "fixed".
#' @param beta Numeric. Shape parameter for exponential power kernel. Ignored unless kernel = "expow".
#' @param dimension Character. \code{"1d"} (default) for the distance along a
#'   line through the point; \code{"2d"} for the radius of the circle holding
#'   \code{prob} of the total weight when the kernel weights cells on a plane.
#' @return Numeric distance in map units at which the cumulative kernel weight
#'   reaches \code{prob} under the selected dimensional definition.
#' @keywords internal
#' @importFrom stats integrate qnorm qgamma uniroot
k_dist <- function(sigma, prob = 0.95, kernel = c("gaussian", "exp", "expow", "fixed"),
                   beta = NULL, dimension = c("1d", "2d")) {
  kernel <- match.arg(kernel)
  dimension <- match.arg(dimension)
  if (prob <= 0 || prob >= 1) stop("prob must be between 0 and 1")
  if (!is.finite(sigma)) return(sigma)
  if (kernel == "expow") {
    if (is.null(beta)) stop("beta must be specified for exponential power kernel")
    if (beta <= 0) stop("beta must be positive")
  }

  if (dimension == "2d") {
    # Planar radius holding `prob` of the total weight. Integrating the weight
    # w(d) over rings of circumference 2 * pi * d gives closed forms:
    #   Gaussian, w = exp(-d^2 / (2 sigma^2)): 1 - exp(-r^2 / (2 sigma^2))
    #   Exponential power, w = exp(-(d / sigma)^beta): the cumulative weight is
    #     a gamma CDF with shape 2 / beta evaluated at (r / sigma)^beta
    #   Negative exponential: exponential power with beta = 1
    #   Fixed buffer, w = 1 for d < sigma: area fraction (r / sigma)^2
    # Each assumes the kernel is not truncated by `max_D`.
    return(switch(kernel,
                  gaussian = sigma * sqrt(-2 * log(1 - prob)),
                  exp      = sigma * qgamma(prob, shape = 2),
                  expow    = sigma * qgamma(prob, shape = 2 / beta)^(1 / beta),
                  fixed    = sigma * sqrt(prob)))
  }

  if (kernel == "gaussian") {
    return(qnorm((1 + prob)/2, mean = 0, sd = sigma))

  } else if (kernel == "exp") {
    return(-sigma * log(1 - prob))

  } else if (kernel == "expow") {
    c_beta <- beta / (2 * sigma * gamma(1 / beta))
    f <- function(x) c_beta * exp(-abs(x / sigma)^beta)
    cdf <- function(x) integrate(f, -x, x, rel.tol = 1e-8)$value
    target_fn <- function(d) cdf(d) - prob
    return(uniroot(target_fn, c(1e-6, 10 * sigma))$root)

  } else if (kernel == "fixed") {
    # Proportion of mass in 1D uniform kernel: step function
    return(sigma * (prob))  # Assume sigma is the full extent
  }
}

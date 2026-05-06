#' @title Scale Distance
#'
#' @description Estimates the distance at which a specified cumulative proportion
#' of the kernel density function is reached. Can be used with a fitted
#' \code{multiScaleR} object to report optimized scale distances, or with
#' manually supplied kernel parameters to explore kernel behavior.
#'
#' @param model A fitted \code{multiScaleR} object from \code{\link{multiScale_optim}}.
#'   When provided, scale distances and 95\% confidence intervals are returned for
#'   all optimized covariates. When omitted, \code{sigma}, \code{kernel}, and
#'   (for \code{expow}) \code{beta} must be supplied via \code{...}.
#' @param prob Numeric between 0 and 1 (exclusive). Cumulative kernel density
#'   threshold used to define the effective distance. Default: \code{0.9}, meaning
#'   the distance enclosing 90\% of the kernel weight.
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
#' covariate and three columns — \code{Mean} (distance at the estimated sigma),
#' \code{Lower} (distance at the lower 95\% CI of sigma), and \code{Upper}
#' (distance at the upper 95\% CI of sigma). Values are rounded to two decimal
#' places.
#'
#' When kernel parameters are supplied directly via \code{...}: a single numeric
#' value giving the distance at which the cumulative kernel density first reaches
#' \code{prob}.
#'
#' @details
#' The effective distance depends on both the kernel type and the scale parameter
#' \code{sigma}:
#' \itemize{
#'   \item \strong{Gaussian}: uses the inverse normal CDF, so the 90\% distance
#'     is approximately 1.65 sigma.
#'   \item \strong{Negative exponential}: uses \code{-sigma * log(1 - prob)}.
#'   \item \strong{Fixed buffer}: returns \code{sigma * prob} (the fraction of
#'     the buffer radius).
#'   \item \strong{Exponential power}: integrates the density numerically; both
#'     \code{sigma} and \code{beta} (shape) must be specified.
#' }
#'
#' Confidence intervals for the fitted \code{model} case are derived from the
#' Hessian-based standard errors (or profile-likelihood intervals when
#' \code{\link{profile_sigma}} has been run and stored on the object).
#' @seealso \code{\link{plot.multiScaleR}}
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
#' kernel_dist(model = opt)
#' kernel_dist(model = opt, prob = 0.95)
#' kernel_dist(sigma = 500, kernel = 'gaussian', prob = 0.95)
#' kernel_dist(sigma = 100, prob = 0.975, kernel = "exp")
#' kernel_dist(sigma = 100, prob = 0.95, kernel = "expow", beta = 1.5)
#' kernel_dist(sigma = 100, kernel = "fixed")
#' }
#'
#' @rdname kernel_dist
#' @export
#' @importFrom insight get_df

kernel_dist <- function(model,
                        prob = 0.9,
                        ...){
  param_list <- list(...)
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
        names <- all.vars(opt_mod@formula)

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
      scale_var_types <- NULL
      if (!is.null(model$kernel_inputs$scale_vars)) {
        scale_var_types <- stats::setNames(
          model$kernel_inputs$scale_vars$type,
          model$kernel_inputs$scale_vars$covariate
        )
      }


      for(i in 1:nrow(ci_)){
        if(!is.nan(ci_[i,2])){
          is_landscape <- !is.null(scale_var_types) &&
            identical(scale_var_types[[rownames(ci_)[[i]]]], "landscape")

          if (isTRUE(is_landscape)) {
            scale_mn <- ci_[i, 1]
            scale_l <- ci_[i, 3]
            scale_u <- ci_[i, 4]
          } else {
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
                             beta = shape_i)
          scale_l <- k_dist(sigma = ci_[i, 3],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = shape_i)
          scale_u <- k_dist(sigma = ci_[i, 4],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = shape_i)
          }
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
                      beta = shp_)

    return(round(scale_d, digits = 2))
  } else {
    stop("Parameters not correctly specified to calculate distance. See Details and try again.")
  }
}


#' @title Distance at Cumulative Kernel Proportion
#' @description Compute the distance at which a given cumulative density is reached for several kernel types.
#' @param sigma Numeric. Scale parameter. For Gaussian and exponential, this is standard deviation or decay rate. For expow, this is the kernel bandwidth.
#' @param prob Numeric. Desired cumulative proportion (e.g., 0.95).
#' @param kernel Character. One of "gaussian", "exp", "expow", or "fixed".
#' @param beta Numeric. Shape parameter for exponential power kernel. Ignored unless kernel = "expow".
#' @return Numeric. Distance at which the cumulative kernel density reaches the specified proportion.
#' @keywords internal
#' @importFrom stats integrate qnorm uniroot
k_dist <- function(sigma, prob = 0.95, kernel = c("gaussian", "exp", "expow", "fixed"), beta = NULL) {
  kernel <- match.arg(kernel)
  if (prob <= 0 || prob >= 1) stop("prob must be between 0 and 1")
  if (!is.finite(sigma)) return(sigma)

  if (kernel == "gaussian") {
    return(qnorm((1 + prob)/2, mean = 0, sd = sigma))

  } else if (kernel == "exp") {
    return(-sigma * log(1 - prob))

  } else if (kernel == "expow") {
    if (is.null(beta)) stop("beta must be specified for exponential power kernel")
    if (beta <= 0) stop("beta must be positive")

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



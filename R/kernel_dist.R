#' @title Scale Distance
#' @description Function to estimate the effective distance encompassing a specified cumulative probability density of the kernel function
#' @param model \code{\link{multiScale_optim}} object of class 'multiScaleR'
#' @param prob Density probability cutoff for calculating distance, Default: 0.9
#' @param ... Parameters to be used if not providing a 'multiScaleR' fitted object. See Details
#' @return Numeric. Distance at which the cumulative kernel density reaches the specified proportion.
#' @details This function is used to determine the distance at which kernel density distributions have influence. If not providing a fitted model, you can plot kernel distributions by specifying (1) sigma, (2) beta (if using exponential power), and (3) the kernel transformation ('exp' = negative exponential, 'gaussian', 'fixed' = fixed buffer, and 'expow' = exponential power)
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

  if(!missing("model")){
    if(!inherits(model, "multiScaleR")){
      stop("Provide a fitted `multiScaleR` model object")
    }
  }

  if((!is.numeric(prob) | prob < 0 | prob > 1)){
    stop("`prob` must be a decimal between 0 and 1")
  }

  if(!missing("model")){
    if(length(param_list) >= 1){
      warning("Calculating fitted scale relationship; Ignoring specified `sigma` and/or `shape` parameters")
    }

    if(!missing("model")){
      # ci_ <- summary(model)$opt_scale

      if(any(class(model$opt_mod) == 'gls')){
        df <- model$opt_mod$dims$N - model$opt_mod$dims$p
        names <- all.vars(formula(model$opt_mod)[-2])

      } else if(any(grepl("^unmarked", class(model$opt_mod)))){
        df <- dim(model$opt_mod@data@y)[1]
        names <- all.vars(model$opt_mod@formula)

      } else {
        df <- get_df(model$opt_mod, type = "residual")
        names <- all.vars(formula(model$opt_mod)[-2])
      }

      ci_ <- ci_func(model$scale_est,
                     df = df,
                     min_D = model$min_D,
                     names = row.names(model$scale_est))

      # browser()

      d <- seq(1, round(max(ci_, na.rm = T)*1000,0), length.out = round(max(ci_, na.rm = T)*1000,0))
      dist_list <- vector('list', nrow(ci_))


      for(i in 1:nrow(ci_)){
        if(!is.nan(ci_[i,2])){
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
                             beta = model$shape_est[i,1])
          scale_l <- k_dist(sigma = ci_[i, 3],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = model$shape_est[i,1])
          scale_u <- k_dist(sigma = ci_[i, 4],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = model$shape_est[i,1])
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
    if(is.null(kern)){
      stop('\nYou must specify `kernel` function; See Details\n')
    }
    if(kern == 'expow' & is.null(shp_)){
      stop('\nBoth a `sigma` and `shape` parameter must be specified when using the `expow` kernel; See Details\n')
    }

    if (kern == "expow") {
      if (is.null(shp_)) stop("beta must be specified for exponential power kernel")
      if (shp_ <= 0) stop("beta must be positive")
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



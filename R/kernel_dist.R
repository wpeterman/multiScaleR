#' @title Scale Distance
#' @description Function to estimate the effective distance encompassing a specified probability density of the kernel density function
#' @param x \code{\link[multiScaleR]{multiScaleR}} object
#' @param prob Density probability cutoff for calculating distance, Default: 0.9
#' @param ... Parameters to be used if not providing a 'multiScaleR' fitted object. See Details
#' @return Distance at which cumulative density of kernel is reached
#' @details This function is used to determine the distance at which kernel density distributions have influence. If not providing a fitted model, you can plot kernel distributions by specifying (1) sigma, (2) shape (if using exponential power), and (3) the kernel transformation ('exp' = negative exponential, 'gaussian', 'fixed' = fixed buffer, and 'expow' = exponential power)
#' @seealso \code{\link[multiScaleR]{plot.multiScaleR}}
#' @examples
#' kernel_dist(x)
#'
#' @rdname kernel_dist
#' @export
kernel_dist <- function(model,
                        prob = 0.9,
                        ...){

  param_list <- list(...)

  if(!missing("model")){
    if(class(model) != 'multiScaleR'){
      stop("Provide a fitted `multiScaleR` model object")
    }
  }


  if((!is.numeric(prob) | prob < 0 | prob > 1)){
    stop("`prob` must be a decimal between 0–1")
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
        df <- insight::get_df(model$opt_mod, type = "residual")
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
          wt_mn <- scale_type(d = d,
                              kernel = model$kernel_inputs$kernel,
                              sigma = ci_[i, 1],
                              shape = model$shape_est[i,1],
                              output = 'wts')


          wt_l <- scale_type(d = d,
                             kernel = model$kernel_inputs$kernel,
                             sigma = ci_[i,3],
                             shape = model$shape_est[i,1],
                             output = 'wts')

          wt_u <- scale_type(d = d,
                             kernel = model$kernel_inputs$kernel,
                             sigma = ci_[i,4],
                             shape = model$shape_est[i,1],
                             output = 'wts')

          scale_mn <- Hmisc::wtd.Ecdf(d, weights = wt_mn)
          scale_mn <- round(scale_mn$x[which(scale_mn$ecdf > prob)[1]], digits = 2)

          scale_l <- Hmisc::wtd.Ecdf(d, weights = wt_l)
          scale_l <- round(scale_l$x[which(scale_l$ecdf > prob)[1]], digits = 2)

          scale_u <- Hmisc::wtd.Ecdf(d, weights = wt_u)
          scale_u <- round(scale_u$x[which(scale_u$ecdf > prob)[1]], digits = 2)
        } else {
          scale_mn <- NaN
          scale_l <- NaN
          scale_u <- NaN
        }


        dist_list[[i]] <- data.frame(mn = scale_mn,
                                     l = scale_l,
                                     u = scale_u)
      }
      dist_out <- do.call(rbind, dist_list)
      rownames(dist_out) <- rownames(ci_)
      colnames(dist_out) <- colnames(ci_)[c(1,3,4)]
    }
    return(dist_out)
  } else if(length(param_list) >= 1){
    sig_ <- param_list$sigma
    shp_ <- param_list$shape
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

    d <- seq(1, round(sig_*1000,0), length.out = round(sig_*1000,0))
    wt <- scale_type(d = d,
                     kernel = kern,
                     sigma = sig_,
                     shape = shp_,
                     output = 'wts')

    mx <- Hmisc::wtd.Ecdf(d, weights = wt)
    mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)

    d <- seq(1, mx, length.out = 100)
    wt <- scale_type(d = d,
                     kernel = kern,
                     sigma = sig_,
                     shape = shp_,
                     output = 'wts')

    scale_d <- round(d[which(Hmisc::wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], 2)
    return(scale_d)
  } else {
    stop("Parameters not correctly specified to calculate distance. See Details and try again.")
  }
}

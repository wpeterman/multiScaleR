#' @title Plot mulitScaleR
#' @description Function to plot fitted multiScaleR objects
#' @param x \code{\link[multiScaleR]{multiScaleR}} object
#' @param prob Density probability cutoff for plotting, Default: 0.9
#' @param scale_dist Logical. If TRUE (Default), The distance at which the specified density probability is achieved
#' @param ... Not used
#' @return Plots of kernel density distributions
#' @details This function is used to visualize kernel density distributions from multiScaleR optimized objects.
#' @examples
#' plot(x)
#'
#' @rdname plot.multiScaleR
#' @export
#' @importFrom Hmisc wtd.Ecdf
#' @importFrom cowplot theme_cowplot
#' @importFrom DescTools Quantile
#' @importFrom ggplot2 aes annotate geom_line geom_vline ggplot xlab ylab theme_light
plot.multiScaleR <- function(x,
                             prob = 0.9,
                             scale_dist = TRUE,
                             ...) {
  if(isTRUE(scale_dist) & (!is.numeric(prob) | prob < 0 | prob > 1)){
    stop("`prob` must be a number between 0–1")
  }

  param_list <- list(...)

  if(!missing("x")){
    if(length(param_list) >= 1){
      warning("Plotting fitted scale relationship; Ignoring specified `sigma` and/or `shape` parameters")
    }
    sig_ <- summary(x)$opt_scale

    if(x$kernel_inputs$kernel == 'expow'){
      shp_ <- summary(x)$opt_shape
    } else {
      shp_ <- NULL
    }

    df_list <- plot_list <- vector('list', nrow(sig_))

    for(i in 1:nrow(sig_)){
      d <- seq(1, 1e6,
               length.out = 1e6)
      wt <- scale_type(d = d,
                       kernel = x$kernel_inputs$kernel,
                       sigma = x$scale_est[[1]][i],
                       shape = x$shape_est[[1]][i],
                       output = 'wts')

      mx <- Hmisc::wtd.Ecdf(d, weights = wt)
      mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)

      # mx <- round(DescTools::Quantile(d,
      #                                 weights = wt,
      #                                 probs = 0.99), digits = -1)

      d <- seq(1, mx, length.out = 100)
      wt <- scale_type(d = d,
                       kernel = x$kernel_inputs$kernel,
                       sigma = x$scale_est[[1]][i],
                       shape = x$shape_est[[1]][i],
                       output = 'wts')

      scale_d <- round(d[which(Hmisc::wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], -1)

      df_list[[i]] <- data.frame(dist = d,
                                 wt = wt)

      if(isTRUE(scale_dist)){
        if(min(prob) >= 0.8){
          ax <- max(d) * 0.075
          ay <- 0.075*max(wt)
        } else {
          ax <- max(d) * 0.8
          ay <- 0.9*max(wt)
        }
      }

      plot_ <- ggplot(data = df_list[[i]], aes(x = dist, y = wt)) +
        geom_line(linewidth = 1.25) +
        {if(isTRUE(scale_dist))
          geom_vline(xintercept = scale_d,
                     linetype = 'dashed',
                     color = 'red')
        } +
        {if(isTRUE(scale_dist))
          annotate('text', x = ax, y = ay,
                   label = paste0(prob*100,"% density \n Distance: ", round(scale_d, 0)))
        } +
        xlab('Distance') +
        ylab('Smoothing weight') +
        cowplot::theme_cowplot()

      plot_list[[i]] <- plot_
    }
    return(plot_list)

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

    d <- seq(1, 1e6, length.out = 1e6)
    wt <- scale_type(d = d,
                     kernel = kern,
                     sigma = sig_,
                     shape = shp_,
                     output = 'wts')

    mx <- round(DescTools::Quantile(d,
                                    weights = wt,
                                    probs = 0.999) * 1.1, digits = -1)

    d <- seq(1, mx, length.out = 100)
    wt <- scale_type(d = d,
                     kernel = kern,
                     sigma = sig_,
                     shape = shp_,
                     output = 'wts')

    scale_d <- DescTools::Quantile(d,
                                   weights = wt,
                                   probs = prob)

    df <- data.frame(dist = d,
                     wt = wt)

    if(isTRUE(scale_dist)){
      if(min(prob) >= 0.8){
        ax <- max(d) * 0.075
        ay <- 0.075*max(wt)
      } else {
        ax <- max(d) * 0.8
        ay <- 0.9*max(wt)
      }
    }
    ggplot(data = df, aes(x = d, y = wt)) +
      geom_line(size = 1.25) +
      {if(isTRUE(scale_dist))
        geom_vline(xintercept = scale_d,
                   linetype = 'dashed',
                   color = 'red')
      } +
      {if(isTRUE(scale_dist))
        annotate('text', x = ax, y = ay,
                 label = paste0(prob*100,"% density \n Distance: ", round(scale_d, 0)))
      } +
      xlab('Distance') +
      ylab('Smoothing weight') +
      cowplot::theme_cowplot()
  } else {
    stop("Parameters not correctly specified to create plot. See Details and try again.")
  }
}

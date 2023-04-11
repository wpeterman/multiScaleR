#' @title Plot kernel densities
#' @description Generic function to plot kernels
#' @param prob Cumulative kernel density to identify scale of effect distance, Default: 0.9
#' @param sigma Value of scaling parameter, sigma
#' @param shape If plotting an Exponential Power kernel. Values between 1-50 are typically valid. (Default = NULL)
#' @param kernel Kernel function to use. Valid functions are c('exp', 'gaussian', fixed', 'expow'). See details
#' @param scale_dist Logical. If TRUE (Default), the distance at which the specified density probability is achieved is added to the plot along with 95\% confidence interval
#' @param add_label Logical. If TRUE (Default), the distance value calculated for 'scale_dist' is added as an annotation to the plot.
#' @param ... Not used
#' @returns ggplot2 objects of kernel density distributions
#' @details
#' This function is used to visualize kernel density distributions without having a fitted multiScaleR optimized object. Requires (1) sigma, (2) shape (if using exponential power), and (3) the kernel transformation ('exp' = negative exponential, 'gaussian', 'fixed' = fixed buffer, and 'expow' = exponential power)
#'
#' @examples
#'#' ## General use of plot method
#' plot_kernel(prob = 0.95,
#'             sigma = 100,
#'             kernel = 'gaussian')
#' plot_kernel(prob = 0.95,
#'             sigma = 100,
#'             kernel = 'exp')
#' plot_kernel(prob = 0.95,
#'             sigma = 100,
#'             kernel = 'fixed')
#' plot_kernel(prob = 0.95,
#'             sigma = 100,
#'             shape = 2.5,
#'             kernel = 'negexp')
#' @rdname plot_kernel
#' @export
#' @importFrom Hmisc wtd.Ecdf
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes annotate geom_line geom_vline ggplot xlab ylab theme_light geom_rect ggtitle
plot_kernel <- function(prob = 0.9,
                        sigma,
                        shape = NULL,
                        kernel,
                        scale_dist = TRUE,
                        add_label = TRUE,
                        ...) {
  if(isTRUE(scale_dist) & (!is.numeric(prob) | prob < 0 | prob > 1)){
    stop("`prob` must be a number between 0–1")
  }

  sig_ <- sigma
  shp_ <- shape
  kern <- kernel

  if(is.null(sig_)){
    stop('\nA value for `sigma` must be provided!\n')
  }
  if(is.null(kern)){
    stop('\nYou must specify `kernel` function; See Details\n')
  }
  if(kern == 'expow' & is.null(shp_)){
    stop('\nBoth a `sigma` and `shape` parameter must be specified when using the `expow` kernel; See Details\n')
  }

  # d <- seq(1, 1e6, length.out = 1e6)
  d <- seq(1, round(max(sig_)*1000,0),
           length.out = round(max(sig_)*1000,0))
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

  scale_d <- round(d[which(Hmisc::wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], -1)

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
    geom_line(linewidth = 1.25) +
    {if(isTRUE(scale_dist))
      geom_vline(xintercept = scale_d,
                 linetype = 'dashed',
                 color = 'red')
    } +
    {if(isTRUE(scale_dist) & isTRUE(add_label))
      annotate('text', x = ax, y = ay,
               label = paste0(prob*100,"% density \n Distance: ", round(scale_d, 0)))
    } +
    xlab('Distance') +
    ylab('Weight') +
    cowplot::theme_cowplot()
}


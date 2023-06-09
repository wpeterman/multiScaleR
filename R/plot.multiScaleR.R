#' @title Plot multiScaleR
#' @description Function to plot fitted multiScaleR objects
#' @param x A `multiScaleR` object from  \code{\link[multiScaleR]{multiScale_optim}} object
#' @param ... Parameters to modify plot. See Details
#' @returns ggplot2 objects of kernel density distributions
#' @details
#' This function is used to visualize kernel density distributions from multiScaleR optimized objects. Plots can be modified by specifying:
#' \tabular{ll}{
#'  \tab * `prob` --> Cumulative kernel density to identify scale of effect distance, Default: 0.9 \cr
#'  \tab * `scale_dist` --> Logical. If TRUE (Default), the distance at which the specified density probability is achieved is added to the plot along with 95\% confidence interval \cr
#'  \tab * `add_label` --> Logical. If TRUE (Default), the distance value calculated for 'scale_dist' is added as an annotation to the plot. \cr
#'  }
#' @seealso \code{\link[multiScaleR]{plot_kernel}}
#' @examples
#' plot(x)
#'
#' plot(x, prob = 0.95)
#'
#' plot(x, scale_dist = FALSE)
#'
#' plot(x, scale_dist = TRUE, add_label = FALSE)
#' @rdname plot.multiScaleR
#' @export
#' @method plot multiScaleR
#' @importFrom Hmisc wtd.Ecdf
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes annotate geom_line geom_vline ggplot xlab ylab theme_light geom_rect ggtitle
#'
plot.multiScaleR <- function(x,
                             ...) {
  param_list <- list(...)

  if(length(param_list) >= 1){
    if('prob' %in% names(param_list)){
      prob <- param_list$prob
    } else {
      prob <- 0.9
    }

    if('scale_dist' %in% names(param_list)){
      scale_dist <- param_list$scale_dist
    } else {
      scale_dist <- TRUE
    }

    if('add_label' %in% names(param_list)){
      add_label <- param_list$add_label
    } else {
      add_label <- TRUE
    }
  } else {
    prob <- 0.9
    scale_dist <- TRUE
    add_label <- TRUE
  }

  if(isTRUE(scale_dist) & (!is.numeric(prob) | prob < 0 | prob > 1)){
    stop("`prob` must be a number between 0–1")
  }

  mod_summary <- summary.multiScaleR(x, prob = prob)
  sig_ <- mod_summary$opt_scale
  shp_ <- mod_summary$opt_shape
  titles <- rownames(sig_)
  dist_tab <- mod_summary$opt_dist

  ## DEBUG
  # browser()

  # if(x$kernel_inputs$kernel == 'expow'){
  #   shp_ <- summary(x)$opt_shape
  # } else {
  #   shp_ <- NULL
  # }

  df_list <- plot_list <- vector('list', sum(!is.nan(sig_$SE)))
  s <- which(!is.nan(sig_$SE))
  for(t in 1:length(s)){
    i <- s[t]
    d <- seq(1, round(max(sig_[i,])*1000,0),
             length.out = round(max(sig_[i,])*1000,0))
    wt <- scale_type(d = d,
                     kernel = x$kernel_inputs$kernel,
                     sigma = x$scale_est[[1]][i],
                     shape = x$shape_est[[1]][i],
                     output = 'wts')

    mx <- Hmisc::wtd.Ecdf(d, weights = wt)
    mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)

    d <- seq(1, mx, length.out = 100)
    wt <- scale_type(d = d,
                     kernel = x$kernel_inputs$kernel,
                     sigma = x$scale_est[[1]][i],
                     shape = x$shape_est[[1]][i],
                     output = 'wts')

    scale_d <- dist_tab[i,1]
    scale_lci <- dist_tab[i,2]
    scale_uci <- dist_tab[i,3]

    df_list[[i]] <- data.frame(dist = d,
                               wt = wt)
    mx_y <- max(wt)

    if(isTRUE(scale_dist)){
      if(min(prob) >= 0.8){
        ax <- max(d) * 0.08
        ay <- 0.08*max(wt)
      } else {
        ax <- max(d) * 0.8
        ay <- 0.9*max(wt)
      }
    }

    # browser()

    plot_ <- ggplot(data = df_list[[i]], aes(x = dist, y = wt)) +
      {if(isTRUE(scale_dist))
        geom_rect(xmin = scale_lci, xmax = scale_uci, ymin = -Inf, ymax = Inf,
                  fill = 'lightgrey', alpha = 0.25)

      } +
      {if(isTRUE(scale_dist))
        geom_vline(xintercept = scale_d,
                   linetype = 'dashed',
                   color = 'red')
      } +
      {if(isTRUE(scale_dist) & isTRUE(add_label))
        annotate('text', x = ax, y = ay,
                 label = paste0(prob*100,"% density \n Distance: ", round(scale_d, 0),
                                "\n  ", " 95% CI: ",round(scale_lci, 0), " - ", round(scale_uci, 0)))
      } +
      geom_line(linewidth = 1.25) +
      ggtitle(titles[i]) +
      xlab('Distance') +
      ylab('Weight') +
      cowplot::theme_cowplot()

    plot_list[[i]] <- plot_
  }
  return(plot_list)

  # } else if(length(param_list) >= 1){
  #   sig_ <- param_list$sigma
  #   shp_ <- param_list$shape
  #   kern <- param_list$kernel
  #
  #   if(is.null(sig_)){
  #     stop('\nA value for `sigma` must be provided!\n')
  #   }
  #   if(is.null(kern)){
  #     stop('\nYou must specify `kernel` function; See Details\n')
  #   }
  #   if(kern == 'expow' & is.null(shp_)){
  #     stop('\nBoth a `sigma` and `shape` parameter must be specified when using the `expow` kernel; See Details\n')
  #   }
  #
  #   # d <- seq(1, 1e6, length.out = 1e6)
  #   d <- seq(1, round(max(sig_)*1000,0),
  #            length.out = round(max(sig_)*1000,0))
  #   wt <- scale_type(d = d,
  #                    kernel = kern,
  #                    sigma = sig_,
  #                    shape = shp_,
  #                    output = 'wts')
  #
  #   mx <- Hmisc::wtd.Ecdf(d, weights = wt)
  #   mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)
  #
  #   d <- seq(1, mx, length.out = 100)
  #   wt <- scale_type(d = d,
  #                    kernel = kern,
  #                    sigma = sig_,
  #                    shape = shp_,
  #                    output = 'wts')
  #
  #   scale_d <- round(d[which(Hmisc::wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], -1)
  #
  #   df <- data.frame(dist = d,
  #                    wt = wt)
  #
  #   if(isTRUE(scale_dist)){
  #     if(min(prob) >= 0.8){
  #       ax <- max(d) * 0.075
  #       ay <- 0.075*max(wt)
  #     } else {
  #       ax <- max(d) * 0.8
  #       ay <- 0.9*max(wt)
  #     }
  #   }
  #   ggplot(data = df, aes(x = d, y = wt)) +
  #     geom_line(linewidth = 1.25) +
  #     {if(isTRUE(scale_dist))
  #       geom_vline(xintercept = scale_d,
  #                  linetype = 'dashed',
  #                  color = 'red')
  #     } +
  #     {if(isTRUE(scale_dist) & isTRUE(add_label))
  #       annotate('text', x = ax, y = ay,
  #                label = paste0(prob*100,"% density \n Distance: ", round(scale_d, 0)))
  #     } +
  #     xlab('Distance') +
  #     ylab('Weight') +
  #     cowplot::theme_cowplot()
  # } else {
  #   stop("Parameters not correctly specified to create plot. See Details and try again.")
  # }
}

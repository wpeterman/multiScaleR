#' Plot method for multiScaleR objects
#'
#' Plot kernel weight distributions from optimized \code{multiScaleR} objects.
#'
#' @param x An object of class \code{multiScaleR}.
#' @param ... Arguments to modify the plot. See Details.
#'
#' @details
#' Supported arguments include:
#' \itemize{
#'   \item \code{prob}: Cumulative weight cutoff for distance scale (default = 0.9).
#'   \item \code{scale_dist}: Logical; add vertical line for distance scale (default = TRUE).
#'   \item \code{add_label}: Logical; annotate scale distance and CI (default = TRUE).
#' }
#'
#' @return A list of ggplot2 objects.
#' @examples
#' \dontrun{
#' plot(x)
#'
#' plot(x, prob = 0.95)
#'
#' plot(x, scale_dist = FALSE)
#'
#' plot(x, scale_dist = TRUE, add_label = FALSE)
#' }
#'
#' @seealso \code{\link[multiScaleR]{plot_kernel}}
#' @export
#' @method plot multiScaleR
#'
#' @importFrom Hmisc wtd.Ecdf
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes annotate geom_line geom_vline ggplot xlab ylab theme_light geom_rect ggtitle
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
    stop("`prob` must be a decimal between 0 and 1")
  }

  mod_summary <- summary.multiScaleR(x, prob = prob)
  sig_ <- mod_summary$opt_scale
  shp_ <- mod_summary$opt_shape
  titles <- rownames(sig_)
  dist_tab <- mod_summary$opt_dist

  df_list <- plot_list <- vector('list', sum(!is.nan(sig_$SE)))
  s <- which(!is.nan(sig_$SE))
  for(t in 1:length(s)){
    i <- s[t]
    d <- seq(1, round(max(sig_[i,])*1000,0),
             length.out = round(max(sig_[i,])*1000,0))
    wt <- scale_type_r(d = d,
                       kernel = x$kernel_inputs$kernel,
                       sigma = x$scale_est[[1]][i],
                       shape = x$shape_est[[1]][i],
                       output = 'wts')

    mx <- wtd.Ecdf(d, weights = wt)
    mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)

    d <- seq(1, mx, length.out = 100)
    wt <- scale_type_r(d = d,
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
      theme_cowplot()

    plot_list[[i]] <- plot_
  }
  invisible(plot_list)
}

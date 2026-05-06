#' Plot method for multiScaleR objects
#'
#' Plots the optimized kernel weight distribution for each spatial covariate in
#' a fitted \code{multiScaleR} object, with optional annotations showing the
#' effective scale of effect and its 95\% confidence interval.
#'
#' @param x A fitted \code{multiScaleR} object returned by
#'   \code{\link{multiScale_optim}}.
#' @param ... Optional named arguments to customize the plot:
#' \describe{
#'   \item{\code{prob}}{Numeric between 0 and 1 (exclusive). Cumulative kernel
#'     density threshold used to mark the effective distance. Default:
#'     \code{0.9} (90\% of the kernel weight falls within the annotated
#'     distance).}
#'   \item{\code{scale_dist}}{Logical. If \code{TRUE} (default), a vertical
#'     dashed line is drawn at the effective distance, and a shaded rectangle
#'     spans the 95\% confidence interval.}
#'   \item{\code{add_label}}{Logical. If \code{TRUE} (default), an annotation
#'     in the upper-right corner reports the cumulative percentage, the
#'     effective distance, and the 95\% CI bounds (in projection units).}
#' }
#'
#' @return A list of \code{ggplot2} objects (one per covariate with a
#'   non-\code{NaN} standard error), returned invisibly. Plots are printed as
#'   a side effect.
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
#' plot(opt)
#'
#' plot(opt, prob = 0.95)
#'
#' plot(opt, scale_dist = FALSE)
#'
#' plot(opt, scale_dist = TRUE, add_label = FALSE)
#' }
#'
#' @seealso \code{\link[multiScaleR]{plot_kernel}}
#' @export
#' @method plot multiScaleR
#'
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
    # d <- seq(1, round(max(sig_[i,])*1000,0),
    #          length.out = round(max(sig_[i,])*1000,0))
    # wt <- scale_type_r(d = d,
    #                    kernel = x$kernel_inputs$kernel,
    #                    sigma = x$scale_est[[1]][i],
    #                    shape = x$shape_est[[1]][i],
    #                    output = 'wts')
    #
    # mx <- wtd.Ecdf(d, weights = wt)
    # mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)

    mx <- round(k_dist(sigma = max(sig_[i,]),
                       prob = 0.9999,
                       kernel = x$kernel_inputs$kernel,
                       beta = x$shape_est[[1]][i]), digits = -1)

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
        annotate('text', x = Inf, y = Inf,
                 hjust = 1.1, vjust = 1.4,
                 label = paste0(prob * 100, "% density\nDistance: ", round(scale_d, 0),
                                "\n95% CI: ", round(scale_lci, 0), " - ", round(scale_uci, 0)))
      } +
      geom_line(linewidth = 1.25) +
      ggtitle(titles[i]) +
      xlab('Distance') +
      ylab('Weight') +
      theme_cowplot()

    plot_list[[i]] <- plot_
  }
  # Print all plots
  lapply(plot_list, print)

  invisible(plot_list)
  }

#' @title Plot kernel densities
#'
#' @description Visualizes the weight (density) of a kernel function as a
#' function of distance, without requiring a fitted \code{multiScaleR} object.
#' Useful for exploring how different kernel types and sigma values translate
#' into spatial influence zones before running \code{\link{multiScale_optim}}.
#' To plot the kernel from a fitted model, use \code{plot(opt)} instead.
#'
#' @param prob Numeric between 0 and 1 (exclusive). Cumulative density threshold
#'   used to mark the effective distance on the plot. Default: \code{0.9} (90\%
#'   of the kernel weight falls within the marked distance).
#' @param sigma Positive numeric. The scale parameter of the kernel, in the
#'   same units as the projection used with \code{\link{kernel_prep}} (e.g.,
#'   metres when using a metric CRS). Controls the width of the kernel: larger
#'   values spread influence over greater distances.
#' @param beta Positive numeric. Shape parameter for the exponential power
#'   kernel (\code{kernel = "expow"}). Controls kernel shape: values near 1
#'   approximate a Laplace (double-exponential) shape; values near 2 approximate
#'   Gaussian; larger values approach a flat-top (platykurtic) shape. Typical
#'   range is 1–50. Ignored for all other kernel types. Default: \code{NULL}.
#' @param kernel Character. The kernel function to visualize. One of:
#'   \describe{
#'     \item{\code{"gaussian"}}{Gaussian (normal) decay; weights follow a
#'       half-normal density with standard deviation \code{sigma}.}
#'     \item{\code{"exp"}}{Negative exponential decay; weights follow
#'       \code{exp(-d / sigma)}.}
#'     \item{\code{"fixed"}}{Fixed-radius step function; all distances up to
#'       \code{sigma} receive equal weight.}
#'     \item{\code{"expow"}}{Exponential power kernel parameterized by
#'       \code{sigma} and \code{beta}.}
#'   }
#' @param scale_dist Logical. If \code{TRUE} (default), a vertical dashed line
#'   is added to the plot at the distance where the cumulative density equals
#'   \code{prob}.
#' @param add_label Logical. If \code{TRUE} (default), an annotation showing
#'   the cumulative percentage and the effective distance value is added to the
#'   plot. Requires \code{scale_dist = TRUE}.
#' @param ... Not used.
#'
#' @return A \code{ggplot2} object showing kernel weight (y-axis) as a function
#'   of distance (x-axis). The object is printed as a side effect and returned
#'   invisibly.
#'
#' @details
#' The x-axis range is set to cover 99.9\% of the cumulative density so that
#' the tails of the distribution are visible. The \code{prob} marker is rounded
#' to the nearest 10 distance units for display purposes.
#'
#' For fitted-model kernel plots (with confidence intervals), use
#' \code{plot(multiScaleR_object)} instead.
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
#'             beta = 2.5,
#'             kernel = 'expow')
#' @rdname plot_kernel
#' @export
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes annotate geom_line geom_vline ggplot xlab ylab theme_light geom_rect ggtitle
plot_kernel <- function(prob = 0.9,
                        sigma,
                        beta = NULL,
                        kernel,
                        scale_dist = TRUE,
                        add_label = TRUE,
                        ...) {
  if(isTRUE(scale_dist) & (!is.numeric(prob) | prob < 0 | prob > 1)){
    stop("`prob` must be a decimal between 0 and 1")
  }

  sig_ <- sigma
  shp_ <- beta
  kern <- kernel

  if(is.null(sig_)){
    stop('\nA value for `sigma` must be provided!\n')
  }
  if(is.null(kern)){
    stop('\nYou must specify `kernel` function; See Details\n')
  }
  if(kern == 'expow' & is.null(shp_)){
    stop('\nBoth a `sigma` and `beta` parameter must be specified when using the `expow` kernel; See Details\n')
  }

  # d <- seq(1, round(max(sig_)*1000,0),
  #          length.out = round(max(sig_)*1000,0))
  # wt <- scale_type_r(d = d,
  #                    kernel = kern,
  #                    sigma = sig_,
  #                    shape = shp_,
  #                    output = 'wts')
  #
  # mx <- wtd.Ecdf(d, weights = wt)
  # mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)

  mx <- round(k_dist(sigma = sig_,
               prob = 0.999,
               kernel = kern,
               beta = shp_), digits = -1)

  d <- seq(1, mx, length.out = 100)
  wt <- scale_type_r(d = d,
                     kernel = kern,
                     sigma = sig_,
                     shape = shp_,
                     output = 'wts')

  # scale_d <- round(d[which(wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], -1)
  scale_d <- round(k_dist(sigma = sig_,
                                prob = prob,
                                kernel = kern,
                                beta = shp_), digits = -1)

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
    theme_cowplot()
}


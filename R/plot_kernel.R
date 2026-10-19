#' @title Plot kernel densities
#'
#' @description Visualizes the weight (density) of a kernel function as a
#' function of distance, without requiring a fitted \code{multiScaleR} object.
#' Useful for exploring how different kernel types and sigma values translate
#' into spatial influence zones before running \code{\link{multiScale_optim}}.
#' To plot the kernel from a fitted model, use \code{plot(opt)} instead.
#'
#' @param prob Numeric between 0 and 1 (exclusive). Proportion of ideal kernel
#'   weight enclosed by the marked distance under the selected \code{dimension}.
#'   Default: \code{0.9}.
#' @param sigma Positive numeric. The scale parameter of the kernel, in the
#'   same units as the projection used with \code{\link{kernel_prep}} (e.g.,
#'   meters when using a metric CRS). Controls the width of the kernel: larger
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
#'   marks the distance enclosing \code{prob} of kernel weight under the
#'   selected dimensional definition.
#' @param add_label Logical. If \code{TRUE} (default), a caption reports the
#'   cumulative percentage, dimension, sigma, and effective distance below the
#'   plot. Requires \code{scale_dist = TRUE}.
#' @param dimension Character. \code{"2d"} (default) marks the radius enclosing
#'   \code{prob} of ideal kernel weight on a plane. \code{"1d"} restores the
#'   original distance enclosing that proportion along a line through the point.
#' @param ... Not used.
#'
#' @return A \code{ggplot2} object showing kernel weight (y-axis) as a function
#'   of distance in map units (x-axis), returned for display or customization.
#'
#' @details
#' The curve shows relative weight assigned to a cell at each distance,
#' normalized over the plotted sample distances. It is not the distribution of
#' total weight in concentric rings. Under the default \code{"2d"} definition,
#' the marker accounts for ring area (proportional to distance), so the area
#' under the displayed curve to its left is not \code{prob}.
#'
#' The x-axis covers at least 99.9\% of weight under the selected definition,
#' or \code{prob} when larger. Marker positions use the unrounded distance;
#' labels round to two decimal places. Sigma remains a kernel parameter, not
#' the marked radius. For a Gaussian kernel with sigma = 300 m, the 90\% marker
#' is about 644 m in 2D and 493 m in 1D.
#'
#' Development versions after 0.7.3 default to 2D; use
#' \code{dimension = "1d"} to reproduce the original distance definition.
#' These are ideal continuous-kernel distances, not site-specific retained
#' weight after \code{max_D}, raster boundaries, missing cells, or cell geometry
#' limit the extraction. See \code{\link{kernel_dist}} for both definitions.
#'
#' For fitted-model kernel plots (with confidence intervals), use
#' \code{plot(multiScaleR_object)} instead.
#'
#' @examples
#' ## Default: radius containing 95% of ideal planar weight
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
#' ## Original line-based definition
#' plot_kernel(sigma = 100, kernel = 'gaussian', dimension = '1d')
#' @seealso \code{\link{kernel_dist}}, \code{\link{plot.multiScaleR}}
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
                        dimension = c("2d", "1d"),
                        ...) {
  dimension <- match.arg(dimension)
  validate_scalar_numeric(prob, "prob", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)

  sig_ <- sigma
  shp_ <- beta
  kern <- kernel

  if(is.null(sig_)){
    stop('\nA value for `sigma` must be provided!\n')
  }
  if(is.null(kern)){
    stop('\nYou must specify `kernel` function; See Details\n')
  }
  validate_scalar_numeric(sig_, "sigma", positive = TRUE)
  kern <- match.arg(kern, c("gaussian", "exp", "expow", "fixed"))
  if(kern == 'expow' & is.null(shp_)){
    stop('\nBoth a `sigma` and `beta` parameter must be specified when using the `expow` kernel; See Details\n')
  }

  if (kern == "expow") validate_scalar_numeric(shp_, "beta", positive = TRUE)
  mx <- k_dist(sigma = sig_, prob = max(0.999, prob), kernel = kern,
                beta = shp_, dimension = dimension)
  if (kern == "fixed") mx <- sig_ * 1.1

  d <- seq(0, mx, length.out = 200)
  wt <- scale_type_r(d = d,
                     kernel = kern,
                     sigma = sig_,
                     shape = shp_,
                     output = 'wts')

  scale_d <- k_dist(sigma = sig_,
                                prob = prob,
                                kernel = kern,
                                beta = shp_, dimension = dimension)

  df <- data.frame(dist = d,
                   wt = wt)

  definition <- if (dimension == "2d") "of ideal planar weight (2D)" else "of line-based weight (1D)"
  ggplot(data = df, aes(x = dist, y = wt)) +
    geom_line(linewidth = 1.25) +
    {if(isTRUE(scale_dist))
      geom_vline(xintercept = scale_d,
                 linetype = 'dashed',
                 color = 'red')
    } +
    ggplot2::labs(caption = if (isTRUE(scale_dist) && isTRUE(add_label))
      paste0(prob * 100, "% ", definition, "\nSigma: ", signif(sig_, 4), "; ",
             if (dimension == "2d") "Radius: " else "Distance: ", round(scale_d, 2))) +
    xlab('Distance (map units)') +
    ylab('Weight') +
    theme_cowplot() +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 10, hjust = 0))
}


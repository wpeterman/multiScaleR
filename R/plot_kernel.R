#' @title Plot a Kernel Weight Curve
#'
#' @description Draws the weight a kernel gives to a raster cell as a function
#' of its distance from a survey point, and marks the effective distance that
#' holds a chosen share of the kernel weight. No fitted model is needed, so
#' use it to explore how kernel type, sigma, and shape translate into spatial
#' reach before you run [multiScale_optim()]. To plot the kernels of a fitted
#' model, with their confidence intervals, use [plot.multiScaleR()]. The
#' kernels are defined in [multiScaleR-kernels].
#'
#' @param prob Numeric between 0 and 1 (exclusive). Share of the ideal kernel
#'   weight held within the marked distance, under the selected `dimension`.
#'   Default: `0.9`.
#' @param sigma Positive numeric scalar. The kernel scale parameter, in the
#'   map units of the projection used with [kernel_prep()] (for example,
#'   meters). Larger values spread the weight over greater distances. For
#'   `kernel = "fixed"` it is the buffer radius. Required.
#' @param shape Positive numeric scalar. Shape of the exponential power kernel
#'   (`kernel = "expow"`), required for that kernel and ignored for the
#'   others. `shape = 1` gives the negative exponential curve, `shape = 2` the
#'   Gaussian form (with standard deviation sigma / sqrt(2)), and larger
#'   values a flat top with a steep edge that approaches a fixed buffer. The
#'   default refit search of [multiScale_optim()] keeps the shape between 0.75
#'   and 50. Default: `NULL`.
#' @param kernel Character. The kernel to draw: `"gaussian"`, `"exp"`
#'   (negative exponential), `"expow"` (exponential power), or `"fixed"`
#'   (fixed-radius buffer). Required. The weight functions are tabulated in
#'   [multiScaleR-kernels].
#' @param scale_dist Logical. If `TRUE` (default), a dashed red vertical line
#'   marks the distance that holds `prob` of the kernel weight under the
#'   selected `dimension`.
#' @param add_label Logical. If `TRUE` (default), a caption below the plot
#'   reports the percentage of weight, the dimensional definition, sigma, and
#'   the marked distance. Shown only when `scale_dist = TRUE`.
#' @param dimension Character. `"2d"` (default) marks the radius of the circle
#'   that holds `prob` of the ideal kernel weight on a plane. `"1d"` marks the
#'   original line-based distance, which holds that share along a transect
#'   through the point. See [kernel_dist()].
#' @param ... Not used.
#' @param beta Deprecated former name of `shape`. It is still accepted, with a
#'   warning, and will be removed in a future release. Supplying both `shape`
#'   and `beta` with different values is an error.
#'
#' @return A `ggplot2` object with distance in map units on the x-axis and
#'   relative kernel weight on the y-axis. It prints when returned at the
#'   console. Add layers or themes to customize it.
#'
#' @details
#' **How to read the plot.** The curve shows the relative weight given to a
#' single cell at each distance. Weights are normalized to sum to one over the
#' 200 plotted distances, so only the shape of the curve and its relative
#' height matter. The curve is not the distribution of total weight in
#' concentric rings. Under the default `"2d"` definition the marker accounts
#' for the area of each ring, which grows with distance, so the area under
#' the plotted curve to the left of the marker is not `prob`. Sigma remains a
#' kernel parameter, not the marked radius. For a Gaussian kernel with
#' sigma = 300 m, the 90% marker sits at about 644 m under `"2d"` and 493 m
#' under `"1d"`.
#'
#' The x-axis extends to the distance holding at least 99.9% of the weight
#' under the selected definition (or `prob`, if larger). For the fixed kernel
#' it extends to 1.1 times sigma. Marker positions use the unrounded distance,
#' and the caption rounds it to two decimal places.
#'
#' Versions 0.7.4 and later mark the `"2d"` radius by default. Set
#' `dimension = "1d"` to reproduce plots from earlier versions. Both are ideal
#' continuous-kernel distances. They do not reflect truncation at `max_D`,
#' raster edges, missing cells, or cell geometry at a particular site.
#'
#' @examples
#' ## Radius holding 95% of the ideal planar weight, sigma = 100 map units
#' plot_kernel(prob = 0.95, sigma = 100, kernel = "gaussian")
#' # The negative exponential reaches farther for the same sigma
#' plot_kernel(prob = 0.95, sigma = 100, kernel = "exp")
#' # Equal weight out to sigma, then zero
#' plot_kernel(prob = 0.95, sigma = 100, kernel = "fixed")
#' # A flatter top and steeper edge than the Gaussian
#' plot_kernel(prob = 0.95, sigma = 100, shape = 5, kernel = "expow")
#'
#' ## The original line-based definition marks a shorter distance
#' plot_kernel(sigma = 100, kernel = "gaussian", dimension = "1d")
#'
#' @seealso [multiScaleR-kernels], [kernel_dist()], [plot.multiScaleR()].
#'   `vignette("v02-user-guide", package = "multiScaleR")` compares the four
#'   kernels side by side.
#' @rdname plot_kernel
#' @export
#' @importFrom cowplot theme_cowplot
#' @importFrom ggplot2 aes annotate geom_line geom_vline ggplot xlab ylab theme_light geom_rect ggtitle
plot_kernel <- function(prob = 0.9,
                        sigma,
                        shape = NULL,
                        kernel,
                        scale_dist = TRUE,
                        add_label = TRUE,
                        dimension = c("2d", "1d"),
                        ...,
                        beta = NULL) {
  dimension <- match.arg(dimension)
  validate_scalar_numeric(prob, "prob", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)

  sig_ <- sigma
  # `shape` keeps the third position that `beta` held, so positional calls
  # are unaffected; `beta` is a deprecated alias (one release).
  shp_arg <- .msr_resolve_shape(shape, beta)
  shp_ <- shp_arg$shape
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
    stop('\nBoth a `sigma` and `shape` (formerly `beta`) parameter must be specified when using the `expow` kernel; See Details\n')
  }

  if (kern == "expow") validate_scalar_numeric(shp_, shp_arg$name, positive = TRUE)
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


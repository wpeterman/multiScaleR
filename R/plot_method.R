#' Plot method for multiScaleR objects
#'
#' Plots the optimized kernel weight distribution for each spatial covariate in
#' a fitted \code{multiScaleR} object, with optional annotations showing the
#' effective scale of effect and its 95\% confidence interval.
#'
#' @param x A fitted \code{multiScaleR} object returned by
#'   \code{\link{multiScale_optim}}.
#' @param dimension Character. \code{"2d"} (default) marks the radius containing
#'   the requested proportion of ideal planar kernel weight. \code{"1d"}
#'   restores the original distance along a line through the point.
#' @param ... Optional named arguments to customize the plot:
#' \describe{
#'   \item{\code{prob}}{Numeric between 0 and 1 (exclusive). Cumulative kernel
#'     weight enclosed by the marked distance under the selected
#'     \code{dimension}. Default: \code{0.9}.}
#'   \item{\code{scale_dist}}{Logical. If \code{TRUE} (default), a vertical
#'     dashed line is drawn at the effective distance, and a shaded rectangle
#'     spans the 95\% confidence interval.}
#'   \item{\code{add_label}}{Logical. If \code{TRUE} (default), a caption
#'     below the plot reports sigma, the dimensional definition,
#'     the effective distance, and 95\% CI bounds (in map units).}
#' }
#'
#' @return A list of \code{ggplot2} objects (one per covariate with a finite
#'   scale and effective distance), returned invisibly. Plots are printed as
#'   a side effect.
#' @details
#' Development versions after 0.7.3 default to 2D reporting, as do
#' \code{\link{kernel_dist}} and \code{\link{plot_kernel}}. Set
#' \code{dimension = "1d"} to reproduce the original line-based distances.
#' This switch changes distance markers and converted confidence limits, not
#' the fitted sigma, regression coefficients, or likelihood.
#'
#' The curve shows cell weight declining with distance, normalized over the
#' plotted sample distances. It is not a radial weight distribution. A planar
#' marker accounts for the increasing area of concentric rings, so the area
#' under the displayed curve to its left is not the stated proportion.
#' Distances refer to an ideal untruncated kernel. Actual retained weight can
#' differ with \code{max_D}, raster edges, missing cells, and cell geometry.
#'
#' Hard-radius covariates show a uniform buffer and mark its full radius,
#' regardless of \code{prob} or \code{dimension}. A shaded interval is drawn
#' only when finite limits are available.
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
#' plot(opt, dimension = "1d")  # Original line-based distances
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
                             dimension = c("2d", "1d"),
                             ...) {
  param_list <- list(...)
  dimension <- match.arg(dimension)

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

  validate_scalar_numeric(prob, "prob", lower = 0, upper = 1,
                          inclusive_lower = FALSE, inclusive_upper = FALSE)

  mod_summary <- summary.multiScaleR(x, prob = prob, dimension = dimension)
  sig_ <- mod_summary$opt_scale
  shp_ <- mod_summary$opt_shape
  titles <- rownames(sig_)
  dist_tab <- mod_summary$opt_dist

  s <- which(is.finite(sig_$Mean) & is.finite(dist_tab[, 1]))
  plot_list <- vector('list', length(s))
  for(t in seq_along(s)){
    i <- s[t]
    hard_radius <- titles[i] %in% mod_summary$hard_radius
    kernel <- if (hard_radius) "fixed" else x$kernel_inputs$kernel
    shape <- if (is.null(shp_)) NULL else shp_[i, 1]
    sigma_max <- max(unlist(sig_[i, c("Mean", "2.5%", "97.5%")]), na.rm = TRUE)
    mx <- if (hard_radius || kernel == "fixed") sigma_max * 1.1 else
      k_dist(sigma = sigma_max, prob = max(0.9999, prob),
              kernel = kernel, beta = shape, dimension = dimension)

    d <- seq(0, mx, length.out = 200)
    wt <- scale_type_r(d = d,
                       kernel = kernel,
                       sigma = sig_$Mean[i],
                       shape = shape,
                       output = 'wts')

    scale_d <- dist_tab[i,1]
    scale_lci <- dist_tab[i,2]
    scale_uci <- dist_tab[i,3]

    df <- data.frame(dist = d, wt = wt)
    definition <- if (dimension == "2d") "of ideal planar weight (2D)" else "of line-based weight (1D)"
    label <- if (hard_radius) {
      paste0("Buffer radius: ", round(scale_d, 2))
    } else {
      paste0(prob * 100, "% ", definition, "\nSigma: ", signif(sig_$Mean[i], 4),
             "; ", if (dimension == "2d") "Radius: " else "Distance: ", round(scale_d, 2))
    }
    if (all(is.finite(c(scale_lci, scale_uci)))) {
      label <- paste0(label, "\n95% CI: ", round(scale_lci, 2), " - ", round(scale_uci, 2))
    }

    plot_ <- ggplot(data = df, aes(x = dist, y = wt)) +
      {if(isTRUE(scale_dist) && all(is.finite(c(scale_lci, scale_uci))))
        geom_rect(xmin = scale_lci, xmax = scale_uci, ymin = -Inf, ymax = Inf,
                  fill = 'lightgrey', alpha = 0.25)

      } +
      {if(isTRUE(scale_dist))
        geom_vline(xintercept = scale_d,
                   linetype = 'dashed',
                   color = 'red')
      } +
      ggplot2::labs(caption = if (isTRUE(scale_dist) && isTRUE(add_label)) label) +
      geom_line(linewidth = 1.25) +
      ggtitle(titles[i]) +
      xlab('Distance (map units)') +
      ylab('Weight') +
      theme_cowplot() +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 10, hjust = 0))

    plot_list[[t]] <- plot_
  }
  # Print all plots
  lapply(plot_list, print)

  invisible(plot_list)
  }

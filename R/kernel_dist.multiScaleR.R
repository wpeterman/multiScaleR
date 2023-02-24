#' @title Scale Distance
#' @description Function to estimate the effective distance encomapssing a specified probability density of the kernel density function
#' @param x \code{\link[multiScaleR]{multiScaleR}} object
#' @param prob Density probability cutoff for calculating distance, Default: 0.9
#' @return Distance
#' @details None
#' @examples
#' kernel_dist(x)
#'
#' @rdname kernel_dist.multiScaleR
#' @export
kernel_dist.multiScaleR <- function(x,
                                    prob = 0.9){
  if(class(x) != 'multiScaleR'){
    stop("x is not of class multiScaleR")
  }
  d <- seq(1, 1e6, length.out = 1e6)
  wt <- scale_type(d = d,
                   kernel = x$kernel_inputs$kernel,
                   sigma = x$scale_est[[1]],
                   shape = x$shape_est[[1]],
                   output = 'wts')

  dist <- DescTools::Quantile(d,
                              weights = wt,
                              probs = prob)
  return(dist)
}

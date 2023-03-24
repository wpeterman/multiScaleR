#' @title Scale Function
#' @description Scaling function to be applied to rasters
#' @param d Vector of distances
#' @param kernel Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow'; Default: 'gaussian')
#' @param sigma Scaling parameter
#' @param shape Shape parameter if using exponential power kernel
#' @param r_stack.df Dataframe values extracted from rasters
#' @param output If NULL, a vector of weights is returned, otherwise a weighted raster values are returned, Default: NULL
#' @return A vector of weights or vector of weighted raster values
#' @details DETAILS
#' @examples
#' ### TO BE COMPLETED ###
#' @rdname scale_type
#' @keywords internal

scale_type <- function(d,
                       kernel = 'gaussian',
                       sigma,
                       shape = NULL,
                       r_stack.df = NULL,
                       output = NULL) {
  if(kernel == 'exp'){
    d1 <- data.frame(exp(-outer(d, sigma, "/")))
    w0 <- mapply("*",d1, (1/(2*pi*sigma^2)))

  } else if(kernel == 'fixed'){
    w0 <- outer(d, sigma, "<") * 1 ## Fixed

  } else if(kernel == 'gaussian'){
    d1 <- data.frame(exp(-outer(d^2, (2*sigma^2), FUN = "/")))
    w0 <- mapply("*", d1, (1/(pi*sigma^2)))

    # browser()
  } else {
    d1 <- data.frame(exp(-(outer(d,shape,"^") / (sigma^shape))))
    w0 <- mapply("*", d1, (shape/((2*pi*sigma^2)*gamma(2/shape))))

  }
  w <- apply(w0, 2, function(x) x/sum(x))

  if(!is.null(output)){
    return(w)
  } else {
    colSums(r_stack.df * w)  ## exact_extract solution
  }
}

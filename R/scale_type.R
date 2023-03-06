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
#' @export
#' @keywords internal

scale_type <- function(d,
                       kernel = 'gaussian',
                       sigma,
                       shape = NULL,
                       r_stack.df = NULL,
                       output = NULL) {
  if(kernel == 'exp'){
    # w0 <- ((1/(2*pi*sigma^2)) * exp(-(d / sigma))) * as.numeric(d!=0)  ## Neg. Exponential
    # w0 <- ((1/(2*pi*sigma^2)) * exp(-outer(d, sigma, "/")))#[,1]

    d1 <- data.frame(exp(-outer(d, sigma, "/")))
    w0 <- mapply("*",d1, (1/(2*pi*sigma^2)))
    # w0 <- d1 %*% diag(1/(2*pi*sigma^2))
    # w0 <- t(t(d1)*(1/(2*pi*sigma^2)))

  } else if(kernel == 'fixed'){
    # w0 <- as.numeric(d != 0 & outer(d, sigma, "<"))  ## Fixed
    w0 <- outer(d, sigma, "<") * 1 ## Fixed
    # w0 <- w0[,1]

  } else if(kernel == 'gaussian'){
    # w0_ <- ((1/(pi*sigma[1]^2)) * exp(-(d^2 / (2*sigma[1]^2)))) #* as.numeric(d!=0)  ## Original

    d1 <- data.frame(exp(-outer(d^2, (2*sigma^2), FUN = "/")))
    w0 <- mapply("*", d1, (1/(pi*sigma^2)))
    # w0 <- d1 %*% diag(1/(pi*sigma^2))
    # w0 <- t(t(d1)*(1/(pi*sigma^2)))
# browser()
  } else {
    # w0 <- ((shape/((2*pi*sigma^2)*gamma(2/shape))) * exp(-(d^shape / (sigma^shape)))) * as.numeric(d!=0)  ## Exponential power
    # w0 <- ((shape/((2*pi*sigma^2)*gamma(2/shape)))  * exp(-(outer(d,shape,"^") / (sigma^shape))))#[,1] ## Exponential power

    d1 <- data.frame(exp(-(outer(d,shape,"^") / (sigma^shape))))
    # w0 <- d1 %c*% (shape/((2*pi*sigma^2)*gamma(2/shape)))
    w0 <- mapply("*", d1, (shape/((2*pi*sigma^2)*gamma(2/shape))))
    # w0 <- t(t(d1)*(shape/((2*pi*sigma^2)*gamma(2/shape))))

  }

  # if(kernel == 'exp'){
  #   # w0 <- ((1/(2*pi*sigma^2)) * exp(-(d / sigma))) * as.numeric(d!=0)  ## Neg. Exponential
  #   w0 <- ((1/(2*pi*sigma^2)) * exp(-outer(d, sigma, "/")))#[,1]
  #
  # } else if(kernel == 'fixed'){
  #   # w0 <- as.numeric(d != 0 & outer(d, sigma, "<"))  ## Fixed
  #   w0 <- outer(d, sigma, "<") * 1 ## Fixed
  #   # w0 <- w0[,1]
  #
  # } else if(kernel == 'gaussian'){
  #   # w0_ <- ((1/(pi*sigma[1]^2)) * exp(-(d^2 / (2*sigma[1]^2)))) #* as.numeric(d!=0)  ## Original
  #   w0 <- ((1/(pi*sigma^2)) * exp(-outer(d^2, (2*sigma^2), FUN = "/")))#[,1] #* as.numeric(d!=0)  ## Original
  #
  # } else {
  #   # w0 <- ((shape/((2*pi*sigma^2)*gamma(2/shape))) * exp(-(d^shape / (sigma^shape)))) * as.numeric(d!=0)  ## Exponential power
  #   w0 <- ((shape/((2*pi*sigma^2)*gamma(2/shape)))  * exp(-(outer(d,shape,"^") / (sigma^shape))))#[,1] ## Exponential power
  # }

  # w <- w0/sum(w0))
  # w <- sapply(w0, function(x) x/sum(x))
  w <- apply(w0, 2, function(x) x/sum(x))

  if(!is.null(output)){
    return(w)
  } else {
    colSums(r_stack.df * w)  ## exact_extract solution
  }
}

#' @title Scale of effect distance
#' @description Calculate and plot the scale of effect distance
#' @param scale_est Estimated scale (sigma)
#' @param proportion Density of kernel weight, Default: 0.95
#' @param plot Logical, Default: FALSE
#' @return Calculated distance
#' @details Only valid for `gaussian`, `exp`, and `fixed`
#' @examples
#' ## TO BE COMPLETED ##
#' @export

scale_dist <- function(scale_est,
                       proportion = 0.95,
                       plot = FALSE){

    d <- 1:1e6
    w0 <- exp(-(d^2 / (2*scale_est^2)))  ## Original
    w0[d==0] <- 0

    ## Calculate weighted average
    wts <- w0/sum(w0)

    ## Find the distance where specified proportion of smooth weight lies
    scale_d <- DescTools::Quantile(d,
                                   weights = wts,
                                   probs = proportion)

    full_d <- DescTools::Quantile(d,
                                  weights = wts,
                                  probs = 0.99)

    if(isTRUE(plot) & length(scale_est) == 1){
      vals <- d <= (full_d * 1.05)
      plot(wts[vals] ~ d[vals], type = 'l', lwd = 2,
           main = 'Scale of Effect', xlab = 'Distance', ylab = 'Smoothing weight')
      abline(v = scale_d, col = 'red', lty = 2)
    }

  return(scale_d)
}



# Test function -----------------------------------------------------------

## Set `plot = T` to generate a plot in addition to returning the distance value
## Default for function is `plot = F`
scale_dist(scale = 50,
           proportion = 0.975,
           plot = T)


## Calculate for mutliple values, like a column of estimates
sapply(c(50,100,500), scale_dist)

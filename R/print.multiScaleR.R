#' @title Print multiScaleR
#' @description Print method for fitted multiScaleR objects
#' @param x Object of class `multiScaleR` created from running  \code{\link[multiScale_optim]{multiScale_optim}}
#' @param ... Additional parameters; not used
#' @return Summary of `multiScaleR` optimized object
#' @details Summary of optimization and model results
#' @examples
#' ## Not Run
#' print(opt)
#' @rdname print
#' @usage
#' print(x)
#' @export
print.multiScaleR <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)

  cat('\n\nKernel used: \n')
  cat(x$kernel_inputs$kernel)

  cat("\n\n***** Optimized Scale of Effect -- Sigma *****\n\n")
  print(x$scale_est)

  if(x$kernel_inputs$kernel == 'expow'){
    cat("\n\n***** Optimized Kernel Shape Parameter *****\n\n")
    print(x$shape_est)
  }

  cat("\n\n***** Optimized Scale of Effect -- Distance *****\n")
  cat("90% Kernel Weight")
  cat("\n\n")
  print(kernel_dist(x))

  cat("\n  ==================================== ")
  cat("\n\n ***** Fitted Model *****\n")
  print(x$opt_mod)
}

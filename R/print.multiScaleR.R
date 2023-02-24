#' @title Print multiScaleR
#' @description Print method for fitted multiScaleR objects
#' @param x Object of class `multiScapeR` created from running  \code{\link[multiScaleR]{multiScale_optim}}
#' @param ... Additional parameters; not used
#' @return OUTPUT_DESCRIPTION
#' @details Summary of optimization and model results
#' @examples
#' ## Not Run
#' print(opt)
#' @rdname print.multiScaleR
#' @usage
#' print(x)
#' @export
print.multiScaleR <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n\n***** Optimized Scale of Effect *****\n\n")
  cat('Kernel used: ')
  cat(x$kernel_inputs$kernel)
  cat("\n\n")
  print(x$scale_est)
  if(x$kernel_inputs$kernel == 'expow'){
    cat("\n\n***** Optimized Kernel Shape Parameter *****\n\n")
    print(x$shape_est)
  }

  cat("\n  ==================================== ")
  cat("\n\n ***** Fitted Model *****\n")
  print(x$opt_mod)
}

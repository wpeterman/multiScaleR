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


# Warning Messages --------------------------------------------------------


  if(1 %in% x$warn_message){
    cat(red("\n WARNING!!!\n",
            "The estimated scale of effect extends beyond the maximum distance specified.\n",
            "Consider increasing " %+% blue$bold("max_D") %+% " in `kernel_prep` to ensure accurate estimation of scale.\n\n"))
  }

  if(2 %in% x$warn_message){
    cat(red("\n WARNING!!!\n",
            "The standard error of one or more `sigma` estimates is >= 50% of the estimated mean value.\n",
            "Carefully assess whether or not this variable is meaningful in your analysis and interpret with caution.\n\n"))
  }

  if(3 %in% x$warn_message){
    cat(red("\n WARNING!!!\n",
            "The standard error of one or more `shape` estimates is >= 50% of the estimated mean value.\n",
            "Carefully assess if the Exponential Power kernel is appropriate, whether or not this variable is meaningful in your analysis, and interpret with caution.\n\n"))
  }
}

#' @title Print multiScaleR summary
#' @description Print summary method for fitted multiScaleR objects
#' @param x Object of class multiScaleR
#' @param ... Not used
#' @return Prints full summary of fitted `multiScaleR` object
#' @examples
#' print(x)
#' @rdname print.summary.multiScaleR
#' @export


print.summary.multiScaleR <- function(x, ...){


  cat("\nCall:\n")
  print(x$call)

  cat("\n\nKernel used:\n")
  cat(x$kernel)

  cat("\n\n***** Optimized Scale of Effect -- Sigma *****\n\n")
  print(x$opt_scale)
  cat("\n\n  ==================================== ")

  if(!is.null(x$opt_shape)){
    cat("\n\n***** Optimized Kernel Shape *****\n\n")
    print(x$opt_shape)
    cat("\n\n  ==================================== ")
  }

  cat("\n\n***** Optimized Scale of Effect -- Distance *****\n")
  dist_print <- (paste0(x$prob*100,"% Kernel Weight"))
  cat(dist_print)
  cat("\n\n")
  print(x$opt_dist)
  cat("\n\n  ==================================== ")


  cat('\n\n *****     Fitted Model Summary     *****')
  print(summary(x$fitted_mod))
}



#' @title multiScaleR summary
#' @description Internal summary function
#' @param object Object of class multiScaleR
#' @param ... Additional parameters. Currently, only `prob` can be specified for calculating the distance scale of effect
#' @export
#' @noRd
#' @keywords internal
summary.multiScaleR <- function(object,...){

  param_list <- list(...)

  # if(length(param_list) >= 1) browser()

  if(any(class(object$opt_mod) == 'gls')){
    df <- object$opt_mod$dims$N - object$opt_mod$dims$p
    names <- all.vars(formula(object$opt_mod)[-2])

  } else if(any(grepl("^unmarked", class(object$opt_mod)))){
    df <- dim(object$opt_mod@data@y)[1]
    names <- all.vars(object$opt_mod@formula)

  } else {
    df <- object$opt_mod$df.residual
    names <- all.vars(formula(object$opt_mod)[-2])
  }

  tab_scale <- ci_func(object$scale_est,
                       df = df,
                       min_D = object$min_D,
                       names = row.names(object$scale_est))


  if(length(param_list) >= 1){
    if('prob' %in% names(param_list)){
      prob <- param_list$prob
    }
  } else {
    prob <- 0.9
  }

  if(!is.null(object$shape_est)){
    tab_shape <- ci_func(object$shape_est,
                         df = object$opt_mod$df.residual,
                         names = row.names(object$shape_est))

    out <- list(opt_scale = tab_scale,
                opt_shape = tab_shape,
                opt_dist = kernel_dist(object, prob = prob),
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                call = object$call)
  } else {
    out <- list(opt_scale = tab_scale,
                opt_shape = NULL,
                opt_dist = kernel_dist(object, prob = prob),
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                call = object$call)
  }


  class(out) <- 'summary.multiScaleR'
  out
}

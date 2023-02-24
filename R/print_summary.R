#' @title Print multiScaleR summary
#' @description Print summary method for fitted multiScaleR objects
#' @param x PARAM_DESCRIPTION
#' @param ... PARAM_DESCRIPTION
#' @return Prints full summary of fitted `multiScaleR` object
#' @examples
#' print(x)
#' @rdname print.summary.multiScaleR
#' @export


print.summary.multiScaleR <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n\n***** Optimized Scale of Effect *****\n\n")
  print(x$opt_scale)
  cat("\n\n  ==================================== ")

  if(!is.null(x$opt_shape)){
    cat("\n\n***** Optimized Kernel Shape *****\n\n")
    print(x$opt_shape)
    cat("\n\n  ==================================== ")
  }

  cat('\n\n *****     Fitted Model Summary     *****')
  print(summary(x$fitted_mod))
}



#' @title multiScaleR summary
#' @description Internal summary function
#' @rdname summary.multiScaleR
#' @export
#' @keywords internal
summary.multiScaleR <- function(object,...){

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
  if(!is.null(object$shape_est)){
    tab_shape <- ci_func(object$shape_est,
                         df = object$opt_mod$df.residual,
                         names = row.names(object$shape_est))

    out <- list(opt_scale = tab_scale,
                opt_shape = tab_shape,
                fitted_mod = object$opt_mod,
                call = object$call)
  } else {
    out <- list(opt_scale = tab_scale,
                opt_shape = NULL,
                fitted_mod = object$opt_mod,
                call = object$call)
  }


  class(out) <- 'summary.multiScaleR'
  out
}

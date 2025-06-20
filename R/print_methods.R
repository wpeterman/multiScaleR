#' @title Print method for summary_multiScaleR
#' @description Print method for objects of class \code{summary_multiScaleR}.
#' @param x A \code{summary_multiScaleR} object
#' @param ... Ignored
#' @export
#' @method print summary_multiScaleR
print.summary_multiScaleR <- function(x, ...){
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


  cat('\n\n *****     Fitted Model Summary     *****\n\n')
  if(any(grepl("^unmarked", class(x$fitted_mod)))) {
    print(x$fitted_mod)
  } else {
    print(summary(x$fitted_mod))
  }

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

#' Summarize multiScaleR objects
#'
#' Summarizes output from \code{multiScale_optim}.
#'
#' @param object An object of class \code{multiScaleR}.
#' @param ... Optional arguments passed to the method (e.g., \code{prob} for cumulative kernel weight threshold).
#'
#' @return An object of class \code{summary_multiScaleR}.
#' @export
#' @method summary multiScaleR
summary.multiScaleR <- function(object,...){

  param_list <- list(...)

  if(any(class(object$opt_mod) == 'gls')){
    df <- object$opt_mod$dims$N - object$opt_mod$dims$p
    names <- all.vars(formula(object$opt_mod)[-2])

  } else if(any(grepl("^unmarked", class(object$opt_mod)))){
    df <- dim(object$opt_mod@data@y)[1]
    names <- all.vars(object$opt_mod@formula)

  } else {
    df <- get_df(object$opt_mod, type = "residual")
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

  ## DEBUG
  # browser()

  if(!is.null(object$shape_est)){
    tab_shape <- ci_func(object$shape_est,
                         df = df,
                         min_D = NULL,
                         names = row.names(object$shape_est))

    out <- list(opt_scale = tab_scale,
                opt_shape = tab_shape,
                opt_dist = kernel_dist(object, prob = prob),
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                warn_message = object$warn_message,
                call = object$call)
  } else {
    out <- list(opt_scale = tab_scale,
                opt_shape = NULL,
                opt_dist = kernel_dist(object, prob = prob),
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                warn_message = object$warn_message,
                call = object$call)
  }


  class(out) <- c('summary_multiScaleR')
  if (sys.nframe() == 1L) print(out)
  out
}

#' @title Print method for multiScaleR
#' @description Print method for objects of class \code{multiScaleR}.
#' @param x A \code{multiScaleR} object
#' @param ... Ignored
#' @export
#' @method print multiScaleR
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

#' @title Print method for multiScaleR_data
#' @description Print method for objects of class \code{multiScaleR_data}.
#' @param x A \code{multiScaleR_data} object
#' @param ... Ignored
#' @export
#' @method print multiScaleR_data
print.multiScaleR_data <- function(x, ...){
  cat("\nThere are ")
  cat(paste0(nrow(x$kernel_dat)," observations at ", ncol(x$kernel_dat), ' spatial covariate(s): \n'))
  cat(colnames(x$kernel_dat))
  cat("\n\nThe specified kernel is:\n")
  cat(x$kernel)
  # cat("\n\nSparse Matrix contains: ")
  # cat(paste0(length(x$D@x), ' elements\n'))
  cat("\n\nNumber of elements: \n")
  cat(paste0(length(x$d_list[[1]])))
  cat("\nMinimum Distance:\n")
  cat(x$min_D)
  cat("\nMaximum Distance:\n")
  cat(x$max_D)
  cat("\nUnit Conversion:\n")
  cat(x$unit_conv)
}

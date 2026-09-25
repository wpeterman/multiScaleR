#' @keywords internal
#' @noRd
format_max_distance_warning <- function(diagnostics) {
  if (is.null(diagnostics) || is.null(diagnostics$max_distance)) {
    return(list(
      headline = "The estimated scale of effect extends beyond the maximum distance specified.",
      suggested_max_D = NA_real_
    ))
  }

  diag <- diagnostics$max_distance

  if (!isTRUE(diag$triggered)) {
    return(list(
      headline = "The estimated scale of effect extends beyond the maximum distance specified.",
      suggested_max_D = diag$suggested_max_D
    ))
  }

  vars <- paste(diag$variables, collapse = ", ")
  list(
    headline = paste0(
      "The estimated scale of effect approaches or exceeds the current search extent",
      if (nzchar(vars)) paste0(" for: ", vars) else "",
      "."
    ),
    suggested_max_D = diag$suggested_max_D
  )
}

.msr_print_next_run <- function(next_run) {
  if (is.null(next_run) || !is.list(next_run)) {
    return(invisible(NULL))
  }

  cat("\n\n***** Follow-up recommendation *****\n")

  if (!is.null(next_run$n_cores)) {
    cat("Conservative recommended `n_cores`: ", next_run$n_cores, "\n", sep = "")
  }

  cat("Recommended `max_D`: ", next_run$max_D, "\n", sep = "")
  cat("Recommended starting `sigma` values:\n")
  print(next_run$start_sigma)

  if (!is.null(next_run$start_shape)) {
    cat("\nRecommended starting `shape` values:\n")
    print(next_run$start_shape)
  }

  if (!is.null(next_run$action)) {
    cat("\nAction:\n")
    cat(next_run$action, "\n")
  }

  invisible(next_run)
}

# Print the distance section shared by `print.multiScaleR()` and
# `print.summary_multiScaleR()`. `dist_2d` is NULL for summary objects created
# before the planar distance was added; those print the one-dimensional table
# alone, as before.
.msr_print_distance <- function(dist_1d, dist_2d, prob, hard_radius) {
  pct <- paste0(prob * 100, "%")
  cat("\n\n***** Optimized Scale of Effect -- Distance *****\n")
  weighted <- setdiff(rownames(dist_1d), hard_radius)

  if (!length(weighted)) {
    cat("Radius (map units) of each hard-radius covariate\n\n")
    print(dist_1d)
    return(invisible(NULL))
  }
  if (is.null(dist_2d)) {
    cat(paste0("1D: ", pct, " kernel weight along a line for weighted covariates; radius for hard-radius covariates"))
    cat("\n\n")
    print(dist_1d)
    return(invisible(NULL))
  }

  cat("Scale distances (map units) at ", pct, ", two definitions.\n", sep = "")
  if (length(hard_radius)) {
    cat("Hard-radius covariates show their radius in both tables.\n")
  }
  cat("\n1D, along a line through the point (`kernel_dist()` default):\n")
  print(dist_1d)
  cat("\n2D, radius of the circle for the ideal untruncated kernel (`dimension = \"2d\"`):\n")
  print(dist_2d)
  invisible(NULL)
}

#' @title Print method for summary_multiScaleR
#' @description Print method for objects of class \code{summary_multiScaleR}.
#' @param x A \code{summary_multiScaleR} object
#' @param ... Ignored
#' @export
#' @method print summary_multiScaleR
#' @return Invisibly returns the input \code{summary_multiScaleR} object
print.summary_multiScaleR <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)

  cat("\n\nKernel used:\n")
  cat(x$kernel)

  cat("\n\n***** Optimized Scale of Effect *****\n\n")
  print(x$opt_scale)
  hard_radius <- x$hard_radius
  if (length(hard_radius)) {
    cat("\nRadius (map units) for: ", paste(hard_radius, collapse = ", "), ".\n", sep = "")
    cat("Hard-radius effects change in steps as cells enter the buffer. ")
    cat("Hessian standard errors are omitted; use `summary(x, profile = TRUE)` ")
    cat("for profile-likelihood limits.\n")
  }
  scale_ci_method <- attr(x$opt_scale, "interval_method")
  if (!is.null(scale_ci_method) && any(scale_ci_method != "wald")) {
    cat("\nProfile-likelihood confidence limits were used for `sigma` where available;\n")
    cat("reported standard errors remain Hessian-based approximations.\n")
  }
  cat("\n\n  ==================================== ")

  if(!is.null(x$opt_shape)){
    cat("\n\n***** Optimized Kernel Shape *****\n\n")
    print(x$opt_shape)
    cat("\n\n  ==================================== ")
  }

  .msr_print_distance(dist_1d = x$opt_dist,
                      dist_2d = x$opt_dist_2d,
                      prob = x$prob,
                      hard_radius = hard_radius)
  cat("\n\n  ==================================== ")


  cat('\n\n *****     Fitted Model Summary     *****\n\n')
  if(any(grepl("^unmarked", class(x$fitted_mod)))) {
    print(x$fitted_mod)
  } else {
    print(summary(x$fitted_mod))
  }

  # Warning Messages --------------------------------------------------------


  if(1 %in% x$warn_message){
    msg <- format_max_distance_warning(x$diagnostics)
    suggest_txt <- if (is.finite(msg$suggested_max_D)) {
      paste0(" to >= ", round(msg$suggested_max_D, 2))
    } else {
      ""
    }
    cat(red("\n WARNING!!!\n",
            msg$headline, "\n",
            "Consider increasing " %+% blue$bold("max_D") %+% " in `kernel_prep`" %+%
              suggest_txt %+% " to ensure accurate estimation of scale.\n\n"))
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
  if (isTRUE(x$diagnostics$sample_size$triggered)) {
    cat("\nSample size: ", x$diagnostics$sample_size$fitted_n,
        " of ", x$diagnostics$sample_size$prepared_n,
        " prepared observations entered the fitted model.\n", sep = "")
  }
  .msr_print_next_run(x$next_run)
  invisible(x)
}

#' Summarize multiScaleR objects
#'
#' Summarizes output from \code{multiScale_optim}.
#'
#' @param object An object of class \code{multiScaleR}.
#' @param profile Logical. If \code{TRUE}, use profile-likelihood confidence
#'   limits for `sigma` when feasible. Defaults to \code{FALSE} so summaries
#'   remain fast; profile results are cached for repeated calls on the same
#'   fitted object during the current R session.
#' @param ... Optional arguments passed to the method (e.g., \code{prob} for cumulative kernel weight threshold).
#'
#' @return An object of class \code{summary_multiScaleR}. Kernel-scale limits
#'   default to Wald intervals. Scale distances, in map units, come in two
#'   forms: \code{opt_dist} (also \code{opt_distance}) holds the
#'   one-dimensional distance from \code{\link{kernel_dist}}, the distance at
#'   which the kernel's cumulative density along a line through the point
#'   reaches \code{prob} (about 1.65 sigma for a Gaussian kernel at 90\%).
#'   \code{opt_dist_2d} holds the two-dimensional distance, the radius of the
#'   circle around the point that holds \code{prob} of the weight of an ideal,
#'   untruncated kernel on a plane (about 2.15 sigma for a Gaussian kernel at
#'   90\%). Both tables have columns
#'   \code{Mean}, \code{2.5\%}, and \code{97.5\%}, and hard-radius covariates
#'   report their radius in both. See \code{\link{kernel_dist}} for when to
#'   use each. For hard-radius landscape and unweighted
#'   surface metrics, Hessian standard errors and Wald limits are shown as
#'   \code{NA}; use \code{profile = TRUE} for profile-likelihood limits.
#'   Profiling errors are reported instead of silently returning invalid limits.
#' @export
#' @method summary multiScaleR
summary.multiScaleR <- function(object, profile = FALSE, ...){

  param_list <- list(...)
  prob <- 0.9

  if (!inherits(object, "multiScaleR")) {
    stop("`object` must inherit from class 'multiScaleR'.")
  }
  if (is.null(object$opt_mod)) {
    stop("`object$opt_mod` is missing; cannot summarize model output.")
  }
  if (is.null(object$scale_est)) {
    stop("`object$scale_est` is missing; cannot summarize optimized scale estimates.")
  }
  if (!is.logical(profile) || length(profile) != 1L || is.na(profile)) {
    stop("`profile` must be TRUE or FALSE.")
  }

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

  tab_scale <- scale_ci_table(object = object,
                              df = df,
                              min_D = object$min_D,
                              names = row.names(object$scale_est),
                              profile = profile)

  if ("prob" %in% names(param_list)) {
    prob <- param_list$prob
  }
  if (!is.numeric(prob) || length(prob) != 1L || is.na(prob) || prob <= 0 || prob >= 1) {
    stop("`prob` must be a single numeric value strictly between 0 and 1.")
  }

  ## DEBUG
  # browser()

  object_profile <- object
  object_profile$profile_scale_est <- tab_scale
  opt_distance <- kernel_dist(object_profile, prob = prob)
  opt_distance_2d <- kernel_dist(object_profile, prob = prob, dimension = "2d")
  hard_radius <- intersect(rownames(tab_scale),
                           .msr_hard_radius_covariates(object$kernel_inputs$scale_vars))
  if (length(hard_radius)) {
    tab_scale[hard_radius, "SE"] <- NA_real_
    interval_method <- attr(tab_scale, "interval_method")
    no_profile <- hard_radius[interval_method[hard_radius] == "wald"]
    if (length(no_profile)) {
      tab_scale[no_profile, c("2.5%", "97.5%")] <- NA_real_
      opt_distance[no_profile, c("2.5%", "97.5%")] <- NA_real_
      opt_distance_2d[no_profile, c("2.5%", "97.5%")] <- NA_real_
    }
  }

  if(!is.null(object$shape_est)){
    tab_shape <- ci_func(object$shape_est,
                         df = df,
                         min_D = NULL,
                         names = row.names(object$shape_est))

    out <- list(opt_scale = tab_scale,
                opt_shape = tab_shape,
                opt_dist = opt_distance,
                opt_distance = opt_distance,
                opt_dist_2d = opt_distance_2d,
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                hard_radius = hard_radius,
                diagnostics = object$diagnostics,
                warn_message = object$warn_message,
                next_run = object$next_run,
                call = object$call)
  } else {
    out <- list(opt_scale = tab_scale,
                opt_shape = NULL,
                opt_dist = opt_distance,
                opt_distance = opt_distance,
                opt_dist_2d = opt_distance_2d,
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                hard_radius = hard_radius,
                diagnostics = object$diagnostics,
                warn_message = object$warn_message,
                next_run = object$next_run,
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
#' @return Invisibly returns the input \code{multiScaleR} object
print.multiScaleR <- function(x, ...){
  cat("\nCall:\n")
  print(x$call)

  cat('\n\nKernel used: \n')
  cat(x$kernel_inputs$kernel)

  hard_radius <- intersect(rownames(x$scale_est),
                           .msr_hard_radius_covariates(x$kernel_inputs$scale_vars))
  display_scale <- x$scale_est
  if (length(hard_radius)) display_scale[hard_radius, "SE"] <- NA_real_
  cat("\n\n***** Optimized Scale of Effect *****\n\n")
  print(display_scale)
  if (length(hard_radius)) {
    cat("\nRadius (map units) for: ", paste(hard_radius, collapse = ", "),
        ". Use `summary(x, profile = TRUE)` for radius intervals.\n", sep = "")
  }

  if(x$kernel_inputs$kernel == 'expow'){
    cat("\n\n***** Optimized Kernel Shape Parameter *****\n\n")
    print(x$shape_est)
  }

  display_distance <- kernel_dist(x)
  display_distance_2d <- kernel_dist(x, dimension = "2d")
  if (length(hard_radius)) {
    display_distance[hard_radius, c("2.5%", "97.5%")] <- NA_real_
    display_distance_2d[hard_radius, c("2.5%", "97.5%")] <- NA_real_
  }
  .msr_print_distance(dist_1d = display_distance,
                      dist_2d = display_distance_2d,
                      prob = 0.9,
                      hard_radius = hard_radius)

  cat("\n  ==================================== ")
  cat("\n\n ***** Fitted Model *****\n")
  print(x$opt_mod)


  # Warning Messages --------------------------------------------------------

  if(1 %in% x$warn_message){
    msg <- format_max_distance_warning(x$diagnostics)
    suggest_txt <- if (is.finite(msg$suggested_max_D)) {
      paste0(" to >= ", round(msg$suggested_max_D, 2))
    } else {
      ""
    }
    cat(red("\n WARNING!!!\n",
            msg$headline, "\n",
            "Consider increasing " %+% blue$bold("max_D") %+% " in `kernel_prep`" %+%
              suggest_txt %+% " to ensure accurate estimation of scale.\n\n"))
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
  if (isTRUE(x$diagnostics$sample_size$triggered)) {
    cat("\nSample size: ", x$diagnostics$sample_size$fitted_n,
        " of ", x$diagnostics$sample_size$prepared_n,
        " prepared observations entered the fitted model.\n", sep = "")
  }
  .msr_print_next_run(x$next_run)
  invisible(x)
}

#' @title Print method for multiScaleR_data
#' @description Print method for objects of class \code{multiScaleR_data}.
#' @param x A \code{multiScaleR_data} object
#' @param ... Ignored
#' @export
#' @method print multiScaleR_data
#' @return Invisibly returns the input \code{multiScaleR_data} object
print.multiScaleR_data <- function(x, ...){
  cat("\nThere are ")
  cat(paste0(nrow(x$kernel_dat)," observations at ", ncol(x$kernel_dat), ' spatial covariate(s): \n'))
  cat(colnames(x$kernel_dat))
  cat("\n\nThe specified kernel is:\n")
  cat(x$kernel)
  # cat("\n\nSparse Matrix contains: ")
  # cat(paste0(length(x$D@x), ' elements\n'))
  cat("\n\nNumber of elements: \n")
  if (!is.null(x$d_list)) {
    cat(paste0(length(x$d_list[[1]])))
  } else {
    cat("cell-level data not stored (lean mode; `store_cell_data = FALSE`)")
  }
  cat("\nDistance binning:\n")
  if (!is.null(x$binned)) {
    cat(paste0("enabled (", x$binned$nbins, " distance bins per covariate)"))
  } else {
    cat("disabled (exact per-cell evaluation)")
  }
  cat("\nMinimum Distance:\n")
  cat(x$min_D)
  cat("\nMaximum Distance:\n")
  cat(x$max_D)
  cat("\nUnit Conversion:\n")
  cat(x$unit_conv)
  invisible(x)
}

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
.msr_print_distance <- function(dist_1d, dist_2d, prob, hard_radius,
                                dimension = "2d") {
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
  tables <- list("2d" = dist_2d, "1d" = dist_1d)
  labels <- c("2d" = "2D, radius of the circle for the ideal untruncated kernel",
              "1d" = "1D, along a line through the point")
  for (dim in c(dimension, setdiff(names(tables), dimension))) {
    cat("\n", labels[[dim]], if (dim == dimension) " (selected)" else " (comparison)",
        ":\n", sep = "")
    print(tables[[dim]])
  }
  invisible(NULL)
}

#' @title Print a summary of a multiScaleR fit
#' @description Prints the object returned by [summary.multiScaleR()]: the
#'   estimated scales with their intervals, the scale distances under both
#'   definitions, the fitted model, any warnings, and the follow-up
#'   recommendation.
#' @param x A `summary_multiScaleR` object from [summary.multiScaleR()].
#' @param dimension Character, `"2d"` or `"1d"`. Selects which distance table
#'   is printed first. The default, `NULL`, uses the definition stored in `x`;
#'   summaries saved before this choice existed print the 1D table first.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @seealso [summary.multiScaleR()], [kernel_dist()] for the two distance
#'   definitions.
#' @examples
#' data(pts)
#' data(count_data)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
#'                         inputs, verbose = FALSE)
#' s <- summary(fit)
#' print(s, dimension = "1d")   # line-based distances first
#' @export
#' @method print summary_multiScaleR
print.summary_multiScaleR <- function(x, dimension = NULL, ...){
  if (is.null(dimension)) {
    dimension <- if (is.null(x$dimension)) "1d" else x$dimension
  }
  dimension <- match.arg(dimension, c("2d", "1d"))
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
    if (identical(x$backend, "conditional_logit")) {
      # Step-selection fits have no calibrated scale intervals to profile.
      cat("Scale intervals are unavailable for conditional_logit fits; ")
      cat("use `profile_sigma(x)` to view the log-likelihood across radii.\n")
    } else if (is.null(x$diagnostics$radius_profile)) {
      cat("Hessian standard errors are omitted; use `summary(x, profile = TRUE)` ")
      cat("for profile-likelihood limits.\n")
    } else {
      cat("Inspect the candidate likelihoods in the radius_profile diagnostics table.\n")
      if (any(x$diagnostics$refit_hessian$active))
        cat("Smooth standard errors and Wald limits condition on the selected hard radii.\n")
    }
  }
  .msr_print_refit_bounds(x$diagnostics)
  scale_ci_method <- attr(x$opt_scale, "interval_method")
  if (!is.null(scale_ci_method) && any(grepl("^profile", scale_ci_method))) {
    cat("\nProfile-likelihood confidence limits were used for `sigma` where available;\n")
    cat("reported standard errors remain Hessian-based approximations.\n")
    if (any(scale_ci_method == "profile-bound-limited"))
      cat("Some profile limits reach search bounds; those limits are not fully resolved.\n")
  } else if (!is.null(scale_ci_method) && any(scale_ci_method == "log_wald")) {
    cat("\nScale confidence limits use joint observed information on the log scale.\n")
    # joint_mixed reports withheld scales below, even when none has a limit.
    if (!identical(x$backend, "joint_mixed")) {
      if (any(scale_ci_method == "unavailable"))
        cat("Scale SEs failing their own checks are withheld; see diagnostics(x)$profile_glm.\n")
      bounded <- x$diagnostics$profile_glm$scale_se_conditional_on
      if (length(bounded))
        cat("Other scale SEs hold ", paste(bounded, collapse = ", "),
            " fixed at the search bound.\n", sep = "")
    }
  }
  if (identical(x$backend, "joint_mixed"))
    .msr_print_mixed_scale_se(x$diagnostics$joint_mixed)
  cat("\n\n  ==================================== ")

  if(!is.null(x$opt_shape)){
    cat("\n\n***** Optimized Kernel Shape *****\n\n")
    print(x$opt_shape)
    cat("\n\n  ==================================== ")
  }

  .msr_print_distance(dist_1d = if (is.null(x$opt_dist_1d)) x$opt_dist else x$opt_dist_1d,
                      dist_2d = x$opt_dist_2d,
                      prob = x$prob,
                      hard_radius = hard_radius,
                      dimension = dimension)
  cat("\n\n  ==================================== ")


  cat('\n\n *****     Fitted Model Summary     *****\n\n')
  if (identical(x$backend, "conditional_logit"))
    cat("Coefficient uncertainty conditions on fitted kernel scales; scale uncertainty is unavailable.\n")
  if (.msr_analytic_backend(x)) {
    cat("Coefficient summary below is conditional on fitted kernel parameters.\n")
    cat("Use vcov(fit) and predict(fit, se.fit = TRUE) for joint uncertainty.\n")
    held <- x$diagnostics$profile_glm$joint_conditional_on
    if (length(held))
      cat("Joint Wald uncertainty holds ", paste(held, collapse = ", "),
          " fixed at the search bound.\n", sep = "")
    else if (length(x$diagnostics$profile_glm$reasons))
      cat("Joint Wald uncertainty unavailable: ",
          paste(x$diagnostics$profile_glm$reasons, collapse = "; "), ".\n", sep = "")
  }
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
  .msr_print_start_stall(x$diagnostics)
  .msr_print_next_run(x$next_run)
  invisible(x)
}

#' Summarize a multiScaleR fit
#'
#' Summarizes a fit from [multiScale_optim()]: the estimated scales of effect
#' with standard errors and intervals, the distances that contain most of each
#' kernel's weight, and the fitted model at those scales.
#'
#' @param object An object of class `multiScaleR` from [multiScale_optim()].
#' @param profile Logical. If `TRUE`, report profile-likelihood confidence
#'   limits for sigma where feasible instead of Wald limits. Default: `FALSE`,
#'   which keeps summaries fast. See "Profile-likelihood intervals" below.
#' @param dimension Character. `"2d"` (default) reports, as the selected
#'   distance, the radius that holds `prob` of the weight of an ideal kernel on
#'   a plane. `"1d"` restores the original distance along a line. Both
#'   definitions remain in the returned object and the printed output.
#' @param ... Optional arguments. `prob`, a number strictly between 0 and 1,
#'   sets the proportion of cumulative kernel weight used for the distances.
#'   Default: `0.9`.
#'
#' @return An object of class `summary_multiScaleR`, a list that prints as a
#'   report. Its elements:
#'   * `opt_scale`: data frame with one row per optimized covariate and columns
#'     `Mean` (estimated sigma, map units), `SE`, `2.5%`, and `97.5%`. Its
#'     `interval_method` attribute records how each row's limits were found,
#'     such as `"wald"`, `"log_wald"`, `"profile"`, or `"unavailable"`.
#'     Profile limits that reach a search bound are labeled `"profile_bound"`
#'     by the default refit backend and `"profile-bound-limited"` by the
#'     analytic and joint mixed backends.
#'   * `opt_shape`: the same table for the shape parameter under
#'     `kernel = "expow"`, otherwise `NULL`.
#'   * `opt_dist_1d`: distances in map units, with columns `Mean`, `2.5%`, and
#'     `97.5%`, at which a kernel's cumulative weight along a line through the
#'     point reaches `prob` (about 1.65 sigma for a Gaussian kernel at 90%).
#'   * `opt_dist_2d`: the radius, in map units and with the same columns, of
#'     the circle around the point that holds `prob` of the weight of an ideal,
#'     untruncated kernel on a plane (about 2.15 sigma for a Gaussian kernel at
#'     90%).
#'   * `opt_dist` and `opt_distance`: two names for the distance table selected
#'     by `dimension` (2D by default, rather than the 1D definition used
#'     through version 0.7.3).
#'   * `dimension`: the selected definition, `"2d"` or `"1d"`.
#'   * `prob`: the cumulative weight proportion used for the distances.
#'   * `fitted_mod`: the model refit at the estimated scales (`opt_mod` of the
#'     fit).
#'   * `kernel`: the kernel name.
#'   * `hard_radius`: names of covariates whose scale is a hard radius.
#'   * `diagnostics`, `warn_message`, `next_run`, `call`: copied from the fit;
#'     see [multiScaleR-class].
#'   * `backend`: the backend name, present for the analytic and
#'     conditional-logistic backends.
#'
#' @details
#' # Read the summary
#'
#' Sigma is the kernel's scale parameter, and it is not itself a distance you
#' can report. The distance tables translate it into the distance that holds
#' 90 percent (by default) of the kernel's weight, which is usually the number
#' to report as the scale of effect. Hard-radius covariates report their
#' radius in both tables. See [kernel_dist()] and [multiScaleR-kernels] for
#' when to use each definition and for the formulas.
#'
#' Kernel-scale limits default to Wald intervals from the Hessian. For
#' hard-radius landscape and unweighted surface metrics, and every radius
#' estimated under `kernel = "fixed"`, Hessian standard errors and Wald limits
#' are shown as `NA`; use `profile = TRUE` for profile-likelihood limits. An
#' error during profiling is reported rather than returning invalid limits.
#'
#' Set `dimension = "1d"` to keep the original `opt_dist` and `opt_distance`
#' values in existing scripts. For extraction that does not depend on the
#' setting, use `opt_dist_1d` or `opt_dist_2d`. Neither choice changes the
#' fitted sigma, coefficients, or likelihood. The 2D radius is for an ideal
#' continuous kernel; `max_D`, raster edges, missing cells, and cell geometry
#' can change the realized weight in the prepared data.
#'
#' # Profile-likelihood intervals
#'
#' The default refit backend caches profile results during the current R
#' session; the analytic and joint mixed backends recompute them. Set
#' `options(multiScaleR.profile_cache_max_bytes = 64 * 1024^2)` to bound the
#' retained cache entries (64 MiB by default), or set it to zero to disable
#' caching and clear retained entries at the next profile request. Oversized
#' entries and custom refits are computed without the cache. A failure to
#' store a cache key warns and computes the same intervals directly, while
#' numerical profile failures still raise errors. The returned refit scale
#' table records cache availability in its `profile_cache` attribute. The
#' option bounds retained entries, not temporary key construction or the
#' memory used by the numerical profile.
#'
#' # Analytic backends
#'
#' For `backend = "profile_glm"` or `"joint_nb2"`, default scale limits use
#' joint information on the log scale and are omitted when the
#' identifiability checks fail (see [multiScaleR-fitting-methods]).
#' `backend = "joint_mixed"` reports the same log-scale Wald limits from its
#' joint Hessian, one scale at a time. A returned
#' Wald interval describes the local curvature around the selected optimum
#' and can extend beyond the searched scale range. Coverage can fall below the
#' nominal level even for an interior fit, so inspect the likelihood over the
#' admissible range when scale inference matters.
#'
#' With `profile = TRUE`, the analytic and `"joint_mixed"` summaries follow
#' the connected likelihood region around the fitted scale, and an endpoint
#' at a search bound is labeled bound-limited. A `"joint_mixed"` profile
#' refits the native `glmmTMB` model at each candidate sigma and reports
#' limits for every scale, including one whose standard error is withheld
#' (see [multiScaleR-fitting-methods]). For a supported one-scale NB2 fit,
#' [scale_inference()] checks the full range and reports the regions it
#' detects. If a full-range profile suggests a competing optimum, compare
#' additional starts, including both bounds. The printed coefficient summary
#' still conditions on the estimated kernels; use [vcov.multiScaleR()] for the
#' joint covariance.
#' @seealso [multiScale_optim()], [print.summary_multiScaleR()],
#'   [profile_sigma()], [kernel_dist()], and
#'   `vignette("v09-estimates-and-intervals", package = "multiScaleR")`.
#' @examples
#' data(pts)
#' data(count_data)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
#'                         inputs, verbose = FALSE)
#' s <- summary(fit)
#' s$opt_scale     # sigma with Wald limits, in map units
#' s$opt_dist_2d   # radius holding 90 percent of the kernel weight
#'
#' # Profile-likelihood limits instead of Wald limits
#' summary(fit, profile = TRUE)$opt_scale
#' @export
#' @method summary multiScaleR
summary.multiScaleR <- function(object, profile = FALSE,
                                dimension = c("2d", "1d"), ...){

  param_list <- list(...)
  prob <- 0.9
  dimension <- match.arg(dimension)

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
  object <- .msr_refit_bound_view(object)
  if (identical(object$backend, "conditional_logit"))
    object$opt_mod <- .analysis_model(object$opt_mod)

  if(any(class(object$opt_mod) == 'gls')){
    df <- object$opt_mod$dims$N - object$opt_mod$dims$p
    names <- all.vars(formula(object$opt_mod)[-2])

  } else if(any(grepl("^unmarked", class(object$opt_mod)))){
    df <- dim(object$opt_mod@data@y)[1]
    names <- .unmarked_model_predictors(object$opt_mod)

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
  opt_distance <- kernel_dist(object_profile, prob = prob, dimension = "1d")
  opt_distance_2d <- kernel_dist(object_profile, prob = prob, dimension = "2d")
  hard_radius <- intersect(rownames(tab_scale),
                           .msr_hard_radius_covariates(object$kernel_inputs$scale_vars,
                             object$kernel_inputs$kernel, rownames(tab_scale)))
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
    if (.msr_analytic_backend(object))
      tab_shape <- .msr_profile_wald_table(object, shape = TRUE)

    out <- list(opt_scale = tab_scale,
                opt_shape = tab_shape,
                opt_dist = if (dimension == "2d") opt_distance_2d else opt_distance,
                opt_distance = if (dimension == "2d") opt_distance_2d else opt_distance,
                opt_dist_1d = opt_distance,
                opt_dist_2d = opt_distance_2d,
                dimension = dimension,
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
                opt_dist = if (dimension == "2d") opt_distance_2d else opt_distance,
                opt_distance = if (dimension == "2d") opt_distance_2d else opt_distance,
                opt_dist_1d = opt_distance,
                opt_dist_2d = opt_distance_2d,
                dimension = dimension,
                fitted_mod = object$opt_mod,
                prob = prob,
                kernel = object$kernel_inputs$kernel,
                hard_radius = hard_radius,
                diagnostics = object$diagnostics,
                warn_message = object$warn_message,
                next_run = object$next_run,
                call = object$call)
  }


  out$backend <- object$backend
  class(out) <- c('summary_multiScaleR')
  if (sys.nframe() == 1L) print(out)
  out
}

#' @title Print a multiScaleR fit
#' @description Prints a fit from [multiScale_optim()]: the kernel, the
#'   estimated scales with their standard errors, the scale distances, the
#'   fitted model, any warnings, and the follow-up recommendation. Use
#'   [summary.multiScaleR()] for intervals.
#' @param x A `multiScaleR` object.
#' @param dimension Character. `"2d"` (default) prints the ideal planar radius
#'   first; `"1d"` prints the original line-based distance first.
#' @param prob Numeric between 0 and 1 (exclusive). Proportion of cumulative
#'   kernel weight used for the distances. Default: `0.9`.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @seealso [multiScaleR-class], [summary.multiScaleR()].
#' @examples
#' data(pts)
#' data(count_data)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
#'                         inputs, verbose = FALSE)
#' fit                      # same as print(fit)
#' print(fit, prob = 0.5)   # distances holding half of each kernel's weight
#' @export
#' @method print multiScaleR
print.multiScaleR <- function(x, dimension = c("2d", "1d"), prob = 0.9, ...){
  dimension <- match.arg(dimension)
  cat("\nCall:\n")
  print(x$call)

  cat('\n\nKernel used: \n')
  cat(x$kernel_inputs$kernel)

  hard_radius <- intersect(rownames(x$scale_est),
                           .msr_hard_radius_covariates(x$kernel_inputs$scale_vars,
                             x$kernel_inputs$kernel, rownames(x$scale_est)))
  display_fit <- .msr_refit_bound_view(x)
  display_scale <- display_fit$scale_est
  if (length(hard_radius)) display_scale[hard_radius, "SE"] <- NA_real_
  cat("\n\n***** Optimized Scale of Effect *****\n\n")
  print(display_scale)
  .msr_print_refit_bounds(display_fit$diagnostics)
  .msr_print_mixed_scale_se(x$diagnostics$joint_mixed)
  if (length(hard_radius)) {
    cat("\nRadius (map units) for: ", paste(hard_radius, collapse = ", "), ". ", sep = "")
    if (identical(x$backend, "conditional_logit")) {
      cat("Use `profile_sigma(x)` to view the log-likelihood across radii.\n")
    } else if (is.null(x$diagnostics$radius_profile)) {
      cat("Use `summary(x, profile = TRUE)` for radius intervals.\n")
    } else {
      cat("Use `diagnostics(x)$radius_profile` for candidate likelihoods.\n")
    }
  }

  if(x$kernel_inputs$kernel == 'expow'){
    cat("\n\n***** Optimized Kernel Shape Parameter *****\n\n")
    print(display_fit$shape_est)
  }

  display_distance <- kernel_dist(x, prob = prob, dimension = "1d")
  display_distance_2d <- kernel_dist(x, prob = prob, dimension = "2d")
  if (length(hard_radius)) {
    display_distance[hard_radius, c("2.5%", "97.5%")] <- NA_real_
    display_distance_2d[hard_radius, c("2.5%", "97.5%")] <- NA_real_
  }
  .msr_print_distance(dist_1d = display_distance,
                      dist_2d = display_distance_2d,
                      prob = prob,
                      hard_radius = hard_radius,
                      dimension = dimension)

  cat("\n  ==================================== ")
  cat("\n\n ***** Fitted Model *****\n")
  if (.msr_analytic_backend(x)) {
    cat("Conditional GLM at the optimized kernel parameters.\n")
    held <- if (isTRUE(x$uncertainty$available)) NULL else x$uncertainty$joint_conditional_on
    if (length(held))
      cat("Joint Wald uncertainty holds ", paste(held, collapse = ", "),
          " fixed at the search bound.\n", sep = "")
    else if (length(x$uncertainty$reasons))
      cat("Joint Wald uncertainty unavailable: ",
          paste(x$uncertainty$reasons, collapse = "; "), ".\n", sep = "")
  }
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
  .msr_print_start_stall(x$diagnostics)
  .msr_print_next_run(x$next_run)
  invisible(x)
}

#' @title Print prepared multiScaleR inputs
#' @description Prints a short report of a [kernel_prep()] object: the number
#'   of points and covariates, the kernel, the binning, and the distance
#'   bounds.
#' @param x A `multiScaleR_data` object.
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @seealso [multiScaleR_data-class], [kernel_prep()].
#' @examples
#' data(pts)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' inputs
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
  if (!is.null(x$d_list)) {
    cat(paste0(length(x$d_list[[1]])))
  } else {
    cat("cell-level data not stored (lean mode; `store_cell_data = FALSE`)")
  }
  cat("\nDistance binning:\n")
  if (!is.null(x$binned) || !is.null(x$surface_binned)) {
    nb <- if (!is.null(x$binned)) x$binned$nbins else x$surface_binned$nbins
    cat(paste0("enabled (", nb, " distance bins per covariate)"))
    if (!is.null(x$surface_binned))
      cat(paste0("\nBinned weighted surface metrics: ",
                 paste(x$surface_binned$covariates, collapse = ", ")))
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

# A refit search that never left its starting values reports them as
# estimates. Repeat the stored explanation wherever the fit is reported.
.msr_print_start_stall <- function(diagnostics) {
  stall <- diagnostics$start_stall
  if (isTRUE(stall$triggered)) {
    cat(red("\n WARNING!!!\n", stall$message, "\n\n"))
  }
  invisible(NULL)
}

# joint_mixed scale SEs: name any withheld scale and its reason. Fits saved
# before per-scale SEs existed have no `scale_se_reasons` and print nothing.
.msr_print_mixed_scale_se <- function(jm) {
  reasons <- jm$scale_se_reasons
  if (is.null(reasons)) return(invisible(NULL))
  withheld <- reasons[!is.na(reasons) & reasons != "available"]
  if (length(withheld))
    cat("\nScale SE withheld for ", paste0(names(withheld), " (", withheld, ")",
        collapse = ", "), ". See diagnostics(x)$joint_mixed.\n", sep = "")
  if (length(jm$scale_se_conditional_on))
    cat("Other scale SEs hold ", paste(jm$scale_se_conditional_on, collapse = ", "),
        " fixed at the search bound.\n", sep = "")
  cat("Coefficient and prediction uncertainty condition on the selected scales.\n")
  invisible(NULL)
}

# Print only on explicit reporting, or when the fit's verbose switch permits it.
.msr_print_refit_bounds <- function(diagnostics) {
  rh <- diagnostics$refit_hessian
  if (isTRUE(rh$bounds$active_boundary)) {
    p <- rh$bounds$parameters
    hit <- p$parameter[p$smooth & (p$at_lower | p$at_upper)]
    cat("\nAt an optimization bound: ", paste(hit, collapse = ", "),
        ". Wald SEs and limits are withheld for ",
        if (length(hit) > 1L) "these parameters" else "this parameter", ".\n", sep = "")
    if (length(rh$conditional_on))
      cat("Wald SEs for the other smooth parameters hold the bounded ones fixed.\n")
    cat("Coefficient uncertainty still conditions on the selected kernels.\n")
  }
  invisible(NULL)
}

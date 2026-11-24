#' @title Convert Sigma into an Effective Distance
#'
#' @description Converts a kernel scale parameter (sigma) into an *effective
#' distance* in map units, the distance that holds a stated proportion of the
#' kernel weight. Use it with a fitted `multiScaleR` object to report the
#' estimated scales of effect as distances, or with kernel parameters you
#' supply to explore how far a kernel reaches before fitting. The kernels, the
#' two distance definitions, and their formulas are described in
#' [multiScaleR-kernels].
#'
#' @param model A fitted `multiScaleR` object from [multiScale_optim()]. When
#'   supplied, distances and 95% confidence limits are returned for every
#'   optimized covariate, and any `sigma`, `kernel`, or `shape` passed through
#'   `...` is ignored with a warning. When omitted, supply `sigma` and
#'   `kernel` (and `shape` for `"expow"`) through `...`.
#' @param prob Numeric between 0 and 1 (exclusive). Proportion of the kernel
#'   weight that the reported distance must hold. Default: `0.9`. Larger
#'   values give longer distances. Under `dimension = "1d"` the proportion is
#'   of the weight along a line through the point. Under `"2d"` it is of the
#'   weight on a plane.
#' @param dimension Character. How the distance is measured. `"2d"` (default)
#'   returns the radius of the circle around the point that holds `prob` of
#'   the ideal kernel weight on an unbounded plane. `"1d"` returns the original
#'   line-based distance, which encloses `prob` of the area under the kernel
#'   curve along a transect through the point. For distance-weighted kernels
#'   the `"2d"` radius is longer than the `"1d"` distance. Hard-radius
#'   covariates in a fitted model report the same radius under both. See
#'   Details.
#' @param ... Kernel parameters used when `model` is not supplied:
#'   * `sigma`: positive numeric scalar. The kernel scale parameter, in the
#'     map units of the points and rasters passed to [kernel_prep()]. For the
#'     Gaussian kernel it is the standard deviation, for the negative
#'     exponential kernel the decay length, and for the fixed kernel the
#'     buffer radius.
#'   * `kernel`: character. One of `"gaussian"`, `"exp"` (negative
#'     exponential), `"expow"` (exponential power), or `"fixed"` (fixed-radius
#'     buffer). Required.
#'   * `shape`: positive numeric scalar. Shape of the exponential power
#'     kernel. Required when `kernel = "expow"` and ignored otherwise.
#'   * `beta`: deprecated former name of `shape`. It is still accepted, with a
#'     warning, and will be removed in a future release. Supplying both
#'     `shape` and `beta` with different values is an error.
#'
#' @return
#' When `model` is supplied, a data frame with one row per optimized covariate
#' (row names are the covariate names) and three columns, all in map units and
#' rounded to two decimal places:
#' * `Mean`: the effective distance at the estimated sigma.
#' * `2.5%`: the distance at the lower 95% limit of sigma.
#' * `97.5%`: the distance at the upper 95% limit of sigma.
#'
#' For a hard-radius covariate (a landscape metric, an unweighted surface
#' metric, or any covariate fitted with `kernel = "fixed"`), the row reports
#' the radius itself under either `dimension`. Limits are `NA` when no
#' interval is available (see Details). A covariate whose sigma estimate is
#' not finite returns `NaN`.
#'
#' When kernel parameters are supplied through `...`, a single numeric value,
#' the distance in map units that holds `prob` of the kernel weight under
#' the chosen `dimension`, rounded to two decimal places.
#'
#' @details
#' # Read the distances
#'
#' Sigma is a kernel parameter, not a radius. An effective distance puts the
#' estimate in map units you can compare across kernels and studies. With
#' sigma = 300 m, a Gaussian kernel holds 90% of its planar weight within
#' about 644 m (`dimension = "2d"`), while the line-based 90% distance is about
#' 493 m (`dimension = "1d"`). Use `"2d"` when you describe the area a kernel
#' draws on, and always state the definition and the value of `prob` you
#' report. The formulas for each kernel are tabulated in
#' [multiScaleR-kernels].
#'
#' Both definitions describe the ideal continuous kernel. They do not account
#' for truncation at `max_D`, raster edges, missing cells, or cell geometry,
#' so the share of weight among the cells actually extracted around a point
#' can differ. If the `"2d"` radius approaches or exceeds `max_D`, increase
#' `max_D` in [kernel_prep()] and refit the model.
#'
#' # Reproduce distances from earlier versions
#'
#' Versions 0.7.4 and later report `dimension = "2d"` distances by default in
#' `kernel_dist()`, `summary()`, `print()`, and kernel plots. Earlier versions
#' reported the line-based distance. The switch changes reporting only. It
#' does not change fitted sigma values, kernel weights, or likelihoods. To
#' reproduce distances reported by earlier versions, set the dimension
#' explicitly:
#'
#' ```
#' kernel_dist(fit, dimension = "1d")
#' summary(fit, dimension = "1d")
#' print(fit, dimension = "1d")
#' plot(fit, dimension = "1d")
#' plot_kernel(sigma = 300, kernel = "gaussian", dimension = "1d")
#' ```
#'
#' Simulation extents in [sim_dat()], the `pct_wt` window in
#' [kernel_scale.raster()], and the `max_D` warning of [multiScale_optim()]
#' keep their line-based definitions.
#'
#' # Confidence limits from a fitted model
#'
#' Each limit of sigma is converted to a distance with the same formula as the
#' estimate. By default the sigma limits are Wald limits from the
#' Hessian-based standard errors. When the fitted object stores its own
#' interval table in `profile_scale_est` (see the
#' [fitted-object help page][multiScaleR-class]), those
#' limits are used instead. `summary(fit, profile = TRUE)` converts
#' profile-likelihood limits to distances in its own output, as described in
#' [summary.multiScaleR()]. For exponential-power kernels, the shape is held
#' at its estimate during this conversion, so the limits do not include
#' uncertainty in the shape. These rules depend on the fitting route:
#'
#' * With `backend = "profile_glm"` or `"joint_nb2"` and an estimated
#'   exponential-power shape, distance limits use joint log-scale delta-method
#'   uncertainty for both sigma and shape. They remain Wald limits even when
#'   sigma limits are profiled in `summary()`, and they are missing when the
#'   fit fails the joint identifiability checks. Point distances remain
#'   available.
#' * With `backend = "conditional_logit"`, point distances are available but
#'   scale uncertainty has not been calibrated, so distance limits are
#'   missing.
#' * With `backend = "joint_mixed"`, distance limits convert the log-scale
#'   Wald limits of sigma and are missing for a scale whose standard error is
#'   withheld. See [multiScaleR-fitting-methods] for their coverage caution.
#' * With the default refit backend, point distances remain available when a
#'   sigma standard error is not finite, but its Wald limits are missing. When
#'   a smooth parameter sits at a search bound, the whole Wald block is
#'   withheld, while explicitly marked profile limits remain available.
#' * A hard-radius covariate reports its radius directly and has no kernel
#'   shape, including in an exponential-power fit. It has no Hessian-based
#'   interval, but stored profile limits are kept. See the section on smooth
#'   kernels and hard radii in [multiScaleR-kernels].
#'
#' Native coefficient standard errors remain conditional on the selected
#' scales. Manually supplied fixed-kernel parameters still return the
#' requested weight quantile under the chosen `dimension`.
#'
#' @seealso [multiScaleR-kernels] for the kernel formulas and both distance
#'   definitions. [summary.multiScaleR()], [plot.multiScaleR()], and
#'   [plot_kernel()] report and draw the same distances. `vignette("v02-user-guide", package = "multiScaleR")`
#'   converts sigma to effective distances on a worked example.
#'
#' @examples
#' ## Distances from kernel parameters, without a fitted model
#' # A Gaussian kernel with sigma = 300 m: planar radius holding 90% of the weight
#' kernel_dist(sigma = 300, kernel = "gaussian")                     # about 644
#' # The original line-based definition gives a shorter distance
#' kernel_dist(sigma = 300, kernel = "gaussian", dimension = "1d")   # about 493
#' # Other kernels and proportions
#' kernel_dist(sigma = 500, kernel = "gaussian", prob = 0.95)
#' kernel_dist(sigma = 100, kernel = "exp", prob = 0.975)
#' kernel_dist(sigma = 100, kernel = "expow", shape = 1.5, prob = 0.95)
#' kernel_dist(sigma = 100, kernel = "fixed")   # radius holding 90% of the buffer area
#'
#' \donttest{
#' ## Distances from a fitted model
#' data("pts")
#' data("count_data")
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#'
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = hab,
#'                              max_D = 250,
#'                              kernel = "gaussian",
#'                              verbose = FALSE)
#' mod <- glm(y ~ hab, family = poisson, data = count_data)
#' opt <- multiScale_optim(fitted_mod = mod,
#'                         kernel_inputs = kernel_inputs,
#'                         verbose = FALSE)
#'
#' # Planar radius (m) holding 90% of the weight, with limits from the
#' # 95% interval of sigma
#' kernel_dist(model = opt)
#' # Line-based distances, as reported before version 0.7.4
#' kernel_dist(model = opt, dimension = "1d")
#' # Radius holding 95% of the weight
#' kernel_dist(model = opt, prob = 0.95)
#' }
#'
#' @rdname kernel_dist
#' @export
#' @importFrom insight get_df

kernel_dist <- function(model,
                        prob = 0.9,
                        dimension = c("2d", "1d"),
                        ...){
  param_list <- list(...)
  dimension <- match.arg(dimension)
  validate_scalar_numeric(prob,
                          "prob",
                          lower = 0,
                          upper = 1,
                          inclusive_lower = FALSE,
                          inclusive_upper = FALSE)

  if(!missing("model")){
    if(!inherits(model, "multiScaleR")){
      stop("Provide a fitted `multiScaleR` model object")
    }
  }

  if(!missing("model")){
    if(length(param_list) >= 1){
      warning("Calculating fitted scale relationship; Ignoring specified `sigma` and/or `shape` parameters")
    }

    if(!missing("model")){
      # ci_ <- summary(model)$opt_scale

      opt_mod <- .analysis_model(model$opt_mod)

      if(any(class(opt_mod) == 'gls')){
        df <- opt_mod$dims$N - opt_mod$dims$p
        names <- all.vars(formula(opt_mod)[-2])

      } else if(any(grepl("^unmarked", class(opt_mod)))){
        df <- dim(opt_mod@data@y)[1]
        names <- .unmarked_model_predictors(opt_mod)

      } else {
        df <- get_df(opt_mod, type = "residual")
        names <- all.vars(formula(opt_mod)[-2])
      }

      model <- .msr_refit_bound_view(model)
      if(!is.null(model$profile_scale_est)){
        ci_ <- model$profile_scale_est
      } else {
        ci_ <- ci_func(model$scale_est,
                       df = df,
                       min_D = model$min_D,
                       names = row.names(model$scale_est))
      }

      # browser()

      dist_list <- vector('list', nrow(ci_))
      hard_radius <- .msr_hard_radius_covariates(model$kernel_inputs$scale_vars,
        model$kernel_inputs$kernel, rownames(ci_))
      if (length(hard_radius)) {
        hard <- intersect(rownames(ci_), hard_radius)
        ci_[hard, "SE"] <- NA_real_
        methods <- attr(ci_, "interval_method")
        no_profile <- if (is.null(model$profile_scale_est)) hard else
          if (is.null(methods)) character() else
            hard[is.na(methods[hard]) | methods[hard] %in% c("wald", "unavailable")]
        ci_[no_profile, c("2.5%", "97.5%")] <- NA_real_
      }
      for(i in 1:nrow(ci_)){
        is_landscape <- rownames(ci_)[[i]] %in% hard_radius
        if (isTRUE(is_landscape) && is.finite(ci_[i, 1])) {
            scale_mn <- ci_[i, 1]
            scale_l <- ci_[i, 3]
            scale_u <- ci_[i, 4]
        } else if (is.finite(ci_[i, 1])) {
            shape_i <- if (!is.null(model$shape_est)) model$shape_est[i,1] else NULL
          # wt_mn <- scale_type_r(d = d,
          #                       kernel = model$kernel_inputs$kernel,
          #                       sigma = ci_[i, 1],
          #                       shape = model$shape_est[i,1],
          #                       output = 'wts')
          #
          #
          # wt_l <- scale_type_r(d = d,
          #                      kernel = model$kernel_inputs$kernel,
          #                      sigma = ci_[i,3],
          #                      shape = model$shape_est[i,1],
          #                      output = 'wts')
          #
          # wt_u <- scale_type_r(d = d,
          #                      kernel = model$kernel_inputs$kernel,
          #                      sigma = ci_[i,4],
          #                      shape = model$shape_est[i,1],
          #                      output = 'wts')
          #
          # scale_mn <- wtd.Ecdf(d, weights = wt_mn)
          # scale_mn <- round(scale_mn$x[which(scale_mn$ecdf > prob)[1]], digits = 2)
          #
          # scale_l <- wtd.Ecdf(d, weights = wt_l)
          # scale_l <- round(scale_l$x[which(scale_l$ecdf > prob)[1]], digits = 2)
          #
          # scale_u <- wtd.Ecdf(d, weights = wt_u)
          # scale_u <- round(scale_u$x[which(scale_u$ecdf > prob)[1]], digits = 2)

          scale_mn <- k_dist(sigma = ci_[i, 1],
                             prob = prob,
                             kernel = model$kernel_inputs$kernel,
                             beta = shape_i,
                             dimension = dimension)
          scale_l <- k_dist(sigma = ci_[i, 3],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = shape_i,
                            dimension = dimension)
          scale_u <- k_dist(sigma = ci_[i, 4],
                            prob = prob,
                            kernel = model$kernel_inputs$kernel,
                            beta = shape_i,
                            dimension = dimension)
          interval_method <- attr(ci_, "interval_method")
          if (!is.finite(ci_[i, 2]) &&
              (is.null(interval_method) || interval_method[i] == "wald"))
            scale_l <- scale_u <- NA_real_
          if (.msr_analytic_backend(model) &&
              model$kernel_inputs$kernel == "expow") {
            # Radius depends on both scale and shape. Mapping a scale interval
            # through a fixed shape omits their joint estimation uncertainty.
            scale_l <- scale_u <- NA_real_
            if (model$uncertainty$available) {
              h <- 1e-5
              g_shape <- (log(k_dist(ci_[i, 1], prob, "expow",
                                     beta = shape_i * exp(h), dimension = dimension)) -
                            log(k_dist(ci_[i, 1], prob, "expow",
                                       beta = shape_i * exp(-h), dimension = dimension))) / (2 * h)
              idx <- c(i, nrow(ci_) + i)
              V <- model$uncertainty$vcov[idx, idx, drop = FALSE]
              se_log <- sqrt(drop(c(1, g_shape) %*% V %*% c(1, g_shape)))
              scale_l <- scale_mn * exp(-stats::qnorm(.975) * se_log)
              scale_u <- scale_mn * exp(stats::qnorm(.975) * se_log)
            }
          }
        } else {
          scale_mn <- NaN
          scale_l <- NaN
          scale_u <- NaN
        }


        dist_list[[i]] <- data.frame(mn = round(scale_mn, digits = 2),
                                     l = round(scale_l, digits = 2),
                                     u = round(scale_u, digits = 2))
      }
      dist_out <- do.call(rbind, dist_list)
      rownames(dist_out) <- rownames(ci_)
      colnames(dist_out) <- colnames(ci_)[c(1,3,4)]
    }
    return(dist_out)
  } else if(length(param_list) >= 1){
    sig_ <- param_list$sigma
    # `shape` is the current name; `beta` is a deprecated alias (one release).
    # Exact names only, so `$` partial matching cannot pick up another argument.
    shp_arg <- .msr_resolve_shape(param_list[["shape"]], param_list[["beta"]])
    shp_ <- shp_arg$shape
    kern <- param_list$kernel

    if(is.null(sig_)){
      stop('\nA value for `sigma` must be provided!\n')
    }
    validate_scalar_numeric(sig_, "sigma", positive = TRUE)
    if(is.null(kern)){
      stop('\nYou must specify `kernel` function; See Details\n')
    }
    kern <- match.arg(kern, c("gaussian", "exp", "expow", "fixed"))
    if(kern == 'expow' & is.null(shp_)){
      stop('\nBoth a `sigma` and `shape` parameter must be specified when using the `expow` kernel; See Details\n')
    }

    if (kern == "expow") {
      # Name the argument the caller actually used.
      validate_scalar_numeric(shp_, shp_arg$name, positive = TRUE)
    }

    # d <- seq(1, round(sig_*1000,0), length.out = round(sig_*1000,0))
    # wt <- scale_type_r(d = d,
    #                    kernel = kern,
    #                    sigma = sig_,
    #                    shape = shp_,
    #                    output = 'wts')
    #
    # mx <- wtd.Ecdf(d, weights = wt)
    # mx <- round(mx$x[which(mx$ecdf > 0.999)[1]], digits = -2)
    #
    # d <- seq(1, mx, length.out = 100)
    # wt <- scale_type_r(d = d,
    #                    kernel = kern,
    #                    sigma = sig_,
    #                    shape = shp_,
    #                    output = 'wts')
    #
    # scale_d <- round(d[which(wtd.Ecdf(d, weights = wt)$ecdf > prob)[1]], 2)

    scale_d <- k_dist(sigma = sig_,
                      prob = prob,
                      kernel = kern,
                      beta = shp_,
                      dimension = dimension)

    return(round(scale_d, digits = 2))
  } else {
    stop("Parameters not correctly specified to calculate distance. See Details and try again.")
  }
}


#' @title Distance at Cumulative Kernel Proportion
#' @description Compute the distance at which a given cumulative proportion of
#'   kernel weight is reached for several kernel types, measured either along a
#'   line through the point (one-dimensional) or as the radius of a circle
#'   around it (two-dimensional).
#' @param sigma Numeric. Scale parameter. For Gaussian and exponential, this is standard deviation or decay length. For expow, this is the kernel bandwidth.
#' @param prob Numeric. Desired cumulative proportion (e.g., 0.95).
#' @param kernel Character. One of "gaussian", "exp", "expow", or "fixed".
#' @param beta Numeric. Shape parameter for exponential power kernel. Ignored unless kernel = "expow".
#' @param dimension Character. \code{"2d"} (default) for the radius of the
#'   circle holding \code{prob} of ideal planar weight; \code{"1d"} for the
#'   original distance along a line through the point.
#' @return Numeric distance in map units at which the cumulative kernel weight
#'   reaches \code{prob} under the selected dimensional definition.
#' @noRd
#' @importFrom stats qnorm qgamma
k_dist <- function(sigma, prob = 0.95, kernel = c("gaussian", "exp", "expow", "fixed"),
                   beta = NULL, dimension = c("2d", "1d")) {
  kernel <- match.arg(kernel)
  dimension <- match.arg(dimension)
  if (prob <= 0 || prob >= 1) stop("prob must be between 0 and 1")
  if (!is.finite(sigma)) return(sigma)
  if (kernel == "expow") {
    if (is.null(beta)) stop("beta must be specified for exponential power kernel")
    if (beta <= 0) stop("beta must be positive")
  }

  if (dimension == "2d") {
    # Planar radius holding `prob` of the total weight. Integrating the weight
    # w(d) over rings of circumference 2 * pi * d gives closed forms:
    #   Gaussian, w = exp(-d^2 / (2 sigma^2)): 1 - exp(-r^2 / (2 sigma^2))
    #   Exponential power, w = exp(-(d / sigma)^beta): the cumulative weight is
    #     a gamma CDF with shape 2 / beta evaluated at (r / sigma)^beta
    #   Negative exponential: exponential power with beta = 1
    #   Fixed buffer, w = 1 for d < sigma: area fraction (r / sigma)^2
    # Each assumes the kernel is not truncated by `max_D`.
    return(switch(kernel,
                  gaussian = sigma * sqrt(-2 * log(1 - prob)),
                  exp      = sigma * qgamma(prob, shape = 2),
                  expow    = sigma * qgamma(prob, shape = 2 / beta)^(1 / beta),
                  fixed    = sigma * sqrt(prob)))
  }

  if (kernel == "gaussian") {
    return(qnorm((1 + prob)/2, mean = 0, sd = sigma))

  } else if (kernel == "exp") {
    return(-sigma * log(1 - prob))

  } else if (kernel == "expow") {
    # Exact line-based quantile, including tails beyond the old 10-sigma bracket.
    return(sigma * qgamma(prob, shape = 1 / beta)^(1 / beta))

  } else if (kernel == "fixed") {
    # Proportion of mass in 1D uniform kernel: step function
    return(sigma * (prob))  # Assume sigma is the full extent
  }
}


#' Resolve the exponential-power shape from `shape` and its deprecated alias
#'
#' `beta` was the public name of the exponential-power shape in
#' `kernel_dist()` and `plot_kernel()`. It collides with the usual symbol for
#' regression coefficients, so `shape` replaces it. `beta` remains accepted
#' for one release.
#'
#' @param shape Value supplied as `shape`, or `NULL`.
#' @param beta Value supplied as the deprecated `beta`, or `NULL`.
#' @return A list with `shape` (the value to use, possibly `NULL`) and `name`
#'   (the argument name the caller used, for error messages).
#' @keywords internal
#' @noRd
.msr_resolve_shape <- function(shape = NULL, beta = NULL) {
  if (is.null(beta)) return(list(shape = shape, name = "shape"))
  if (!is.null(shape) && !isTRUE(all.equal(shape, beta))) {
    stop("`shape` and the deprecated `beta` were both supplied with different ",
         "values. Supply only `shape`.", call. = FALSE)
  }
  warning("`beta` is deprecated in favor of `shape` and will be removed in a ",
          "future release. Use `shape` instead.", call. = FALSE)
  list(shape = beta, name = if (is.null(shape)) "beta" else "shape")
}

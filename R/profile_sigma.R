#' Profile Model Fit Across Sigma Parameter Space
#'
#' Evaluates the model log-likelihood and AICc at a series of sigma values
#' spanning the optimization range for each spatial covariate. This provides
#' a diagnostic view of the likelihood surface and helps assess whether the
#' optimized sigma is at a clear minimum or on a flat plateau.
#'
#' @param x A fitted \code{multiScaleR} object returned by
#'   \code{\link{multiScale_optim}}.
#' @param n_pts Positive integer (>= 3). Number of sigma values to evaluate
#'   for each covariate along the profile grid. More points give a smoother
#'   profile at the cost of additional model refits. Default: \code{10}.
#'   Ignored when \code{sigma_values} is supplied.
#' @param metric Character. Metric to display on the y-axis of the profile.
#'   One of \code{"AICc"} (default) or \code{"LL"} (log-likelihood). The
#'   optimized sigma minimizes negative log-likelihood, so the \code{"LL"}
#'   profile should peak at the optimized value and \code{"AICc"} should
#'   show a minimum.
#' @param verbose Logical. Print per-covariate progress messages. Default:
#'   \code{TRUE}.
#' @param n_cores Optional positive integer. Number of CPU cores to use when
#'   profiling covariates in parallel. Parallel profiling is applied across
#'   covariates with a PSOCK cluster. Default: \code{NULL} (serial profiling).
#' @param spacing Character. Spacing of the automatically generated sigma grid.
#'   \code{"log"} (default) concentrates evaluation points at small sigma
#'   values where the likelihood surface typically changes more rapidly.
#'   \code{"linear"} provides equal spacing across the range.
#' @param sigma_values Optional positive numeric vector of sigma values (in the
#'   original projection units) to evaluate directly. When supplied,
#'   \code{spacing}, \code{sigma_range}, and \code{n_pts} are ignored. Duplicate
#'   values are removed and remaining values are sorted before profiling.
#'   Must contain at least 3 unique values.
#' @param sigma_range Optional positive numeric vector of length 2 giving the
#'   minimum and maximum sigma values (in projection units) for the generated
#'   profile grid. Defaults to the optimization range stored in \code{x}
#'   (\code{min_D} to \code{max_D}). Ignored when \code{sigma_values} is
#'   supplied.
#'
#' @return A list of class \code{sigma_profile} containing:
#' \describe{
#'   \item{profiles}{A data frame with columns \code{variable}, \code{sigma},
#'     \code{LL}, and \code{AICc}.}
#'   \item{opt_sigma}{A named numeric vector of the optimized sigma for each
#'     covariate.}
#'   \item{metric}{The metric used for profiling.}
#'   \item{spacing}{The profile grid type: \code{"log"}, \code{"linear"}, or
#'     \code{"custom"}.}
#'   \item{sigma_grid}{The sigma values evaluated for each covariate.}
#' }
#'
#' @details
#' For each spatial covariate, sigma is varied across a sequence of candidate
#' values, while all other sigma values are held at their optimized values. By
#' default this is a log-spaced sequence from the minimum to maximum distance
#' considered during optimization. Users can instead request linear spacing with
#' \code{spacing = "linear"} or supply exact values with \code{sigma_values}.
#' At each evaluation point the model is refit and the log-likelihood extracted.
#' AICc is computed from the log-likelihood, the number of regression parameters
#' (including sigma), and the number of observations.
#'
#' Log-spacing concentrates evaluation points at smaller sigma values where the
#' likelihood surface often changes more rapidly, and spaces them out at larger
#' sigma values where the surface tends to be flatter.
#'
#' Linear spacing can be useful when the profile needs equal resolution across a
#' specific sigma interval. User-supplied \code{sigma_values} are sorted and
#' duplicate values are removed before profiling to avoid redundant refits.
#'
#' @seealso \code{\link{plot.sigma_profile}}, \code{\link{multiScale_optim}}
#'
#' @examples
#' \donttest{
#' data('pts')
#' data('count_data')
#' hab <- terra::rast(system.file('extdata',
#'                    'hab.tif', package = 'multiScaleR'))
#'
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = hab,
#'                              max_D = 250,
#'                              kernel = 'gaussian')
#'
#' mod <- glm(y ~ hab,
#'            family = poisson,
#'            data = count_data)
#'
#' opt <- multiScale_optim(fitted_mod = mod,
#'                         kernel_inputs = kernel_inputs)
#'
#' ## Profile sigma
#' prof <- profile_sigma(opt)
#' plot(prof)
#'
#' ## More evaluation points
#' prof <- profile_sigma(opt, n_pts = 20)
#' plot(prof)
#'
#' ## Linearly spaced values between explicit limits
#' prof <- profile_sigma(opt, n_pts = 8, spacing = "linear",
#'                       sigma_range = c(25, 250))
#' plot(prof)
#'
#' ## User-supplied sigma values
#' prof <- profile_sigma(opt, sigma_values = c(25, 50, 100, 150, 250))
#' plot(prof)
#'
#' ## Profile log-likelihood instead of AICc
#' prof <- profile_sigma(opt, metric = "LL")
#' plot(prof)
#' }
#'
#' @export
#' @importFrom insight get_parameters n_obs
#' @importFrom stats logLik
profile_sigma <- function(x,
                          n_pts = 10,
                          metric = c("AICc", "LL"),
                          verbose = TRUE,
                          n_cores = NULL,
                          spacing = c("log", "linear"),
                          sigma_values = NULL,
                          sigma_range = NULL) {

  if (!inherits(x, "multiScaleR")) {
    stop("`x` must be a fitted `multiScaleR` object.", call. = FALSE)
  }
  validate_scalar_logical(verbose, "verbose")
  if (!is.null(n_cores)) {
    validate_scalar_numeric(n_cores, "n_cores", integerish = TRUE, positive = TRUE)
  }
  metric <- match.arg(metric)
  spacing <- match.arg(spacing)

  # Extract components from the multiScaleR object
  kernel_inputs <- x$kernel_inputs
  opt_context   <- x$opt_context
  join_by       <- x$join_by
  kernel        <- kernel_inputs$kernel
  unit_conv     <- kernel_inputs$unit_conv

  # Optimized sigma values (on the scaled parameter space used by optim)
  scale_est <- x$scale_est
  covs      <- rownames(scale_est)
  n_covs    <- length(covs)
  opt_par   <- scale_est$Mean / unit_conv  # back to scaled space

  if (n_covs < 1) {
    stop("`x` does not contain any optimized sigma values to profile.", call. = FALSE)
  }

  # Shape parameters if expow kernel
  if (kernel == "expow" && !is.null(x$shape_est)) {
    opt_shape <- x$shape_est$Mean
  } else {
    opt_shape <- NULL
  }

  # Sigma search range (on original scale)
  min_D <- x$min_D
  max_D <- x$max_D

  # Compute K (number of model parameters) for AICc
  mod <- x$opt_mod
  k_base <- .msr_parameter_count(mod)
  n <- .msr_model_nobs(mod)
  # Add sigma parameters (and shape if expow)
  k_total <- k_base + n_covs
  if (kernel == "expow") k_total <- k_total + n_covs

  if (!is.null(sigma_values)) {
    validate_numeric_vector(sigma_values, "sigma_values", positive = TRUE)
    sigma_seq <- sort(unique(sigma_values))

    if (length(sigma_seq) < 3) {
      stop("`sigma_values` must contain at least 3 unique values.", call. = FALSE)
    }

    grid_type <- "custom"
  } else {
    validate_scalar_numeric(n_pts, "n_pts", integerish = TRUE, positive = TRUE)

    if (n_pts < 3) {
      stop("`n_pts` must be >= 3 for a meaningful profile.", call. = FALSE)
    }

    if (is.null(sigma_range)) {
      sigma_range <- c(min_D, max_D)
      sigma_range[1] <- max(sigma_range[1], 1e-3)
    } else {
      validate_numeric_vector(sigma_range, "sigma_range", length_ = 2,
                              positive = TRUE)
      sigma_range <- sort(sigma_range)
    }

    if (sigma_range[1] >= sigma_range[2]) {
      stop("`sigma_range` must contain two distinct values.", call. = FALSE)
    }

    sigma_seq <- if (spacing == "log") {
      exp(seq(log(sigma_range[1]), log(sigma_range[2]), length.out = n_pts))
    } else {
      seq(sigma_range[1], sigma_range[2], length.out = n_pts)
    }

    grid_type <- spacing
  }
  n_grid <- length(sigma_seq)

  # Only the profiled covariate's sigma changes across the grid; every other
  # covariate sits at its optimum and produces an identical kernel-weighted
  # column on every grid point. Compute those optimum columns once and reuse
  # them, recomputing just the profiled column per point. This avoids the
  # dominant cost (the per-point covariate weighting) for the held covariates.
  # If the baseline cannot be built (e.g. a lean object missing data a covariate
  # needs), `base_cov_w` is NULL and each point falls back to a full evaluation.
  base_par <- if (is.null(opt_shape)) opt_par else c(opt_par, opt_shape)
  base_cov_w <- tryCatch(
    .msr_kernel_cov_w(par = base_par,
                      d_list = kernel_inputs$d_list,
                      cov_df = kernel_inputs$raw_cov,
                      kernel = kernel,
                      opt_context = opt_context),
    error = function(e) NULL
  )

  profile_one_covariate <- function(j) {
    ll_vec   <- numeric(n_grid)
    aicc_vec <- numeric(n_grid)
    cov_name <- covs[j]

    for (i in seq_len(n_grid)) {
      # Build parameter vector: replace j-th sigma, keep others at optimum
      par_i <- opt_par
      par_i[j] <- sigma_seq[i] / unit_conv  # scale to optim space

      # Append shape parameters if expow
      if (!is.null(opt_shape)) {
        par_i <- c(par_i, opt_shape)
      }

      # Reuse the cached optimum columns and recompute only the profiled
      # covariate. A NULL baseline (or any failure here) means `cov_w_i` stays
      # NULL and `kernel_scale_fn()` recomputes every covariate as before.
      cov_w_i <- NULL
      if (!is.null(base_cov_w)) {
        cov_w_i <- tryCatch({
          col_j <- .msr_kernel_cov_w(
            par = par_i,
            d_list = kernel_inputs$d_list,
            cov_df = kernel_inputs$raw_cov,
            kernel = kernel,
            opt_context = opt_context,
            covariates = cov_name
          )
          cw <- base_cov_w
          cw[, cov_name] <- col_j[, cov_name]
          cw
        }, error = function(e) NULL)
      }

      # Evaluate negative log-likelihood
      neg_ll <- tryCatch(
        kernel_scale_fn(
          par         = par_i,
          d_list      = kernel_inputs$d_list,
          cov_df      = kernel_inputs$raw_cov,
          kernel      = kernel,
          fitted_mod  = opt_context$fitted_mod,
          join_by     = join_by,
          mod_return  = NULL,
          opt_context = opt_context,
          cov_w       = cov_w_i
        ),
        error = function(e) NA_real_
      )

      # kernel_scale_fn returns negative LL; convert
      ll_val <- -neg_ll

      # AICc: -2*LL + 2*K + 2*K*(K+1)/(n-K-1)
      aic_val <- -2 * ll_val + 2 * k_total
      aicc_val <- if ((n - k_total - 1) > 0) {
        aic_val + (2 * k_total * (k_total + 1)) / (n - k_total - 1)
      } else {
        NA_real_
      }

      ll_vec[i]   <- ll_val
      aicc_vec[i] <- aicc_val
    }

    data.frame(
      variable = covs[j],
      sigma    = sigma_seq,
      LL       = ll_vec,
      AICc     = aicc_vec,
      stringsAsFactors = FALSE
    )
  }

  # Profile each covariate
  results <- vector("list", n_covs)
  use_parallel <- !is.null(n_cores) && n_cores > 1 && n_covs > 1

  if (use_parallel) {
    if (isTRUE(verbose)) {
      cat(sprintf("Profiling %d covariates in parallel using %d cores...\n",
                  n_covs, n_cores))
    }

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cluster_prep(mod, cl)

    results_try <- try(
      parallel::parLapply(cl = cl,
                          X = seq_len(n_covs),
                          fun = profile_one_covariate),
      silent = TRUE
    )

    if (inherits(results_try, "try-error")) {
      err_msg <- attr(results_try, "condition")
      if (!is.null(err_msg)) {
        stop("Parallel sigma profiling failed: ",
             conditionMessage(err_msg),
             call. = FALSE)
      }
      stop("Parallel sigma profiling failed due to an unknown error.",
           call. = FALSE)
    }

    results <- results_try
  } else {
    if (!is.null(n_cores) && n_cores > 1 && n_covs == 1 && isTRUE(verbose)) {
      cat("Only one covariate was available; using serial profiling.\n")
    }
    for (j in seq_len(n_covs)) {
      if (isTRUE(verbose)) {
        cat(sprintf("Profiling '%s' (%d of %d)...\n", covs[j], j, n_covs))
      }
      results[[j]] <- profile_one_covariate(j)
    }
  }

  profiles <- do.call(rbind, results)
  rownames(profiles) <- NULL

  out <- list(
    profiles   = profiles,
    opt_sigma  = stats::setNames(scale_est$Mean, covs),
    metric     = metric,
    spacing    = grid_type,
    sigma_grid = sigma_seq
  )
  class(out) <- "sigma_profile"

  if (isTRUE(verbose)) cat("Done.\n")

  return(out)
}


#' Plot Sigma Profile
#'
#' Plots the profiled log-likelihood or AICc across sigma values for each
#' spatial covariate. The optimized sigma is marked with a vertical red line.
#'
#' @param x An object of class \code{sigma_profile} from \code{\link{profile_sigma}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return A list of \code{ggplot} objects (one per covariate), returned
#'   invisibly. Plots are printed as a side effect.
#'
#' @seealso \code{\link{profile_sigma}}
#'
#' @export
#' @method plot sigma_profile
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline geom_point xlab ylab ggtitle
#' @importFrom cowplot theme_cowplot
plot.sigma_profile <- function(x, ...) {

  if (!inherits(x, "sigma_profile")) {
    stop("`x` must be a `sigma_profile` object.", call. = FALSE)
  }

  profiles  <- x$profiles
  opt_sigma <- x$opt_sigma
  metric    <- x$metric
  covs      <- names(opt_sigma)

  y_col  <- if (metric == "LL") "LL" else "AICc"
  y_label <- if (metric == "LL") "Log-likelihood" else "AICc"

  plot_list <- lapply(covs, function(v) {
    df <- profiles[profiles$variable == v, ]

    y_vals <- df[[y_col]]
    df$y_val <- y_vals

    p <- ggplot(df, aes(x = sigma, y = y_val)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_vline(xintercept = opt_sigma[v],
                 linetype = "dashed",
                 color = "red",
                 linewidth = 0.8) +
      xlab(expression(sigma)) +
      ylab(y_label) +
      ggtitle(v) +
      theme_cowplot()

    p
  })

  names(plot_list) <- covs
  lapply(plot_list, print)
  invisible(plot_list)
}

utils::globalVariables(c("sigma", "y_val"))

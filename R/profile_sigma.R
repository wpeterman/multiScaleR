#' Profile Model Fit Across Sigma Parameter Space
#'
#' Evaluates the model log-likelihood and AICc at a series of sigma values
#' spanning the optimization range for each spatial covariate. This provides
#' a diagnostic view of the likelihood surface and helps assess whether the
#' optimized sigma is at a clear minimum or on a flat plateau.
#'
#' @param x A fitted \code{multiScaleR} object.
#' @param n_pts Integer. Number of sigma values to evaluate for each covariate.
#'   Default is 10.
#' @param metric Character. Which metric to profile: \code{"AICc"} (default)
#'   or \code{"LL"} (log-likelihood).
#' @param verbose Logical. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list of class \code{sigma_profile} containing:
#' \describe{
#'   \item{profiles}{A data frame with columns \code{variable}, \code{sigma},
#'     \code{LL}, and \code{AICc}.}
#'   \item{opt_sigma}{A named numeric vector of the optimized sigma for each
#'     covariate.}
#'   \item{metric}{The metric used for profiling.}
#' }
#'
#' @details
#' For each spatial covariate, sigma is varied across a log-spaced sequence from
#' the minimum to maximum distance considered during optimization, while all
#' other sigma values are held at their optimized values. At each evaluation
#' point the model is refit and the log-likelihood extracted. AICc is computed
#' from the log-likelihood, the number of regression parameters (including
#' sigma), and the number of observations.
#'
#' Log-spacing concentrates evaluation points at smaller sigma values where the
#' likelihood surface often changes more rapidly, and spaces them out at larger
#' sigma values where the surface tends to be flatter.
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
                          verbose = TRUE) {

  if (!inherits(x, "multiScaleR")) {
    stop("`x` must be a fitted `multiScaleR` object.", call. = FALSE)
  }
  validate_scalar_numeric(n_pts, "n_pts", integerish = TRUE, positive = TRUE)
  validate_scalar_logical(verbose, "verbose")
  metric <- match.arg(metric)

  if (n_pts < 3) {
    stop("`n_pts` must be >= 3 for a meaningful profile.", call. = FALSE)
  }

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
  if (any(grepl("^unmarked", class(mod)))) {
    k_base <- length(all.vars(formula(mod@formula))) + 1
    n <- dim(mod@data@siteCovs)[1]
  } else {
    k_base <- nrow(get_parameters(mod)) + 1
    n <- n_obs(mod)
  }
  # Add sigma parameters (and shape if expow)
  k_total <- k_base + n_covs
  if (kernel == "expow") k_total <- k_total + n_covs

  # Log-spaced sigma sequence from min_D to max_D
  # Use a small offset above zero to avoid log(0)
  log_min <- log(max(min_D, 1e-3))
  log_max <- log(max_D)
  sigma_seq <- exp(seq(log_min, log_max, length.out = n_pts))

  # Profile each covariate
  results <- vector("list", n_covs)

  for (j in seq_len(n_covs)) {
    if (isTRUE(verbose)) {
      cat(sprintf("Profiling '%s' (%d of %d)...\n", covs[j], j, n_covs))
    }

    ll_vec   <- numeric(n_pts)
    aicc_vec <- numeric(n_pts)

    for (i in seq_len(n_pts)) {
      # Build parameter vector: replace j-th sigma, keep others at optimum
      par_i <- opt_par
      par_i[j] <- sigma_seq[i] / unit_conv  # scale to optim space

      # Append shape parameters if expow
      if (!is.null(opt_shape)) {
        par_i <- c(par_i, opt_shape)
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
          opt_context = opt_context
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

    results[[j]] <- data.frame(
      variable = covs[j],
      sigma    = sigma_seq,
      LL       = ll_vec,
      AICc     = aicc_vec,
      stringsAsFactors = FALSE
    )
  }

  profiles <- do.call(rbind, results)
  rownames(profiles) <- NULL

  out <- list(
    profiles  = profiles,
    opt_sigma = stats::setNames(scale_est$Mean, covs),
    metric    = metric
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

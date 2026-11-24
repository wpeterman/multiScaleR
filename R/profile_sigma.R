#' Profile Model Fit Across Sigma
#'
#' Evaluates the log-likelihood and AICc of a fitted model at a grid of sigma
#' values (or hard radii) for each optimized covariate. The resulting curve
#' shows whether an estimate sits at a clear optimum or on a flat plateau,
#' which tells you how well the data identify the scale of effect. Plot the
#' result with [plot.sigma_profile()].
#'
#' @param x A fitted `multiScaleR` object returned by [multiScale_optim()].
#' @param n_pts Positive integer, at least 3. Number of sigma values in the
#'   generated grid for each covariate. More values give a smoother curve at
#'   the cost of more model refits. Ignored when `sigma_values` is supplied.
#'   Default: `10`.
#' @param metric Character. Which quantity [plot.sigma_profile()] draws:
#'   `"AICc"` (default) or `"LL"` (log-likelihood). Both are always computed.
#'   The log-likelihood is highest, and AICc lowest, at the best-supported
#'   sigma.
#' @param verbose Logical. Print a progress message for each covariate.
#'   Default: `TRUE`.
#' @param n_cores Optional positive integer. Number of CPU cores used to
#'   profile covariates in parallel, one covariate per worker (a PSOCK
#'   cluster). Default: `NULL` (serial). The analytic, conditional-logistic,
#'   and joint mixed backends always profile serially.
#' @param spacing Character. Spacing of the generated grid. `"log"` (default)
#'   places more values at small sigma, where the likelihood usually changes
#'   fastest. `"linear"` spaces them evenly. Ignored when `sigma_values` is
#'   supplied.
#' @param sigma_values Optional positive numeric vector of sigma values, in
#'   map units, to evaluate directly. Duplicates are removed and the values
#'   are sorted, and at least 3 unique values are required. When supplied,
#'   `spacing`, `sigma_range`, and `n_pts` are ignored.
#' @param sigma_range Optional positive numeric vector of length 2 giving the
#'   smallest and largest sigma, in map units, of the generated grid. Default:
#'   the search range of the fit, `min_D` to `max_D`. Ignored when
#'   `sigma_values` is supplied.
#'
#' @return A list of class `sigma_profile` with these elements:
#' * `profiles`: data frame with one row per covariate and grid value, and
#'   columns `variable` (covariate name), `sigma` (map units), `LL`
#'   (log-likelihood of the refit model), and `AICc`. AICc is `NA` when the
#'   number of observations `n` is not larger than `K + 1`, where `K` is the
#'   parameter count.
#' * `opt_sigma`: named numeric vector with the estimated sigma of each
#'   covariate, in map units.
#' * `metric`: the metric chosen for plotting.
#' * `spacing`: the grid type, `"log"`, `"linear"`, or `"custom"` (for
#'   `sigma_values`).
#' * `sigma_grid`: the sigma values evaluated, shared by all covariates.
#'
#' @details
#' # How the profile is computed
#'
#' For each covariate in turn, sigma is set to each grid value while every
#' other sigma (and any shape) stays at its estimate. The model is refit at
#' each value, so its coefficients adjust, and the log-likelihood is
#' recorded. Because the other scales are not re-optimized, a model with
#' several covariates gives a slice through the likelihood surface rather than
#' a full profile likelihood. AICc counts all free likelihood parameters
#' (including ancillary, `unmarked` detection, and variance parameters) plus
#' the estimated sigma and shape parameters, as in [aic_tab()], and uses the
#' model's observation count.
#'
#' **How to read the profile.** A sharp peak in log-likelihood (or trough in
#' AICc) at the estimate means the data identify the scale well. A flat curve
#' over a wide range means many scales fit about equally well, so treat the
#' estimate with caution. A curve still improving at `max_D` suggests the
#' scale may lie beyond the search range. In that case, increase `max_D` in
#' [kernel_prep()] and refit. For a profile-likelihood interval for sigma,
#' use `summary(fit, profile = TRUE)`. See
#' `vignette("v09-estimates-and-intervals", package = "multiScaleR")`.
#'
#' If every refit fails and returns the optimizer's penalty value, profiling
#' stops with an error instead of returning an artificially flat profile.
#' Keep the model formula and data available to the stored model call, or
#' give [multiScale_optim()] a working `refit_fn`.
#'
#' # Profile a hard radius
#'
#' Landscape metrics, unweighted surface metrics, and every covariate fitted
#' with `kernel = "fixed"` have a hard radius, so a cell either counts fully
#' or not at all. The covariate changes only when the radius passes a cell
#' center, so the likelihood changes in steps. The Hessian standard error is
#' therefore not a useful measure of uncertainty, and the default summary
#' shows `NA` for the standard error and Wald interval. Use this profile to
#' inspect the likelihood, and report the estimate as a radius in map units.
#' The plot joins grid values with straight lines, so it does not show the
#' steps between them, and an interval read from a profile cutoff does not by
#' itself establish coverage for a stepped likelihood. Candidate scales that
#' make a covariate undefined at a modeled observation are rejected during
#' fitting, so all likelihoods compare the same observations. Hard-radius
#' covariates can be profiled in exponential-power fits, and their unused
#' shape positions do not affect the profile or the parameter count. See
#' [multiScaleR-kernels].
#'
#' # Analytic, conditional-logistic, and joint mixed fits
#'
#' Fits from `backend = "profile_glm"`, `"joint_nb2"`, `"conditional_logit"`,
#' and `"joint_mixed"` profile serially, and every sigma value must lie
#' within the fitted bounds (`min_D` to `max_D`). For a `"conditional_logit"`
#' fit, each grid value is the log-likelihood of a native
#' `survival::coxph()` refit at that sigma, with the coefficients refit and
#' every other scale held at its estimate. That backend withholds scale
#' uncertainty, so `summary(fit, profile = TRUE)` returns no profile interval
#' for it. The grid describes the likelihood surface only.
#'
#' For a `"joint_mixed"` fit, each grid value is the log-likelihood of a new
#' native `glmmTMB` fit at that sigma, the same refit the backend uses to
#' check each start. The coefficients, the random-intercept standard
#' deviation, and any NB2 dispersion are re-estimated, the random intercept is
#' integrated out with the Laplace approximation, and every other scale stays
#' at its estimate. A refit that fails gives `NA` at that value, and
#' profiling stops if every refit for a covariate fails.
#'
#' @seealso [plot.sigma_profile()], [multiScale_optim()],
#'   [summary.multiScaleR()], [multiScaleR-kernels].
#'
#' @examples
#' \donttest{
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
#' mod <- glm(y ~ hab, family = poisson,
#'            data = data.frame(count_data["y"], kernel_inputs$kernel_dat))
#' opt <- multiScale_optim(fitted_mod = mod,
#'                         kernel_inputs = kernel_inputs,
#'                         verbose = FALSE)
#'
#' ## Default: 10 log-spaced values from min_D to max_D
#' prof <- profile_sigma(opt, verbose = FALSE)
#' # Log-likelihood and AICc at each sigma (m)
#' prof$profiles
#' plot(prof)   # AICc is lowest near the estimated sigma (dashed line)
#'
#' ## Linearly spaced values between explicit limits
#' prof <- profile_sigma(opt, n_pts = 8, spacing = "linear",
#'                       sigma_range = c(25, 250), verbose = FALSE)
#'
#' ## Your own sigma values, drawn as log-likelihood
#' prof <- profile_sigma(opt, sigma_values = c(25, 50, 75, 100, 150, 250),
#'                       metric = "LL", verbose = FALSE)
#' plot(prof)   # log-likelihood peaks near the estimate
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
    opt_shape <- .msr_shape_values(x)
  } else {
    opt_shape <- NULL
  }

  # Sigma search range (on original scale)
  min_D <- x$min_D
  max_D <- x$max_D

  # Compute K (number of model parameters) for AICc
  mod <- x$opt_mod
  n <- .msr_model_nobs(mod)
  k_total <- .msr_selection_parameter_count(x)

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
    # exp(log(x)) can round just outside x (exp(log(20)) < 20). Pin the
    # endpoints so the default grid stays inside the fitted bounds.
    sigma_seq[c(1L, n_pts)] <- sigma_range

    grid_type <- spacing
  }
  n_grid <- length(sigma_seq)

  clogit <- identical(x$backend, "conditional_logit")
  mixed <- identical(x$backend, "joint_mixed")
  if (.msr_analytic_backend(x) || clogit || mixed) {
    if (!is.null(n_cores) && n_cores > 1L)
      stop(x$backend, " likelihood profiling is serial.", call. = FALSE)
    if (any(sigma_seq < min_D | sigma_seq > max_D))
      stop(sprintf(paste0(
        "%s profiles need sigma values within the fitted bounds, %s to %s ",
        "map units. Supply `sigma_values` or `sigma_range` inside them, or ",
        "omit both to profile the whole fitted range."),
        x$backend, format(min_D), format(max_D)), call. = FALSE)
    # A conditional_logit fit has no refit context. Each grid value is a
    # native clogit refit, with the other scales held at their estimates.
    # A joint_mixed grid value is a cold native glmmTMB refit, again with the
    # other scales at their estimates. A refit that fails gives NA there.
    profile_ll <- if (clogit) .msr_clogit_profile_loglik(x) else if (mixed) {
      mixed_profile <- .msr_mixed_profiler(x)
      function(j, z) tryCatch(-mixed_profile$at(j, z, reoptimize = FALSE),
                              error = function(e) NA_real_)
    } else
      function(j, z) -.msr_profile_at(x, j, z, reoptimize = FALSE)
    profiles <- do.call(rbind, lapply(seq_len(n_covs), function(j) {
      if (verbose && (clogit || mixed)) cat("Profiling sigma for", covs[j], "\n")
      ll <- vapply(log(sigma_seq / unit_conv), function(z) profile_ll(j, z),
        numeric(1))
      if (mixed && !any(is.finite(ll)))
        stop(sprintf("Could not profile `%s`: every native refit failed.", covs[j]),
             call. = FALSE)
      data.frame(variable = covs[j], sigma = sigma_seq, LL = ll,
        AICc = if (n > k_total + 1) -2 * ll + 2 * k_total +
          2 * k_total * (k_total + 1) / (n - k_total - 1) else NA_real_)
    }))
    return(structure(list(profiles = profiles,
      opt_sigma = stats::setNames(scale_est$Mean, covs), metric = metric,
      spacing = grid_type, sigma_grid = sigma_seq), class = "sigma_profile"))
  }

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

    if (all(!is.finite(ll_vec) | ll_vec <= -(1e6^10))) {
      first_par <- opt_par
      first_par[j] <- sigma_seq[1] / unit_conv
      if (!is.null(opt_shape)) first_par <- c(first_par, opt_shape)
      detail <- tryCatch({
        kernel_scale_fn(
          par = first_par, d_list = kernel_inputs$d_list,
          cov_df = kernel_inputs$raw_cov, kernel = kernel,
          fitted_mod = opt_context$fitted_mod, join_by = join_by,
          mod_return = TRUE, opt_context = opt_context
        )
        "All likelihood evaluations returned the optimizer penalty."
      }, error = function(e) conditionMessage(e))
      stop(sprintf("Could not profile `%s`: %s", cov_name, detail),
           call. = FALSE)
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


#' Plot a Sigma Profile
#'
#' Draws one panel per covariate from a [profile_sigma()] result. Each panel
#' plots AICc or the log-likelihood (whichever `metric` was chosen in
#' [profile_sigma()]) against sigma, with points at the evaluated values and a
#' dashed red vertical line at the estimated sigma.
#'
#' @param x A `sigma_profile` object from [profile_sigma()].
#' @param ... Not used.
#'
#' @return A named list of `ggplot` objects, one per covariate, returned
#'   invisibly. Each plot is printed as a side effect. Modify a stored plot
#'   and print it again to customize it.
#'
#' @details
#' **How to read the plot.** The x-axis is sigma in map units, and the
#' y-axis is AICc (lower is better) or the log-likelihood (higher is better).
#' A clear trough in AICc, or peak in log-likelihood, at the dashed line means
#' the data identify the scale of effect well. A flat curve means many
#' scales fit about equally well. Points are joined with straight lines, so
#' a hard-radius covariate's stepped likelihood is not shown between grid
#' values. See [profile_sigma()] for how the values are computed.
#'
#' @seealso [profile_sigma()]
#'
#' @examples
#' ## A stored profile from the Quick Start fit (two covariates, 15 values each)
#' prof <- readRDS(system.file("extdata", "cache-quickstart.rds",
#'                             package = "multiScaleR"))$prof
#' # One AICc panel each for land1 and land2. AICc is lowest near the
#' # estimated sigma (dashed line)
#' p <- plot(prof)
#'
#' # The stored ggplot objects can be modified and drawn again
#' p$land1 + ggplot2::labs(x = "sigma (m)")
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


# Build a profile function for a conditional_logit fit. The fitted result
# keeps only prediction state, so rebuild the conditional-logit contract from
# the original model and inputs. The returned function gives the native
# coxph() log-likelihood, with coefficients refit, at log(sigma / unit_conv)
# `z` for scale `j` and every other scale at its estimate.
.msr_clogit_profile_loglik <- function(x) {
  inputs <- x$kernel_inputs
  inputs$min_D <- x$min_D
  inputs$max_D <- x$max_D
  engine <- x$profile_fit$control$kernel_engine
  if (is.null(engine)) engine <- "reference"
  contract <- .msr_clogit_contract(x$fitted_mod_original, inputs, engine)
  kp_hat <- unname(x$profile_fit$par)
  function(j, z) {
    kp <- kp_hat
    kp[j] <- z
    as.numeric(stats::logLik(.msr_clogit_native(contract, kp)))
  }
}

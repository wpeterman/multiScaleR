.empty_diagnostics <- function() {
  list(
    max_distance = NULL,
    sigma_precision = NULL,
    shape_precision = NULL
  )
}


.max_distance_diagnostic <- function(scale_D, max_D) {
  est_D <- scale_D[, 1]
  ratio <- max_D / est_D
  triggered <- ratio < 2

  list(
    code = "max_distance",
    triggered = any(triggered, na.rm = TRUE),
    threshold_ratio = 2,
    variables = names(est_D)[which(triggered)],
    effective_distance = est_D,
    max_D = max_D,
    ratio = ratio,
    suggested_max_D = max(est_D * 2, na.rm = TRUE)
  )
}


.precision_diagnostic <- function(est, code) {
  ratio <- est[, 2] / est[, 1]
  triggered <- ratio >= 0.5

  list(
    code = code,
    triggered = any(triggered, na.rm = TRUE),
    threshold_ratio = 0.5,
    variables = row.names(est)[which(triggered)],
    estimate = est[, 1],
    se = est[, 2],
    se_to_mean = ratio
  )
}


.next_run_recommendation <- function(scale_est,
                                     shape_est,
                                     max_D,
                                     diagnostics) {
  start_sigma <- scale_est[, "Mean"]
  names(start_sigma) <- row.names(scale_est)

  recommended_max_D <- max_D
  max_distance_diag <- diagnostics$max_distance
  if (isTRUE(max_distance_diag$triggered) &&
      is.finite(max_distance_diag$suggested_max_D) &&
      max_distance_diag$suggested_max_D > recommended_max_D) {
    recommended_max_D <- max_distance_diag$suggested_max_D
  }

  start_par_sigma <- start_sigma / recommended_max_D
  names(start_par_sigma) <- names(start_sigma)

  start_shape <- NULL
  start_par <- start_par_sigma
  if (!is.null(shape_est)) {
    start_shape <- shape_est[, "Mean"]
    names(start_shape) <- row.names(shape_est)
    start_shape[attr(shape_est, "not_applicable")] <- 1
    start_par <- c(start_par_sigma, start_shape)
    names(start_par) <- c(
      paste0(names(start_par_sigma), "_sigma"),
      paste0(names(start_shape), "_shape")
    )
  }

  flags <- vapply(diagnostics,
                  function(x) is.list(x) && isTRUE(x$triggered),
                  logical(1))
  flags <- names(flags)[flags]

  action <- if (length(flags) == 0) {
    "Use the optimized parameters as efficient starting values if refitting."
  } else {
    "Refit from the optimized parameters and address the flagged diagnostics."
  }

  # The RAM-based parallel-worker recommendation is intentionally NOT computed
  # here. Querying system RAM shells out to an external process (PowerShell on
  # Windows) and rebuilds the optimization context, which adds latency to every
  # `multiScale_optim()` call and can block under load. Call
  # `estimate_multiscale_ram()` directly when a worker-count recommendation is
  # wanted. `n_cores` is kept (as NULL) for backward compatibility.
  list(
    max_D = recommended_max_D,
    start_sigma = start_sigma,
    start_shape = start_shape,
    start_par = start_par,
    n_cores = NULL,
    flags = flags,
    action = action
  )
}


.msr_log_sequence <- function(lower, upper, n) {
  if (n <= 1 || isTRUE(all.equal(lower, upper))) {
    return(rep(lower, max(1, n)))
  }

  exp(seq(log(lower), log(upper), length.out = n))
}


.msr_objective_value <- function(par,
                                 fitted_mod,
                                 kernel_inputs,
                                 join_by,
                                 opt_context) {
  value <- try(
    kernel_scale_fn(
      par = par,
      d_list = kernel_inputs$d_list,
      cov_df = kernel_inputs$raw_cov,
      kernel = kernel_inputs$kernel,
      fitted_mod = fitted_mod,
      join_by = join_by,
      opt_context = opt_context
    ),
    silent = TRUE
  )

  if (inherits(value, "try-error") || length(value) != 1 || !is.finite(value)) {
    return(Inf)
  }

  as.numeric(value)
}


.msr_screen_run <- function(par,
                            lwr,
                            uppr,
                            fitted_mod,
                            kernel_inputs,
                            join_by,
                            opt_context,
                            screen_maxit) {
  out <- try(
    optim(
      par = par,
      fn = kernel_scale_fn,
      hessian = FALSE,
      lower = lwr,
      upper = uppr,
      method = "L-BFGS-B",
      fitted_mod = fitted_mod,
      d_list = kernel_inputs$d_list,
      cov_df = kernel_inputs$raw_cov,
      control = list(maxit = screen_maxit),
      kernel = kernel_inputs$kernel,
      join_by = join_by,
      opt_context = opt_context
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    return(list(
      par = par,
      value = .msr_objective_value(
        par = par,
        fitted_mod = fitted_mod,
        kernel_inputs = kernel_inputs,
        join_by = join_by,
        opt_context = opt_context
      ),
      convergence = 1L
    ))
  }

  out
}


.msr_prescreen_start <- function(par,
                                 n_covs,
                                 lwr,
                                 uppr,
                                 fitted_mod,
                                 kernel_inputs,
                                 join_by,
                                 opt_context,
                                 screen_n_sigma,
                                 screen_n_jitter,
                                 screen_maxit,
                                 screen_jitter_sd,
                                 n_cores,
                                 PSOCK,
                                 verbose) {
  if (isTRUE(verbose)) {
    cat("Prescreening starting values with log-spaced sigma scans.\n")
  }

  sigma_lwr <- lwr[seq_len(n_covs)]
  sigma_uppr <- uppr[seq_len(n_covs)]
  marginal_par <- par

  for (j in seq_len(n_covs)) {
    sigma_grid <- .msr_log_sequence(sigma_lwr[j], sigma_uppr[j], screen_n_sigma)
    obj_vals <- vapply(
      sigma_grid,
      function(candidate) {
        par_j <- marginal_par
        par_j[j] <- candidate
        .msr_objective_value(
          par = par_j,
          fitted_mod = fitted_mod,
          kernel_inputs = kernel_inputs,
          join_by = join_by,
          opt_context = opt_context
        )
      },
      numeric(1)
    )

    finite_idx <- which(is.finite(obj_vals))
    if (length(finite_idx) > 0) {
      marginal_par[j] <- sigma_grid[finite_idx[which.min(obj_vals[finite_idx])]]
    }
  }

  if (screen_n_jitter < 1 || screen_maxit < 1) {
    return(marginal_par)
  }

  if (isTRUE(verbose)) {
    cat("Running short screening optimizations around the prescreened start.\n")
  }

  candidates <- vector("list", screen_n_jitter + 1L)
  candidates[[1]] <- marginal_par

  for (i in seq_len(screen_n_jitter)) {
    candidate <- marginal_par
    candidate[seq_len(n_covs)] <- exp(
      log(marginal_par[seq_len(n_covs)]) +
        stats::rnorm(n_covs, mean = 0, sd = screen_jitter_sd)
    )
    candidate[seq_len(n_covs)] <- pmin(
      pmax(candidate[seq_len(n_covs)], sigma_lwr),
      sigma_uppr
    )
    candidates[[i + 1L]] <- candidate
  }

  candidate_mat <- unique(as.data.frame(do.call(rbind, candidates)))
  candidates <- lapply(
    seq_len(nrow(candidate_mat)),
    function(i) as.numeric(candidate_mat[i, ])
  )

  screen_cores <- if (is.null(n_cores)) 1L else as.integer(n_cores)
  use_parallel_screen <- screen_cores > 1L && length(candidates) > 1L

  run_screen_serial <- function() {
    lapply(
      candidates,
      .msr_screen_run,
      lwr = lwr,
      uppr = uppr,
      fitted_mod = fitted_mod,
      kernel_inputs = kernel_inputs,
      join_by = join_by,
      opt_context = opt_context,
      screen_maxit = screen_maxit
    )
  }

  if (use_parallel_screen) {
    n_screen_cores <- min(screen_cores, length(candidates))
    if (isTRUE(verbose)) {
      cat(sprintf("Running %d screening attempts in parallel using %d cores.\n",
                  length(candidates), n_screen_cores))
    }

    cl <- if (.Platform$OS.type == "unix" && isFALSE(PSOCK)) {
      makeForkCluster(n_screen_cores)
    } else {
      makeCluster(n_screen_cores)
    }
    on.exit(stopCluster(cl), add = TRUE)
    cluster_prep(fitted_mod, cl)

    screened_try <- try(
      parLapply(
        cl = cl,
        X = candidates,
        fun = .msr_screen_run,
        lwr = lwr,
        uppr = uppr,
        fitted_mod = fitted_mod,
        kernel_inputs = kernel_inputs,
        join_by = join_by,
        opt_context = opt_context,
        screen_maxit = screen_maxit
      ),
      silent = TRUE
    )

    if (inherits(screened_try, "try-error")) {
      if (isTRUE(verbose)) {
        cat("Parallel screening failed; retrying screening attempts serially.\n")
      }
      screened <- run_screen_serial()
    } else {
      screened <- screened_try
    }
  } else {
    screened <- run_screen_serial()
  }

  screen_values <- vapply(
    screened,
    function(x) {
      value <- x$value[1]
      if (length(value) != 1 || !is.finite(value)) Inf else as.numeric(value)
    },
    numeric(1)
  )

  if (!any(is.finite(screen_values))) {
    return(marginal_par)
  }

  screened[[which.min(screen_values)]]$par
}

#' @title Multiscale optimization
#'
#' @description Identifies the kernel scale of effect (sigma) for each raster
#' covariate by maximizing the log-likelihood of a fitted statistical model.
#' Repeatedly replaces kernel-weighted covariate values at different scales and
#' refits the model, using \code{optim} (or \code{optimParallel} for parallel
#' execution) with the L-BFGS-B algorithm.
#'
#' @param fitted_mod A fitted model object whose covariates include the
#'   kernel-weighted variables defined in \code{kernel_inputs}. Supported
#'   classes include \code{lm}, \code{glm}, \code{gls} (nlme), and
#'   \code{unmarkedFit} (unmarked). Many other model classes are also supported
#'   via \code{stats::update()} or a custom \code{refit_fn}. A predictor cannot
#'   be named \code{"x"}, which is reserved for internal coordinate data.
#' @param kernel_inputs A list of class \code{"multiScaleR_data"} created by
#'   \code{\link{kernel_prep}}. Must contain elements \code{raw_cov},
#'   \code{d_list}, \code{min_D}, \code{max_D}, \code{unit_conv}, and
#'   \code{kernel}.
#' @param join_by Default: \code{NULL}. For \code{unmarked} models, a data
#'   frame with one unique key row per spatial site in \code{kernel_inputs},
#'   in the same order as those sites. Its named key columns must also occur in
#'   the fitted model's \code{siteCovs}. Use it when model rows repeat sites,
#'   such as one row per site and year. See Details.
#' @param par Optional numeric vector of starting values for the optimizer.
#'   Values must be divided by \code{max_D} to match the internal scaled
#'   parameter space. Length must equal the number of optimized covariates
#'   (or twice that for \code{kernel = "expow"}, where shape parameters follow
#'   sigma parameters). Default: \code{NULL}; starting values are chosen
#'   automatically.
#'   The expow native-refit interface retains one shape position per optimized
#'   variable, including hard radii. A hard radius ignores that position;
#'   provide a finite positive placeholder when supplying \code{par}.
#' @param start_strategy Character. How to choose starting values when
#'   \code{par = NULL}: \code{"single"} (default) uses one shared default
#'   start, while \code{"screen"} performs a low-cost prescreen of sigma
#'   values and short screening optimizations before launching one full
#'   optimization.
#' @param screen_n_sigma Positive integer. Number of log-spaced sigma values
#'   evaluated per covariate during the prescreen when
#'   \code{start_strategy = "screen"}. Default: \code{5}.
#' @param screen_n_jitter Non-negative integer. Number of jittered candidate
#'   starts evaluated with short screening optimizations after the marginal
#'   prescreen. Default: \code{6}. Set to \code{0} to skip the short jittered
#'   runs.
#' @param screen_maxit Positive integer. Maximum iterations used for each short
#'   screening optimization. Default: \code{8}.
#' @param screen_jitter_sd Non-negative numeric scalar. Standard deviation of
#'   multiplicative sigma jitter on the log scale during screening. Default:
#'   \code{0.5}.
#' @param n_cores Positive integer. Number of cores for parallel optimization
#'   via \code{optimParallel}. Default: \code{NULL} (single-threaded). When
#'   \code{start_strategy = "screen"}, the same \code{n_cores} setting is also
#'   used for short screening attempts before the full optimization.
#'   Parallel optimization is beneficial for models with many covariates or
#'   slow log-likelihood evaluations. When \code{n_cores > 1},
#'   \code{multiScale_optim()} checks \code{\link{estimate_multiscale_ram}} and
#'   stops before screened starts or worker creation if the requested workers
#'   exceed the conservative RAM-based limit for the current inputs. Use
#'   \code{n_cores = 1} for serial execution, or rebuild inputs with
#'   \code{kernel_prep(..., bin = TRUE, store_cell_data = FALSE)} for lean mode
#'   when every covariate is kernel-weighted. Landscape and surface metrics
#'   cannot use lean mode. An unavailable RAM estimate does not block the run;
#'   a passing estimate does not guarantee sufficient memory. No worker count
#'   is reduced automatically.
#' @param PSOCK Logical. On Windows, a PSOCK cluster is always used. On Unix,
#'   a FORK cluster is used by default (faster). Set \code{TRUE} to force a
#'   PSOCK cluster on Unix. Default: \code{FALSE}.
#' @param verbose Logical. Print optimization status and diagnostic blocks to the
#'   console. Default: \code{TRUE}.
#' @param refit_fn Optional function for refitting \code{fitted_mod} during
#'   optimization. When \code{NULL} (default), \code{stats::update()} (or the
#'   unmarked equivalent) is used. Provide this only when the default refit
#'   path is insufficient. See Details.
#' @param backend \code{"refit"} (default) preserves the existing model-refit
#'   optimizer. \code{"profile_glm"} opts into the analytic-profile backend for
#'   independent Gaussian, Poisson, binomial, and NB2 observations. See the
#'   formula, data, and uncertainty contracts in Details. \code{"joint_nb2"}
#'   jointly estimates kernel parameters, coefficients, and NB2 dispersion,
#'   including the exact Poisson boundary. It requires an NB2 log-link model.
#'   \code{"joint_mixed"} opts into joint Poisson/NB2 random-intercept fitting
#'   with optional zero inflation; see its restricted contracts in Details.
#'   \code{"conditional_logit"} jointly estimates smooth kernel scales and
#'   conditional-logistic coefficients for one used step per matched choice set.
#'   It accepts retained survival models or amt wrappers; see Details.
#' @param profile_control Named list used by the analytic backends:
#'   For \code{conditional_logit}, accepted entries are \code{n_starts} (6),
#'   \code{starts}, \code{iter.max} (500), \code{gradient_tolerance} (0.001),
#'   and \code{kernel_engine} (\code{"reference"} or \code{"stream"}).
#'   For \code{joint_mixed}, only \code{n_starts} (default 6) and
#'   \code{starts} are accepted. Each row of \code{starts} contains positive
#'   kernel scales in internal distance units (map distance divided by
#'   \code{unit_conv}), in model predictor order. All starts are retained.
#'   For the analytic GLM backends, available controls are
#'   \code{n_starts} (6), \code{starts} (optional positive numeric matrix),
#'   \code{theta_bounds} (\code{c(1e-4, 1e5)}), \code{theta_start} (optional),
#'   \code{gradient_tolerance} (0.001), \code{iter.max} (500), and
#'   \code{eval.max} (1000). Set \code{kernel_engine = "stream"} to opt into
#'   contiguous matrix traversal for shared kernel values and derivatives,
#'   reduced copying of prepared matrices, and cached response-only NB2 terms;
#'   \code{"reference"} is the default. This retains the same starts,
#'   likelihood, and uncertainty checks. The stream
#'   engine uses additional temporary vectors proportional to the number of
#'   observations and can reduce traversal cost for larger binned matrices.
#'   For smooth weighted SQ predictors in \code{conditional_logit}, it also
#'   computes values and scale derivatives in compiled code, reusing work
#'   buffers across sites. Temporary storage is proportional to the largest
#'   retained site buffer. The reference engine and default remain unchanged.
#'   For NB2 models in either analytic backend, set
#'   \code{coefficient_start = "warm"} to reuse
#'   conditional regression coefficients within each search. The default is
#'   \code{"cold"}. Warm continuation resets at every independent start and
#'   falls back to cold initialization if its score or rank checks fail. Final
#'   likelihood checks always use independently initialized coefficients.
#'   Counts and fallback reasons are stored in
#'   \code{profile_fit$coefficient_audit}. All independent starts are retained.
#'   For \code{joint_nb2} only, \code{worker_pool} accepts a reusable PSOCK
#'   cluster created with \code{parallel::makePSOCKcluster()}. Its default,
#'   \code{NULL}, keeps fitting serial. Keep \code{n_cores = 1}; the pool size
#'   determines the number of workers for independent starts. Workers must
#'   load the same installed package code as the controller. Each start keeps
#'   its cold terminal audit, and results retain their original start order.
#'   You own the pool: reuse it across calls and close it with
#'   \code{parallel::stopCluster()} when finished, including after errors.
#'   Pool connections are not stored in the fitted object. Worker coefficient
#'   audits are stored in \code{profile_fit$worker_coefficient_audits};
#'   \code{coefficient_audit} then covers controller-side final evaluation.
#'   Subsequent profile and bootstrap inference remains serial. Separate
#'   processes duplicate data and package memory; inspect
#'   \code{estimate_multiscale_ram()} before creating a pool. Pool startup
#'   also checks the binned payload against 75 percent of physical RAM, allowing
#'   four payload copies and 1 GiB per process for the controller and workers.
#'   This planning estimate is not a guarantee of peak memory or available RAM.
#'   Startup
#'   can outweigh fitting gains for short jobs. Set math-library thread limits
#'   before creating workers to avoid nested parallel computation.
#'   Each row of \code{starts} contains kernel scales
#'   divided by \code{unit_conv}, followed by shapes for \code{"expow"}, in
#'   model predictor order. Explicit starts replace the deterministic scale
#'   starts; \code{par}, if supplied, is prepended. NB2 dispersion starts cycle
#'   through the initial theta, 0.3, 2, and 20 within the specified bounds.
#'   Theta bounds must be positive, increasing, and at most 1e6.
#'   For \code{"joint_nb2"}, the default theta bounds are \code{c(1e-4, Inf)};
#'   the upper bound must remain infinite, and the lower bound must be at least
#'   \code{1e-8}. Each kernel start is tried at the exact Poisson boundary and
#'   at an interior dispersion, giving 12 fits with the default six kernel
#'   starts. \code{iter.max} limits each L-BFGS-B search; \code{eval.max} applies
#'   to subsequent profile/inference refits, which use conditional fitting.
#'
#' @return A list of class \code{"multiScaleR"} containing:
#' \describe{
#'   \item{\code{scale_est}}{Data frame with one row per optimized covariate and
#'     two columns: \code{Mean} (optimized sigma on the original projection
#'     scale) and \code{SE} (Hessian-based standard error). \code{SE} is
#'     \code{NA} for hard-radius landscape and unweighted surface metrics,
#'     whose likelihood changes in steps. Row names are covariate names.}
#'   \item{\code{shape_est}}{Data frame of the same structure as
#'     \code{scale_est} for the shape parameter, or \code{NULL} when
#'     \code{kernel != "expow"}. Hard-radius variables have NA shape and shape
#'     SE because shape is inapplicable; their ignored shape positions are
#'     excluded from likelihood parameter counts.}
#'   \item{\code{optim_results}}{The raw list returned by \code{stats::optim}
#'     or \code{optimParallel::optimParallel}, including \code{par},
#'     \code{value} (negative log-likelihood), \code{hessian}, and convergence
#'     codes. Also contains \code{par_unscale} (sigma values back on the
#'     projection scale) and \code{hessian_unscale}.}
#'   \item{\code{opt_mod}}{The refitted model at the optimized sigma values.
#'     This is the model object to use for inference, prediction, and
#'     \code{\link{plot_marginal_effects}}.}
#'   \item{\code{fitted_mod_original}}{The original \code{fitted_mod} passed
#'     by the user, stored for reference.}
#'   \item{\code{min_D}}{Numeric. Lower bound for sigma during optimization.}
#'   \item{\code{max_D}}{Numeric. Upper bound for sigma during optimization.}
#'   \item{\code{kernel_inputs}}{The \code{kernel_inputs} list (minus
#'     \code{min_D} and \code{max_D}, which are stored separately).}
#'   \item{\code{scl_params}}{Named list with \code{mean} and \code{sd} vectors
#'     for each covariate: the centering and scaling parameters from the
#'     optimized kernel data. Used by \code{\link{kernel_scale.raster}} when
#'     \code{scale_center = TRUE}.}
#'   \item{\code{join_by}}{The \code{join_by} data frame, or \code{NULL}.}
#'   \item{\code{opt_context}}{Internal optimization context object storing the
#'     model-class-specific refit logic. Retained for use by
#'     \code{\link{profile_sigma}}.}
#'   \item{\code{profile_scale_est}}{Profile-likelihood CI data frame, or
#'     \code{NULL} until \code{\link{profile_sigma}} has been run.}
#'   \item{\code{diagnostics}}{List of diagnostic objects; see
#'     \code{\link{diagnostics}}.}
#'   \item{\code{next_run}}{List of recommended values for a follow-up fit.
#'     Includes \code{max_D} for the next \code{\link{kernel_prep}} call,
#'     optimized \code{start_sigma} values in map units, optional
#'     \code{start_shape} values for \code{kernel = "expow"}, and
#'     \code{start_par} on the internal scaled parameter space expected by
#'     \code{multiScale_optim}. \code{n_cores} is \code{NULL}; a conservative
#'     parallel worker-count suggestion is available on demand from
#'     \code{\link{estimate_multiscale_ram}} (it is no longer computed
#'     automatically, as that query can be slow and shells out to an external
#'     process).}
#'   \item{\code{warn_message}}{Integer vector of triggered warning codes:
#'     1 = max-distance, 2 = sigma precision, 3 = shape precision.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{profile_fit}, \code{uncertainty}, \code{dispersion_est}}{Present
#'     for both analytic GLM backends. These retain transformed parameter estimates,
#'     starts and convergence scores, the profile and joint information
#'     matrices, identifiability diagnostics, and NB2 theta. The ordinary
#'     \code{opt_mod} and its coefficient summary condition on optimized kernel
#'     parameters. Use \code{\link{vcov.multiScaleR}} and
#'     \code{\link{predict.multiScaleR}} for joint uncertainty.}
#'   \item{\code{mixed_fit}}{For \code{joint_mixed}, retained kernel engine,
#'     log-scale estimates and bounds in internal distance units, independent
#'     starts, every start outcome and warning, selected-start index, and
#'     required prediction variables. \code{uncertainty$available} is
#'     \code{FALSE}; \code{next_run} and \code{opt_context} are \code{NULL}.}
#' }
#' @details
#' \strong{Opt-in analytic GLM contracts}
#'
#' The opt-in \code{backend = "conditional_logit"} uses the exact finite-choice-set
#' likelihood with one used and at least one available step per set. Use a
#' retained \code{survival::clogit} model (\code{model = TRUE}) or an amt
#' \code{fit_clogit} wrapper, Gaussian or exponential inputs, and one
#' bare \code{strata(choice_id)} term. IDs must identify sets uniquely across
#' animals and bursts. Row names must match prepared point IDs. Explicitly
#' prepare complete sets before fitting; omitted missing rows, subsets, weights,
#' robust clustering and special Cox terms are rejected. Numeric kernel predictors
#' may enter native R interaction, crossing and nesting terms, including crosses
#' with factors or other optimized predictors. Fixed transformations and offsets
#' retain their native meaning and saved training values. Optimized
#' kernel-weighted surface SQ predictors are also supported when the matching
#' raw buffers are retained. These use the weighted root-mean-square deviation,
#' not hard-radius sample SD. Their scales remain independent, and prediction
#' requires the same source and metric specifications and training normalization.
#' Fixed-radius metrics can be supplied as fixed predictors. Their original
#' prepared numeric columns must be retained in the native model frame.
#' New predictions require these columns explicitly unless
#' \code{include_fixed = TRUE} verifies and reuses their original normalization.
#' Optimization of hard-radius
#' landscape/surface metrics, learned dynamic bases, arbitrary or
#' nondifferentiable functions, fixed context inside a dynamic expression, and
#' kernel-dependent offsets remain unsupported. Kernel-only scalar arithmetic
#' and supported smooth base functions use the shared response contract.
#' A quoted choice-set ID is still one bare variable; special-call text inside
#' a quoted predictor name is not interpreted as a Cox special term.
#'
#' Conditional-logistic fitting uses serial independent L-BFGS-B starts, with
#' cold coefficient initialization and a cold native survival refit at every
#' terminal point. Eligibility requires convergence, a small projected score,
#' agreement in likelihood and coefficients, well-conditioned coefficient
#' information and agreement with native covariance. Failed attempts and warnings
#' remain explicit. After the original joint point passes the native agreement
#' audits, cold native coefficients may refine its terminal score if the objective
#' does not increase and the score improves. Both scores and the refinement flag
#' are retained; the eligibility threshold is unchanged. Failed attempts
#' remain in \code{profile_fit$start_results}; if all fail, the error condition
#' contains \code{start_results}. \code{profile_fit$starts} stores log-scale
#' starting values and \code{winner} identifies the selected record. Final
#' \code{opt_mod} retains the native conditional fit and amt metadata when
#' supplied. If the original fitting call no longer exposes its raw data,
#' retained model-frame values support fixed transformations, factors and
#' compound offsets without attempting to invert them. Cold native refits
#' retain the original prediction transformations, including learned polynomial
#' bases. New-data prediction still requires the formula's raw inputs.
#' The original amt proposal distributions remain unchanged; apply amt's
#' movement-distribution updates to \code{opt_mod} when appropriate for the
#' movement terms included in the model. Its internal prediction contract omits duplicate fitting data to
#' reduce serialized size; the original model remains in
#' \code{fitted_mod_original}. Updating the conditional fit requires explicit data with already
#' standardized kernel covariates; otherwise rerun the original multiscale fit.
#' Scale SEs and joint prediction uncertainty are withheld. Native
#' coefficient covariance conditions on the selected scales and does not correct
#' temporal dependence. Movement terms must match the available-step proposal;
#' the backend does not add them, validate a sampling interval, or calibrate
#' population inference. Keep available steps fixed throughout optimization.
#'
#' The opt-in \code{backend = "joint_mixed"} requires optional packages RTMB,
#' TMB, and glmmTMB. It jointly estimates kernel scales and model parameters,
#' integrating random intercepts by the Laplace approximation. It supports
#' log-link Poisson and NB2 glmmTMB fits with one conditional random intercept,
#' constant dispersion, and optional fixed-effects zero inflation. Use
#' Gaussian or negative exponential kernels, native R crossings, nesting or
#' factor interactions of numeric spatial predictors, retained model row names
#' matching prepared point IDs, and
#' \code{n_cores = 1}. Random slopes, multiple grouping terms, spatial terms in
#' zero-inflation or dispersion formulas, non-unit weights, priors, parameter
#' maps, and REML are unsupported. Raw formula variables must be recoverable
#' from the fitted model; use formula \code{offset()} terms. Each start uses a
#' fresh native fit, a joint Hessian check, and an independent native terminal
#' refit with likelihood agreement within \code{1e-5}. Failed starts remain in
#' \code{mixed_fit$start_results}; if all fail, the error condition contains
#' \code{start_results}. No Hessian checks are bypassed.
#'
#' For \code{joint_mixed}, \code{opt_mod} is the final native glmmTMB fit.
#' \code{mixed_fit} records log scales, bounds, starts, all outcomes, and the
#' selected start. Narrow raw prediction data retain fixed-transform inputs
#' for effect curves after temporary native-refit data are no longer available.
#' \code{diagnostics$joint_mixed} reports eligibility and scale
#' boundary flags. Scale uncertainty has not been calibrated: scale SEs and
#' legacy scale Hessians are \code{NA}, and joint intervals are unavailable.
#' Native coefficient and prediction covariance condition on the selected
#' scales. Use \code{predict(fit, se.fit = TRUE, uncertainty = "conditional")}
#' or \code{predict(fit$opt_mod, ...)}. This backend does not extend the
#' single-predictor NB2 GLM contract of \code{scale_inference()}.
#'
#' The analytic optimization backends support optimized kernel means and kernel-weighted
#' surface SQ (root mean square deviation), including combinations of these
#' predictors. Weighted SQ uses retained raw cell buffers with Gaussian or
#' exponential kernels; estimation of other surface metrics' radii requires
#' \code{backend = "refit"}. Declared fixed-radius metrics can enter the native
#' numeric design without adding a scale parameter. Their original prepared
#' values and normalization must remain recoverable from the native model frame.
#' Kernel means require \code{bin = TRUE}; surface-only models need explicit
#' preparation row IDs and raw buffers, without binned summaries.
#' The two analytic GLM backends require a smooth kernel,
#' and a \code{glm} or \code{MASS::glm.nb} object retaining its model frame.
#' For \code{"profile_glm"}, supported links are Gaussian identity, Poisson log, binomial logit, and
#' NB2 log. NB2 has variance \eqn{\mu + \mu^2/\theta}; theta is estimated
#' alongside kernel parameters, even when the initial GLM used a fixed theta.
#' This backend is experimental. Development simulations found scale-interval
#' undercoverage even among fits passing curvature checks. For single-kernel
#' NB2 models, \code{\link{scale_inference}} examines the full likelihood range
#' and optionally supplies bootstrap calibration; neither method guarantees
#' reliable inference under weak identification.
#' For \code{"profile_glm"}, regression coefficients are solved by IRLS at each candidate. A full-rank
#' NB2 fit whose iterations fail receives a coefficient-likelihood refinement
#' using the analytic score and observed information. Gaussian residual
#' standard deviation is profiled at its maximum likelihood estimate.
#'
#' \code{"joint_nb2"} uses L-BFGS-B to estimate kernel parameters, coefficients,
#' and \eqn{\rho = \log(1 + \alpha)} jointly, where \eqn{\alpha = 1/\theta}.
#' At \code{rho = 0}, the likelihood is exactly Poisson. Each candidate must
#' also agree with an independently optimized conditional coefficient fit
#' within \code{1e-5} negative log-likelihood units. Start records retain
#' failures, projected gradients, conditional likelihood gaps, and warnings.
#' \code{dispersion_est} adds \code{alpha} and \code{exact_poisson}; a boundary
#' estimate has \code{theta = Inf}. This is still an NB2 model, so the fitted
#' model's likelihood parameter count includes dispersion at the boundary.
#' Its conditional coefficient summary uses the Poisson GLM summary at that
#' boundary; inspect \code{dispersion_est} for the NB2 boundary metadata.
#' The \code{profile_fit} parameter vector uses log kernel parameters followed
#' by \code{rho}; \code{joint_par} also includes regression coefficients.
#' Its marginal covariance and Hessian are missing when joint uncertainty is
#' unavailable. Interior joint covariance returned by \code{vcov()} retains
#' log(theta) coordinates, as for \code{"profile_glm"}.
#'
#' Numeric kernel predictors may enter native R product, crossing and nesting
#' terms, including factor interactions and products of optimized predictors.
#' Native terms and fitted contrast coding define the coefficient layout;
#' all kernel, coefficient and cross derivatives enter joint uncertainty.
#' GAMs retain their smooth and penalty targets through the native-refit backend;
#' they are not treated as ordinary unpenalized GLMs here. Kernel-only scalar
#' arithmetic inside \code{I()}, and the base functions \code{exp}, \code{log},
#' \code{sqrt}, \code{sin}, \code{cos}, and \code{tan}, retain their native
#' coefficient basis and analytic derivatives. Their inputs must be optimized
#' kernel values or literal constants. Nonfinite domains are errors.
#' Fixed inputs within a dynamic expression, learned dynamic bases,
#' nondifferentiable or arbitrary functions, kernel-dependent offsets,
#' hard-radius metrics, random effects,
#' and zero inflation are unsupported. Put offsets in the formula with
#' \code{offset()}, rather than in the separate model argument.
#'
#' The retained model-frame rows define the fitting sample. Their row names
#' must match unique prepared point IDs; preserve these through subsetting.
#' Kernel means and standard deviations use only these rows. No observations
#' are dropped during optimization. Binomial responses must be binary or
#' \code{cbind(successes, failures)} with positive trial totals. Positive integer
#' frequency weights are supported for count and binomial models; Gaussian
#' weights must all equal one. Proportion responses and precision weights are
#' unsupported.
#'
#' These backends default to serial full multi-start fits; use \code{n_cores = 1},
#' \code{start_strategy = "single"}, and \code{profile_control}. It does not
#' accept \code{refit_fn} or \code{join_by}. Every winning start must pass both
#' optimizer convergence and projected-gradient checks.
#' The optional \code{profile_control$worker_pool} distributes independent
#' \code{joint_nb2} starts across a caller-owned PSOCK pool. It does not reduce
#' the number of starts or change uncertainty criteria. For repeated fitting,
#' create a pool once, pass \code{profile_control = list(worker_pool = pool)}
#' to each call, and close the pool afterward. Use \code{on.exit()} inside
#' your calling function to ensure cleanup if a fit fails.
#'
#' Joint observed information includes kernel, coefficient, and dispersion
#' cross derivatives. Default scale limits are normal approximations on the
#' log scale. Wald uncertainty is withheld at bounds, with nonpositive
#' information or scaled-information condition number at least 1e8, log-kernel
#' standard errors exceeding one, or competing starts more than 0.5 log units
#' apart within a 95 percent likelihood threshold. These are diagnostic
#' screening rules, not a proof of identifiability. A separate near-Poisson
#' diagnostic flags \code{max(mu/theta) < 0.01}. A finite upper theta bound
#' does not establish a finite dispersion maximum.
#' Joint Wald results are also withheld when a Poisson fit at the selected
#' scales is within the same likelihood threshold of the NB2 fit. This
#' conservative dispersion diagnostic does not constitute a formal model test;
#' reoptimizing Poisson scales could reveal additional compatibility.
#' An exact Poisson boundary fit never receives joint Wald uncertainty.
#' Explicit \code{uncertainty = "conditional"} predictions at that boundary
#' condition on both the fitted kernel parameters and the Poisson dispersion
#' boundary. They do not quantify uncertainty about the boundary or scales.
#' The conditional \code{opt_mod} retains standard GLM methods, but
#' \code{update(opt_mod)} is disabled because its original call references
#' pre-optimization kernel values. Refit through \code{multiScale_optim()}.
#'
#' \code{summary(x, profile = TRUE)} reoptimizes other kernel parameters and
#' NB2 dispersion while profiling each sigma. It reports the connected
#' likelihood region around the selected solution; limits that reach a search
#' bound are labeled \code{"profile-bound-limited"}, not interpreted as fully
#' resolved confidence limits. \code{profile_sigma()} instead holds other
#' kernel parameters fixed, while reestimating coefficients and NB2 dispersion.
#' Neither operation guarantees discovery of unvisited likelihood modes.
#'
#' \strong{Optimization approach}
#'
#' The default refit optimizer uses the L-BFGS-B algorithm (bounded quasi-Newton) to
#' maximize the log-likelihood over sigma (and shape for \code{"expow"}) while
#' holding all regression coefficients at their fitted values for each candidate
#' scale. Standard errors are Hessian-based approximations from the outer
#' optimization. Nonfinite, singular, or non-positive-definite active Hessians
#' produce missing SEs and a retained diagnostic; the raw Hessian is preserved.
#' Summary methods report profile-likelihood confidence intervals
#' for sigma when \code{\link{profile_sigma}} has been run on the object;
#' otherwise they fall back to Hessian-based intervals.
#' Candidate scales that make a covariate undefined at a modeled observation
#' are rejected, so the optimizer compares likelihoods for the same sample.
#' The final \code{sample_size} diagnostic reports prepared and fitted row
#' counts; for repeated-site \code{unmarked} data it also reports the number
#' of distinct prepared spatial sites. For hard-radius landscape and
#' unweighted surface metrics, Hessian uncertainty is not meaningful; use
#' \code{summary(x, profile = TRUE)} for profile-likelihood intervals.
#'
#' \strong{Binned acceleration}
#'
#' When \code{kernel_inputs} was created by \code{\link{kernel_prep}} with
#' \code{bin = TRUE} (the default), kernel-weighted covariates are evaluated
#' from precomputed distance-binned summaries rather than by iterating every
#' buffer cell on each optimizer evaluation. This makes the per-evaluation cost
#' independent of \code{max_D} and typically speeds up optimization by an order
#' of magnitude or more for large buffers, with a negligible binning
#' approximation. Inputs created with \code{store_cell_data = FALSE} ("lean"
#' mode) carry only the binned summaries; these still optimize, profile, and
#' refit normally. Landscape-metric covariates always use the exact per-cell
#' path and therefore require cell-level data.
#'
#' \strong{Parallel optimization}
#'
#' To ensure that fitted model function calls are properly serialized for
#' parallel workers, fit models using fully namespace-qualified functions. For
#' example, fit a negative binomial model as
#' \code{fitted_mod <- MASS::glm.nb(y ~ x, data = df)} rather than loading the
#' namespace implicitly.
#'
#' \strong{Joining unmarked multi-year data}
#'
#' When \code{unmarked} model rows repeat spatial sites, provide a
#' \code{join_by} data frame with one row per prepared site. For example, if
#' the spatial sites have IDs \code{c("A", "B")} and model rows have site IDs
#' \code{c("A", "B", "A", "B")} for two years, use
#' \code{join_by = data.frame(site = c("A", "B"))} and include the four-row
#' \code{site} column in the model's \code{siteCovs}. Multiple key columns
#' are allowed. Prepared keys must be unique and every model row must match.
#' The model's response and site-covariate row order is preserved. Without
#' \code{join_by}, the number and order of model rows must equal the prepared
#' spatial sites. When both objects retain the same row IDs, their order is
#' checked. Otherwise, the order must be supplied correctly by the caller.
#' Optimized covariates must be static \code{siteCovs};
#' \code{obsCovs} and \code{yearlySiteCovs} cannot be optimized through this
#' site-level workflow.
#'
#' \strong{Custom refit functions}
#'
#' By default, \code{multiScale_optim()} refits the model using
#' \code{stats::update()} (standard models) or the \code{unmarked} equivalent.
#' When the default path is insufficient, supply \code{refit_fn}. The function
#' must accept arguments \code{model}, \code{data}, and \code{context}, and
#' return a fitted model whose log-likelihood can be extracted by
#' \code{stats::logLik()} or \code{insight::get_loglikelihood()}.
#'
#' A minimal custom refit function:
#' \preformatted{
#' refit_fn <- function(model, data, context) \{
#'   stats::update(model, data = data)
#' \}
#' }
#'
#' For models requiring explicit call reconstruction:
#' \preformatted{
#' refit_fn <- function(model, data, context) \{
#'   call <- model$call
#'   call$data <- quote(data)
#'   eval(call, envir = list(data = data), enclos = parent.frame())
#' \}
#' }
#'
#' With PSOCK clusters, \code{refit_fn} must serialize cleanly. Prefer
#' namespace-qualified calls (e.g., \code{stats::update()}) and avoid
#' closures over large local objects.
#'
#' The callback receives \code{context$return_model = FALSE} for scalar
#' objective evaluations and \code{TRUE} when a fitted model is requested,
#' including the final selected-model refit. For supported unmarked models,
#' you can omit intermediate Hessian calculations while retaining the original
#' model's uncertainty settings whenever a fitted model is returned:
#' \preformatted{
#' refit_unmarked <- function(model, data, context) \{
#'   if (identical(context$return_model, FALSE))
#'     unmarked::update(model, data = data, se = FALSE)
#'   else
#'     unmarked::update(model, data = data)
#' \}
#' }
#' Supply this function as \code{refit_fn} to opt in. This changes neither
#' coefficient starts nor the model's likelihood. It does not add uncertainty
#' for kernel scales or override an original \code{se = FALSE} setting.
#' Verify this contract for your model class before using a custom refitter.
#'
#' \strong{Wrapper model objects}
#'
#' Some packages return wrapper objects around another fitted model (e.g.,
#' \code{amt::fit_clogit()} wraps a \code{survival::clogit()} fit in its
#' \code{model} component). When possible, \code{multiScale_optim()} unwraps
#' these automatically. For \code{amt::fit_clogit()}, pass \code{model = TRUE}
#' when fitting so the nested \code{clogit} model retains its model frame.
#'
#' When `start_strategy = "screen"` and `par` is left as `NULL`,
#' `multiScale_optim()` first scouts the sigma space with one-dimensional
#' log-spaced scans, then optionally tests a few jittered starts using short,
#' Hessian-free optimizations. The marginal sigma scans are serial, but the
#' short screening attempts use `n_cores` when parallel optimization is
#' requested. These screening steps are used only to choose one starting vector
#' for the single full optimization. On Windows, those parallel screening
#' attempts use PSOCK workers, so the same PSOCK serialization guidance applies
#' as in the full optimization. For reproducible screened starts, call
#' `set.seed()` before `multiScale_optim()`.
#'
#' @seealso \code{\link[multiScaleR]{kernel_dist}}
#' @examples
#' \donttest{
#' set.seed(555)
#'
#' points <- vect(cbind(c(5,7,9,11,13),
#'                      c(13,11,9,7,5)))
#'
#' mat_list <- list(r1 = rast(matrix(rnorm(20^2),
#'                                   nrow = 20)),
#'                  r2 = rast(matrix(rnorm(20^2),
#'                                   nrow = 20)))
#' rast_stack <- rast(mat_list)
#' kernel_inputs <- kernel_prep(pts = points,
#'                              raster_stack = rast_stack,
#'                              max_D = 5,
#'                              kernel = 'gaussian',
#'                              sigma = NULL)
#' ## Example response data
#' y <- rnorm(5)
#'
#' ## Create data frame with raster variables
#' dat <- data.frame(y = y,
#'                   kernel_inputs$kernel_dat)
#' mod1 <- glm(y ~ r1 + r2,
#'             data = dat)
#'
#' ## NOTE: This code is only for demonstration
#' ## Optimization results will have no meaning
#' opt_mod <- multiScale_optim(fitted_mod = mod1,
#'                             kernel_inputs = kernel_inputs,
#'                             par = NULL,
#'                             n_cores = NULL)
#'
#' ## Using package data
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
#' ## Optimize scale
#' opt <- multiScale_optim(fitted_mod = mod,
#'                         kernel_inputs = kernel_inputs)
#'
#' ## Summary of fitted model
#' summary(opt)
#'
#' ## 'True' parameter values data were simulated from:
#' # hab scale = 75
#' # Intercept = 0.5,
#' # hab slope estimate = 0.75
#'
#' ## Plot and visualize kernel density
#' plot(opt)
#'
#'
#' ## Apply optimized kernel to the environmental raster
#' opt_hab <- kernel_scale.raster(hab, multiScaleR = opt)
#'
#' plot(c(hab, opt_hab))
#'
#' ## Project model; scale & center
#' opt_hab.s_c <- kernel_scale.raster(raster_stack = hab,
#'                                    multiScaleR = opt,
#'                                    scale_center = TRUE)
#'
#' mod_pred <- predict(opt_hab.s_c, opt$opt_mod, type = 'response')
#' plot(mod_pred)
#'
#' ## Custom refit hook for model classes that need explicit control.
#' ## This example still uses glm(), but the same pattern can be used for
#' ## classes whose default update path is not sufficient.
#' refit_glm <- function(model, data, context) {
#'   stats::update(model, data = data)
#' }
#'
#' opt_custom <- multiScale_optim(fitted_mod = mod,
#'                                kernel_inputs = kernel_inputs,
#'                                refit_fn = refit_glm)
#'
#' ## Opt-in NB2 fit: theta and kernel scales are estimated together.
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   ki_nb <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                        kernel = "gaussian", bin = TRUE)
#'   nb_data <- data.frame(ki_nb$kernel_dat)
#'   nb_data$y <- rnbinom(nrow(nb_data),
#'                       mu = exp(0.5 + 0.8 * nb_data$hab), size = 2)
#'   nb_initial <- glm(y ~ hab, data = nb_data,
#'                      family = MASS::negative.binomial(2))
#'   nb_fit <- multiScale_optim(nb_initial, ki_nb, backend = "profile_glm",
#'                              n_cores = 1,
#'                              profile_control = list(n_starts = 6))
#'   diagnostics(nb_fit)$profile_glm
#'   predict(nb_fit)  # Expected counts for the retained fitting rows.
#'   if (nb_fit$uncertainty$available) {
#'     predict(nb_fit, interval = "confidence")
#'     vcov(nb_fit)  # Includes kernel/coefficient/dispersion covariance.
#'   }
#' }
#'}
#' @rdname multiScale_optim
#' @export
#' @importFrom optimParallel optimParallel
#' @importFrom stats optim
#' @importFrom parallel clusterEvalQ makeCluster setDefaultCluster stopCluster makeForkCluster clusterExport parLapply
#' @importFrom crayon %+% green red bold blue
#' @importFrom pscl zeroinfl
multiScale_optim <- function(fitted_mod,
                             kernel_inputs,
                             join_by = NULL,
                             par = NULL,
                             start_strategy = c("single", "screen"),
                             screen_n_sigma = 5,
                             screen_n_jitter = 6,
                             screen_maxit = 8,
                             screen_jitter_sd = 0.5,
                             n_cores = NULL,
                             PSOCK = FALSE,
                             verbose = TRUE,
                             refit_fn = NULL,
                             backend = c("refit", "profile_glm", "joint_nb2", "joint_mixed", "conditional_logit"),
                             profile_control = list()){
  backend <- match.arg(backend)
  if(is.null(fitted_mod)){
    stop("`fitted_mod` must be a fitted model object.")
  }
  user_supplied_par <- !is.null(par)
  start_strategy <- match.arg(start_strategy)
  validate_scalar_logical(PSOCK, "PSOCK")
  validate_scalar_logical(verbose, "verbose")
  validate_scalar_numeric(screen_n_sigma,
                          "screen_n_sigma",
                          lower = 3,
                          integerish = TRUE)
  validate_scalar_numeric(screen_n_jitter,
                          "screen_n_jitter",
                          lower = 0,
                          integerish = TRUE)
  validate_scalar_numeric(screen_maxit,
                          "screen_maxit",
                          positive = TRUE,
                          integerish = TRUE)
  validate_scalar_numeric(screen_jitter_sd,
                          "screen_jitter_sd",
                          lower = 0)
  if(!is.null(refit_fn) && !is.function(refit_fn)){
    stop("`refit_fn` must be a function if provided.", call. = FALSE)
  }

  # Check fitted_mod class
  # if (!inherits(fitted_mod, c("glm", "lm", "gls", "unmarked"))) {
  #   stop("fitted_mod must be of class 'glm', 'lm', 'gls', or 'unmarked'.")
  # }

  # Check kernel_inputs structure
  if (!is.list(kernel_inputs) || !all(c("raw_cov", "min_D", "max_D", "unit_conv", "kernel", "d_list") %in% names(kernel_inputs))) {
    stop("kernel_inputs must be a list with required elements: 'raw_cov', 'min_D', 'max_D', 'unit_conv', 'kernel', and 'd_list'.")
  }

  # Check join_by if provided
  if (!is.null(join_by) && !is.data.frame(join_by)) {
    stop("join_by must be a data frame if provided.")
  }

  # Check par values if provided
  if (!is.null(par) && !is.numeric(par)) {
    stop("par must be numeric if provided.")
  }

  # Check n_cores
  if (!is.null(n_cores) && (!is.numeric(n_cores) || n_cores < 1)) {
    stop("n_cores must be a positive integer if provided.")
  }
  if(!is.null(n_cores) && (length(n_cores) != 1 || is.na(n_cores) || n_cores != as.integer(n_cores))){
    stop("n_cores must be a positive integer if provided.")
  }
  if (!is.null(kernel_inputs$raw_cov) || !is.null(kernel_inputs$d_list)) {
    if(!is.list(kernel_inputs$raw_cov) || !is.list(kernel_inputs$d_list) ||
       length(kernel_inputs$raw_cov) == 0 || length(kernel_inputs$d_list) == 0 ||
       length(kernel_inputs$raw_cov) != length(kernel_inputs$d_list)){
      stop("`kernel_inputs$raw_cov` and `kernel_inputs$d_list` must be non-empty lists of equal length.")
    }
  } else if (is.null(kernel_inputs$binned)) {
    stop("`kernel_inputs` must contain either cell-level data (`raw_cov`/`d_list`) or precomputed `binned` summaries.")
  }
  validate_scalar_numeric(kernel_inputs$min_D, "kernel_inputs$min_D", positive = TRUE)
  validate_scalar_numeric(kernel_inputs$max_D, "kernel_inputs$max_D", positive = TRUE)
  validate_scalar_numeric(kernel_inputs$unit_conv, "kernel_inputs$unit_conv", positive = TRUE)
  kernel_inputs$kernel <- match.arg(kernel_inputs$kernel, c("gaussian", "exp", "expow", "fixed"))
  if (backend == "conditional_logit") {
    return(.msr_fit_clogit_backend(fitted_mod, kernel_inputs, join_by, par,
      start_strategy, n_cores, refit_fn, profile_control, verbose, match.call()))
  }
  if (backend == "joint_mixed") {
    return(.msr_fit_mixed_backend(fitted_mod, kernel_inputs, join_by, par,
      start_strategy, n_cores, refit_fn, profile_control, verbose, match.call()))
  }
  if (backend %in% c("profile_glm", "joint_nb2")) {
    return(.msr_fit_profile_backend(
      fitted_mod, kernel_inputs, join_by, par, start_strategy, n_cores,
      refit_fn, profile_control, verbose, match.call(), backend))
  }
  if (length(profile_control))
    stop("`profile_control` requires backend = 'profile_glm' or 'joint_nb2'.", call. = FALSE)
  analysis_mod <- .analysis_model(fitted_mod)

  # Extract variables from fitted model
  if (any(class(fitted_mod) == 'gls')) {
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- .model_predictors(analysis_mod)
  } else if (any(grepl("^unmarked", class(fitted_mod)))) {
    mod_vars <- .unmarked_model_predictors(fitted_mod)
  } else {
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- .model_predictors(analysis_mod)
  }
  validate_reserved_variable_names(mod_vars, "model predictor")

  if (!is.null(kernel_inputs$scale_vars)) {
    available_sources <- if (!is.null(kernel_inputs$raw_cov)) {
      colnames(kernel_inputs$raw_cov[[1]])
    } else if (!is.null(kernel_inputs$binned)) {
      unname(kernel_inputs$binned$sources)
    } else {
      character(0)
    }
    missing_sources <- setdiff(kernel_inputs$scale_vars$source,
                               available_sources)
    if (length(missing_sources) > 0) {
      stop("The raster surfaces provided do not match the variables used in your fitted model. Ensure names of surfaces match model variable names.",
           call. = FALSE)
    }
  }

  optimized_covariates <- .msr_optimized_covariates(
    kernel_inputs$scale_vars,
    cov_df = kernel_inputs$raw_cov
  )

  # Ensure model variables are in kernel_inputs
  r_vars <- mod_vars[mod_vars %in% optimized_covariates]
  if (length(r_vars) == 0) {
    stop("The raster surfaces provided do not match the variables used in your fitted model. Ensure names of surfaces match model variable names.")
  }

  n_covs <- length(r_vars)

  # Validate par length
  if (!is.null(par)) {
    expected_length <- if (kernel_inputs$kernel == 'expow') n_covs * 2 else n_covs
    if (length(par) != expected_length) {
      stop("The length of par does not match the expected number of covariates.")
    }
  }

  # Determine if parallel optimization should be used
  use_parallel <- !is.null(n_cores) && n_cores > 1

  if(any(class(fitted_mod) == 'gls')){
    mod_class <- 'gls'
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- .model_predictors(analysis_mod)
    r_vars <- mod_vars[which(mod_vars %in% optimized_covariates)]
    n_covs <- length(r_vars)
  } else if(any(grepl("^unmarked", class(fitted_mod)))) {
    mod_class <- 'unmarked'
    mod_vars <- .unmarked_model_predictors(fitted_mod)
    r_vars <- mod_vars[which(mod_vars %in% optimized_covariates)]
    n_covs <- length(r_vars)
    fitType <- fitted_mod@fitType
  } else {
    mod_class <- 'other'
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- .model_predictors(analysis_mod)
    r_vars <- mod_vars[which(mod_vars %in% optimized_covariates)]
    n_covs <- length(r_vars)
  }

  cnt <- 0
  opt_results <- data.frame()
  class(opt_results) <- 'try-error'
  par_starts <- seq((kernel_inputs$min_D/kernel_inputs$unit_conv) * 5,
                    (kernel_inputs$max_D/kernel_inputs$unit_conv) * 0.8,
                    length = 5)

  lwr <- rep(kernel_inputs$min_D, n_covs) / kernel_inputs$unit_conv
  uppr <- rep(kernel_inputs$max_D, n_covs) / kernel_inputs$unit_conv

  ## Modify to only confirm used variables in formula are present in stack

  if(!(any(mod_vars %in% optimized_covariates))){
    stop("\nThe raster surfaces provided do not match the variables used in your fitted model!\n\nMake sure raster surfaces being used in fitted model are included,\nand make sure names of surfaces match the names of variables in the model.\n")
  }


  if(!is.null(par) && kernel_inputs$kernel != 'expow' && length(par) != n_covs){
    stop("\nYou have specified starting values for `par`, the scale of effect. The number of starting values provided does not match the number of covariates in the model. Please correct and run again. \n")
  }

  if(!is.null(par) && kernel_inputs$kernel == 'expow' && length(par) != n_covs*2){
    stop("\nYou have specified an Exponential Power kernel, which has starting values (`par`), for both the scale of effect and shape. The number of starting values should be 2x the number of covariates in the model. Please correct and run again. \n")
  }

  # Reject an over-budget request before main setup, screening, or workers.
  .msr_stop_unsafe_parallel_ram(
    kernel_inputs = kernel_inputs,
    fitted_mod = fitted_mod,
    join_by = join_by,
    refit_fn = refit_fn,
    n_cores = n_cores,
    PSOCK = (.Platform$OS.type != "unix" || isTRUE(PSOCK))
  )

  opt_context <- build_opt_context(fitted_mod = fitted_mod,
                                   cov_df = kernel_inputs$raw_cov,
                                   join_by = join_by,
                                   refit_fn = refit_fn,
                                   scale_vars = kernel_inputs$scale_vars,
                                   unit_conv = kernel_inputs$unit_conv,
                                   resolution = kernel_inputs$resolution,
                                   n_cols = kernel_inputs$n_cols,
                                   binned = kernel_inputs$binned,
                                   d_list = kernel_inputs$d_list)

  if(kernel_inputs$kernel == 'expow'){
    lwr <- c(lwr, rep(0.75, n_covs))
    uppr <- c(uppr, rep(50, n_covs))
  }
  # if(kernel_inputs$kernel == 'expow'){
  #   lwr <- sqrt(c(lwr, rep(0.75, n_covs)))
  #   uppr <- sqrt(c(uppr, rep(50, n_covs)))
  # }

  if(is.null(par) & kernel_inputs$kernel != 'expow'){
    # par <- runif(n_covs, lwr, uppr)
    par <- rep(par_starts[2], n_covs)
    # par <- exp(rep(par_starts[2], n_covs))
  }

  if(is.null(par) & kernel_inputs$kernel == 'expow'){
    # par <- runif(n_covs, lwr[1:n_covs], uppr[1:n_covs])
    # par <- runif(n_covs, lwr[(n_covs+1):(n_covs*2)], uppr[(n_covs+1):(n_covs*2)])
    par <- rep(par_starts[2], n_covs)
    par <- c(par, rep(2, n_covs))
    # par <- exp(rep(par_starts[2], n_covs))
    # par <- c(par, sqrt(rep(2, n_covs)))
  }

  if (!user_supplied_par && identical(start_strategy, "screen")) {
    par <- .msr_prescreen_start(
      par = par,
      n_covs = n_covs,
      lwr = lwr,
      uppr = uppr,
      fitted_mod = fitted_mod,
      kernel_inputs = kernel_inputs,
      join_by = join_by,
      opt_context = opt_context,
      screen_n_sigma = screen_n_sigma,
      screen_n_jitter = screen_n_jitter,
      screen_maxit = screen_maxit,
      screen_jitter_sd = screen_jitter_sd,
      n_cores = n_cores,
      PSOCK = PSOCK,
      verbose = verbose
    )
  }

  opt_results <- data.frame()
  class(opt_results) <- 'try-error'

  if (use_parallel) {
    ## Initiate parallel cluster
    if (.Platform$OS.type == "unix" & isFALSE(PSOCK)) {
      cl <- makeForkCluster(n_cores)
    } else {
      cl <- makeCluster(n_cores)
    }
    setDefaultCluster(cl = cl)
    cluster_prep(fitted_mod, cl)

    ## Run parallel optimization
    opt_results <- try(
      optimParallel(
        par = par,
        fn = kernel_scale_fn,
        hessian = TRUE,
        lower = lwr,
        upper = uppr,
        fitted_mod = fitted_mod,
        d_list = kernel_inputs$d_list,
        cov_df = kernel_inputs$raw_cov,
        kernel = kernel_inputs$kernel,
        join_by = join_by,
        opt_context = opt_context,
        control = list(maxit = 1000),
        parallel = list(forward = FALSE, loginfo = TRUE)
      ),
      silent = TRUE
    )

    ## Stop parallel cluster
    setDefaultCluster(cl = NULL)
    stopCluster(cl)

    ## If parallel optimization failed
    if (inherits(opt_results, "try-error")) {
      err_msg <- attr(opt_results, "condition")
      if (verbose) {
        cat("\n\n*** Parallel optimization failed ***\n")
        if (!is.null(err_msg)) {
          cat("Error message:\n", conditionMessage(err_msg), "\n", sep = "")
        } else {
          cat("Unknown error occurred during parallel optimization.\n")
        }
        cat("Parallel optimization aborted.\n")
        cat("Try running with `n_cores = 1` to use standard optimization.\n\n")
      }
      if (!is.null(err_msg)) {
        stop("Parallel optimization failed: ", conditionMessage(err_msg), call. = FALSE)
      }
      stop("Parallel optimization failed due to an unknown error.", call. = FALSE)
    }

  } else {
    ## Standard (non-parallel) optimization
    opt_results <- try(
      optim(
        par = par,
        fn = kernel_scale_fn,
        hessian = TRUE,
        lower = lwr,
        upper = uppr,
        method = "L-BFGS-B",
        fitted_mod = fitted_mod,
        d_list = kernel_inputs$d_list,
        cov_df = kernel_inputs$raw_cov,
        control = list(maxit = 1000),
        kernel = kernel_inputs$kernel,
        join_by = join_by,
        opt_context = opt_context
      ),
      silent = TRUE
    )

    ## Handle standard optimization errors
    if (inherits(opt_results, "try-error")) {
      err_msg <- attr(opt_results, "condition")
      if (verbose) {
        cat("\n\n*** Standard optimization failed ***\n")
        if (!is.null(err_msg)) {
          cat("Error message:\n", conditionMessage(err_msg), "\n", sep = "")
        } else {
          cat("Unknown error occurred during standard optimization.\n")
        }
      }
      if (!is.null(err_msg)) {
        stop("Standard optimization failed: ", conditionMessage(err_msg), call. = FALSE)
      }
      stop("Standard optimization failed due to an unknown error.", call. = FALSE)
    }
  }


  if(inherits(opt_results, "try-error")){


    warning("\nOptimization has failed!\nConsider providing more informative starting values for optimization.")
    return(opt_results)
  } else {

    if (isTRUE(verbose)) cat('\n\nOptimization complete\n')

    opt_results$par_unscale <- c(opt_results$par[1:n_covs] * kernel_inputs$unit_conv,
                                 opt_results$par[(n_covs + 1):(n_covs * 2)])
    opt_results$hessian_unscale <- opt_results$hessian #* kernel_inputs$unit_conv
    # opt_results$par_unscale <- c(log(opt_results$par[1:n_covs] * kernel_inputs$unit_conv),
    #                              opt_results$par[(n_covs + 1):(n_covs * 2)]^2)
    # opt_results$hessian_unscale <- opt_results$hessian #* kernel_inputs$unit_conv

    hard_radius <- .msr_hard_radius_covariates(kernel_inputs$scale_vars)
    n_se <- n_covs * if (kernel_inputs$kernel == "expow") 2L else 1L
    active <- rep(!r_vars %in% hard_radius, length.out = n_se)
    uncertainty <- .msr_refit_hessian_se(opt_results$hessian_unscale, active,
      c(rep(kernel_inputs$unit_conv, n_covs),
        if (kernel_inputs$kernel == "expow") rep(1, n_covs)))
    res <- uncertainty$SE

    scale_est <- data.frame(Mean = opt_results$par_unscale[1:n_covs],
                            SE = res[1:n_covs])
    rownames(scale_est) <- r_vars
    hard_radius <- .msr_hard_radius_covariates(kernel_inputs$scale_vars)
    scale_est[intersect(rownames(scale_est), hard_radius), "SE"] <- NA_real_

    if(kernel_inputs$kernel == 'expow') {
      shape_est <- data.frame(Mean = opt_results$par_unscale[(n_covs + 1):(n_covs * 2)],
                              SE = res[(n_covs + 1):(n_covs * 2)])
      rownames(shape_est) <- r_vars
      inactive_shape <- intersect(r_vars, hard_radius)
      if (length(inactive_shape)) {
        shape_est[inactive_shape, c("Mean", "SE")] <- NA_real_
        attr(shape_est, "not_applicable") <- inactive_shape
      }

    } else {
      shape_est <- NULL
    }

    # browser()

    final_mod <- kernel_scale_fn(par = c(opt_results$par),
                                 kernel = kernel_inputs$kernel,
                                 d_list = kernel_inputs$d_list,
                                 cov_df = kernel_inputs$raw_cov,
                                 fitted_mod = fitted_mod,
                                 join_by = join_by,
                                 mod_return = TRUE,
                                 opt_context = opt_context)

    scl_params <- .msr_merge_scl_params(primary = final_mod$scl_params,
                                        fallback = kernel_inputs$scl_params)

    out <- list(scale_est = scale_est,
                shape_est = shape_est,
                optim_results = opt_results,
                opt_mod = final_mod$mod,
                fitted_mod_original = fitted_mod,
                min_D = kernel_inputs$min_D,
                max_D = kernel_inputs$max_D,
                kernel_inputs = kernel_inputs[
                  setdiff(names(kernel_inputs), c("min_D", "max_D"))
                ],
                scl_params = scl_params,
                join_by = join_by,
                opt_context = opt_context,
                profile_scale_est = NULL,
                diagnostics = .empty_diagnostics(),
                next_run = NULL,
                warn_message = 0,
                call = match.call())

    out$diagnostics$refit_hessian <- uncertainty$diagnostic
    class(out) <- 'multiScaleR'

    # Keep the existing warning threshold independent of distance reporting.
    scale_D <- kernel_dist(out, dimension = "1d")
    max_dist_diag <- .max_distance_diagnostic(scale_D = scale_D,
                                              max_D = kernel_inputs$max_D)
    out$diagnostics$max_distance <- max_dist_diag

    fitted_n <- .msr_model_nobs(out$opt_mod)
    prepared_sites <- nrow(kernel_inputs$kernel_dat)
    # A repeated-site unmarked frame has more model rows than spatial sites.
    # Its data slots also retain rows dropped from the likelihood.
    prepared_n <- if (inherits(fitted_mod, "unmarkedFit")) {
      fitted_n <- fitted_n - length(out$opt_mod@sitesRemoved)
      nrow(fitted_mod@data@siteCovs)
    } else {
      prepared_sites
    }
    sample_size <- list(
      code = "sample_size",
      triggered = fitted_n < prepared_n,
      fitted_n = fitted_n,
      prepared_n = prepared_n
    )
    if (inherits(fitted_mod, "unmarkedFit")) {
      sample_size$prepared_sites <- prepared_sites
    }
    out$diagnostics$sample_size <- sample_size
    if (isTRUE(out$diagnostics$sample_size$triggered)) {
      warning(sprintf(
        "The fitted model uses %d of %d prepared observations. Check the initial model's missing-data handling before comparing fits.",
        fitted_n, prepared_n
      ), call. = FALSE)
    }

    if (isTRUE(max_dist_diag$triggered)) {
      out$warn_message <- c(out$warn_message, 1)
      if (isTRUE(verbose)) cat(red("\n WARNING!!!\n",
              "The estimated scale of effect approaches or exceeds the current search extent.\n",
              "Consider increasing " %+% blue$bold("max_D") %+% " in `kernel_prep` to >="  %+% green$bold(max_dist_diag$suggested_max_D) %+% " to ensure accurate estimation of scale.\n\n"))
    }

    smooth_scale <- scale_est[!rownames(scale_est) %in% hard_radius, , drop = FALSE]
    if(nrow(smooth_scale) && any((smooth_scale[,1] / smooth_scale[,2]) < 2,
                                 na.rm = TRUE)){
      out$diagnostics$sigma_precision <- .precision_diagnostic(smooth_scale,
                                                                "sigma_precision")
      out$warn_message <- c(out$warn_message, 2)
      if (isTRUE(verbose)) cat(red("\n WARNING!!!\n",
              "The standard error of one or more `sigma` estimates is >= 50% of the estimated mean value.\n",
              "Carefully assess whether or not this variable is meaningful in your analysis and interpret with caution.\n\n"))
    }

    if(!is.null(shape_est) && any((shape_est[,1] / shape_est[,2]) < 2, na.rm = T)){
      out$diagnostics$shape_precision <- .precision_diagnostic(shape_est, "shape_precision")
      out$warn_message <- c(out$warn_message, 3)
      if (isTRUE(verbose)) cat(red("\n WARNING!!!\n",
              "The standard error of one or more `shape` estimates is >= 50% of the estimated mean value.\n",
              "Carefully assess if the Exponential Power kernel is appropriate, whether or not this variable is meaningful in your analysis, and interpret with caution.\n\n"))
    }

    out$next_run <- .next_run_recommendation(
      scale_est = scale_est,
      shape_est = shape_est,
      max_D = kernel_inputs$max_D,
      diagnostics = out$diagnostics
    )

    return(out)
  }
}


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
#'   via \code{stats::update()} or a custom \code{refit_fn}.
#' @param kernel_inputs A list of class \code{"multiScaleR_data"} created by
#'   \code{\link{kernel_prep}}. Must contain elements \code{raw_cov},
#'   \code{d_list}, \code{min_D}, \code{max_D}, \code{unit_conv}, and
#'   \code{kernel}.
#' @param join_by Default: \code{NULL}. A data frame used to join site-level
#'   spatial covariates to repeated observations for \code{unmarked} models
#'   where sites are surveyed across multiple years. The column name in
#'   \code{join_by} must match a column in the data used to fit the
#'   \code{unmarked} model. See Details.
#' @param par Optional numeric vector of starting values for the optimizer.
#'   Values must be divided by \code{max_D} to match the internal scaled
#'   parameter space. Length must equal the number of optimized covariates
#'   (or twice that for \code{kernel = "expow"}, where shape parameters follow
#'   sigma parameters). Default: \code{NULL}; starting values are chosen
#'   automatically.
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
#'   warns if the requested workers exceed the conservative RAM budget for the
#'   current \code{kernel_inputs} payload.
#' @param PSOCK Logical. On Windows, a PSOCK cluster is always used. On Unix,
#'   a FORK cluster is used by default (faster). Set \code{TRUE} to force a
#'   PSOCK cluster on Unix. Default: \code{FALSE}.
#' @param verbose Logical. Print optimization status and warnings to the
#'   console. Default: \code{TRUE}.
#' @param refit_fn Optional function for refitting \code{fitted_mod} during
#'   optimization. When \code{NULL} (default), \code{stats::update()} (or the
#'   unmarked equivalent) is used. Provide this only when the default refit
#'   path is insufficient. See Details.
#'
#' @return A list of class \code{"multiScaleR"} containing:
#' \describe{
#'   \item{\code{scale_est}}{Data frame with one row per optimized covariate and
#'     two columns: \code{Mean} (optimized sigma on the original projection
#'     scale) and \code{SE} (Hessian-based standard error). Row names are
#'     covariate names.}
#'   \item{\code{shape_est}}{Data frame of the same structure as
#'     \code{scale_est} for the shape parameter, or \code{NULL} when
#'     \code{kernel != "expow"}.}
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
#' }
#' @details
#' \strong{Optimization approach}
#'
#' The optimizer uses the L-BFGS-B algorithm (bounded quasi-Newton) to
#' maximize the log-likelihood over sigma (and shape for \code{"expow"}) while
#' holding all regression coefficients at their fitted values for each candidate
#' scale. Standard errors are Hessian-based approximations from the outer
#' optimization. Summary methods report profile-likelihood confidence intervals
#' for sigma when \code{\link{profile_sigma}} has been run on the object;
#' otherwise they fall back to Hessian-based intervals.
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
#' When using \code{unmarked} models where sites are surveyed across multiple
#' years but the spatial covariates are constant across years, provide a
#' \code{join_by} data frame to match each site's kernel-weighted covariate
#' value to its observations. The column name in \code{join_by} must match a
#' column in the data used to fit the \code{unmarked} model.
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
                             refit_fn = NULL){
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
  analysis_mod <- .analysis_model(fitted_mod)

  # Extract variables from fitted model
  if (any(class(fitted_mod) == 'gls')) {
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- .model_predictors(analysis_mod)
  } else if (any(grepl("^unmarked", class(fitted_mod)))) {
    mod_vars <- all.vars(formula(fitted_mod@formula))
  } else {
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- .model_predictors(analysis_mod)
  }

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
    mod_vars <- all.vars(formula(fitted_mod@formula))
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

  opt_context <- build_opt_context(fitted_mod = fitted_mod,
                                   cov_df = kernel_inputs$raw_cov,
                                   join_by = join_by,
                                   refit_fn = refit_fn,
                                   scale_vars = kernel_inputs$scale_vars,
                                   unit_conv = kernel_inputs$unit_conv,
                                   resolution = kernel_inputs$resolution,
                                   n_cols = kernel_inputs$n_cols,
                                   binned = kernel_inputs$binned)

  .msr_warn_unsafe_parallel_ram(
    kernel_inputs = kernel_inputs,
    fitted_mod = fitted_mod,
    join_by = join_by,
    refit_fn = refit_fn,
    n_cores = n_cores,
    PSOCK = (.Platform$OS.type != "unix" || isTRUE(PSOCK))
  )

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

    cat('\n\nOptimization complete\n')

    opt_results$par_unscale <- c(opt_results$par[1:n_covs] * kernel_inputs$unit_conv,
                                 opt_results$par[(n_covs + 1):(n_covs * 2)])
    opt_results$hessian_unscale <- opt_results$hessian #* kernel_inputs$unit_conv
    # opt_results$par_unscale <- c(log(opt_results$par[1:n_covs] * kernel_inputs$unit_conv),
    #                              opt_results$par[(n_covs + 1):(n_covs * 2)]^2)
    # opt_results$hessian_unscale <- opt_results$hessian #* kernel_inputs$unit_conv

    i_hess <- try(solve(opt_results$hessian_unscale), silent = TRUE)

    if(!inherits(i_hess, "try-error")){
      res <- c(sqrt(diag(i_hess)[1:n_covs]) * kernel_inputs$unit_conv,
               sqrt(diag(i_hess)[(n_covs + 1):(n_covs * 2)]))

    } else {
      n_se <- if (kernel_inputs$kernel == "expow") n_covs * 2 else n_covs
      res <- rep(Inf, n_se)

    }

    scale_est <- data.frame(Mean = opt_results$par_unscale[1:n_covs],
                            SE = res[1:n_covs])
    rownames(scale_est) <- r_vars

    if(kernel_inputs$kernel == 'expow') {
      shape_est <- data.frame(Mean = opt_results$par_unscale[(n_covs + 1):(n_covs * 2)],
                              SE = res[(n_covs + 1):(n_covs * 2)])
      rownames(shape_est) <- r_vars

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

    class(out) <- 'multiScaleR'

    scale_D <- kernel_dist(out)
    max_dist_diag <- .max_distance_diagnostic(scale_D = scale_D,
                                              max_D = kernel_inputs$max_D)
    out$diagnostics$max_distance <- max_dist_diag

    if (isTRUE(max_dist_diag$triggered)) {
      out$warn_message <- c(out$warn_message, 1)
      cat(red("\n WARNING!!!\n",
              "The estimated scale of effect approaches or exceeds the current search extent.\n",
              "Consider increasing " %+% blue$bold("max_D") %+% " in `kernel_prep` to >="  %+% green$bold(max_dist_diag$suggested_max_D) %+% " to ensure accurate estimation of scale.\n\n"))
    }

    if(any((scale_est[,1] / scale_est[,2]) < 2, na.rm = T)){
      out$diagnostics$sigma_precision <- .precision_diagnostic(scale_est, "sigma_precision")
      out$warn_message <- c(out$warn_message, 2)
      cat(red("\n WARNING!!!\n",
              "The standard error of one or more `sigma` estimates is >= 50% of the estimated mean value.\n",
              "Carefully assess whether or not this variable is meaningful in your analysis and interpret with caution.\n\n"))
    }

    if(!is.null(shape_est) && any((shape_est[,1] / shape_est[,2]) < 2, na.rm = T)){
      out$diagnostics$shape_precision <- .precision_diagnostic(shape_est, "shape_precision")
      out$warn_message <- c(out$warn_message, 3)
      cat(red("\n WARNING!!!\n",
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


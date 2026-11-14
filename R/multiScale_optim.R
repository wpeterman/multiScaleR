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


# Flag a refit search that never left its starting values.
#
# L-BFGS-B stops after one evaluation when its finite-difference gradient is
# exactly zero. That happens when the objective is flat around the start:
# every optimized scale is a hard radius (the likelihood changes only when a
# whole cell enters or leaves the buffer), or every refit at the start failed
# and returned the penalty value (for example, glmmTMB reports an NA
# log-likelihood when its Hessian is not positive definite). A search can
# also try a few points, find no improvement, and return its start, as with
# an irregular refit objective. Either way the fit reports its starting values
# as estimates, so the user must be told.
# A start on a bound is excluded, because a zero projected gradient there is
# a legitimate stop. A user-supplied or screened start that survives a line
# search may already be an optimum, so only a zero-step stop flags it. The
# fixed default start is never a plausible optimum. Returns NULL when the
# search moved.
.msr_start_stall_diagnostic <- function(opt_results, start, lower, upper,
                                        unit_conv, covariates, hard_only,
                                        screened, default_start = FALSE,
                                        penalty = 1e6^10) {
  evals <- opt_results$counts[[1L]]
  stalled <- identical(as.numeric(opt_results$par), as.numeric(start)) &&
    all(start > lower & start < upper) &&
    (isTRUE(evals <= 1L) || isTRUE(default_start))
  if (!stalled) return(NULL)

  n_covs <- length(covariates)
  start_sigma <- stats::setNames(start[seq_len(n_covs)] * unit_conv, covariates)
  reason <- if (isTRUE(opt_results$value >= penalty)) {
    "refit_failed"
  } else if (isTRUE(hard_only)) {
    "hard_radius"
  } else {
    "no_progress"
  }
  message <- if (reason == "refit_failed") {
    paste0(
      "The optimizer returned its starting values because no refitted model ",
      "gave a finite log-likelihood at the starting scales. Check that the ",
      "starting model converged. A glmmTMB model with a non-positive-definite ",
      "Hessian has no log-likelihood, so refit it with explicit `start` values ",
      "and run multiScale_optim() again.")
  } else if (screened) {
    paste0(
      "The final search did not move from the screened start, so the ",
      "estimates are the best screened candidate rather than a local optimum. ",
      if (reason == "hard_radius") paste0(
        "Every optimized scale is a hard radius, whose likelihood changes in ",
        "steps. "),
      "Increase `screen_n_sigma` or compare fits from several `par` values.")
  } else if (reason == "hard_radius") {
    paste0(
      "The optimizer returned its starting values. Every optimized scale is a ",
      "hard radius, whose likelihood changes in steps, so a gradient search ",
      "cannot move from its start. Use start_strategy = \"screen\" or supply `par`.")
  } else {
    paste0(
      "The optimizer returned its starting values because it found no ",
      "improvement near them. The likelihood there may be flat or too ",
      "irregular for a gradient search. Supply `par` or use ",
      "start_strategy = \"screen\".")
  }
  list(code = "start_stall",
       triggered = TRUE,
       reason = reason,
       screened = isTRUE(screened),
       start = start_sigma,
       function_evaluations = evals,
       message = message)
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

#' @title Estimate scales of effect by maximum likelihood
#'
#' @description Estimates the scale of effect (sigma) of each kernel-weighted
#' covariate in a fitted model by maximizing the model's log-likelihood. The
#' *scale of effect* is the distance over which a landscape variable influences
#' the response. By default, `multiScale_optim()` recomputes the kernel-weighted
#' covariates at each candidate scale, refits your model, and searches with the
#' L-BFGS-B algorithm through [stats::optim()] (or
#' [optimParallel::optimParallel()] for parallel searches). Other backends
#' estimate the scales and the model parameters in one likelihood; see
#' [multiScaleR-fitting-methods].
#'
#' @param fitted_mod A fitted model object whose covariates include the
#'   kernel-weighted variables defined in `kernel_inputs`. Supported classes
#'   include `lm`, `glm`, `gls` (nlme), and `unmarkedFit` (unmarked). Many
#'   other model classes also work through `stats::update()` or a custom
#'   `refit_fn`. A predictor cannot be named `"x"`, which is reserved for
#'   internal coordinate data.
#' @param kernel_inputs A list of class `"multiScaleR_data"` created by
#'   [kernel_prep()]. It must contain the elements `raw_cov`, `d_list`, `min_D`,
#'   `max_D`, `unit_conv`, and `kernel`. See [multiScaleR_data-class].
#' @param join_by Default: `NULL`. For `unmarked` models whose rows repeat
#'   sites (for example one row per site and year), a data frame with one
#'   unique key row per spatial site in `kernel_inputs`, in the same order as
#'   those sites. Its named key columns must also occur in the fitted model's
#'   `siteCovs`. See "Fit unmarked models with repeated sites" below.
#' @param par Optional numeric vector of starting values. Default: `NULL`,
#'   which chooses starts automatically (see `start_strategy`). For the refit
#'   backend, divide each starting sigma by `max_D` to match the internal scaled
#'   parameter space. For the analytic backends, give scales in internal
#'   distance units (map distance divided by `unit_conv`), as for the rows of
#'   `profile_control$starts`. The length must equal the number of optimized
#'   covariates, or twice that for `kernel = "expow"`, where the shape values
#'   follow the sigma values. The `"expow"` layout keeps one shape position per
#'   optimized covariate, including hard radii, which ignore it; give a finite
#'   positive placeholder there.
#' @param start_strategy Character. How to choose starting values when
#'   `par = NULL`. `"single"` uses one shared default start. `"screen"` runs a
#'   cheap prescreen of sigma values and a few short searches before the one
#'   full search (see "Choose starting values" below). The default, `NULL`,
#'   uses `"screen"` when `kernel_inputs` was prepared with `kernel = "fixed"`
#'   and the backend is `"refit"` without a `radius_grid`, and `"single"`
#'   otherwise. Under a fixed kernel every optimized scale is a hard radius,
#'   and a gradient search cannot leave a single start on that stepped
#'   likelihood. A refit search that returns its starting values gives a
#'   warning and records `start_stall` in [diagnostics()].
#' @param screen_n_sigma Positive integer. Number of log-spaced sigma values
#'   evaluated per covariate during the prescreen when
#'   `start_strategy = "screen"`. Default: `5`.
#' @param screen_n_jitter Non-negative integer. Number of jittered candidate
#'   starts tried with short searches after the prescreen. Default: `6`. Set
#'   to `0` to skip the short jittered runs.
#' @param screen_maxit Positive integer. Maximum iterations for each short
#'   screening search. Default: `8`.
#' @param screen_jitter_sd Non-negative numeric scalar. Standard deviation of
#'   the multiplicative sigma jitter, on the log scale, used during screening.
#'   Default: `0.5`.
#' @param n_cores Positive integer. Number of cores for a parallel search
#'   through `optimParallel`. Default: `NULL` (serial). With
#'   `start_strategy = "screen"`, the short screening searches also use these
#'   cores. Parallel searches help most with many covariates or slow model
#'   fits. When `n_cores > 1`, `multiScale_optim()` checks
#'   [estimate_multiscale_ram()] and stops before screening or creating
#'   workers if the requested workers exceed its conservative RAM-based limit.
#'   Use `n_cores = 1` for a serial search, or rebuild the inputs with
#'   `kernel_prep(..., bin = TRUE, store_cell_data = FALSE)` (lean mode) when
#'   every covariate is kernel-weighted. Landscape and surface metrics cannot
#'   use lean mode. An unavailable RAM estimate does not block the run, a
#'   passing estimate does not guarantee enough memory, and the worker count is
#'   never reduced automatically.
#' @param PSOCK Logical. On Windows, a PSOCK cluster is always used. On Unix,
#'   a FORK cluster is used by default because it is faster. Set `TRUE` to
#'   force a PSOCK cluster on Unix. Default: `FALSE`.
#' @param verbose Logical. Print the search status and diagnostic messages to
#'   the console. Default: `TRUE`.
#' @param refit_fn Optional function for refitting `fitted_mod` during the
#'   search. When `NULL` (default), `stats::update()` (or the `unmarked`
#'   equivalent) is used. Supply one only when the default refit is not enough;
#'   see "Supply a custom refit function" below.
#' @param covariate_transform Character. How optimized covariates enter the
#'   model. `"candidate"` (default) centers and scales the kernel or metric
#'   values at each candidate scale, so a coefficient describes a change of one
#'   standard deviation. `"raw"` keeps the values in their own units (for a
#'   binary layer, the kernel-weighted proportion) and records centers of zero
#'   and scales of one for prediction. Fit the starting model in those raw
#'   units, including any fixed spline knots. `"raw"` requires
#'   `backend = "refit"`. See [multiScaleR-fitting-methods].
#' @param backend Character. The fitting method. `"refit"` (default) refits
#'   `fitted_mod` at each candidate scale and works with any model class that
#'   can be refit. `"profile_glm"` (Gaussian, Poisson, binomial, and NB2 GLMs),
#'   `"joint_nb2"` (NB2 models, including the exact Poisson boundary),
#'   `"joint_mixed"` (Poisson and NB2 `glmmTMB` models with one random
#'   intercept), and `"conditional_logit"` (step-selection models) estimate
#'   the scales and the model parameters in one likelihood. Fit a Gaussian
#'   model with `glm(family = gaussian())` for the analytic backends, because
#'   `lm()` fits are rejected. [multiScaleR-fitting-methods] gives each
#'   backend's requirements and the uncertainty it reports.
#' @param profile_control Named list of settings for the analytic backends,
#'   such as the number of starts (`n_starts`, default 6), a matrix of explicit
#'   `starts`, iteration limits, NB2 theta bounds, and a PSOCK `worker_pool`
#'   for `"joint_nb2"`. Default: `list()`, which uses every default. The refit
#'   backend accepts no entries. The "Control settings" section of
#'   [multiScaleR-fitting-methods] lists every entry, the backends that accept
#'   it, and its default.
#' @param radius_grid Experimental, so its behavior may change in a later
#'   release without a deprecation period. With `backend = "refit"`, an
#'   optional named list of candidate radius vectors, one for each optimized
#'   hard-radius predictor. For a single hard radius, supply the vector
#'   directly or name it in a list. Each vector needs at least two distinct
#'   radii in map units, between `kernel_inputs$min_D` and
#'   `kernel_inputs$max_D`. The fit evaluates every combination and
#'   reoptimizes any smooth kernel bandwidths and shapes at each one. The
#'   default `NULL` uses the backend's usual search. See "Profile hard radii
#'   over a grid" in [multiScaleR-fitting-methods].
#'
#' @return An object of class `"multiScaleR"`, a list. The elements you will
#'   use most are:
#'   * `scale_est`: the estimated sigma of each optimized covariate (`Mean`)
#'     and its standard error (`SE`), in map units;
#'   * `shape_est`: the same for the shape parameter under `kernel = "expow"`,
#'     otherwise `NULL`;
#'   * `opt_mod`: your model refit at the estimated scales, to use for
#'     coefficients, prediction, and effect plots;
#'   * `diagnostics`: search and uncertainty diagnostics, read with
#'     [diagnostics()];
#'   * `next_run`: suggested `max_D` and starting values for a follow-up fit.
#'
#'   [multiScaleR-class] describes every element, including those added by the
#'   analytic backends, and lists the methods that read the object.
#'
#' @details
#' # How the search works
#'
#' A *kernel* is a distance-weighting function that converts a raster layer
#' into one covariate value per survey point. Its scale, sigma, sets how fast
#' the weights decline with distance (see [multiScaleR-kernels]). The default
#' refit backend proposes sigma values (and shape values for `"expow"`),
#' recomputes the kernel-weighted covariates, refits `fitted_mod`, and keeps
#' the values that maximize the log-likelihood. Standard errors come from the
#' Hessian of this outer search. Hard-radius covariates have no Hessian
#' standard error because their likelihood changes in steps. Profile
#' likelihood intervals are available with `summary(fit, profile = TRUE)` or
#' [profile_sigma()]. [multiScaleR-fitting-methods] gives the rules for each
#' backend.
#'
#' # Read the warnings
#'
#' `multiScale_optim()` prints a warning, and records a code in
#' `warn_message`, in three situations:
#' * **The scale approaches `max_D`.** The estimate is near or beyond the
#'   search extent, so the true scale may be larger. Rebuild the inputs with
#'   [kernel_prep()] using at least the suggested `max_D`, which is also stored
#'   in `fit$next_run$max_D`, and refit.
#' * **A sigma standard error is at least half its estimate.** The scale is
#'   poorly determined. Assess whether the covariate is meaningful in your
#'   analysis and interpret it with caution.
#' * **A shape standard error is at least half its estimate.** Assess whether
#'   the exponential power kernel suits your data.
#'
#' A search that never leaves its starting values also warns and records
#' `start_stall` in [diagnostics()]. Supply different starts with `par`, or use
#' `start_strategy = "screen"`.
#'
#' # Choose starting values
#'
#' With `start_strategy = "screen"` and `par = NULL`, `multiScale_optim()`
#' first scouts the sigma range with one-dimensional, log-spaced scans of each
#' covariate. It then tests a few jittered starts with short searches that skip
#' the Hessian. The marginal scans are serial, but the short searches use
#' `n_cores` when you request a parallel search. Screening only chooses one
#' starting vector for the single full search. On Windows those parallel
#' screening searches use PSOCK workers, so the serialization advice under
#' "Speed up the search" applies to them too. Call `set.seed()` before
#' `multiScale_optim()` for reproducible screened starts.
#'
#' To continue from an earlier fit, rebuild the inputs with `fit$next_run$max_D`
#' and pass `fit$next_run$start_par` as `par`. Those starts are scaled to that
#' `max_D`.
#'
#' # Speed up the search
#'
#' When [kernel_prep()] was run with `bin = TRUE` (the default),
#' kernel-weighted covariates are computed from precomputed distance-binned
#' summaries instead of every buffer cell at each evaluation. The cost of an
#' evaluation then no longer depends on `max_D`, which usually speeds the
#' search by an order of magnitude or more for large buffers, with a
#' negligible binning approximation. Inputs created with
#' `store_cell_data = FALSE` (lean mode) carry only the binned summaries; they
#' still optimize, profile, and refit normally. Landscape-metric covariates
#' always use the exact per-cell path and so need the cell-level data.
#'
#' For a parallel search, fit the starting model with fully
#' namespace-qualified functions so that its call serializes cleanly to the
#' workers. For example, fit a negative binomial model as
#' `fitted_mod <- MASS::glm.nb(y ~ x, data = df)` rather than relying on an
#' attached package.
#'
#' # Fit unmarked models with repeated sites
#'
#' When the rows of an `unmarked` model repeat spatial sites, for example a
#' stacked layout with one row per site and year, the preparation object
#' still has one row per unique site. Add a site ID to the model's `siteCovs`
#' and pass one matching ID per prepared site through `join_by`. For spatial
#' sites `c("A", "B")` and model rows `c("A", "B", "A", "B")` for two years,
#' use `join_by = data.frame(site = c("A", "B"))` and include the four-row
#' `site` column in `siteCovs`. Multiple key columns are allowed. Prepared keys
#' must be unique, every model row must match one, and the IDs in the model
#' may repeat in any order. The model's response and site-covariate row order
#' is preserved.
#'
#' The layout looks like this, where `counts_by_year` is a list of count
#' matrices, each with one row per site (in the order of
#' `kernel_inputs$kernel_dat`) and one column per visit:
#'
#' ```r
#' site_id <- seq_len(nrow(kernel_inputs$kernel_dat))
#' y_stacked <- do.call(rbind, counts_by_year)
#' site_row <- rep(site_id, times = length(counts_by_year))
#' site_covs <- kernel_inputs$kernel_dat[site_row, , drop = FALSE]
#' site_covs$site_id <- site_id[site_row]
#' site_covs$Year <- factor(rep(seq_along(counts_by_year),
#'                              each = length(site_id)))
#' umf <- unmarked::unmarkedFramePCount(y = y_stacked, siteCovs = site_covs)
#' mod <- unmarked::pcount(~Year ~ Year + hab, data = umf, K = 100)
#' fit <- multiScale_optim(mod, kernel_inputs,
#'                         join_by = data.frame(site_id = site_id))
#' ```
#'
#' Check that the rows of `y_stacked` and `site_covs` refer to the same
#' site-year records before fitting. Without `join_by`, model rows must match
#' the prepared sites one for one, in the same order. When both objects keep
#' the same row IDs their order is checked; otherwise you are responsible for
#' it. An open-population `unmarked::pcountOpen()` model instead stores one
#' row per site with the years across count columns, so it uses the same
#' static site covariates without `join_by`. Optimized covariates must be
#' static `siteCovs`; `obsCovs` and `yearlySiteCovs` cannot be optimized
#' through this site-level workflow. Sparse counts can leave scales weakly
#' identified, so examine the scale profile and [diagnostics()].
#'
#' # Supply a custom refit function
#'
#' By default, `multiScale_optim()` refits the model with `stats::update()`
#' (standard models) or the `unmarked` equivalent. When that is not enough,
#' supply `refit_fn`. The function must accept the arguments `model`, `data`,
#' and `context`, and return a fitted model whose log-likelihood
#' `stats::logLik()` or `insight::get_loglikelihood()` can extract. A minimal
#' refit function:
#'
#' ```r
#' refit_fn <- function(model, data, context) {
#'   stats::update(model, data = data)
#' }
#' ```
#'
#' For models that need their call rebuilt explicitly:
#'
#' ```r
#' refit_fn <- function(model, data, context) {
#'   call <- model$call
#'   call$data <- quote(data)
#'   eval(call, envir = list(data = data), enclos = parent.frame())
#' }
#' ```
#'
#' With PSOCK clusters, `refit_fn` must serialize cleanly. Prefer
#' namespace-qualified calls such as `stats::update()` and avoid closures over
#' large local objects.
#'
#' The function receives `context$return_model = FALSE` for scalar objective
#' evaluations and `TRUE` when a fitted model is wanted, including the final
#' refit at the selected scales. For supported `unmarked` models, you can skip
#' the intermediate Hessian calculations while keeping the original model's
#' uncertainty settings whenever a fitted model is returned:
#'
#' ```r
#' refit_unmarked <- function(model, data, context) {
#'   if (identical(context$return_model, FALSE))
#'     unmarked::update(model, data = data, se = FALSE)
#'   else
#'     unmarked::update(model, data = data)
#' }
#' ```
#'
#' This changes neither the coefficient starts nor the model's likelihood. It
#' does not add uncertainty for the kernel scales or override an original
#' `se = FALSE` setting. Verify this behavior for your model class before you
#' rely on a custom refit function.
#'
#' # Fit wrapper model objects
#'
#' Some packages return a wrapper around another fitted model. For example,
#' `amt::fit_clogit()` wraps a `survival::clogit()` fit in its `model`
#' component. `multiScale_optim()` unwraps these automatically when it can.
#' For `amt::fit_clogit()`, pass `model = TRUE` when fitting so the nested
#' `clogit` model keeps its model frame. For complete choice sets with learned
#' predictors such as `base::scale()`, use `na.action = stats::na.fail` to keep
#' the fitted transformation attributes, because some `survival` versions drop
#' them when applying omission, even with complete data.
#'
#' @seealso [multiScaleR-fitting-methods] for the backends,
#'   [multiScaleR-class] for the returned object, [kernel_prep()] to prepare
#'   the inputs, [profile_sigma()] for profile-likelihood intervals,
#'   [kernel_scale.raster()] to map the fitted covariates, [kernel_dist()] for
#'   effective distances, and
#'   `vignette("v01-quickstart", package = "multiScaleR")` and
#'   `vignette("v02-user-guide", package = "multiScaleR")` for worked
#'   examples.
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
multiScale_optim <- function(fitted_mod,
                             kernel_inputs,
                             join_by = NULL,
                             par = NULL,
                             start_strategy = NULL,
                             screen_n_sigma = 5,
                             screen_n_jitter = 6,
                             screen_maxit = 8,
                             screen_jitter_sd = 0.5,
                             n_cores = NULL,
                             PSOCK = FALSE,
                             verbose = TRUE,
                             refit_fn = NULL,
                             backend = c("refit", "profile_glm", "joint_nb2", "joint_mixed", "conditional_logit"),
                             profile_control = list(),
                             covariate_transform = c("candidate", "raw"),
                             radius_grid = NULL){
  backend <- match.arg(backend)
  if (!is.null(radius_grid) && backend != "refit")
    stop("radius_grid requires backend = 'refit'.", call. = FALSE)
  covariate_transform <- match.arg(covariate_transform)
  if (methods::is(fitted_mod, "unmarkedFitGDRField")) {
    if (!identical(backend, "refit") || !identical(PSOCK, FALSE) ||
        (!is.null(n_cores) && !(is.numeric(n_cores) && is.null(dim(n_cores)) &&
          length(n_cores) == 1L && is.finite(n_cores) && n_cores == 1)))
      stop("Native GDR response-field scale fitting requires backend = 'refit', n_cores = 1 and PSOCK = FALSE; native fit/tape worker transport is not qualified.", call. = FALSE)
    n_cores <- 1L
  }
  if (covariate_transform == "raw" && backend != "refit")
    stop("Raw optimized covariates currently require backend = 'refit'; analytic backends retain candidate standardization.", call. = FALSE)
  if(is.null(fitted_mod)){
    stop("`fitted_mod` must be a fitted model object.")
  }
  user_supplied_par <- !is.null(par)
  # A fixed kernel makes every optimized scale a hard radius. Its likelihood
  # changes in steps, so a gradient search cannot leave a single start. Screen
  # starting values by default in that case. An explicit choice always wins,
  # and a radius grid keeps its own exhaustive search.
  if (is.null(start_strategy))
    start_strategy <- if (identical(kernel_inputs$kernel, "fixed") &&
                          backend == "refit" && is.null(radius_grid))
      "screen" else "single"
  start_strategy <- match.arg(start_strategy, c("single", "screen"))
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
  if (!is.null(radius_grid))
    radius_grid <- .msr_refit_radius_grid(radius_grid, r_vars, kernel_inputs,
      fitted_mod, join_by, refit_fn, n_cores, PSOCK, start_strategy)

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
  if (covariate_transform == "raw") opt_context$covariate_transform <- "raw"

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

  grid_final_mod <- NULL
  if (!is.null(radius_grid)) {
    if (!identical(unname(r_vars), opt_context$covs))
      stop("The native model and radius grid disagree on fitted spatial predictor order.", call. = FALSE)
    opt_results <- .msr_refit_radius_profile(radius_grid, kernel_inputs, opt_context,
      fitted_mod, par, lwr, uppr, user_supplied_par)
    grid_final_mod <- opt_results$final_mod
    radius_profile <- opt_results$radius_profile
    opt_results[c("final_mod", "radius_profile")] <- NULL
  } else if (use_parallel) {
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

    hard_radius <- .msr_hard_radius_covariates(kernel_inputs$scale_vars,
      kernel_inputs$kernel, r_vars)
    n_se <- n_covs * if (kernel_inputs$kernel == "expow") 2L else 1L
    active <- rep(!r_vars %in% hard_radius, length.out = n_se)
    uncertainty <- .msr_refit_hessian_se(opt_results$hessian_unscale, active,
      c(rep(kernel_inputs$unit_conv, n_covs),
        if (kernel_inputs$kernel == "expow") rep(1, n_covs)))
    if (!is.null(radius_grid) && any(active))
      uncertainty$diagnostic$uncertainty <- "outer_refit_wald_conditional_on_selected_hard_radii"
    if (identical(uncertainty$diagnostic$status, "hessian_includes_failed_refit"))
      warning("A model refit near the estimated scales failed while the Hessian was ",
              "computed, so Wald standard errors for the scales are withheld. The ",
              "estimates stand, but the likelihood near them is fragile. See ",
              "diagnostics(fit)$refit_hessian and use summary(fit, profile = TRUE) ",
              "for profile-likelihood limits.", call. = FALSE)
    bounds <- .msr_refit_bounds(opt_results$par, lwr, uppr, active,
      c(rep(kernel_inputs$unit_conv, n_covs),
        if (kernel_inputs$kernel == "expow") rep(1, n_covs)),
      c(paste0("sigma:", r_vars),
        if (kernel_inputs$kernel == "expow") paste0("shape:", r_vars)))
    uncertainty <- .msr_refit_bound_uncertainty(uncertainty, bounds)
    res <- uncertainty$SE

    scale_est <- data.frame(Mean = opt_results$par_unscale[1:n_covs],
                            SE = res[1:n_covs])
    rownames(scale_est) <- r_vars
    hard_radius <- .msr_hard_radius_covariates(kernel_inputs$scale_vars,
      kernel_inputs$kernel, r_vars)
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

    final_mod <- grid_final_mod
    if (is.null(final_mod)) final_mod <- kernel_scale_fn(par = c(opt_results$par),
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
                opt_mod = .msr_native_fixed_metadata(final_mod$mod),
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
    if (!is.null(radius_grid)) out$diagnostics$radius_profile <- radius_profile
    if (isTRUE(verbose)) .msr_print_refit_bounds(out$diagnostics)
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

    # A search that never left its start reports the starts as estimates.
    start_stall <- .msr_start_stall_diagnostic(
      opt_results = opt_results,
      start = par,
      lower = lwr,
      upper = uppr,
      unit_conv = kernel_inputs$unit_conv,
      covariates = r_vars,
      hard_only = all(r_vars %in% hard_radius),
      screened = !user_supplied_par && identical(start_strategy, "screen"),
      default_start = !user_supplied_par && identical(start_strategy, "single")
    )
    if (!is.null(start_stall)) {
      out$diagnostics$start_stall <- start_stall
      warning(start_stall$message, call. = FALSE)
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


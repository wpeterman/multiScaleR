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
#'   sigma parameters). Default: \code{NULL} — starting values are chosen
#'   automatically.
#' @param n_cores Positive integer. Number of cores for parallel optimization
#'   via \code{optimParallel}. Default: \code{NULL} (single-threaded). Parallel
#'   optimization is beneficial for models with many covariates or slow
#'   log-likelihood evaluations.
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
#'     for each covariate — the centering and scaling parameters from the
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
#' @importFrom parallel clusterEvalQ makeCluster setDefaultCluster stopCluster makeForkCluster clusterExport
#' @importFrom crayon %+% green red bold blue
#' @importFrom pscl zeroinfl
multiScale_optim <- function(fitted_mod,
                             kernel_inputs,
                             join_by = NULL,
                             par = NULL,
                             n_cores = NULL,
                             PSOCK = FALSE,
                             verbose = TRUE,
                             refit_fn = NULL){
  if(is.null(fitted_mod)){
    stop("`fitted_mod` must be a fitted model object.")
  }
  validate_scalar_logical(PSOCK, "PSOCK")
  validate_scalar_logical(verbose, "verbose")
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
  if(!is.list(kernel_inputs$raw_cov) || !is.list(kernel_inputs$d_list) ||
     length(kernel_inputs$raw_cov) == 0 || length(kernel_inputs$d_list) == 0 ||
     length(kernel_inputs$raw_cov) != length(kernel_inputs$d_list)){
    stop("`kernel_inputs$raw_cov` and `kernel_inputs$d_list` must be non-empty lists of equal length.")
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
    missing_sources <- setdiff(kernel_inputs$scale_vars$source,
                               colnames(kernel_inputs$raw_cov[[1]]))
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
                                   n_cols = kernel_inputs$n_cols)

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
              "The estimated scale of effect extends beyond the maximum distance specified.\n",
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
    return(out)
  }
}


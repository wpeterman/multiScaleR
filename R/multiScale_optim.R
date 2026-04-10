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
#' @description Function to conduct multiscale optimization
#' @param fitted_mod Model object of class glm, lm, gls, or unmarked
#' @param kernel_inputs Object created from running \code{\link[multiScaleR]{kernel_prep}}
#' @param join_by Default: NULL. A data frame containing the variable used to join spatial point data with observation data (see Details)
#' @param par Optional starting values for parameter estimation. If provided, should be divided by the `max_D` value to be appropriately scaled. Default: NULL
#' @param n_cores If attempting to optimize in parallel, the number of cores to use. Default: NULL
#' @param PSOCK Logical. If attempting to optimize in parallel on a Windows machine, a PSOCK cluster will be created. If using a Unix OS a FORK cluster will be created. You can force a Unix system to create a PSOCK cluster by setting to TRUE. Default: FALSE
#' @param verbose Logical. Print status of optimization to the console. Default: TRUE
#' @return Returns a list of class `multiScaleR` containing scale estimates, shape estimates (if using kernel = 'expow'), optimization results, and the final optimized model.
#' @details Identifies the kernel scale, and uncertainty of that scale, for each raster within the context of the fitted model provided. Summary methods use profile-likelihood confidence intervals for `sigma` when feasible, while reported standard errors remain Hessian-based approximations from the outer optimization.
#'
#' To ensure that fitted model function calls are properly parallelized, fit models directly from the packages. For example, fit a negative binomial distribution from the MASS package as `fitted_mod <- MASS::glm.nb(y ~ x, data = df)`
#'
#' There may situations when using `unmarked` where sites are sampled across multiple years, but spatial environmental values are relevant for all years. In this situation, you want to join the scaled landscape variables from each site to each observation at a site. This can be achieved by providing a data frame object containing the values (e.g. site names) that will be used to join spatial data to sites. The name of the column in the `join_by` data frame must match a column name in the data used to fit your `unmarked` model.
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
                             verbose = TRUE){
  if(is.null(fitted_mod)){
    stop("`fitted_mod` must be a fitted model object.")
  }
  validate_scalar_logical(PSOCK, "PSOCK")
  validate_scalar_logical(verbose, "verbose")

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

  # Extract variables from fitted model
  if (any(class(fitted_mod) == 'gls')) {
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- unlist(find_predictors(fitted_mod))
  } else if (any(grepl("^unmarked", class(fitted_mod)))) {
    mod_vars <- all.vars(formula(fitted_mod@formula))
  } else {
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- unlist(find_predictors(fitted_mod))
  }

  # Ensure model variables are in kernel_inputs
  r_vars <- mod_vars[mod_vars %in% colnames(kernel_inputs$raw_cov[[1]])]
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
    mod_vars <- unlist(find_predictors(fitted_mod))
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
  } else if(any(grepl("^unmarked", class(fitted_mod)))) {
    mod_class <- 'unmarked'
    mod_vars <- all.vars(formula(fitted_mod@formula))
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
    fitType <- fitted_mod@fitType
  } else {
    mod_class <- 'other'
    # mod_vars <- find_predictors(fitted_mod)[[1]]
    mod_vars <- unlist(find_predictors(fitted_mod))
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
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

  if(!(any(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]])))){
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
                                   join_by = join_by)

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

    i_hess <- try(solve(opt_results$hessian_unscale))

    if(class(i_hess)[1] != 'try-error'){
      res <- c(sqrt(diag(i_hess)[1:n_covs]) * kernel_inputs$unit_conv,
               sqrt(diag(i_hess)[(n_covs + 1):(n_covs * 2)]))

    } else {
      res <- rep('Inf', n_covs)

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
                scl_params = final_mod$scl_params,
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

    if(any((shape_est[,1] / shape_est[,2]) < 2, na.rm = T)){
      out$diagnostics$shape_precision <- .precision_diagnostic(shape_est, "shape_precision")
      out$warn_message <- c(out$warn_message, 3)
      cat(red("\n WARNING!!!\n",
              "The standard error of one or more `shape` estimates is >= 50% of the estimated mean value.\n",
              "Carefully assess if the Exponential Power kernel is appropriate, whether or not this variable is meaningful in your analysis, and interpret with caution.\n\n"))
    }
    return(out)
  }
}


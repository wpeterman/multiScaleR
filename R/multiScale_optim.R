#' @title Multiscale optimization
#' @description Function to conduct multiscale optimization
#' @param fitted_mod Model object of class glm, lm, gls, or unmarked
#' @param kernel_inputs Object created from running \code{\link[multiScaleR]{kernel_prep}}
#' @param method Optimizer to be used. Default: 'L-BFGS-B', which is the only optimization method supported by \code{\link[optimParallel]{optimParallel}}
#' @param par Optional starting values for parameter estimation. If provided, must be scaled by the `unit_conv` used during `kernel_prep` (which defaults to 1e5 if not specified), Default: NULL
#' @param n_cores If attempting to optimize in parallel, the number of cores to use. Default: NULL
#' @return Returns a list of class `multiScaleR` containing scale estimates, shape estimates (if using kernel = 'expow'), optimization results, and the final optimized model.
#' @details Identifies the kernel scale, and uncertainty of that scale, for each raster within the context of the fitted model provided.
#' @examples
#' ## NOT RUN
#' pts <- vect(cbind(c(3,5,7),
#'                   c(7,5,3)))
#'
#' mat_list <- list(r1 = rast(matrix(rnorm(100),
#'                                   nrow = 10)),
#'                  r2 = rast(matrix(rnorm(100),
#'                                   nrow = 10)))
#' rast_stack <- rast(mat_list)
#' kernel_inputs <- kernel_prep(pts = pts,
#'                              raster_stack = rast_stack,
#'                              max_D = 2,
#'                              kernel = 'gaussian',
#'                              sigma = NULL)
#' ## Example response data
#' y <- rnorm(3)
#'
#' ## Create data frame with raster variables
#' df <- data.frame(y = y,
#'                  kernel_inputs$kernel_dat)
#' mod <- glm(y ~ r1 + r2,
#'            data = df)
#'
#' ## NOTE: This code is only for demonstration
#' ## Optimization results will have no meaning
#' opt <- multiScale_optim(fitted_mod = mod,
#'                         kernel_inputs = kernel_inputs,
#'                         method ='L-BFGS-B',
#'                         par = NULL,
#'                         opt_parallel = FALSE,
#'                         n_cores = NULL)
#'
#' @usage
#' multiScale_optim(fitted_mod,
#'                  kernel_inputs,
#'                  method ='L-BFGS-B',
#'                  par = NULL,
#'                  opt_parallel = FALSE,
#'                  n_cores = NULL)
#' @seealso
#'  \code{\link[optimParallel]{optimParallel}}
#' @rdname multiScale_optim
#' @export
#' @importFrom optimParallel optimParallel
#' @importFrom stats optim
#' @importFrom parallel clusterEvalQ makeCluster setDefaultCluster stopCluster

multiScale_optim <- function(fitted_mod,
                             kernel_inputs,
                             method ='L-BFGS-B',
                             par = NULL,
                             opt_parallel = FALSE,
                             n_cores = NULL){

  ## Need to add line to selectively identify the formula variables that are spatial raster layers

  if(any(class(fitted_mod) == 'gls')){
    mod_class <- 'gls'
    mod_vars <- insight::find_predictors(fitted_mod)[[1]]
    # mod_vars <- all.vars(formula(fitted_mod)[-2])
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
  } else if(any(grepl("^unmarked", class(fitted_mod)))) {
    mod_class <- 'unmarked'
    mod_vars <- all.vars(formula(fitted_mod@formula))
    r_vars <- mod_vars[which(mod_vars %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(r_vars)
    fitType <- fitted_mod@fitType
  } else {
    mod_class <- 'glm'
    # mod_vars <- all.vars(formula(fitted_mod)[-2])
    mod_vars <- insight::find_predictors(fitted_mod)[[1]]
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


  while(class(opt_results) == 'try-error' & cnt < (length(par_starts))){
    if(is.numeric(n_cores) & isTRUE(opt_parallel)){
      ## Initiate parallel cluster
      cl <- makeCluster(n_cores)     # set the number of processor cores
      setDefaultCluster(cl = cl)     # set'cl'as default cluster
      if(mod_class == 'unmarked'){
        clusterEvalQ(cl, library('unmarked'))
      }

      opt_results <- try(optimParallel::optimParallel(par = par,
                                                      fn = kernel_scale_fn,
                                                      hessian = TRUE,
                                                      lower = lwr,
                                                      upper = uppr,
                                                      method = method,
                                                      fitted_mod = fitted_mod,
                                                      kernel_inputs = kernel_inputs,
                                                      control = list(maxit = 1000),
                                                      parallel = list(forward = F,
                                                                      loginfo = T)), silent = T)
      ## Stop parallel cluster
      setDefaultCluster(cl = NULL)
      on.exit(stopCluster(cl))

      if(class(opt_results) == 'try-error'){
        cat('\n\nParallel optimization failed; attempting standard optimization ')

        opt_results <- try(optim(par = par,
                                 fn = kernel_scale_fn,
                                 hessian = TRUE,
                                 lower = lwr,
                                 upper = uppr,
                                 method = method,
                                 fitted_mod = fitted_mod,
                                 control = list(maxit = 1000),
                                 kernel_inputs = kernel_inputs),
                           silent = T)
      }

      cnt <- cnt + 1

      if(class(opt_results) == 'try-error'){
        cat('\n\nAttempt ')
        cat(cnt)
        cat(" failed\n")
        cat('raw par =  ')
        cat(par)
        cat('\nunscaled par =  ')
        cat(par * kernel_inputs$unit_conv)
        cat('\n')

        if(kernel_inputs$kernel != 'expow'){
          par <- rep(par_starts[cnt], n_covs)
          # par <- exp(rep(par_starts[cnt], n_covs))
        }

        if(kernel_inputs$kernel == 'expow'){
          par <- rep(par_starts[cnt], n_covs)
          par <- c(par, rep(2, n_covs))
          # par <- exp(rep(par_starts[cnt], n_covs))
          # par <- c(par, sqrt(rep(2, n_covs)))
        }
      }
    } else {
      opt_results <- try(optim(par = par,
                               fn = kernel_scale_fn,
                               hessian = TRUE,
                               lower = lwr,
                               upper = uppr,
                               method = method,
                               fitted_mod = fitted_mod,
                               control = list(maxit = 1000),
                               kernel_inputs = kernel_inputs),
                         silent = T)

      cnt <- cnt + 1

      if(class(opt_results) == 'try-error'){
        cat('\n\nAttempt ')
        cat(cnt)
        cat(" failed\n")
        cat('raw par =  ')
        cat(par)
        cat('\nunscaled par =  ')
        cat(par * kernel_inputs$unit_conv)
        cat('\n')

        if(kernel_inputs$kernel != 'expow'){
          par <- rep(par_starts[cnt], n_covs)
          # par <- exp(rep(par_starts[cnt], n_covs))
        }

        if(kernel_inputs$kernel == 'expow'){
          par <- rep(par_starts[cnt], n_covs)
          par <- c(par, rep(2, n_covs))
          # par <- exp(rep(par_starts[cnt], n_covs))
          # par <- c(par, sqrt(rep(2, n_covs)))
        }
      }
    } # End if else
  } # End while

  if(class(opt_results) == 'try-error'){


    warning("\nOptimization has failed!\nConsider providing more informative starting values for optimization.")
    return(opt_results)
  } else {

    cat('\n\nOptimization complete\n')
# browser()
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

    final_mod <- kernel_scale_fn(par = c(opt_results$par),
                                 kernel_inputs = kernel_inputs,
                                 fitted_mod = fitted_mod,
                                 mod_return = TRUE)

    out <- list(scale_est = scale_est,
                shape_est = shape_est,
                optim_results = opt_results,
                opt_mod = final_mod$mod,
                min_D = kernel_inputs$min_D,
                kernel_inputs = kernel_inputs,
                scl_params = final_mod$scl_params,
                call = match.call())

    class(out) <- 'multiScaleR'
    return(out)
  }
}

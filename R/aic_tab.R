# AIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title multiScaleR model selection 
#' @description Function to create AIC(c) table of fitted models
#' @param mod_list List containing fitted `multiScaleR` objects
#' @param AICc Use second order AIC in ranking models (Default = TRUE). See Details
#' @param mod_names Optional. Specify names for fitted model objects. By default, the right hand side of the fitted `multiScaleR` model, in combination with the kernel, will be used as the model name.
#' @param verbose (Default = FALSE) Should the table be printed to the console
#' @param ... Additional arguments (Not used)
#' 
#' @return Data frame of class `aictab` with AIC summary table for provided models
#' @export
#' 
#' @details
#' aic_tab creates a model selection table using \code{\link[AICcmodavg]{aictabCustom}} from the `AICcmodavg` package
#' 
#' @usage 
#' aic_tab(mod_list,
#'         AICc = TRUE,
#'         mod_names = NULL,
#'         verbose = FALSE,
#'         ...)
#' @rdname aic_tab
#' @importFrom insight find_formula get_loglikelihood get_parameters n_obs
#' @importFrom AICcmodavg aictabCustom bictabCustom
#' 
#' @examples 
#' ## Simulate data
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
#' mod1 <- glm(y ~ r1,
#'             data = dat)
#' mod2 <- glm(y ~ r2,
#'             data = dat)
#' mod3 <- glm(y ~ r1 + r2,
#'             data = dat)             
#' 
#' ## NOTE: This code is only for demonstration
#' ## Optimization results will have no meaning
#' 
#' opt_mod1 <- multiScale_optim(fitted_mod = mod1,
#'                              kernel_inputs = kernel_inputs,
#'                              method ='L-BFGS-B',
#'                              par = NULL,
#'                              opt_parallel = FALSE,
#'                              n_cores = NULL)
#'                              
#' opt_mod2 <- multiScale_optim(fitted_mod = mod2,
#'                              kernel_inputs = kernel_inputs,
#'                              method ='L-BFGS-B',
#'                              par = NULL,
#'                              opt_parallel = FALSE,
#'                              n_cores = NULL)
#' opt_mod3 <- multiScale_optim(fitted_mod = mod3,
#'                              kernel_inputs = kernel_inputs,
#'                              method ='L-BFGS-B',
#'                              par = NULL,
#'                              opt_parallel = FALSE,
#'                              n_cores = NULL)   
#'                              
#'  ## AIC table
#'  mod_list <- list(opt_mod1, opt_mod2, opt_mod3)
#'  
#'  aic_tab(mod_list = mod_list,
#'          AICc = FALSE)
#'          
#'  ## AICc table with specified names           
#'  aic_tab(mod_list = mod_list,
#'          AICc = TRUE,
#'          mod_names = c('mod1', 'mod2', 'mod3'))                                                
aic_tab <- function(mod_list,
                    AICc = TRUE,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){
  
  p <- list(...)
  
  ## All models comparable
  mod_dims <- as.vector(lapply(mod_list, function(x) insight::n_obs(x$opt_mod)))
  if(length(unique.default(mod_dims)) != 1L) {
    stop(cat("\n You are attempting to compare models with different number of sample locations. These are not valid comparisons. \n"))
  }
  
  if(length(unique(sapply(mod_list, function(x) class(x)))) != 1L){
    stop(cat("\n All objects in list must be of class `multiScaleR` \n"))
  }
  opt_list <- lapply(mod_list, function(x) x$opt_mod)
  
  mod_eq <- as.vector(sapply(opt_list, function(x) (insight::find_formula(x)$conditional)[-2]))
  mod_kernel <- as.vector(sapply(mod_list, function(x) x$kernel_inputs$kernel))
  k <- as.vector(sapply(opt_list, function(x) nrow(insight::get_parameters(x))[1] + 1))
  k[(mod_kernel %in% 'expow')] <- k[(mod_kernel %in% 'expow')] + 1
  
  if(is.null(mod_names)){
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }
  
  mod_loglik <- as.vector(sapply(opt_list, function(x) insight::get_loglikelihood(x)))
  mod_df <- as.vector(sapply(opt_list, function(x) insight::n_obs(x)))
  # mod_AIC <- as.vector(sapply(opt_list, function(x) x$aic))
  
  
  # AIC ---------------------------------------------------------------------
  if(AICc == FALSE) {
    
    tab <- AICcmodavg::aictabCustom(logL = mod_loglik,
                                    K = k,
                                    modnames = mod_names,
                                    second.ord = FALSE,
                                    nobs = mod_df,
                                    sort = TRUE)
    
    if(isTRUE(verbose)){
      return(print(tab))
    } else {
      return(tab)
    }
    
    # AICc --------------------------------------------------------------------
    
  } else {
    
    tab <- AICcmodavg::aictabCustom(logL = mod_loglik,
                                    K = k,
                                    modnames = mod_names,
                                    second.ord = TRUE,
                                    nobs = mod_df,
                                    sort = TRUE)
    
    if(isTRUE(verbose)){
      return(print(tab))
    } else {
      return(tab)
    }
  } 
} ## End function


# BIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title multiScaleR model selection 
#' @description Function to create BIC table of fitted models
#' @param mod_list List containing fitted `multiScaleR` objects
#' @param mod_names Optional. Specify names for fitted model objects. By default, the right hand side of the fitted `multiScaleR` model, in combination with the kernel, will be used as the model name.
#' @param verbose (Default = FALSE) Should the table be printed to the console
#' @param ... Additional arguments (Not used)
#' 
#' @return Data frame of class `bictab` with BIC summary table for provided models
#' @export
#' 
#' @details
#' bic_tab creates a model selection table using \code{\link[AICcmodavg]{bictabCustom}} from the `AICcmodavg` package
#' 
#' @usage 
#' bic_tab(mod_list,
#'         mod_names = NULL,
#'         verbose = FALSE,
#'         ...)
#' @rdname bic_tab
#' 
#' @examples 
#' ## Simulate data
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
#' mod1 <- glm(y ~ r1,
#'             data = dat)
#' mod2 <- glm(y ~ r2,
#'             data = dat)
#' mod3 <- glm(y ~ r1 + r2,
#'             data = dat)             
#' 
#' ## NOTE: This code is only for demonstration
#' ## Optimization results will have no meaning
#' 
#' opt_mod1 <- multiScale_optim(fitted_mod = mod1,
#'                              kernel_inputs = kernel_inputs,
#'                              method ='L-BFGS-B',
#'                              par = NULL,
#'                              opt_parallel = FALSE,
#'                              n_cores = NULL)
#'                              
#' opt_mod2 <- multiScale_optim(fitted_mod = mod2,
#'                              kernel_inputs = kernel_inputs,
#'                              method ='L-BFGS-B',
#'                              par = NULL,
#'                              opt_parallel = FALSE,
#'                              n_cores = NULL)
#' opt_mod3 <- multiScale_optim(fitted_mod = mod3,
#'                              kernel_inputs = kernel_inputs,
#'                              method ='L-BFGS-B',
#'                              par = NULL,
#'                              opt_parallel = FALSE,
#'                              n_cores = NULL)   
#'                              
#'  ## BIC table
#'  mod_list <- list(opt_mod1, opt_mod2, opt_mod3)
#'  
#'  bic_tab(mod_list = mod_list)
#'          
#'  ## BIC table with specified names        
#'  bic_tab(mod_list = mod_list,
#'          mod_names = c('mod1', 'mod2', 'mod3'))

bic_tab <- function(mod_list,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){
  
  p <- list(...)
  
  ## All models comparable
  mod_dims <- as.vector(lapply(mod_list, function(x) insight::n_obs(x$opt_mod)))
  if(length(unique.default(mod_dims)) != 1L) {
    stop(cat("\n You are attempting to compare models with different number of sample locations. These are not valid comparisons. \n"))
  }
  
  if(length(unique(sapply(mod_list, function(x) class(x)))) != 1L){
    stop(cat("\n All objects in list must be of class `multiScaleR` \n"))
  }
  opt_list <- lapply(mod_list, function(x) x$opt_mod)
  
  mod_eq <- as.vector(sapply(opt_list, function(x) (insight::find_formula(x)$conditional)[-2]))
  mod_kernel <- as.vector(sapply(mod_list, function(x) x$kernel_inputs$kernel))
  k <- as.vector(sapply(opt_list, function(x) nrow(insight::get_parameters(x))[1] + 1))
  k[(mod_kernel %in% 'expow')] <- k[(mod_kernel %in% 'expow')] + 1
  
  if(is.null(mod_names)){
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }
  
  mod_loglik <- as.vector(sapply(opt_list, function(x) insight::get_loglikelihood(x)))
  mod_df <- as.vector(sapply(opt_list, function(x) insight::n_obs(x)))
  # mod_AIC <- as.vector(sapply(opt_list, function(x) x$aic))
  
  
  # * BIC ---------------------------------------------------------------------
  
  tab <- AICcmodavg::bictabCustom(logL = mod_loglik,
                                  K = k,
                                  modnames = mod_names,
                                  nobs = mod_df,
                                  sort = TRUE)
  
  if(isTRUE(verbose)){
    return(print(tab))
  } else {
    return(tab)
  }
}  ## End function

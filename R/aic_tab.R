# AIC table ---------------------------------------------------------------
#' @author Bill Peterman
#' @title multiScaleR model selection
#'
#' @description Creates a model selection table ranked by AIC or AICc for a
#' list of fitted models. Mixed lists containing both \code{multiScaleR}
#' objects and plain model objects (e.g., \code{glm}) are supported, allowing
#' comparison of scale-optimized models against fixed-scale alternatives.
#'
#' @param mod_list List containing fitted \code{multiScaleR} objects (from
#'   \code{\link{multiScale_optim}}) and/or plain fitted model objects (e.g.,
#'   \code{glm}, \code{lm}). All models must have been fit to the same set of
#'   observations; an error is raised if sample sizes differ.
#' @param AICc Logical. If \code{TRUE} (default), the second-order corrected
#'   AIC (AICc) is used. Set to \code{FALSE} for standard AIC. AICc is
#'   recommended when the ratio of observations to parameters is small
#'   (n / K < 40 is a common rule of thumb).
#' @param mod_names Optional character vector of model names. Length must equal
#'   \code{length(mod_list)}. By default, names are constructed as
#'   \code{[kernel]right-hand-side-formula} for \code{multiScaleR} objects.
#' @param verbose Logical. If \code{TRUE}, the table is also printed to the
#'   console. Default: \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A data frame of class \code{"aictab"} (from the \code{AICcmodavg}
#' package) sorted by ascending AIC(c), containing:
#' \describe{
#'   \item{\code{Modnames}}{Model name (from \code{mod_names} or auto-generated).}
#'   \item{\code{K}}{Number of estimated parameters (regression coefficients +
#'     sigma, + shape for \code{"expow"} kernel).}
#'   \item{\code{AICc} or \code{AIC}}{Information criterion value.}
#'   \item{\code{Delta_AICc} or \code{Delta_AIC}}{Difference from the top model.}
#'   \item{\code{ModelLik}}{Relative likelihood (\code{exp(-0.5 * Delta)}).}
#'   \item{\code{AICcWt} or \code{AICWt}}{Akaike weight.}
#'   \item{\code{LL}}{Log-likelihood.}
#'   \item{\code{Cum.Wt}}{Cumulative Akaike weight.}
#' }
#'
#' @details
#' The table is constructed using \code{\link[AICcmodavg]{aictabCustom}} from
#' the \code{AICcmodavg} package. For \code{multiScaleR} objects, the optimized
#' (refitted) model stored in \code{opt_mod} is used for log-likelihood and
#' parameter count. The sigma parameter (and shape for \code{"expow"}) is added
#' to K automatically. For plain model objects in a mixed list, sigma is not
#' added to K.
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
#'                              par = NULL,
#'                              n_cores = NULL)
#'
#' opt_mod2 <- multiScale_optim(fitted_mod = mod2,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
#'                              n_cores = NULL)
#' opt_mod3 <- multiScale_optim(fitted_mod = mod3,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
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

  class_list <- lapply(mod_list, function(x) class(x))
  msclr <- which(class_list == 'multiScaleR')
  opt_list <- lapply(mod_list[msclr], function(x) x$opt_mod)

  if(length(msclr) != length(mod_list)){
    opt_list <- c(opt_list, mod_list[-msclr])
  }

  ## All models comparable
  mod_dims <- as.vector(lapply(opt_list, function(x) n_obs(x)))
  if(length(unique.default(mod_dims)) != 1L) {
    stop("\nYou are attempting to compare models with different number of sample locations. These are not valid comparisons.\n")
  }

  mod_eq <- as.vector(sapply(opt_list, function(x) (find_formula(x)$conditional)[-2]))
  mod_kernel <- as.vector(sapply(mod_list[msclr], function(x) x$kernel_inputs$kernel))
  if(length(msclr) != length(mod_list)){
    mod_kernel <- c(mod_kernel, rep('NA', length(mod_list) - length(msclr)))
  }
  if(any(sapply(opt_list, function(x) any(grepl("^unmarked", class(x)))))){
    k <- as.vector(sapply(opt_list[msclr], function(x) length(all.vars(formula(x@formula)))+1))
    mod_df <- as.vector(sapply(opt_list, function(x) dim(x@data@siteCovs)[1]))
    if(length(msclr) != length(mod_list)){
      k2 <- as.vector(sapply(opt_list[-msclr], function(x) length(all.vars(formula(x@formula)))))
      k <- c(k, k2)

    }
  } else {
    k <- as.vector(sapply(opt_list[msclr], function(x) nrow(get_parameters(x))[1])+1)
    mod_df <- as.vector(sapply(opt_list, function(x) n_obs(x)))
    if(length(msclr) != length(mod_list)){
      k2 <- as.vector(sapply(opt_list[-msclr], function(x) nrow(get_parameters(x))[1]))
      k <- c(k, k2)

    }
  }

  k[(mod_kernel %in% 'expow')] <- k[(mod_kernel %in% 'expow')] + 1

  if(is.null(mod_names)){
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }

  mod_loglik <- as.vector(sapply(opt_list, function(x) logLik(x)))

  # AIC ---------------------------------------------------------------------
  if(AICc == FALSE) {

    tab <- aictabCustom(logL = mod_loglik,
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

    tab <- aictabCustom(logL = mod_loglik,
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
#'
#' @description Creates a model selection table ranked by BIC for a list of
#' fitted models. Mixed lists containing both \code{multiScaleR} objects and
#' plain model objects are supported.
#'
#' @param mod_list List containing fitted \code{multiScaleR} objects (from
#'   \code{\link{multiScale_optim}}) and/or plain fitted model objects. All
#'   models must have been fit to the same number of observations.
#' @param mod_names Optional character vector of model names. Length must equal
#'   \code{length(mod_list)}. By default, names are constructed as
#'   \code{[kernel]right-hand-side-formula} for \code{multiScaleR} objects.
#' @param verbose Logical. If \code{TRUE}, the table is also printed to the
#'   console. Default: \code{FALSE}.
#' @param ... Additional arguments (not used).
#'
#' @return A data frame of class \code{"bictab"} (from the \code{AICcmodavg}
#' package) sorted by ascending BIC, containing:
#' \describe{
#'   \item{\code{Modnames}}{Model name.}
#'   \item{\code{K}}{Number of estimated parameters (regression coefficients +
#'     sigma, + shape for \code{"expow"} kernel).}
#'   \item{\code{BIC}}{Bayesian information criterion value.}
#'   \item{\code{Delta_BIC}}{Difference from the top model.}
#'   \item{\code{BICWt}}{BIC weight.}
#'   \item{\code{Cum.Wt}}{Cumulative BIC weight.}
#'   \item{\code{LL}}{Log-likelihood.}
#' }
#'
#' @details
#' The table is constructed using \code{\link[AICcmodavg]{bictabCustom}} from
#' the \code{AICcmodavg} package. BIC penalizes model complexity more heavily
#' than AIC as sample size grows, and is often preferred when the goal is
#' identifying the single best-supported model rather than model averaging.
#' Sigma (and shape for \code{"expow"}) is counted in K.
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
#'                              par = NULL,
#'                              n_cores = NULL)
#'
#' opt_mod2 <- multiScale_optim(fitted_mod = mod2,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
#'                              n_cores = NULL)
#' opt_mod3 <- multiScale_optim(fitted_mod = mod3,
#'                              kernel_inputs = kernel_inputs,
#'                              par = NULL,
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

  class_list <- lapply(mod_list, function(x) class(x))
  msclr <- which(class_list == 'multiScaleR')
  opt_list <- lapply(mod_list[msclr], function(x) x$opt_mod)

  if(length(msclr) != length(mod_list)){
    opt_list <- c(opt_list, mod_list[-msclr])
  }

  ## All models comparable
  mod_dims <- as.vector(lapply(opt_list, function(x) n_obs(x)))
  if(length(unique.default(mod_dims)) != 1L) {
    stop("\nYou are attempting to compare models with different number of sample locations. These are not valid comparisons.\n")
  }

  mod_eq <- as.vector(sapply(opt_list, function(x) (find_formula(x)$conditional)[-2]))
  mod_kernel <- as.vector(sapply(mod_list[msclr], function(x) x$kernel_inputs$kernel))
  if(length(msclr) != length(mod_list)){
    mod_kernel <- c(mod_kernel, rep('NA', length(mod_list) - length(msclr)))
  }


  if(any(sapply(opt_list, function(x) any(grepl("^unmarked", class(x)))))){
    k <- as.vector(sapply(opt_list[msclr], function(x) length(all.vars(formula(x@formula)))+1))
    mod_df <- as.vector(sapply(opt_list, function(x) dim(x@data@siteCovs)[1]))
    if(length(msclr) != length(mod_list)){
      k2 <- as.vector(sapply(opt_list[-msclr], function(x) length(all.vars(formula(x@formula)))))
      k <- c(k, k2)

    }
  } else {
    k <- as.vector(sapply(opt_list[msclr], function(x) nrow(get_parameters(x))[1])+1)
    mod_df <- as.vector(sapply(opt_list, function(x) n_obs(x)))
    if(length(msclr) != length(mod_list)){
      k2 <- as.vector(sapply(opt_list[-msclr], function(x) nrow(get_parameters(x))[1]))
      k <- c(k, k2)

    }
  }
  k[(mod_kernel %in% 'expow')] <- k[(mod_kernel %in% 'expow')] + 1

  if(is.null(mod_names)){
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }

  mod_loglik <- as.vector(sapply(opt_list, function(x) logLik(x)))
  mod_df <- as.vector(sapply(opt_list, function(x) n_obs(x)))

  # * BIC ---------------------------------------------------------------------

  tab <- bictabCustom(logL = mod_loglik,
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

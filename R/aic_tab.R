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
.msr_selection_model <- function(model) {
  if (inherits(model, "multiScaleR")) {
    return(model$opt_mod)
  }

  model
}


.msr_kernel_label <- function(model) {
  if (inherits(model, "multiScaleR")) {
    return(model$kernel_inputs$kernel)
  }

  "NA"
}


.msr_model_nobs <- function(model) {
  analysis_mod <- .analysis_model(model)

  if (any(grepl("^unmarked", class(analysis_mod)))) {
    return(dim(analysis_mod@data@siteCovs)[1])
  }

  as.integer(n_obs(analysis_mod)[1])
}


.msr_parameter_count <- function(model) {
  analysis_mod <- .analysis_model(model)

  params <- tryCatch(
    get_parameters(analysis_mod),
    error = function(e) NULL
  )
  if (is.data.frame(params) && nrow(params) > 0) {
    return(nrow(params))
  }

  coefs <- tryCatch(
    stats::coef(analysis_mod),
    error = function(e) NULL
  )
  if (length(coefs) > 0) {
    return(length(coefs))
  }

  loglik_df <- tryCatch(
    attr(stats::logLik(analysis_mod), "df"),
    error = function(e) NULL
  )
  if (length(loglik_df) == 1 && is.finite(loglik_df)) {
    return(as.integer(loglik_df))
  }

  if (any(grepl("^unmarked", class(analysis_mod)))) {
    return(length(all.vars(formula(analysis_mod@formula))))
  }

  stop("Could not determine the number of fitted parameters for one or more models.",
       call. = FALSE)
}


.msr_selection_parameter_count <- function(model) {
  if (!inherits(model, "multiScaleR")) {
    return(.msr_parameter_count(model))
  }

  extra_k <- 0L
  if (!is.null(model$scale_est)) {
    extra_k <- extra_k + nrow(model$scale_est)
  }
  if (!is.null(model$shape_est)) {
    extra_k <- extra_k + nrow(model$shape_est)
  }

  .msr_parameter_count(model$opt_mod) + extra_k
}


.msr_valid_observation_ids <- function(ids, n_obs) {
  length(ids) == n_obs &&
    !anyNA(ids) &&
    all(nzchar(ids)) &&
    !anyDuplicated(ids)
}


.msr_observation_ids <- function(model) {
  analysis_mod <- .analysis_model(model)
  mod_n <- .msr_model_nobs(model)

  if (any(grepl("^unmarked", class(analysis_mod)))) {
    dat <- analysis_mod@data@siteCovs
    ids <- row.names(dat)
    if (.msr_valid_observation_ids(ids, mod_n)) {
      return(as.character(ids))
    }
    return(as.character(seq_len(mod_n)))
  }

  dat <- .model_data(analysis_mod, effects = "all")
  if (is.null(dat)) {
    dat <- extract_model_data(model)
  }

  if (is.data.frame(dat)) {
    ids <- row.names(dat)
    if (.msr_valid_observation_ids(ids, nrow(dat))) {
      return(as.character(ids))
    }
  }

  residual_ids <- tryCatch(
    names(stats::residuals(analysis_mod)),
    error = function(e) NULL
  )
  if (.msr_valid_observation_ids(residual_ids, mod_n)) {
    return(as.character(residual_ids))
  }

  fitted_ids <- tryCatch(
    names(stats::fitted(analysis_mod)),
    error = function(e) NULL
  )
  if (.msr_valid_observation_ids(fitted_ids, mod_n)) {
    return(as.character(fitted_ids))
  }

  NULL
}


.msr_validate_comparable_models <- function(mod_list) {
  mod_dims <- vapply(mod_list, .msr_model_nobs, integer(1))
  if (length(unique.default(mod_dims)) != 1L) {
    stop("\nYou are attempting to compare models with different number of sample locations. These are not valid comparisons.\n")
  }

  obs_ids <- lapply(mod_list, .msr_observation_ids)
  if (all(vapply(obs_ids, Negate(is.null), logical(1)))) {
    ref_ids <- sort(obs_ids[[1]])
    same_obs <- vapply(
      obs_ids[-1],
      function(x) identical(sort(x), ref_ids),
      logical(1)
    )
    if (length(same_obs) > 0 && !all(same_obs)) {
      stop("\nYou are attempting to compare models fit to different observation sets. These are not valid comparisons.\n")
    }
  }

  mod_dims
}


aic_tab <- function(mod_list,
                    AICc = TRUE,
                    mod_names = NULL,
                    verbose = FALSE,
                    ...){

  p <- list(...)

  opt_list <- lapply(mod_list, .msr_selection_model)
  mod_df <- .msr_validate_comparable_models(opt_list)
  mod_eq <- as.vector(sapply(opt_list, function(x) (find_formula(x)$conditional)[-2]))
  mod_kernel <- vapply(mod_list, .msr_kernel_label, character(1))
  k <- vapply(mod_list, .msr_selection_parameter_count, integer(1))

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

  opt_list <- lapply(mod_list, .msr_selection_model)
  mod_df <- .msr_validate_comparable_models(opt_list)
  mod_eq <- as.vector(sapply(opt_list, function(x) (find_formula(x)$conditional)[-2]))
  mod_kernel <- vapply(mod_list, .msr_kernel_label, character(1))
  k <- vapply(mod_list, .msr_selection_parameter_count, integer(1))

  if(is.null(mod_names)){
    mod_names <- paste0("[" ,mod_kernel, "]",mod_eq)
  }

  mod_loglik <- as.vector(sapply(opt_list, function(x) logLik(x)))
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

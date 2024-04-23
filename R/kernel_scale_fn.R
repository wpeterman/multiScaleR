#' @title Kernel scaling function
#' @description Function for internal use with optim
#' @param par list of parameters
#' @param kernel_inputs object created from `kernel_prep`
#' @param fitted_mod fitted model object
#' @param join_by Data frame to join unmarked frame during optimization
#' @param mod_return Default: NULL
#' @return Either estimated parameters or the fitted model using provided parameters
#' @details For internal use
#' @rdname kernel_scale_fn
#' @keywords internal
#' @importFrom insight get_data get_loglikelihood
#' @importFrom stats formula logLik
#' @importFrom unmarked logLik update

kernel_scale_fn <- function(par,
                            kernel_inputs,
                            fitted_mod,
                            join_by = NULL,
                            mod_return = NULL){


  # D <- kernel_inputs$D
  d_list <- kernel_inputs$d_list

  cov_df <- kernel_inputs$raw_cov
  kernel <- kernel_inputs$kernel
  n_ind <- length(d_list)
  mod <- fitted_mod

  # browser()

  if(any(class(mod) == 'gls')){
    mod_class <- 'gls'
    covs <- insight::find_predictors(mod)$conditional
    dat <- insight::get_data(mod)
    # dat <- nlme::getData(mod)
    # covs <- all.vars(formula(mod)[-2])
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  } else if(any(grepl("^unmarked", class(mod)))) {
    mod_class <- 'unmarked'
    dat <- mod@data@siteCovs
    covs <- all.vars(mod@formula)
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  } else if(any(class(mod) == 'glm')) {
    mod_class <- 'glm'
    covs <- insight::find_predictors(mod)$conditional
    dat <- insight::get_data(mod)
    dat0 <- insight::get_data(mod)
    # dat <- fitted_mod$data
    # covs <- all.vars(formula(mod)[-2])
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  } else {
    mod_class <- 'other'
    covs <- insight::find_predictors(mod)$conditional
    dat <- insight::get_data(mod)
    covs <- covs[which(covs %in% colnames(kernel_inputs$raw_cov[[1]]))]
    n_covs <- length(covs)
  }

  if(is.null(dat)){
    stop("Data from original model not saved to data frame. Try using `glm`.\n UDPATE FUNCTION TO GENERALIZE!!!")
  }


  cov.w <- vector('list', n_ind)
  sigma <- par[1:n_covs]
  if(kernel == 'expow'){
    shape <- par[(n_covs + 1):(n_covs * 2)]
  } else {
    shape <- NULL
  }

  if(any(sigma < 0)){
    obj <- 1e6^10
    return(obj)
  }

  for(i in 1:n_ind){

    cov.w[[i]] <-
      scale_type(d_list[[i]],
                 kernel = kernel,
                 sigma = sigma,
                 shape = shape,
                 r_stack.df = cov_df[[i]][,covs])
  } ## End for loop

  # browser()

  df <- data.frame(do.call(rbind, cov.w))
  colnames(df) <- covs
  if(mod_class == 'unmarked'){
    umf <- mod@data
    scl_df <- scale(df)
    # scl_df <- (df)
    if(!is.null(join_by)){
      drop_cols <- which(colnames(umf@siteCovs) %in% covs)
      umf@siteCovs <- umf@siteCovs[,-drop_cols]
      scl_df <- data.frame(scl_df, join_by)
      umf@siteCovs <- merge(umf@siteCovs, scl_df, by = colnames(join_by))
      mod_u <- update(mod, data = umf)
    } else {
      umf@siteCovs[,covs] <- scl_df
      mod_u <- update(mod, data = umf)
    }

  } else {
    scl_df <- scale(df)
    # scl_df <- (df)
    dat[,covs] <- as.data.frame(scl_df)
    mod_u <- update(mod, data = dat)

    ## For DEBUGGING
    # mod_u <- try(update(mod, data = dat))
    # if(class(mod_u)[[1]] == 'try-error'){
    #   browser()
    # }
    # browser()

  }

  # browser()

  if(is.null(mod_return)){
    obj <- data.frame()
    class(obj) <- 'try-error'

    if(mod_class == 'unmarked'){
      obj <- try(mod_u@negLogLike)
    }

    if(class(obj) == 'try-error'){
      obj <- try(logLik(mod_u)[1] * -1)
    }

    if(class(obj) == 'try-error'){
      obj <- insight::get_loglikelihood(mod_u)[1] * -1
    }

  } else {
    obj <- list(mod = mod_u,
                scl_params = list(mean = attr(scl_df, "scaled:center"),
                                  sd = attr(scl_df, "scaled:scale")))
    # obj <- list(mod = mod,
    #             scl_params = NULL)
  }
  return(obj)

}

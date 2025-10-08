#' @title Kernel scaling function
#' @description Function for internal use with optim
#' @param par list of parameters
#' @param d_list List of distance vectors
#' @param cov_df List of data frames with values extracted from rasters
#' @param kernel Kernel used
#' @param fitted_mod fitted model object
#' @param join_by Data frame to join unmarked frame during optimization
#' @param mod_return Default: NULL
#' @return Either estimated parameters or the fitted model using provided parameters
#' @details For internal use
#' @rdname kernel_scale_fn
#' @keywords internal
#' @importFrom insight get_data get_loglikelihood find_predictors
#' @importFrom stats formula logLik model.frame
#' @importFrom unmarked logLik update
#' @useDynLib multiScaleR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods is
kernel_scale_fn <- function(par,
                            d_list,
                            cov_df,
                            kernel,
                            fitted_mod,
                            join_by = NULL,
                            mod_return = NULL){

  n_ind <- length(d_list)
  mod <- fitted_mod

  # browser()

  if(any(class(mod) == 'gls')){
    mod_class <- 'gls'
    # covs <- find_predictors(mod)$conditional
    covs <- unlist(find_predictors(mod))
    dat <- get_data(mod)
    covs <- covs[which(covs %in% colnames(cov_df[[1]]))]
    n_covs <- length(covs)
  } else if(any(grepl("^unmarked", class(mod)))) {
    mod_class <- 'unmarked'
    dat <- mod@data@siteCovs
    covs <- all.vars(mod@formula)
    covs <- covs[which(covs %in% colnames(cov_df[[1]]))]
    n_covs <- length(covs)
  } else if(any(class(mod) == 'glm')) {
    mod_class <- 'glm'
    # covs <- find_predictors(mod)$conditional
    covs <- unlist(find_predictors(mod))
    dat <- get_data(mod)
    dat0 <- get_data(mod)
    covs <- covs[which(covs %in% colnames(cov_df[[1]]))]
    n_covs <- length(covs)
  } else {
    mod_class <- 'other'
    dat <- get_data(mod, effects = 'all')
    if(is.null(dat)){
      dat <- extract_model_data(mod)
    }
    # covs <- find_predictors(mod)$conditional
    covs <- unlist(find_predictors(mod))
    covs <- covs[which(covs %in% colnames(cov_df[[1]]))]
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
    if(n_covs == 1){
      cov.w[[i]] <-
        scale_type(d_list[[i]],
                   kernel = kernel,
                   sigma = sigma,
                   shape = shape,
                   r_stack.df = cov_df[[i]])
    } else {
      cov.w[[i]] <-
        scale_type(d_list[[i]],
                   kernel = kernel,
                   sigma = sigma,
                   shape = shape,
                   r_stack.df = cov_df[[i]][,covs])
    }
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

    if(inherits(obj, "try-error")){
      obj <- try(logLik(mod_u)[1] * -1)
    }

    if(inherits(obj, "try-error")){
      obj <- get_loglikelihood(mod_u)[1] * -1
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


## Custom data extraction
extract_model_data <- function(model) {
  # Try common extraction methods
  data <- tryCatch({
    # Method 1: model.frame() (works for most stats models)
    mf <- model.frame(model)
    if (!is.null(mf)) return(mf)

    # Method 2: Check for @frame (lme4, glmmTMB)
    if (is(model, "merMod") || is(model, "glmmTMB")) {
      return(model@frame)
    }

    # Method 3: Check for $data (some packages store it here)
    if (!is.null(model$data)) {
      return(model$data)
    }

    # Method 4: Try eval(model$call$data) (if data was passed in the call)
    if (!is.null(model$call$data)) {
      data_name <- as.character(model$call$data)
      if (exists(data_name, envir = parent.frame())) {
        return(get(data_name, envir = parent.frame()))
      }
    }

    # Method 5: For spaMM models (specific checks)
    if (inherits(model, "HLfit")) {
      if (!is.null(model$data)) return(model$data)
      # Alternative for spaMM: model$fr
      if (!is.null(model$fr)) return(model$fr)
    }

    # Fallback: Return NULL if no method worked
    warning("Could not extract data from model object.")
    NULL
  }, error = function(e) {
    warning("Failed to extract data: ", e$message)
    NULL
  })

  return(data)
}

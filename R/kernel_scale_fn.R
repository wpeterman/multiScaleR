.refit_model <- function(model, data, opt_context) {
  if (!is.null(opt_context$refit_fn)) {
    refit <- opt_context$refit_fn(model = model,
                                  data = data,
                                  context = opt_context)
    if (is.null(refit)) {
      stop("`refit_fn` must return a fitted model object, not NULL.", call. = FALSE)
    }
    return(refit)
  }

  if (identical(opt_context$mod_class, "unmarked")) {
    return(unmarked::update(model, data = data))
  }

  stats::update(model, data = data)
}


.finite_first_numeric <- function(value) {
  value <- suppressWarnings(as.numeric(value)[1])
  length(value) == 1 && is.finite(value)
}


.neg_loglik_model <- function(model, mod_class) {
  if (identical(mod_class, "unmarked")) {
    value <- try(model@negLogLike, silent = TRUE)
    if (!inherits(value, "try-error") && .finite_first_numeric(value)) {
      return(as.numeric(value)[1])
    }

    value <- try(unmarked::logLik(model)[1] * -1, silent = TRUE)
    if (!inherits(value, "try-error") && .finite_first_numeric(value)) {
      return(as.numeric(value)[1])
    }
  }

  value <- try(stats::logLik(model)[1] * -1, silent = TRUE)
  if (!inherits(value, "try-error") && .finite_first_numeric(value)) {
    return(as.numeric(value)[1])
  }

  value <- try(insight::get_loglikelihood(model)[1] * -1, silent = TRUE)
  if (!inherits(value, "try-error") && .finite_first_numeric(value)) {
    return(as.numeric(value)[1])
  }

  stop("Could not extract a finite log-likelihood from the refitted model.", call. = FALSE)
}


#' @title Kernel scaling function
#' @description Function for internal use with optim
#' @param par list of parameters
#' @param d_list List of distance vectors
#' @param cov_df List of data frames with values extracted from rasters
#' @param kernel Kernel used
#' @param fitted_mod fitted model object
#' @param join_by Data frame to join unmarked frame during optimization
#' @param mod_return Default: NULL
#' @param opt_context Cached optimization context created internally by `multiScale_optim()`
#' @return Either estimated parameters or the fitted model using provided parameters
#' @details For internal use
#' @rdname kernel_scale_fn
#' @keywords internal
#' @importFrom insight get_data get_loglikelihood find_predictors
#' @importFrom stats formula logLik model.frame
#' @useDynLib multiScaleR, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods is
kernel_scale_fn <- function(par,
                            d_list,
                            cov_df,
                            kernel,
                            fitted_mod,
                            join_by = NULL,
                            mod_return = NULL,
                            opt_context = NULL){

  if(is.null(opt_context)){
    opt_context <- build_opt_context(fitted_mod = fitted_mod,
                                     cov_df = cov_df,
                                     join_by = join_by)
  }

  if (!is.null(opt_context$complete_idx)) {
    d_list <- d_list[opt_context$complete_idx]
    cov_df <- cov_df[opt_context$complete_idx]
  }

  n_ind <- length(d_list)

  mod <- opt_context$fitted_mod
  mod_class <- opt_context$mod_class
  covs <- opt_context$covs
  n_covs <- opt_context$n_covs

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

  cov.w <- matrix(NA_real_, nrow = n_ind, ncol = n_covs)
  colnames(cov.w) <- covs

  for(i in seq_len(n_ind)){
    cov.w[i, ] <-
      scale_type(d_list[[i]],
                 kernel = kernel,
                 sigma = sigma,
                 shape = shape,
                 r_stack.df = cov_df[[i]][,covs,drop = FALSE])
  } ## End for loop

  scl_df <- scale(cov.w)
  refit_error <- NULL

  if(mod_class == 'unmarked'){
    umf <- opt_context$umf_template

    mod_u <- tryCatch({
      if(!is.null(join_by)){
        scl_df_join <- data.frame(scl_df, join_by, check.names = FALSE)
        umf@siteCovs <- merge(umf@siteCovs, scl_df_join, by = opt_context$join_cols)
        .refit_model(mod, data = umf, opt_context = opt_context)
      } else {
        umf@siteCovs[opt_context$covs] <- as.data.frame(scl_df)
        .refit_model(mod, data = umf, opt_context = opt_context)
      }
    }, error = function(e) {
      refit_error <<- conditionMessage(e)
      NULL
    })

  } else {
    dat <- opt_context$data_template
    mod_u <- tryCatch({
      dat[opt_context$covs] <- as.data.frame(scl_df)
      .refit_model(mod, data = dat, opt_context = opt_context)
    }, error = function(e) {
      refit_error <<- conditionMessage(e)
      NULL
    })

    ## For DEBUGGING
    # mod_u <- try(update(mod, data = dat))
    # if(class(mod_u)[[1]] == 'try-error'){
    #   browser()
    # }
    # browser()

  }

  if (is.null(mod_u)) {
    if (is.null(mod_return)) {
      return(1e6^10)
    }

    stop(
      paste0(
        "Failed to refit the model with the optimized covariates. ",
        "Make sure the fitted model can be updated with a modified data argument. ",
        "Original error: ", refit_error
      ),
      call. = FALSE
    )
  }

  if(is.null(mod_return)){
    obj <- try(.neg_loglik_model(mod_u, mod_class), silent = TRUE)
    if(inherits(obj, "try-error") || !is.finite(obj)){
      return(1e6^10)
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


build_opt_context <- function(fitted_mod,
                              cov_df,
                              join_by = NULL,
                              refit_fn = NULL) {
  mod <- fitted_mod
  cov_names <- colnames(cov_df[[1]])

  if(any(class(mod) == 'gls')){
    mod_class <- 'gls'
    mod_vars <- unname(unlist(find_predictors(mod)))
    dat <- get_data(mod)
  } else if(any(grepl("^unmarked", class(mod)))) {
    mod_class <- 'unmarked'
    mod_vars <- all.vars(mod@formula)
    dat <- mod@data@siteCovs
  } else if(any(class(mod) == 'glm')) {
    mod_class <- 'glm'
    mod_vars <- unname(unlist(find_predictors(mod)))
    dat <- get_data(mod)
  } else {
    mod_class <- 'other'
    mod_vars <- unname(unlist(find_predictors(mod)))
    dat <- get_data(mod, effects = 'all')
    if(is.null(dat)){
      dat <- extract_model_data(mod)
    }
  }

  if(is.null(dat)){
    stop(
      paste0(
        "Could not recover the original model data needed to refit the model during optimization. ",
        "Refit the model with an explicit `data = ...` argument and try again."
      ),
      call. = FALSE
    )
  }

  covs <- unname(mod_vars[which(mod_vars %in% cov_names)])
  n_covs <- length(covs)
  if (n_covs == 0) {
    stop(
      "No raster covariates from `cov_df` were found in the fitted model. Check that raster layer names match the model terms.",
      call. = FALSE
    )
  }

  out <- list(fitted_mod = mod,
              mod_class = mod_class,
              covs = covs,
              n_covs = n_covs,
              refit_fn = refit_fn)

  complete_idx <- which(stats::complete.cases(dat))

  if(mod_class == 'unmarked'){
    umf_template <- mod@data
    if(!is.null(join_by)){
      if(!all(colnames(join_by) %in% colnames(umf_template@siteCovs))){
        stop(
          "Columns in `join_by` must also be present in the unmarked site covariates.",
          call. = FALSE
        )
      }
      drop_cols <- which(colnames(umf_template@siteCovs) %in% covs)
      if(length(drop_cols) > 0){
        umf_template@siteCovs <- umf_template@siteCovs[,-drop_cols,drop = FALSE]
      }
      out$join_cols <- colnames(join_by)
    } else {
      out$site_cov_idx <- match(covs, colnames(umf_template@siteCovs))
    }
    out$complete_idx <- complete_idx
    out$umf_template <- umf_template
  } else {
    out$data_template <- dat[complete_idx, , drop = FALSE]
    out$complete_idx <- complete_idx
    out$cov_idx <- match(covs, colnames(dat))
  }

  return(out)
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
    warning(
      "Could not extract data from the model object; optimization may fail unless the model was fit with an explicit `data = ...` argument."
    )
    NULL
  }, error = function(e) {
    warning(
      paste0(
        "Failed to extract model data needed for refitting: ",
        e$message
      )
    )
    NULL
  })

  return(data)
}

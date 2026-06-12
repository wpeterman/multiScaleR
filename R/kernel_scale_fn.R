.analysis_model <- function(model) {
  if (inherits(model, "fit_clogit") &&
      is.list(model) &&
      !is.null(model$model)) {
    return(model$model)
  }

  model
}


.has_nested_analysis_model <- function(model) {
  !identical(.analysis_model(model), model)
}


.model_predictors <- function(model) {
  vars <- tryCatch(
    unname(unlist(find_predictors(model))),
    error = function(e) character()
  )

  if (length(vars) > 0) {
    return(vars)
  }

  formula <- try(stats::formula(model), silent = TRUE)
  if (inherits(formula, "try-error") || length(formula) < 2) {
    return(character())
  }

  response_vars <- if (length(formula) >= 3) {
    all.vars(formula[[2]])
  } else {
    character()
  }
  setdiff(all.vars(formula), response_vars)
}


.model_data <- function(model, ...) {
  suppressWarnings(tryCatch(
    get_data(model, ...),
    error = function(e) NULL
  ))
}


.refit_from_call <- function(model, data) {
  refit_call <- try(stats::getCall(model), silent = TRUE)
  if (inherits(refit_call, "try-error") || is.null(refit_call)) {
    refit_call <- model$call
  }

  if (is.null(refit_call) || !is.call(refit_call)) {
    stop("The model does not contain a reusable fitting call.", call. = FALSE)
  }

  if (!is.null(refit_call$data)) {
    refit_call$data <- quote(data)
  }

  fit_fun <- refit_call[[1]]
  if (is.symbol(fit_fun)) {
    fun_name <- as.character(fit_fun)
    pkg <- extract_call_function_package(model)
    if (is.null(pkg) && inherits(model, c("coxph", "clogit"))) {
      pkg <- "survival"
    }

    if (!is.null(pkg) &&
        requireNamespace(pkg, quietly = TRUE) &&
        exists(fun_name, envir = asNamespace(pkg), mode = "function",
               inherits = FALSE)) {
      refit_call[[1]] <- as.call(
        list(as.name("::"), as.name(pkg), as.name(fun_name))
      )
    }
  }

  eval(refit_call, envir = list(data = data), enclos = parent.frame())
}


.update_model_data <- function(model, data) {
  refit <- try(stats::update(model, data = data), silent = TRUE)
  if (!inherits(refit, "try-error")) {
    return(refit)
  }

  call_refit <- try(.refit_from_call(model, data), silent = TRUE)
  if (!inherits(call_refit, "try-error")) {
    return(call_refit)
  }

  stop(conditionMessage(attr(refit, "condition")), call. = FALSE)
}


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

  if (.has_nested_analysis_model(model)) {
    model$model <- .update_model_data(.analysis_model(model), data = data)
    return(model)
  }

  .update_model_data(model, data = data)
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

  if (.has_nested_analysis_model(model)) {
    value <- try(stats::logLik(.analysis_model(model))[1] * -1, silent = TRUE)
    if (!inherits(value, "try-error") && .finite_first_numeric(value)) {
      return(as.numeric(value)[1])
    }

    value <- try(insight::get_loglikelihood(.analysis_model(model))[1] * -1,
                 silent = TRUE)
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


# Build the raw (unscaled) kernel-weighted covariate matrix that `kernel_scale_fn()`
# refits on. Extracted from `kernel_scale_fn()` so that callers which vary one
# covariate at a time (notably `profile_sigma()`) can recompute a single column and
# reuse the rest. `covariates = NULL` computes every model covariate, reproducing
# the original in-line behavior exactly; a character subset computes only those
# columns (same rows, same values).
.msr_kernel_cov_w <- function(par,
                              d_list,
                              cov_df,
                              kernel,
                              opt_context,
                              covariates = NULL) {
  binned <- opt_context$binned

  n_total <- if (!is.null(d_list)) {
    length(d_list)
  } else if (!is.null(binned)) {
    binned$n_points
  } else {
    0L
  }

  eval_idx <- NULL
  if (!is.null(opt_context$complete_idx)) {
    eval_idx <- opt_context$complete_idx
    if (length(eval_idx) == 0 ||
        anyNA(eval_idx) ||
        any(eval_idx < 1) ||
        any(eval_idx > n_total)) {
      eval_idx <- opt_context$data_idx
    }

    if (!is.null(eval_idx) && length(eval_idx) > 0 &&
        !anyNA(eval_idx) &&
        all(eval_idx >= 1) &&
        all(eval_idx <= n_total)) {
      if (!is.null(d_list)) d_list <- d_list[eval_idx]
      if (!is.null(cov_df)) cov_df <- cov_df[eval_idx]
    } else {
      eval_idx <- NULL
    }
  }

  n_ind <- if (!is.null(d_list)) {
    length(d_list)
  } else if (!is.null(eval_idx)) {
    length(eval_idx)
  } else {
    n_total
  }

  covs <- opt_context$covs
  n_covs <- opt_context$n_covs
  scale_vars <- opt_context$scale_vars

  sigma <- par[1:n_covs]
  if (kernel == 'expow') {
    shape <- par[(n_covs + 1):(n_covs * 2)]
  } else {
    shape <- NULL
  }

  target_covs <- if (is.null(covariates)) covs else covariates

  cov.w <- matrix(NA_real_, nrow = n_ind, ncol = length(target_covs))
  colnames(cov.w) <- target_covs

  # Resolve point/row identifiers, preferring the cell-level names and falling
  # back to the binned point ids (or the complete-case indices) in lean mode.
  row_ids <- if (!is.null(d_list)) {
    names(d_list)
  } else if (!is.null(binned) && !is.null(binned$point_ids)) {
    if (!is.null(eval_idx)) binned$point_ids[eval_idx] else binned$point_ids
  } else {
    NULL
  }
  if (is.null(row_ids) || length(row_ids) != n_ind ||
      anyNA(row_ids) || any(!nzchar(row_ids)) || anyDuplicated(row_ids)) {
    if (!is.null(eval_idx) && length(eval_idx) == n_ind) {
      row_ids <- as.character(eval_idx)
    } else {
      row_ids <- as.character(seq_len(n_ind))
    }
  }
  rownames(cov.w) <- row_ids

  # Kernel covariates that have precomputed binned summaries use the O(nbins)
  # fast path; everything else (landscape metrics, or legacy objects without
  # binned data) uses the exact per-cell evaluation.
  binned_covs <- if (!is.null(binned)) intersect(target_covs, binned$covariates) else character(0)
  opt_vars <- .msr_optimized_scale_vars(scale_vars)
  param_covariates <- if (is.null(opt_vars)) covs else covs[covs %in% opt_vars$covariate]

  if (length(binned_covs) == length(target_covs) && length(target_covs) > 0) {
    # Every requested covariate has a binned summary: pure fast path.
    cov.w[, binned_covs] <- .msr_binned_kernel_means(
      binned = binned,
      covs = binned_covs,
      param_covariates = param_covariates,
      sigma = sigma,
      shape = shape,
      kernel = kernel,
      eval_idx = eval_idx
    )
  } else {
    if (is.null(d_list) || is.null(cov_df)) {
      stop(
        "Cell-level data are required to evaluate one or more covariates but were not stored. Re-run `kernel_prep()` with `store_cell_data = TRUE`.",
        call. = FALSE
      )
    }

    # Exact per-cell evaluation over the requested covariates. `sig_idx` maps the
    # requested columns back to their positions in the full sigma/shape vectors so
    # a subset is weighted with exactly the same parameters as the full matrix.
    sig_idx <- match(target_covs, covs)
    for (i in seq_len(n_ind)) {
      if (is.null(scale_vars)) {
        cov.w[i, ] <-
          scale_type(d_list[[i]],
                     kernel = kernel,
                     sigma = sigma[sig_idx],
                     shape = if (is.null(shape)) NULL else shape[sig_idx],
                     r_stack.df = cov_df[[i]][, target_covs, drop = FALSE])
      } else {
        # `.msr_eval_scale_vars()` indexes `sigma`/`shape` by each covariate's
        # position within the `covariates` it is handed, so a subset must receive
        # the subset-ordered parameters (`sig_idx` maps back to the full vector).
        cov.w[i, ] <- .msr_eval_scale_vars(
          d = d_list[[i]],
          cov_df = cov_df[[i]],
          scale_vars = scale_vars,
          sigma = sigma[sig_idx],
          shape = if (is.null(shape)) NULL else shape[sig_idx],
          kernel = kernel,
          unit_conv = opt_context$unit_conv,
          resolution = opt_context$resolution,
          n_cols = opt_context$n_cols,
          covariates = target_covs,
          validate = FALSE
        )
      }
    } ## End for loop

    # Overwrite the kernel covariate columns with their (faster, equivalent)
    # binned values so mixed kernel + landscape models still benefit.
    if (length(binned_covs) > 0 && !is.null(scale_vars)) {
      cov.w[, binned_covs] <- .msr_binned_kernel_means(
        binned = binned,
        covs = binned_covs,
        param_covariates = param_covariates,
        sigma = sigma,
        shape = shape,
        kernel = kernel,
        eval_idx = eval_idx
      )
    }
  }

  cov.w
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
#' @param cov_w Optional precomputed (unscaled) kernel-weighted covariate matrix.
#'   When supplied, the per-covariate weighting step is skipped and this matrix is
#'   scaled and refit directly. Used by `profile_sigma()` to reuse the columns of
#'   covariates that are held at their optimum. Default: NULL (compute internally).
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
                            opt_context = NULL,
                            cov_w = NULL){

  if(is.null(opt_context)){
    opt_context <- build_opt_context(fitted_mod = fitted_mod,
                                     cov_df = cov_df,
                                     join_by = join_by)
  }

  mod <- opt_context$fitted_mod
  mod_class <- opt_context$mod_class
  n_covs <- opt_context$n_covs

  sigma <- par[1:n_covs]
  if (any(sigma < 0)) {
    return(1e6^10)
  }

  if (is.null(cov_w)) {
    cov.w <- .msr_kernel_cov_w(par = par,
                               d_list = d_list,
                               cov_df = cov_df,
                               kernel = kernel,
                               opt_context = opt_context)
  } else {
    cov.w <- cov_w
  }

  scale_validation_error <- tryCatch(
    {
      validate_covariates_before_scale(cov.w, context = "optimized covariates")
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(scale_validation_error)) {
    if (is.null(mod_return)) {
      return(1e6^10)
    }
    stop(scale_validation_error, call. = FALSE)
  }

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
    scl_df_dat <- align_scaled_covariates(scl_df = scl_df,
                                          dat = dat,
                                          opt_context = opt_context)
    mod_u <- tryCatch({
      dat[opt_context$covs] <- as.data.frame(scl_df_dat)
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


align_scaled_covariates <- function(scl_df, dat, opt_context) {
  dat_rows <- row.names(dat)
  scl_rows <- row.names(scl_df)

  if (length(dat_rows) == nrow(dat) && length(scl_rows) == nrow(scl_df) &&
      all(dat_rows %in% scl_rows)) {
    return(scl_df[dat_rows, , drop = FALSE])
  }

  idx <- opt_context$complete_idx
  if (!is.null(idx) && length(idx) == nrow(dat) &&
      all(!is.na(idx)) && all(idx >= 1) && all(idx <= nrow(scl_df))) {
    return(scl_df[idx, , drop = FALSE])
  }

  if (nrow(scl_df) == nrow(dat)) {
    return(scl_df)
  }

  stop(
    paste0(
      "Could not align kernel-weighted covariates with the fitted model data. ",
      "Make sure the model data rows correspond to the rows of `kernel_inputs`."
    ),
    call. = FALSE
  )
}


build_opt_context <- function(fitted_mod,
                              cov_df,
                              join_by = NULL,
                              refit_fn = NULL,
                              scale_vars = NULL,
                              unit_conv = 1,
                              resolution = NULL,
                              n_cols = NULL,
                              binned = NULL) {
  mod <- fitted_mod
  analysis_mod <- .analysis_model(mod)
  cov_names <- .msr_optimized_covariates(scale_vars, cov_df = cov_df)

  if(any(class(mod) == 'gls')){
    mod_class <- 'gls'
    mod_vars <- .model_predictors(analysis_mod)
    dat <- .model_data(analysis_mod)
  } else if(any(grepl("^unmarked", class(mod)))) {
    mod_class <- 'unmarked'
    mod_vars <- all.vars(mod@formula)
    dat <- mod@data@siteCovs
  } else if(any(class(mod) == 'glm')) {
    mod_class <- 'glm'
    mod_vars <- .model_predictors(analysis_mod)
    dat <- .model_data(analysis_mod)
  } else {
    mod_class <- 'other'
    mod_vars <- .model_predictors(analysis_mod)
    dat <- .model_data(analysis_mod, effects = 'all')
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
              refit_fn = refit_fn,
              scale_vars = scale_vars,
              unit_conv = unit_conv,
              resolution = resolution,
              n_cols = n_cols,
              binned = .msr_subset_binned(binned, covs))

  n_sites <- if (!is.null(cov_df)) {
    length(cov_df)
  } else if (!is.null(binned)) {
    binned$n_points
  } else {
    nrow(dat)
  }
  complete_cases <- .model_complete_indices(dat = dat,
                                            n_sites = n_sites)
  data_idx <- complete_cases$data_idx
  complete_idx <- complete_cases$complete_idx

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
    out$data_template <- dat[data_idx, , drop = FALSE]
    out$complete_idx <- complete_idx
    out$data_idx <- data_idx
    out$cov_idx <- match(covs, colnames(dat))
  }

  return(out)
}


.model_complete_indices <- function(dat, n_sites) {
  data_idx <- which(stats::complete.cases(dat))
  complete_idx <- data_idx

  row_ids <- row.names(dat)
  has_numeric_row_ids <- length(row_ids) == nrow(dat) &&
    all(grepl("^[0-9]+$", row_ids))

  if (isTRUE(has_numeric_row_ids)) {
    row_idx <- suppressWarnings(as.integer(row_ids))
    valid_row_idx <- !anyNA(row_idx) &&
      !anyDuplicated(row_idx) &&
      all(row_idx >= 1) &&
      all(row_idx <= n_sites)

    if (isTRUE(valid_row_idx)) {
      complete_idx <- row_idx[data_idx]
    }
  }

  list(data_idx = data_idx, complete_idx = complete_idx)
}


.recover_surv_data <- function(response, response_vars) {
  if (!inherits(response, "Surv") || length(response_vars) == 0) {
    return(list())
  }

  response <- as.matrix(response)
  response_data <- list()

  if (length(response_vars) == 1) {
    response_data[[response_vars[[1]]]] <- response[, ncol(response)]
  } else if (length(response_vars) == 2) {
    response_data[[response_vars[[1]]]] <- response[, 1]
    response_data[[response_vars[[2]]]] <- response[, ncol(response)]
  } else if (length(response_vars) >= 3 && ncol(response) >= 3) {
    response_data[[response_vars[[1]]]] <- response[, 1]
    response_data[[response_vars[[2]]]] <- response[, 2]
    response_data[[response_vars[[3]]]] <- response[, ncol(response)]
  }

  response_data
}


.model_frame_to_data <- function(model, mf) {
  if (is.null(mf) || !is.data.frame(mf)) {
    return(NULL)
  }

  formula <- try(stats::formula(model), silent = TRUE)
  if (inherits(formula, "try-error") || length(formula) < 2) {
    return(mf)
  }

  vars <- all.vars(formula)
  if (length(vars) == 0) {
    return(mf)
  }

  response_vars <- character()
  response_data <- list()
  terms <- try(stats::terms(model), silent = TRUE)
  response_idx <- 0
  if (!inherits(terms, "try-error")) {
    response_idx <- attr(terms, "response")
  }

  if (length(formula) >= 3 && isTRUE(response_idx > 0)) {
    response_vars <- all.vars(formula[[2]])
    response <- mf[[response_idx]]
    response_data <- .recover_surv_data(response, response_vars)
  }

  out <- data.frame(row.names = row.names(mf), check.names = FALSE)
  for (var in vars) {
    if (var %in% names(mf)) {
      out[[var]] <- mf[[var]]
    } else if (var %in% names(response_data)) {
      out[[var]] <- response_data[[var]]
    } else {
      strata_col <- paste0("strata(", var, ")")
      if (strata_col %in% names(mf)) {
        out[[var]] <- mf[[strata_col]]
      }
    }
  }

  if (all(vars %in% names(out))) {
    return(out)
  }

  mf
}


## Custom data extraction
extract_model_data <- function(model) {
  # Try common extraction methods
  data <- tryCatch({
    analysis_mod <- .analysis_model(model)

    # Method 1: model.frame() (works for most stats models)
    mf_error <- NULL
    mf <- tryCatch(
      model.frame(analysis_mod),
      error = function(e) {
        mf_error <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(mf)) {
      return(.model_frame_to_data(analysis_mod, mf))
    }

    # Some wrapped survival fits retain the original model frame here when
    # model.frame() cannot be reconstructed cleanly outside the fitting call.
    if (is.data.frame(analysis_mod$model)) {
      return(.model_frame_to_data(analysis_mod, analysis_mod$model))
    }

    # Method 2: Check for @frame (lme4, glmmTMB)
    if (is(analysis_mod, "merMod") && !is.null(analysis_mod@frame)) {
      return(analysis_mod@frame)
    }
    if (inherits(analysis_mod, "glmmTMB") && !is.null(analysis_mod$frame)) {
      return(analysis_mod$frame)
    }

    # Method 3: Check for $data (some packages store it here)
    if (!is.null(analysis_mod$data)) {
      return(analysis_mod$data)
    }

    # Method 4: Try eval(model$call$data) (if data was passed in the call)
    if (!is.null(analysis_mod$call$data)) {
      data_name <- as.character(analysis_mod$call$data)
      if (exists(data_name, envir = parent.frame())) {
        return(get(data_name, envir = parent.frame()))
      }
    }

    # Method 5: For spaMM models (specific checks)
    if (inherits(analysis_mod, "HLfit")) {
      if (!is.null(analysis_mod$data)) return(analysis_mod$data)
      # Alternative for spaMM: model$fr
      if (!is.null(analysis_mod$fr)) return(analysis_mod$fr)
    }

    # Fallback: Return NULL if no method worked
    if (!is.null(mf_error)) {
      warning(
        paste0(
          "Failed to extract model data needed for refitting: ",
          mf_error
        )
      )
    } else {
      warning(
        "Could not extract data from the model object; optimization may fail unless the model was fit with an explicit `data = ...` argument."
      )
    }
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

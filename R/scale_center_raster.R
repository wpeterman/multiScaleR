#' Scale and Center Raster Layers Using Model Parameters
#'
#' Standardizes raster covariates in a `terra::SpatRaster` using mean and standard deviation values extracted from a fitted model `multiScaleR` model.
#'
#' @param r A `terra::SpatRaster` containing covariate layers to be scaled. All layer names must match those used in the `multiScaleR` model.
#' @param multiScaleR A `multiScaleR` object created using `kernel_prep` or `multiScale_optim`.
#' @param clamp Logical. If `TRUE`, scaled values are clamped to the covariate range in the model data.
#' @param pct_mx Numeric greater than -0.5 and no larger than 0.99. If `clamp` is `TRUE`, each clamping limit moves outward by `pct_mx` times the observed range of that covariate at the survey points: `lower = min - pct_mx * (max - min)` and `upper = max + pct_mx * (max - min)`. For example, a covariate observed from -1 to 3 is clamped to -1.4 and 3.4 with `pct_mx = 0.1`, and to -0.6 and 2.6 with `pct_mx = -0.1`. Default = 0 (the observed range).
#'
#' @return A `terra::SpatRaster` with each layer scaled and optionally clamped.
#' @noRd
scale_center_raster <- function(r,
                                multiScaleR,
                                clamp = FALSE,
                                pct_mx = 0){
  if(!inherits(r, "SpatRaster")){
    stop("`r` must be a `SpatRaster` object.")
  }
  validate_multiScaleR_input(multiScaleR)
  validate_scalar_logical(clamp, "clamp")
  validate_pct_mx(pct_mx)

  covs <- names(r)
  if(!is.list(multiScaleR$scl_params) ||
     !all(c("mean", "sd") %in% names(multiScaleR$scl_params))){
    stop("`multiScaleR` must contain `scl_params` with `mean` and `sd` elements.")
  }

  if(inherits(multiScaleR, "multiScaleR_data")){

    mns <- multiScaleR$scl_params$mean[covs]
    sds <- multiScaleR$scl_params$sd[covs]
    dat <- multiScaleR$kernel_dat
    r_ <- vector('list', length(covs))
    names(r_) <- covs
  } else {
    mns <- multiScaleR$scl_params$mean[covs]
    sds <- multiScaleR$scl_params$sd[covs]
    r_ <- vector('list', length(covs))
    names(r_) <- covs
    mod <- multiScaleR$opt_mod

    if(isFALSE(all(names(mns) %in% names(r)))){
      stop('Rasters scaled in model absent from input r.')
    }

    if(any(class(mod) == 'gls')){
      mod_class <- 'gls'
      dat <- get_data(mod, effects = 'all')

    } else if(any(grepl("^unmarked", class(mod)))) {
      mod_class <- 'unmarked'
      dat <- mod@data@siteCovs

    } else {
      mod_class <- 'other'
      dat <- extract_model_data(mod)
      if(is.null(dat)){
        dat <- get_data(mod, effects = 'all')
      }
    }
  }

  if(any(is.na(mns)) || any(is.na(sds))){
    stop("Scaling parameters are missing for one or more raster layers in `r`.")
  }
  if(any(sds == 0)){
    stop("One or more scaling standard deviations are zero; cannot scale raster values.")
  }

  if(isTRUE(clamp)){
    if(!all(covs %in% names(dat))){
      stop("Original model data needed for clamping do not contain all raster layers in `r`.")
    }
    # Nonspatial site keys or factors must not coerce numeric raster columns.
    min_vals <- vapply(dat[covs], min, numeric(1))
    max_vals <- vapply(dat[covs], max, numeric(1))
    # Move each limit by a fraction of the observed range. Multiplying the
    # limits themselves by (1 + pct_mx) moved them the wrong way for negative
    # values and by unequal amounts when the range was not anchored at zero.
    pad <- pct_mx * (max_vals - min_vals)
    min_vals <- min_vals - pad
    max_vals <- max_vals + pad

    for(i in 1:length(covs)){
      lyr <- as.character(covs[[i]])
      r_s <- (r[[lyr]] - mns[[lyr]]) / sds[[lyr]]
      r_[[lyr]] <- clamp(r_s,
                         lower = min_vals[[lyr]],
                         upper = max_vals[[lyr]],
                         values = TRUE)
    }
  } else {
    for(i in 1:length(covs)){
      lyr <- as.character(covs[i])
      r_[[lyr]] <- (r[[lyr]] - mns[[lyr]]) / sds[[lyr]]

    }
  }
  r_ <- rast(r_)
  return(r_)
}

#' Validate the clamping expansion fraction
#'
#' `pct_mx` moves each clamping limit by `pct_mx` times the observed covariate
#' range, so values of -0.5 or less would make the lower limit meet or pass the
#' upper limit.
#'
#' @param pct_mx Value to validate.
#' @return `NULL`, invisibly. Called for its error.
#' @keywords internal
#' @noRd
validate_pct_mx <- function(pct_mx) {
  validate_scalar_numeric(pct_mx, "pct_mx")
  if (pct_mx <= -0.5 || pct_mx > 0.99) {
    stop("`pct_mx` must be greater than -0.5 and no larger than 0.99. ",
         "Each clamping limit moves by `pct_mx` times the observed range, ",
         "so -0.5 or less leaves no interval between the lower and upper limits.",
         call. = FALSE)
  }
  invisible(NULL)
}

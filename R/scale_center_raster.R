#' Scale and Center Raster Layers Using Model Parameters
#'
#' Standardizes raster covariates in a `terra::SpatRaster` using mean and standard deviation values extracted from a fitted model `multiScaleR` model.
#'
#' @param r A `terra::SpatRaster` containing covariate layers to be scaled. All layer names must match those used in the `multiScaleR` model.
#' @param multiScaleR A `multiScaleR` object created using `kernel_prep` or `multiScale_optim`.
#' @param clamp Logical. If `TRUE`, scaled values are clamped to the covariate range in the model data.
#' @param pct_mx Numeric. If `clamp` is `TRUE`, this value specifies the amount (percentage; positive or negative) by which to expand/contract the min/max range when clamping. Default = 0.
#'
#' @return A `terra::SpatRaster` with each layer scaled and optionally clamped.
#' @keywords internal
scale_center_raster <- function(r,
                                multiScaleR,
                                clamp = FALSE,
                                pct_mx = 0){
  if(!inherits(r, "SpatRaster")){
    stop("`r` must be a `SpatRaster` object.")
  }
  validate_multiScaleR_input(multiScaleR)
  validate_scalar_logical(clamp, "clamp")
  validate_scalar_numeric(pct_mx,
                          "pct_mx",
                          lower = -0.99,
                          upper = 0.99)

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
    min_vals <- apply(dat, 2, min) * (pct_mx + 1)
    max_vals <- apply(dat, 2, max) * (pct_mx + 1)

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

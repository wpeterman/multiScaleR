#' Scale and Center Raster Layers Using Model Parameters
#'
#' Standardizes raster covariates in a `terra::SpatRaster` using mean and standard deviation values extracted from a fitted model `multiScaleR` model.
#'
#' @param r A `terra::SpatRaster` containing covariate layers to be scaled. All layer names must match those used in the `multiScaleR` model.
#' @param multiScaleR A `multiScaleR` object created using `kernel_prep` or `multiScale_optim`.
#' @param clamp Logical. If `TRUE`, scaled values are clamped to the covariate range in the model data, expanded by `pct_mx`.
#' @param pct_mx Numeric. Proportion by which to expand the min/max range when clamping. Default is `0.15`.
#'
#' @return A `terra::SpatRaster` with each layer scaled and optionally clamped.
#' @keywords internal
#' @examples
#' ## NOT RUN
#'
#' # scaled_r <- scale_center_raster(r, opt_mod, clamp = TRUE)
#'
#' @usage
#' scale_center_raster(r,
#'                     multiScaleR,
#'                     clamp = FALSE,
#'                     pct_mx = 0.15)
#' @rdname scale_center_raster
scale_center_raster <- function(r,
                                multiScaleR,
                                clamp = FALSE,
                                pct_mx = 0.15){
  covs <- names(r)

  if(inherits(multiScaleR, "multiScaleR_data")){
    # covs_ <- names(multiScaleR$scl_params$mean)
    # covs <- intersect(covs, covs_)

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

    } else if(any(class(mod) == 'glm')){
      mod_class <- 'glm'
      dat <- mod$data

    } else {
      mod_class <- 'other'
      dat <- extract_model_data(mod)
      if(is.null(dat)){
        dat <- get_data(mod, effects = 'all')
      }
    }
  }

  if(isTRUE(clamp)){
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

#' @title Create scaled rasters
#' @description Function to create scaled rasters
#' @param raster_stack Stack of combined `SpatRaster` layers
#' @param sigma Vector of parameters listed in order to scale each raster
#' @param multiScaleR If scale optimization with `multiScale_optim` has been completed,
#' provide the `multiScaleR` object here. You can also pass an object of class
#' `"multiScaleR_data"` created using `kernel_prep`. Default: NULL
#' @param shape Vector of parameters listed in order to scale each raster if using
#' `'expow'` kernel. Default: NULL
#' @param kernel Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow';
#' Default: 'gaussian')
#' @param pct_wt The percentage of the weighted density to include when applying the
#' kernel smoothing function, Default: 0.975
#' @param fft Logical. If TRUE (Default), a fast Fourier transformation will be used
#' to smooth the raster surface. See details.
#' @param scale_center Logical. If `TRUE`, raster values are scaled and centered
#' according to the data used to fit the model. Necessary when predicting model
#' results across the landscape.
#' @param clamp Logical. If `TRUE`, scaled values are clamped to the covariate range
#' in the model data.
#' @param pct_mx Numeric. If `clamp` is `TRUE`, this value specifies the amount
#' (percentage; positive or negative) by which to expand or contract the min/max
#' range when clamping. Can range from -0.99–0.99 (Default = 0).
#' @param na.rm Logical. If TRUE (Default), NA values are removed from the weighted
#' mean calculation.
#' @param verbose Logical. Print status of raster scaling to the console. Default: TRUE
#' @param ... Not used
#'
#' @return A `SpatRaster` object containing scaled rasters. If a `multiScaleR` object
#' is provided and the fitted model includes additional site-level covariates that are
#' not represented in `raster_stack`, constant ("dummy") raster layers will be added
#' for those covariates to facilitate spatial prediction.
#'
#' @details
#' The fast Fourier transformation is substantially faster when scaling large raster
#' surfaces with large kernel areas. There will be some edge effects on the outer boundaries.
#'
#' When a fitted `multiScaleR` object is supplied, `kernel_scale.raster` inspects the
#' underlying model to identify all predictors used during model fitting. If the model
#' includes site-level covariates (i.e., predictors not associated with raster layers),
#' these variables are not spatially explicit and therefore cannot be directly projected.
#'
#' In such cases, the function automatically creates constant raster layers for these
#' covariates. By default, these layers are filled with a value of 0, which corresponds
#' to the reference value for centered and scaled predictors. This ensures compatibility
#' with the fitted model during prediction (e.g., when using `terra::predict`).
#'
#' Users should be aware that these dummy rasters do not represent spatial variation in
#' the covariate. Instead, they define a fixed projection scenario (e.g., predictions
#' conditional on the covariate being equal to 0). For non-centered variables or
#' alternative projection scenarios, users should manually modify or replace these layers.
#'
#' Categorical (factor) covariates are not automatically converted to raster layers.
#' If present in the model, a warning is issued and no dummy raster is created.
#'
#' @examples
#' ## Not Run
#' r1 <- rast(matrix(rnorm(25^2),
#'                   nrow = 25))
#'
#' r1_s <- kernel_scale.raster(r1,
#'                             sigma = 4,
#'                             kernel = 'gaussian')
#' plot(c(r1, r1_s))
#'
#' @rdname kernel_scale.raster
#' @export
#' @importFrom terra crs rast subset cellFromRowCol crds focal res xyFromCell nlyr clamp nrow ncol values as.matrix

kernel_scale.raster <- function(raster_stack,
                                sigma = NULL,
                                multiScaleR = NULL,
                                shape = NULL,
                                kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                                pct_wt = 0.975,
                                fft = TRUE,
                                scale_center = FALSE,
                                clamp = FALSE,
                                pct_mx = 0,
                                na.rm = TRUE,
                                verbose = TRUE,
                                ...){

  args <- list(...)

  if ("scale_opt" %in% names(args)) {
    warning("Argument 'scale_opt' is deprecated. Use 'multiScaleR' instead.", call. = FALSE)
    if (is.null(multiScaleR) & !is.null(args$scale_opt)) {
      multiScaleR <- args$scale_opt
    }
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")){
    covs <- rownames(multiScaleR$scale_est)
    sigma <- multiScaleR$scale_est[,1]
    shape <- multiScaleR$shape_est[,1]
    kernel <- multiScaleR$kernel_inputs$kernel

    if(!any(covs %in% names(raster_stack))){
      stop('optimized covariate(s) not present in the provided SpatRaster!')
    } else {
      var <- intersect(covs, names(raster_stack))
      c <- which(covs %in% var)
      raster_stack <- subset(raster_stack, var)
      sigma <- sigma[c]
      shape <- shape[c]
    }
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")){
    covs <- colnames(multiScaleR$kernel_dat)
    var <- intersect(covs, names(raster_stack))
    c <- which(var %in% covs)

    sigma <- multiScaleR$sigma[c] * multiScaleR$unit_conv
    shape <- multiScaleR$shape[c]
    kernel <- multiScaleR$kernel

    if(isFALSE(covs %in% names(raster_stack))){
      stop('optimized covariate is not present in the provided SpatRaster!')
    } else {
      raster_stack <- subset(raster_stack, covs)
    }
  } else {
    covs <- var <- names(raster_stack)
  }

  kernel <- match.arg(kernel)

  if(is.null(sigma)){
    stop("sigma values must be specified\n")
  }

  if(!inherits(raster_stack, "SpatRaster")){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }

  if(!is.logical(na.rm)){
    stop("`na.rm` must be TRUE or FALSE")
  }

  if(!is.logical(fft)){
    stop("`fft` must be TRUE or FALSE")
  }

  if(length(sigma) != nlyr(raster_stack)){
    warning("Number of sigma values must equal the number of raster layers!!!  \n  All raster layers will be smoothed using the same sigma value")
    sigma <- rep(sigma[1], nlyr(raster_stack))
  }


  smooth_list <- wt_list <-  vector('list', length(sigma))
  out <- rast(raster_stack[[1]])

  for(i in 1:length(sigma)){
    # lyr <- covs[i]
    lyr <- var[i]

    mx <- kernel_dist(kernel = kernel,
                      sigma = sigma[i],
                      shape = shape[i],
                      prob = pct_wt)

    r_res <- res(raster_stack)[1]
    focal_d <- ceiling(mx / r_res) * 2

    if((focal_d %% 2) == 0) {
      focal_d <- focal_d + 1
    }

    mat <- matrix(rnorm(focal_d^2), focal_d, focal_d)
    r_wt <- rast(mat)
    crs(r_wt) <- crs(raster_stack)
    cntr_crd <- cellFromRowCol(r_wt, focal_d/2, focal_d/2)
    cntr_crd <- xyFromCell(r_wt, ceiling(length(mat)/2))
    cell_crds <- crds(r_wt)
    r_wt[] <- rdist(cntr_crd, cell_crds)[1,] * r_res
    r_wt[] <- scale_type_r(d = as.vector(r_wt),
                           kernel = kernel,
                           sigma = sigma[i],
                           shape = shape[i],
                           output = 'wts')


    wt_mat <- as.matrix(r_wt, wide = T)

    if(verbose){
      cat(paste0("\nSmoothing spatRaster ",i, " of ", length(sigma), ": ",lyr," at sigma = ",floor(sigma[i]),"\n"))
    }
    if(isTRUE(fft)){
      mat <- as.matrix(raster_stack[[lyr]], wide = T)

      fft_mat <- fft_convolution(mat,
                                 wt_mat,
                                 fun = "mean",
                                 na.rm = na.rm)
      values(out) <- as.vector(fft_mat)

      smooth_list[[i]] <- out
    } else {
      smooth_list[[i]] <- focal(raster_stack[[lyr]],
                                w = wt_mat,
                                fun = "mean",
                                na.rm = na.rm,
                                expand = F)
    }

  }
  smooth_stack <- rast(smooth_list)
  names(smooth_stack) <- names(raster_stack)

  if (isTRUE(scale_center) &&
      (inherits(multiScaleR, "multiScaleR_data") || inherits(multiScaleR, "multiScaleR"))) {
    smooth_stack <- scale_center_raster(r = smooth_stack,
                                        multiScaleR = multiScaleR,
                                        clamp = clamp,
                                        pct_mx = pct_mx)
  }

  # add dummy rasters for site-level covariates in fitted model
  if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
    smooth_stack <- .add_site_covariate_rasters(
      smooth_stack = smooth_stack,
      multiScaleR = multiScaleR,
      raster_covs = var
    )
  }

  return(smooth_stack)
}


## Helper function
# helper: add constant rasters for non-raster covariates used by fitted model
.add_site_covariate_rasters <- function(smooth_stack, multiScaleR, raster_covs) {

  if (is.null(multiScaleR) || !inherits(multiScaleR, "multiScaleR")) {
    return(smooth_stack)
  }

  mod <- multiScaleR$opt_mod
  if (is.null(mod)) {
    return(smooth_stack)
  }

  mf <- try(stats::model.frame(mod), silent = TRUE)
  if (inherits(mf, "try-error")) {
    warning("Could not inspect fitted model frame; site-level covariates were not added.",
            call. = FALSE)
    return(smooth_stack)
  }

  # response name
  resp <- names(mf)[1]

  # raw variables in the fitted model frame, excluding response
  model_vars <- setdiff(names(mf), resp)

  # variables already represented by smoothed rasters
  site_vars <- setdiff(model_vars, raster_covs)

  if (length(site_vars) == 0) {
    return(smooth_stack)
  }

  # template raster
  tmpl <- smooth_stack[[1]]

  for (v in site_vars) {
    x <- mf[[v]]

    # factor / character / complex terms are not safe to add as dummy rasters
    if (is.factor(x) || is.character(x)) {
      warning(
        paste0(
          "Site-level covariate '", v, "' is categorical. ",
          "A dummy raster cannot be added safely without explicitly defining a projection value."
        ),
        call. = FALSE
      )
      next
    }

    # only allow atomic numeric/logical/integer covariates
    if (!(is.numeric(x) || is.integer(x) || is.logical(x))) {
      warning(
        paste0(
          "Site-level covariate '", v, "' is not numeric/integer/logical. ",
          "Skipping dummy raster creation."
        ),
        call. = FALSE
      )
      next
    }

    # create constant zero raster
    rr <- tmpl * 0
    names(rr) <- v
    smooth_stack <- c(smooth_stack, rr)
  }

  smooth_stack
}

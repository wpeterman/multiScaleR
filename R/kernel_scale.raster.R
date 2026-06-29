#' @title Apply kernel smoothing to raster layers
#'
#' @description Applies a kernel smoothing function to one or more raster
#' layers, producing a spatially weighted mean (or landscape metric) at the
#' scale identified by \code{\link{multiScale_optim}}. Primarily used to
#' generate prediction rasters for spatial model projection.
#'
#' @param raster_stack A \code{SpatRaster} object containing the source raster
#'   layer(s) to be smoothed. Layer names must match the covariate names used
#'   during model fitting (or the \code{source} names in \code{scale_vars}).
#' @param sigma Numeric vector of kernel scale parameter(s) in the same units
#'   as the raster projection (e.g., metres). One value per raster layer.
#'   Ignored when \code{multiScaleR} is provided; sigma values are extracted
#'   from the fitted object automatically.
#' @param multiScaleR A fitted object of class \code{"multiScaleR"} (from
#'   \code{\link{multiScale_optim}}) or class \code{"multiScaleR_data"} (from
#'   \code{\link{kernel_prep}}). When provided, sigma, shape, kernel, and
#'   scale_vars are all extracted automatically. Default: \code{NULL}.
#' @param shape Numeric vector of shape parameters for the exponential power
#'   kernel (\code{kernel = "expow"}), one per raster layer. Ignored when
#'   \code{multiScaleR} is provided. Default: \code{NULL}.
#' @param kernel Character. Kernel function used for smoothing. One of
#'   \code{"gaussian"} (default), \code{"exp"}, \code{"fixed"}, or
#'   \code{"expow"}. Ignored when \code{multiScaleR} is provided.
#' @param scale_vars Optional variable specifications created with
#'   \code{\link{msr_vars}}. Use when projecting landscape metrics or
#'   explicitly defined covariates without passing a fitted \code{multiScaleR}
#'   object. When \code{multiScaleR} is provided, \code{scale_vars} is
#'   extracted automatically.
#' @param pct_wt Numeric between 0 and 1 (exclusive). Cumulative kernel density
#'   cutoff used to determine the focal window size for smoothing. A larger
#'   value (e.g., \code{0.99}) captures more of the kernel tail but increases
#'   computation time. Default: \code{0.975}.
#' @param fft Logical. If \code{TRUE} (default), smoothing is performed via
#'   Fast Fourier Transform (FFT) convolution, which is substantially faster
#'   for large rasters and wide kernels. Some edge effects may occur at raster
#'   boundaries. Set to \code{FALSE} to use \code{terra::focal}, which avoids
#'   edge effects but is slower.
#' @param scale_center Logical. If \code{TRUE}, the smoothed raster values are
#'   centered and scaled using the mean and standard deviation from the model
#'   fitting data (extracted from \code{multiScaleR$scl_params}). Required
#'   when using the output with \code{terra::predict} and a fitted model that
#'   was trained on scaled covariates. Requires \code{multiScaleR} to be
#'   provided. Default: \code{FALSE}.
#' @param clamp Logical. If \code{TRUE}, scaled raster values are clamped to
#'   the observed covariate range from the model fitting data, preventing
#'   extrapolation beyond the training range. Only active when
#'   \code{scale_center = TRUE}. Default: \code{FALSE}.
#' @param pct_mx Numeric between -0.99 and 0.99. When \code{clamp = TRUE},
#'   expands (\code{> 0}) or contracts (\code{< 0}) the clamping range
#'   relative to the observed min/max by this proportion. Default: \code{0}
#'   (exact training range).
#' @param na.rm Logical. If \code{TRUE} (default), \code{NA} cells are excluded
#'   from the kernel-weighted mean calculation. If \code{FALSE}, any window
#'   containing a \code{NA} cell will produce a \code{NA} output.
#' @param verbose Logical. Print layer-level progress messages. Default:
#'   \code{TRUE}.
#' @param ... Reserved for deprecated arguments. Currently only \code{scale_opt}
#'   (deprecated alias for \code{multiScaleR}) is handled.
#'
#' @return A \code{SpatRaster} with one layer per covariate defined in
#' \code{scale_vars} (or per raster layer when \code{scale_vars} is not used).
#' Layer names match the covariate names from the model. When a fitted
#' \code{multiScaleR} object is provided and the model contains site-level
#' predictors (i.e., predictors not derived from raster layers), constant
#' ("dummy") raster layers filled with zeros are appended to make the output
#' compatible with \code{terra::predict}.
#'
#' @details
#' \strong{Typical usage} after running \code{\link{multiScale_optim}}:
#' \preformatted{
#' opt_hab <- kernel_scale.raster(raster_stack = hab, multiScaleR = opt)
#' plot(c(hab, opt_hab))
#'
#' ## Scale and center for prediction
#' opt_hab_sc <- kernel_scale.raster(raster_stack = hab,
#'                                   multiScaleR = opt,
#'                                   scale_center = TRUE)
#' preds <- terra::predict(opt_hab_sc, opt$opt_mod, type = "response")
#' }
#'
#' \strong{FFT vs. focal smoothing}
#'
#' The FFT convolution (\code{fft = TRUE}) is substantially faster for large
#' rasters or wide kernels. It produces minor edge artefacts near raster
#' boundaries, typically within one kernel-width of the edge. For analyses
#' focused on interior areas this is usually negligible. Use
#' \code{fft = FALSE} for exact focal smoothing when edge accuracy matters.
#'
#' \strong{Dummy rasters for site-level covariates}
#'
#' When a fitted \code{multiScaleR} object is supplied, the function inspects
#' the model frame to find any predictors that are not represented by raster
#' layers (e.g., survey effort, habitat type measured in the field). These
#' cannot be projected spatially and are therefore assigned constant zero
#' rasters, which correspond to the reference value for centered and scaled
#' covariates. These dummy layers are required for \code{terra::predict} to
#' work but do not represent real spatial variation. Replace them manually for
#' non-zero projection scenarios. Categorical predictors are skipped with a
#' warning.
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
                                scale_vars = NULL,
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

  if(!inherits(raster_stack, "SpatRaster")){
    stop('Raster layers must be provided as a `SpatRaster` object from `terra`')
  }
  validate_scalar_logical(fft, "fft")
  validate_scalar_logical(scale_center, "scale_center")
  validate_scalar_logical(clamp, "clamp")
  validate_scalar_logical(na.rm, "na.rm")
  validate_scalar_logical(verbose, "verbose")
  validate_scalar_numeric(pct_wt,
                          "pct_wt",
                          lower = 0,
                          upper = 1,
                          inclusive_lower = FALSE,
                          inclusive_upper = FALSE)
  validate_scalar_numeric(pct_mx,
                          "pct_mx",
                          lower = -0.99,
                          upper = 0.99)

  if(!is.null(multiScaleR)){
    validate_multiScaleR_input(multiScaleR)
  }

  kernel <- match.arg(kernel)

  if (is.null(scale_vars) && !is.null(multiScaleR)) {
    if (inherits(multiScaleR, "multiScaleR") &&
        !is.null(multiScaleR$kernel_inputs$scale_vars)) {
      scale_vars <- multiScaleR$kernel_inputs$scale_vars
    } else if (inherits(multiScaleR, "multiScaleR_data") &&
               !is.null(multiScaleR$scale_vars)) {
      scale_vars <- multiScaleR$scale_vars
    }
  }

  if (!is.null(scale_vars)) {
    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      kernel <- multiScaleR$kernel_inputs$kernel
    } else if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")) {
      kernel <- multiScaleR$kernel
    }

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      scale_vars <- .msr_scale_vars_in_model(scale_vars = scale_vars,
                                             multiScaleR = multiScaleR)
    }

    scale_vars <- .msr_validate_scale_vars(scale_vars = scale_vars,
                                           raster_stack = raster_stack,
                                           kernel = kernel)
    opt_vars <- .msr_optimized_scale_vars(scale_vars)
    n_optimized <- nrow(opt_vars)

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      sigma <- multiScaleR$scale_est[opt_vars$covariate, "Mean"]
      shape <- if (!is.null(multiScaleR$shape_est)) {
        multiScaleR$shape_est[opt_vars$covariate, "Mean"]
      } else {
        NULL
      }
    } else if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR_data")) {
      sigma <- multiScaleR$sigma * multiScaleR$unit_conv
      shape <- multiScaleR$shape
    } else if (n_optimized == 0 && is.null(sigma)) {
      sigma <- numeric(0)
    }

    if (n_optimized > 0 && is.null(sigma)) {
      stop("sigma values must be specified for optimized `scale_vars`.\n",
           call. = FALSE)
    }
    if (n_optimized > 0) {
      validate_numeric_vector(sigma,
                              "sigma",
                              length_ = n_optimized,
                              positive = TRUE)
    }
    if (kernel == "expow" && n_optimized > 0) {
      validate_numeric_vector(shape,
                              "shape",
                              length_ = n_optimized,
                              positive = TRUE)
    }

    if(isTRUE(scale_center) && is.null(multiScaleR)){
      warning("`scale_center = TRUE` was ignored because no `multiScaleR` object was provided.",
              call. = FALSE)
    }

    if(isTRUE(clamp) && !isTRUE(scale_center)){
      warning("`clamp = TRUE` was ignored because `scale_center` is FALSE.",
              call. = FALSE)
    }

    smooth_stack <- .msr_scale_vars_raster(
      raster_stack = raster_stack,
      scale_vars = scale_vars,
      sigma = sigma,
      shape = shape,
      kernel = kernel,
      pct_wt = pct_wt,
      fft = fft,
      na.rm = na.rm,
      verbose = verbose
    )

    if (isTRUE(scale_center) &&
        (inherits(multiScaleR, "multiScaleR_data") || inherits(multiScaleR, "multiScaleR"))) {
      smooth_stack <- scale_center_raster(r = smooth_stack,
                                          multiScaleR = multiScaleR,
                                          clamp = clamp,
                                          pct_mx = pct_mx)
    }

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      smooth_stack <- .add_site_covariate_rasters(
        smooth_stack = smooth_stack,
        multiScaleR = multiScaleR,
        raster_covs = names(smooth_stack)
      )
    }

    return(smooth_stack)
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

    if(!all(covs %in% names(raster_stack))){
      stop('optimized covariate is not present in the provided SpatRaster!')
    } else {
      raster_stack <- subset(raster_stack, covs)
    }
  } else {
    covs <- var <- names(raster_stack)
  }

  if(is.null(sigma)){
    stop("sigma values must be specified\n")
  }
  validate_numeric_vector(sigma,
                          "sigma",
                          positive = TRUE)

  if(length(sigma) != nlyr(raster_stack)){
    warning("Number of sigma values must equal the number of raster layers!!!  \n  All raster layers will be smoothed using the same sigma value")
    sigma <- rep(sigma[1], nlyr(raster_stack))
  }

  if(kernel == "expow"){
    validate_numeric_vector(shape,
                            "shape",
                            length_ = length(sigma),
                            positive = TRUE)
  }

  if(isTRUE(scale_center) && is.null(multiScaleR)){
    warning("`scale_center = TRUE` was ignored because no `multiScaleR` object was provided.",
            call. = FALSE)
  }

  if(isTRUE(clamp) && !isTRUE(scale_center)){
    warning("`clamp = TRUE` was ignored because `scale_center` is FALSE.",
            call. = FALSE)
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

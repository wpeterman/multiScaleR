#' @title Apply kernel smoothing to raster layers
#'
#' @description Applies a kernel smoothing function to one or more raster
#' layers, producing a spatially weighted mean, surface metric, or landscape metric at the
#' scale identified by \code{\link{multiScale_optim}}. Primarily used to
#' generate prediction rasters for spatial model projection.
#'
#' @param raster_stack A \code{SpatRaster} object containing the source raster
#'   layer(s) to be smoothed. Layer names must match the covariate names used
#'   during model fitting (or the \code{source} names in \code{scale_vars}).
#' @param sigma Numeric vector of kernel scale parameter(s) in the same units
#'   as the raster projection (e.g., meters). One value per optimized covariate
#'   in declaration order, including separate derived variables from one source.
#'   Ignored when \code{multiScaleR} is provided; sigma values are extracted
#'   from the fitted object automatically.
#' @param multiScaleR A fitted object of class \code{"multiScaleR"} (from
#'   \code{\link{multiScale_optim}}) or class \code{"multiScaleR_data"} (from
#'   \code{\link{kernel_prep}}). When provided, sigma, shape, kernel, and
#'   scale_vars are all extracted automatically. Default: \code{NULL}.
#' @param shape Numeric vector of shape parameters for the exponential power
#'   kernel (\code{kernel = "expow"}), one per optimized covariate. A
#'   hard-radius variable ignores its shape position. Ignored when
#'   \code{multiScaleR} is provided. Default: \code{NULL}.
#' @param kernel Character. Kernel function used for smoothing. One of
#'   \code{"gaussian"} (default), \code{"exp"}, \code{"fixed"}, or
#'   \code{"expow"}. Ignored when \code{multiScaleR} is provided.
#' @param scale_vars Optional variable specifications created with
#'   \code{\link{msr_vars}}. Use when projecting landscape metrics or
#'   explicitly defined covariates without passing a fitted \code{multiScaleR}
#'   object. When \code{multiScaleR} is provided, \code{scale_vars} is
#'   extracted automatically.
#' @param pct_wt Numeric between 0 and 1 (exclusive). Cumulative kernel weight
#'   along a line through the focal cell used to determine the smoothing window
#'   radius. A larger value (e.g., \code{0.99}) includes more of the kernel tail
#'   but increases computation time. This is the historical one-dimensional
#'   cutoff, not the proportion of total weight within the planar window.
#'   Default: \code{0.975}.
#' @param fft Logical. If \code{TRUE} (default), smoothing is performed via
#'   Fast Fourier Transform (FFT) convolution, which is substantially faster
#'   for large rasters and wide kernels. Set to \code{FALSE} to use direct
#'   focal calculations. Both engines support kernel means and the declared
#'   surface and landscape metrics, with the same available-support and
#'   missing-value rules. Direct focal calculations can be slower.
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
#' @param offset_values Named numeric vector of values for variables
#'   inside model offsets, such as \code{c(n_weeks = 4)} for
#'   \code{offset(log(n_weeks))}. These are raw input values, not transformed
#'   offsets. When omitted, a log offset uses 1 (zero on the log scale) and
#'   a direct offset uses 0. Set this argument for a meaningful projection
#'   scenario.
#' @param include_fixed Logical. For a fitted native model with verifiable
#'   retained numeric predictors and training row IDs,
#'   set \code{TRUE}
#'   with \code{scale_center = TRUE} to standardize modeled fixed-radius metrics
#'   using their verified original training preparation. The default
#'   \code{FALSE} uses only the supplied object's scaling parameters. The opt-in
#'   requires retained training buffers, matching fixed metric definitions and
#'   raster resolution, and unchanged native training predictor values.
#'   Retain the native model frame (for GLM and coxph, \code{model = TRUE}) or
#'   an equivalent recoverable training frame. Fixed inputs must be retained
#'   under their original names, with row IDs matching prepared sites. See
#'   \code{\link{prepare_unmarked_data}} for the verification contract.
#' @param match_buffer Logical. The default \code{FALSE} retains the historical
#'   \code{pct_wt} projection window. Set \code{TRUE} to use the retained
#'   preparation's \code{max_D} polygon-intersection support and finite-cell
#'   normalization, including partially intersecting boundary cells.
#'   Requires a fitted or prepared \code{multiScaleR} object, matching projected
#'   square-cell resolution and \code{na.rm = TRUE}. Both smoothing engines
#'   support this option.
#'   Metric radii larger than the prepared buffer raise an error. Central source
#'   NA cells remain masked. Binning, rasterization, and cell-center differences
#'   remain separate approximations.
#' @param site_values Optional named list of scalar nonspatial projection
#'   values in the native predictor units, for example
#'   \code{list(season = "spring", effort = 2, surveyed = TRUE)}. Requires a
#'   fitted object. Factor and character values must be fitted levels.
#'   Numeric predictors retain their historical zero default when omitted;
#'   categorical predictors need an explicit value. Unknown names and attempts
#'   to overwrite computed spatial layers raise an error. Use
#'   \code{offset_values} for variables occurring only inside offsets.
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
#' Offset variables are appended under their input names so transformed offset
#' expressions can be evaluated by prediction. With a fitted object, the
#' returned R object also retains a preparation descriptor for checked
#' \code{\link{predict.multiScaleR}} raster prediction. Ordinary raster file
#' formats need not retain this descriptor.
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
#'                                   scale_center = TRUE,
#'                                   match_buffer = TRUE)
#' preds <- predict(opt, opt_hab_sc, type = "response")
#' }
#'
#' \strong{FFT vs. focal smoothing}
#'
#' FFT and direct focal calculations evaluate the same declared neighborhood
#' and finite-cell normalization. Weighted means renormalize over available
#' mapped support, including raster edges. Central source NA cells remain
#' masked. With \code{na.rm = FALSE}, incomplete footprints, including missing
#' cells beyond the map extent, return NA. Small floating-point differences
#' between engines can occur. Constant weighted SQ support returns zero;
#' standardized higher moments remain undefined when spread is zero.
#'
#' \strong{Dummy rasters for site-level covariates}
#'
#' When a fitted \code{multiScaleR} object is supplied, the function inspects
#' the model frame to find any predictors that are not represented by raster
#' layers (e.g., survey effort, habitat type measured in the field). These
#' cannot be projected spatially and are therefore assigned constant zero
#' rasters, which correspond to the reference value for centered and scaled
#' covariates. These dummy layers are required for \code{terra::predict} to
#' work but do not represent real spatial variation. Set \code{site_values}
#' for other scenarios. Categorical predictors are skipped with a warning
#' unless a fitted level is explicitly supplied. Their raster category labels
#' preserve the fitted levels. Checked prediction also restores logical types.
#' Native unmarked models require explicit nonspatial prediction scenarios.
#' When every native predictor is represented by a spatial layer, no dummy
#' layers are needed. Supply remaining variables through native prediction data
#' or \code{terra::predict(..., const = scenario)}.
#' Offset terms use their underlying input variables as layer names. For
#' \code{offset(log(n_weeks))}, the default layer is \code{n_weeks = 1},
#' representing one week of effort. Set \code{offset_values = c(n_weeks = 4)}
#' to predict for four weeks. Direct offsets default to zero.
#' Compound offsets, such as \code{offset(log(effort * duration))}, require
#' explicit raw values for each underlying input. Literal constant offsets
#' need no layer. Ambiguous or conflicting default scenarios raise an error.
#'
#' \strong{Verified fixed-radius metrics in native projection}
#'
#' Set \code{include_fixed = TRUE} with \code{scale_center = TRUE} when the
#' fitted model contains declared fixed-radius metrics whose scaling parameters
#' are retained only in its original preparation. The original fixed transforms
#' and training normalization must reproduce the native model covariates
#' exactly. Changed definitions or training inputs raise an error. This opt-in
#' does not alter the caller's fitted object, metric definitions, or native
#' uncertainty interpretation. Point-buffer preparation and raster smoothing
#' retain their respective truncation and approximation rules.
#' Native projection requires a recoverable training frame and matching
#' original prepared row IDs for this verification. Transformations that do not
#' retain the fixed predictor under its original name are not reconstructed.
#' For conditional-logistic step selection, native \code{type = "risk"} with
#' \code{reference = "zero"} gives a relative selection score. Supply required
#' stratum variables explicitly through \code{const}. These scores are not
#' choice probabilities or occurrence probabilities; choice probabilities require
#' normalization within explicit sets of alternatives.
#'
#' \strong{Matching the prepared neighborhood}
#'
#' \code{match_buffer = TRUE} uses the same projected \code{sf::st_buffer()}
#' and polygon-intersection cell selection as point preparation. Partial boundary
#' cells receive their usual full-cell kernel weights, matching the fitting
#' convention rather than adding coverage-fraction weights. Finite cells within
#' mapped support are renormalized at raster edges; central source NA cells remain
#' masked. Constant weighted SQ support has zero population spread in this mode.
#' The historical \code{pct_wt} window remains the default when
#' \code{match_buffer = FALSE}; both window choices use the finite-support,
#' source-mask, and missing-value corrections described above.
#'
#' This option aligns the extraction footprint at raster cell centers. It does
#' not undo distance binning, reconstruct raster contents or original CRS metadata,
#' replace missing source coverage, or remove differences between an arbitrary
#' sampled point and a raster cell center. Use the same intended sources and map
#' units. A larger prepared buffer can increase projection memory and runtime.
#' Hard-neighborhood metric radii above \code{max_D} are refused because their
#' truncated cell-intersection footprint is not equivalent to a smaller radius.
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
                                offset_values = NULL,
                                include_fixed = FALSE,
                                match_buffer = FALSE,
                                site_values = NULL,
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
  validate_scalar_logical(include_fixed, "include_fixed")
  validate_scalar_logical(match_buffer, "match_buffer")
  site_values <- .msr_site_values(site_values)
  if (length(site_values) && !inherits(multiScaleR, "multiScaleR"))
    stop("site_values requires a fitted multiScaleR object.", call. = FALSE)
  if (!is.null(offset_values)) {
    if (!is.numeric(offset_values) || is.null(names(offset_values)) ||
        any(!nzchar(names(offset_values))) || anyDuplicated(names(offset_values)) ||
        any(!is.finite(offset_values))) {
      stop("`offset_values` must be a named finite numeric vector.", call. = FALSE)
    }
  }
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
  buffer_radius <- NULL
  if (match_buffer) {
    if (!(inherits(multiScaleR, "multiScaleR") ||
          inherits(multiScaleR, "multiScaleR_data")))
      stop("match_buffer requires a fitted or prepared multiScaleR object.", call. = FALSE)
    if (!isTRUE(na.rm))
      stop("match_buffer requires na.rm = TRUE.", call. = FALSE)
    prepared <- if (inherits(multiScaleR, "multiScaleR"))
      multiScaleR$kernel_inputs else multiScaleR
    buffer_radius <- multiScaleR$max_D
    validate_scalar_numeric(buffer_radius, "prepared max_D", positive = TRUE)
    validate_scalar_numeric(prepared$resolution, "prepared resolution", positive = TRUE)
    if (!nzchar(terra::crs(raster_stack)) || isTRUE(terra::is.lonlat(raster_stack)) ||
        !isTRUE(all.equal(terra::res(raster_stack), rep(prepared$resolution, 2L))))
      stop("match_buffer requires projected square cells at the prepared resolution.", call. = FALSE)
  }
  fixed_contract <- NULL
  if (include_fixed) {
    if (!isTRUE(scale_center) || !inherits(multiScaleR, "multiScaleR"))
      stop("include_fixed requires a supported fitted native model and scale_center = TRUE.", call. = FALSE)
    if (!isTRUE(all.equal(terra::res(raster_stack),
          rep(multiScaleR$kernel_inputs$resolution, 2L))))
      stop("Prediction raster resolution must match the fit.", call. = FALSE)
    fixed_contract <- .msr_fixed_training_contract(multiScaleR)
    fixed_specs <- multiScaleR$kernel_inputs$scale_vars
    fixed_specs <- fixed_specs[fixed_specs$covariate %in% fixed_contract$covariates, , drop = FALSE]
    if (any(fixed_specs$radius > multiScaleR$max_D))
      stop("A fixed radius exceeded the fitting buffer; reprepare with an adequate buffer or supply explicitly verified projection layers.", call. = FALSE)
    # This local list copy leaves the caller's fitted object unchanged.
    multiScaleR$scl_params$mean[fixed_contract$covariates] <- fixed_contract$center
    multiScaleR$scl_params$sd[fixed_contract$covariates] <- fixed_contract$spread
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
    if (include_fixed && length(fixed_contract$covariates)) {
      fixed <- fixed_contract$covariates
      expected <- multiScaleR$kernel_inputs$scale_vars
      expected <- expected[match(fixed, expected$covariate), , drop = FALSE]
      actual <- scale_vars[match(fixed, scale_vars$covariate), , drop = FALSE]
      if (anyNA(actual$covariate) ||
          !isTRUE(all.equal(as.data.frame(expected), as.data.frame(actual),
                            check.attributes = FALSE)))
        stop("Fixed-radius source and metric definitions must match the fit.", call. = FALSE)
    }
    opt_vars <- .msr_optimized_scale_vars(scale_vars)
    n_optimized <- nrow(opt_vars)

    if (!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")) {
      sigma <- multiScaleR$scale_est[opt_vars$covariate, "Mean"]
      shape <- if (!is.null(multiScaleR$shape_est)) {
        .msr_shape_values(multiScaleR)[opt_vars$covariate]
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
      verbose = verbose,
      buffer_radius = buffer_radius
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
        raster_covs = names(smooth_stack),
        offset_values = offset_values, site_values = site_values
      )
    }

    return(.msr_stamp_projection(smooth_stack, multiScaleR, scale_center,
      match_buffer, scale_vars, include_fixed, site_values,
      list(fft = fft, clamp = clamp, pct_mx = pct_mx, na.rm = na.rm,
           pct_wt = pct_wt, offset_values = offset_values)))
  }

  if(!is.null(multiScaleR) && inherits(multiScaleR, "multiScaleR")){
    covs <- rownames(multiScaleR$scale_est)
    sigma <- multiScaleR$scale_est[,1]
    shape <- .msr_shape_values(multiScaleR)
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

    if (match_buffer) {
      smooth_list[[i]] <- .msr_kernel_raster_one(raster_stack[[lyr]], sigma[i],
        shape[i], kernel, pct_wt, fft, na.rm, buffer_radius = buffer_radius)
      next
    }

    mx <- kernel_dist(kernel = kernel,
                      dimension = "1d",
                      sigma = sigma[i],
                      beta = shape[i],
                      prob = pct_wt)

    r_res <- res(raster_stack)[1]
    focal_d <- ceiling(mx / r_res) * 2

    if((focal_d %% 2) == 0) {
      focal_d <- focal_d + 1
    }

    mat <- matrix(0, focal_d, focal_d)
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
                                na.policy = "omit",
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
      raster_covs = var,
      offset_values = offset_values, site_values = site_values
    )
  }

  return(.msr_stamp_projection(smooth_stack, multiScaleR, scale_center,
    match_buffer, scale_vars, include_fixed, site_values,
    list(fft = fft, clamp = clamp, pct_mx = pct_mx, na.rm = na.rm,
         pct_wt = pct_wt, offset_values = offset_values)))
}

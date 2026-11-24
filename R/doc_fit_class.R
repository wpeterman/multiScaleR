#' Fitted multiScaleR objects
#'
#' @description
#' [multiScale_optim()] returns an object of class `"multiScaleR"`, a list
#' that holds the estimated scales of effect, your model refit at those
#' scales, and the record of the search. This page describes every element
#' and the functions that read the object.
#'
#' @details
#' # Elements
#'
#' Every fit has these elements:
#' * `scale_est`: data frame with one row per optimized covariate (row names
#'   are covariate names) and two columns. `Mean` is the estimated sigma in
#'   map units, and `SE` is its Hessian-based standard error. `SE` is `NA` for
#'   hard-radius landscape and unweighted surface metrics and for every radius
#'   estimated under `kernel = "fixed"`, because their likelihood changes in
#'   steps.
#' * `shape_est`: data frame with the same layout for the shape parameter
#'   when `kernel = "expow"`, otherwise `NULL`. Hard-radius covariates have
#'   `NA` shape and shape SE because a shape does not apply to them; their
#'   unused shape positions are left out of the likelihood parameter count.
#' * `optim_results`: the raw list returned by [stats::optim()] or
#'   [optimParallel::optimParallel()], including `par`, `value` (the negative
#'   log-likelihood), `hessian`, and convergence codes. It also holds
#'   `par_unscale` (sigma values in map units) and `hessian_unscale`. For a
#'   radius grid, this list describes the selected combination: its Hessian
#'   covers only the smooth parameters at those fixed radii, with `NA` for
#'   hard-radius rows and columns and any inactive shape. See
#'   `diagnostics(fit)$radius_profile` for every candidate.
#' * `opt_mod`: your model refit at the estimated scales. Use it for
#'   coefficients, inference, prediction, and [plot_marginal_effects()]. Its
#'   coefficients are on the model's link scale (for example log for a Poisson
#'   model), and by default each optimized covariate is standardized, so a
#'   coefficient describes the change in the linear predictor for one standard
#'   deviation of the covariate.
#' * `fitted_mod_original`: the `fitted_mod` you supplied, kept for
#'   reference.
#' * `min_D`, `max_D`: the lower and upper bounds for sigma during the
#'   search, in map units.
#' * `kernel_inputs`: the `kernel_inputs` list (see [multiScaleR_data-class]),
#'   without `min_D` and `max_D`, which are stored separately.
#' * `scl_params`: named list with `mean` and `sd` vectors, the centering and
#'   scaling values of each optimized covariate at the estimated scales.
#'   [kernel_scale.raster()] uses them when `scale_center = TRUE` so that maps
#'   match the model's covariates. With `covariate_transform = "raw"`, the
#'   centers are zero and the scales one; prediction never relearns these
#'   values from new data.
#' * `join_by`: the `join_by` data frame, or `NULL`.
#' * `opt_context`: internal record of the model-class-specific refit logic,
#'   kept for [profile_sigma()]. With `covariate_transform = "raw"` it also
#'   records that setting. `NULL` for the analytic backends.
#' * `profile_scale_est`: data frame of profile-likelihood intervals, or
#'   `NULL` until [profile_sigma()] has been run.
#' * `diagnostics`: list of diagnostic records; read it with [diagnostics()].
#' * `next_run`: list of suggested settings for a follow-up fit: `max_D` for
#'   the next [kernel_prep()] call (raised when an estimate approaches the
#'   current `max_D`), the estimated `start_sigma` values in map units,
#'   `start_shape` for `kernel = "expow"`, and `start_par`, the starts on the
#'   internal scale that `multiScale_optim()` expects, scaled to that `max_D`.
#'   `n_cores` is `NULL`; ask [estimate_multiscale_ram()] for a parallel
#'   worker count. `NULL` for some analytic backends.
#' * `warn_message`: integer vector of warning codes. It starts with `0` and
#'   adds `1` when an estimate approaches `max_D`, `2` when a sigma standard
#'   error is at least half its estimate, and `3` when a shape standard error
#'   is at least half its estimate. See "Read the warnings" in
#'   [multiScale_optim()].
#' * `call`: the matched call.
#'
#' # Elements added by the analytic backends
#'
#' * `profile_fit`, `uncertainty`, `dispersion_est`: present for
#'   `"profile_glm"` and `"joint_nb2"`. They keep the transformed parameter
#'   estimates, the starts and their convergence scores, the profile and joint
#'   information matrices, identifiability diagnostics, and NB2 theta.
#'   `uncertainty$available` says whether joint uncertainty passed its checks,
#'   and `uncertainty$reasons` says why not.
#'   `uncertainty$joint_available` is also `TRUE` when joint uncertainty holds
#'   bounded kernel parameters fixed, and `uncertainty$joint_conditional_on`
#'   names them. `opt_mod` and its coefficient
#'   summary condition on the estimated kernel parameters. Use
#'   [vcov.multiScaleR()] and [predict.multiScaleR()] for joint uncertainty.
#' * `profile_fit` and `uncertainty` for `"conditional_logit"`, which keep
#'   the starts, every start outcome, and the selected start. Scale
#'   uncertainty is not available.
#' * `mixed_fit` for `"joint_mixed"`: the kernel engine, log-scale estimates
#'   and bounds in internal distance units, the independent starts, every start
#'   outcome and warning, the selected start, and the variables that
#'   prediction and profile refits need. `scale_est$SE` and the log-scale
#'   Wald limits in `profile_scale_est` come from the selected start's joint Hessian, one
#'   scale at a time. `uncertainty$available` is `FALSE` because coefficient
#'   and prediction uncertainty still condition on the scales. `next_run` and
#'   the `max_D` warning follow the analytic GLM backends, and `opt_context` is
#'   `NULL`.
#'
#' [multiScaleR-fitting-methods] explains what each backend records and why
#' uncertainty may be withheld.
#'
#' # Where each quantity lives
#'
#' | Quantity | Where to find it | Scale |
#' |:---------|:-----------------|:------|
#' | Scale of effect (sigma) | `fit$scale_est$Mean`, `summary(fit)` | Map units |
#' | Kernel shape | `fit$shape_est$Mean` | Shape units of the exponential power kernel |
#' | Effective distance | `summary(fit)`, [kernel_dist()] | Map units |
#' | Regression coefficients | `coef(fit$opt_mod)`, `summary(fit)` | Link scale, per standard deviation of each covariate by default |
#' | NB2 dispersion | `fit$dispersion_est$theta` (analytic backends) | Natural scale |
#' | Profile-likelihood intervals | `fit$profile_scale_est`, `summary(fit, profile = TRUE)` | Map units |
#'
#' The *effective distance* converts sigma into the distance that contains a
#' set proportion (by default 90 percent) of the kernel's total weight, which
#' is usually the easiest number to report. See [multiScaleR-kernels].
#'
#' # Functions that read the object
#'
#' | Function | What it gives you |
#' |:---------|:------------------|
#' | `print()`, [summary.multiScaleR()] | Scales, effective distances, intervals, and coefficients |
#' | [plot.multiScaleR()] | The fitted kernel weight curve for each covariate |
#' | [diagnostics()] | Search, sample-size, and uncertainty diagnostics |
#' | [profile_sigma()] | Profile-likelihood intervals for sigma |
#' | [predict.multiScaleR()], [vcov.multiScaleR()] | Predictions and covariance, with joint uncertainty where available |
#' | [kernel_scale.raster()] | Maps of the covariates at the estimated scales |
#' | [plot_marginal_effects()] | Effect curves for each covariate |
#' | [aic_tab()], [bic_tab()] | Model comparison by AIC or BIC |
#' | [cross_validate()] | Out-of-sample prediction scores |
#' | [scale_inference()] | Likelihood-range and bootstrap assessment of single-kernel NB2 scales |
#'
#' @seealso [multiScale_optim()], [multiScaleR-fitting-methods],
#'   [multiScaleR-prediction], and
#'   `vignette("v01-quickstart", package = "multiScaleR")`.
#' @examples
#' data(pts)
#' data(count_data)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' mod <- glm(y ~ hab, family = poisson, data = count_data)
#' fit <- multiScale_optim(mod, inputs, verbose = FALSE)
#'
#' names(fit)
#' fit$scale_est            # sigma (Mean) and its SE, in map units
#' coef(fit$opt_mod)        # coefficients on the log (link) scale
#' fit$warn_message         # 0 plus a code for each warning raised
#' fit$next_run$max_D       # max_D suggested for a follow-up fit
#' @name multiScaleR-class
NULL

#' Prepared inputs for scale estimation
#'
#' @description
#' [kernel_prep()] and [kernel_prep_by_group()] return an object of class
#' `"multiScaleR_data"`, a list that holds everything the fitting functions
#' need about the survey points and the raster cells around them. Pass it as
#' `kernel_inputs` to [multiScale_optim()], [multiScale_bayes()],
#' [cross_validate()], and related functions. This page describes every
#' element.
#'
#' @details
#' # Elements
#'
#' * `kernel_dat`: data frame of the initial kernel-weighted covariate values,
#'   centered and scaled to unit variance, with one row per point and one
#'   column per covariate. Row names match the row names of `pts` or are
#'   sequential integers. Fit your starting model to these columns.
#' * `d_list`: named list, one element per point, of distances from the point
#'   to every raster cell within the buffer. `NULL` in lean mode.
#' * `raw_cov`: named list, one element per point, of matrices holding the raw
#'   raster cell values within the buffer, aligned with `d_list`. Each is
#'   stored as a dense or sparse matrix, whichever is smaller. `NULL` in lean
#'   mode.
#' * `kernel`: the kernel name (`"gaussian"`, `"exp"`, `"expow"`, or
#'   `"fixed"`).
#' * `sigma`: initial sigma values on the internal scale (map units divided by
#'   `unit_conv`).
#' * `shape`: initial shape values for `"expow"`, otherwise `NULL`.
#' * `min_D`: approximate minimum distance between cells, used as the lower
#'   bound for sigma during the search, in map units.
#' * `max_D`: the `max_D` you supplied, the upper bound for sigma, in map
#'   units.
#' * `n_covs`: number of covariates whose scale is estimated.
#' * `unit_conv`: internal distance scaling factor, equal to `max_D`.
#' * `scale_vars`: data frame of class `"multiScaleR_vars"` describing every
#'   covariate specification (see [msr_vars()]).
#' * `resolution`: raster cell resolution in map units.
#' * `cell_size`: cell width and height in map units. Cross-validation uses
#'   them to keep intersecting cells whose centers lie beyond `max_D`.
#' * `n_cols`: number of raster columns, used by adjacency-based landscape
#'   metrics.
#' * `binned`: when `bin = TRUE`, a list of distance-binned summaries
#'   (representative bin distances and value and count matrices for each
#'   covariate) that speed up the search. `NULL` when `bin = FALSE` or when no
#'   covariate is a kernel-weighted mean.
#' * `surface_binned`: when `bin = TRUE`, linearly binned power sums for
#'   kernel-weighted surface metrics (`surface_var(weighted = TRUE)` with
#'   `"sq"`, `"ssk"`, or `"sku"`) that speed up the refit search. Each cell's
#'   centered value powers are split between the two nearest distance nodes.
#'   `NULL` when `bin = FALSE`, when `kernel = "fixed"`, or when no covariate is
#'   such a metric. The cell-level data stay stored.
#' * `cell_data_stored`: `FALSE` only in lean mode
#'   (`store_cell_data = FALSE`), when `d_list` and `raw_cov` were dropped.
#' * `scl_params`: named list with `mean` and `sd`, the centering and scaling
#'   values applied to `kernel_dat`.
#' * `raster_group`: from [kernel_prep_by_group()] only, the raster group of
#'   each point.
#'
#' Lean mode keeps only the binned summaries. It works only when every
#' covariate is a kernel-weighted mean, because landscape and surface metrics
#' need the cell-level data. Printing the object reports the number of points
#' and covariates, the kernel, the binning, and the distance bounds.
#'
#' @seealso [kernel_prep()], [kernel_prep_by_group()], [msr_vars()],
#'   [multiScaleR-class], and
#'   `vignette("v02-user-guide", package = "multiScaleR")`.
#' @examples
#' data(pts)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' inputs                    # short summary of the prepared inputs
#' head(inputs$kernel_dat)   # standardized covariate values at the start
#' inputs$max_D              # upper bound for sigma, in map units
#' @name multiScaleR_data-class
NULL

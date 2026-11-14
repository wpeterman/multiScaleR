#' multiScaleR: estimate scales of effect for spatial covariates
#'
#' @description
#' The *scale of effect* is the distance over which a landscape variable
#' influences an ecological response, such as how far from a pond the
#' surrounding forest still affects salamander counts. multiScaleR estimates it
#' directly. It weights the raster cells around each survey point with a
#' distance *kernel*, and it chooses the kernel scale (sigma) that maximizes
#' the likelihood of your regression model. You fit your model once at
#' starting values, and multiScaleR does the search, reports each scale with
#' its uncertainty, and helps you compare, check, and map the result.
#'
#' @details
#' \if{html}{
#'   \figure{logo.png}{options: width=150 alt='Package logo'}
#' }
#'
#' # Start here
#'
#' Work through `vignette("v01-quickstart", package = "multiScaleR")`, then
#' `vignette("v02-user-guide", package = "multiScaleR")`. Run
#' `browseVignettes("multiScaleR")` to see all 15 vignettes, which cover
#' mapping, landscape and surface metrics, year-specific rasters, fitting
#' methods, cross-validation, intervals, step selection, spatial fields,
#' distance-removal models, and Bayesian models. [multiScaleR-glossary]
#' defines the terms they use.
#'
#' # The basic workflow
#'
#' | Step | Functions |
#' |:-----|:----------|
#' | 1. Describe the covariates | [kernel_var()], [landscape_var()], [surface_var()], [msr_vars()] |
#' | 2. Prepare the points and rasters | [kernel_prep()], [kernel_prep_by_group()], [estimate_multiscale_ram()] |
#' | 3. Fit your model at starting values | Your usual function, such as `glm()` or `unmarked::pcount()` |
#' | 4. Estimate the scales | [multiScale_optim()] or [multiScale_bayes()] |
#' | 5. Read and check the fit | `summary()`, `plot()`, [diagnostics()], [profile_sigma()], [kernel_dist()], [plot_kernel()] |
#' | 6. Compare models | [aic_tab()], [bic_tab()], [cross_validate()], [cv_tab()], [select_cv()] |
#' | 7. Predict and map | `predict()`, [kernel_scale.raster()], [plot_marginal_effects()], [prepare_prediction()] |
#'
#' # Specialized models and tools
#'
#' * **Choosing a fitting method**: [multiScaleR-fitting-methods] describes
#'   each backend of [multiScale_optim()].
#' * **Cross-validation for other model types**: [cross_validate_glm()],
#'   [cross_validate_clogit()], [cross_validate_field()],
#'   [cross_validate_mixed()], [cross_validate_bayes()], with fold builders
#'   [spatial_cv_folds()] and [grouped_cv_folds()].
#' * **Occupancy and abundance models from `unmarked`**:
#'   [prepare_unmarked_data()], [validate_unmarked_data()].
#' * **Distance-removal models**: [prepare_gdr_histories()],
#'   [prepare_gdr_intervals()], [expand_gdr_intervals()], [gdr_start_values()],
#'   [make_gdr_refitter()], [fit_gdr_response_field()],
#'   [response_field_geometry()].
#' * **Bayesian models**, in which a scale may vary among years or with a
#'   covariate: [multiScale_bayes()], [check_bayes_setup()],
#'   [simulate_bayes_prior()], [posterior_scales()], [predict_scales()],
#'   [diagnose_bayes_scales()], [diagnose_bayes_probabilities()],
#'   [diagnose_bayes_laplace()].
#' * **Assessing NB2 scales over the full likelihood range**:
#'   [scale_inference()].
#' * **Simulating data with known scales of effect**: [sim_rast()],
#'   [sim_dat()], [sim_dat_unmarked()], [sim_dat_gdr()].
#'
#' # Reference topics
#'
#' | Topic | What it covers |
#' |:------|:---------------|
#' | [multiScaleR-kernels] | Kernel families, sigma and shape, hard radii, effective distance |
#' | [multiScaleR-fitting-methods] | Backends, their requirements, and the uncertainty each reports |
#' | [multiScaleR-prediction] | Prediction quantities and the uncertainty each route carries |
#' | [multiScaleR-glossary] | Definitions of the terms used in the documentation |
#' | [multiScaleR-class] | The fitted object from [multiScale_optim()] |
#' | [multiScaleR_data-class] | The prepared inputs from [kernel_prep()] |
#' | [multiScaleR_bayes-class] | The fitted object from [multiScale_bayes()] |
#' | [multiScaleR_cv-class] | Cross-validation results |
#'
#' # Experimental tools
#'
#' These tools are experimental. Their arguments and results may change in a
#' later release without a deprecation period: [cross_validate_bayes()],
#' [cross_validate_mixed()], the count, binary, and positive families of
#' [cross_validate_field()], [diagnose_bayes_laplace()], and the
#' `radius_grid` argument of [multiScale_optim()].
#'
#' # Citation
#'
#' Run `citation("multiScaleR")` for the paper describing the method.
#'
#' @examples
#' # Estimate the scale of effect of habitat on simulated counts
#' data(pts)
#' data(count_data)
#' hab <- terra::rast(system.file("extdata", "hab.tif",
#'                                package = "multiScaleR"))
#' inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
#'                       kernel = "gaussian", verbose = FALSE)
#' start_model <- glm(y ~ hab, family = poisson, data = count_data)
#' fit <- multiScale_optim(start_model, inputs, verbose = FALSE)
#' summary(fit)   # sigma near the simulated value of 75 map units
"_PACKAGE"

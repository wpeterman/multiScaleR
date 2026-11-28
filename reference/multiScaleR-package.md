# multiScaleR: estimate scales of effect for spatial covariates

The *scale of effect* is the distance over which a landscape variable
influences an ecological response, such as how far from a pond the
surrounding forest still affects salamander counts. multiScaleR
estimates it directly. It weights the raster cells around each survey
point with a distance *kernel*, and it chooses the kernel scale (sigma)
that maximizes the likelihood of your regression model. You fit your
model once at starting values, and multiScaleR does the search, reports
each scale with its uncertainty, and helps you compare, check, and map
the result.

## Details

![Package logo](figures/logo.png)

## Start here

Work through
[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md),
then
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md).
Run `browseVignettes("multiScaleR")` to see all 15 vignettes, which
cover mapping, landscape and surface metrics, year-specific rasters,
fitting methods, cross-validation, intervals, step selection, spatial
fields, distance-removal models, and Bayesian models.
[multiScaleR-glossary](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md)
defines the terms they use.

## The basic workflow

|  |  |
|----|----|
| Step | Functions |
| 1\. Describe the covariates | [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md), [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md), [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md), [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md) |
| 2\. Prepare the points and rasters | [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md), [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md), [`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md) |
| 3\. Fit your model at starting values | Your usual function, such as [`glm()`](https://rdrr.io/r/stats/glm.html) or [`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html) |
| 4\. Estimate the scales | [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md) or [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) |
| 5\. Read and check the fit | [`summary()`](https://rspatial.github.io/terra/reference/summary.html), [`plot()`](https://rspatial.github.io/terra/reference/plot.html), [`diagnostics()`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md), [`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md), [`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md), [`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md) |
| 6\. Compare models | [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md), [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md), [`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md), [`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md), [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md) |
| 7\. Predict and map | [`predict()`](https://rspatial.github.io/terra/reference/predict.html), [`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md), [`plot_marginal_effects()`](https://wpeterman.github.io/multiScaleR/reference/plot_marginal_effects.md), [`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md) |

## Specialized models and tools

- **Choosing a fitting method**:
  [multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md)
  describes each backend of
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- **Cross-validation for other model types**:
  [`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
  [`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md),
  [`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md),
  [`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md),
  [`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md),
  with fold builders
  [`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md)
  and
  [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md).

- **Occupancy and abundance models from `unmarked`**:
  [`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md),
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md).

- **Distance-removal models**:
  [`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md),
  [`prepare_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_intervals.md),
  [`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md),
  [`gdr_start_values()`](https://wpeterman.github.io/multiScaleR/reference/gdr_start_values.md),
  [`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md),
  [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
  [`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md).

- **Bayesian models**, in which a scale may vary among years or with a
  covariate:
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
  [`check_bayes_setup()`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md),
  [`simulate_bayes_prior()`](https://wpeterman.github.io/multiScaleR/reference/simulate_bayes_prior.md),
  [`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
  [`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
  [`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md),
  [`diagnose_bayes_probabilities()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_probabilities.md),
  [`diagnose_bayes_laplace()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_laplace.md).

- **Assessing NB2 scales over the full likelihood range**:
  [`scale_inference()`](https://wpeterman.github.io/multiScaleR/reference/scale_inference.md).

- **Simulating data with known scales of effect**:
  [`sim_rast()`](https://wpeterman.github.io/multiScaleR/reference/sim_rast.md),
  [`sim_dat()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md),
  [`sim_dat_unmarked()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_unmarked.md),
  [`sim_dat_gdr()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md).

## Reference topics

|  |  |
|----|----|
| Topic | What it covers |
| [multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md) | Kernel families, sigma and shape, hard radii, effective distance |
| [multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md) | Backends, their requirements, and the uncertainty each reports |
| [multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md) | Prediction quantities and the uncertainty each route carries |
| [multiScaleR-glossary](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md) | Definitions of the terms used in the documentation |
| [multiScaleR](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md) | The fitted object from [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md) |
| [multiScaleR_data](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_data-class.md) | The prepared inputs from [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md) |
| [multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md) | The fitted object from [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) |
| [multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md) | Cross-validation results |

## Experimental tools

These tools are experimental. Their arguments and results may change in
a later release without a deprecation period:
[`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md),
[`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md),
the count, binary, and positive families of
[`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md),
[`diagnose_bayes_laplace()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_laplace.md),
and the `radius_grid` argument of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

## Citation

Run `citation("multiScaleR")` for the paper describing the method.

## See also

Useful links:

- <https://wpeterman.github.io/multiScaleR/>

- <https://github.com/wpeterman/multiScaleR>

- Report bugs at <https://github.com/wpeterman/multiScaleR/issues>

## Author

**Maintainer**: Bill Peterman <Peterman.73@osu.edu>
([ORCID](https://orcid.org/0000-0001-5229-9268))

Authors:

- Bill Peterman <Peterman.73@osu.edu>
  ([ORCID](https://orcid.org/0000-0001-5229-9268))

## Examples

``` r
# Estimate the scale of effect of habitat on simulated counts
data(pts)
data(count_data)
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian", verbose = FALSE)
start_model <- glm(y ~ hab, family = poisson, data = count_data)
fit <- multiScale_optim(start_model, inputs, verbose = FALSE)
summary(fit)   # sigma near the simulated value of 75 map units
#> 
#> Call:
#> multiScale_optim(fitted_mod = start_model, kernel_inputs = inputs, 
#>     verbose = FALSE)
#> 
#> 
#> Kernel used:
#> gaussian
#> 
#> ***** Optimized Scale of Effect *****
#> 
#>         Mean       SE     2.5%    97.5%
#> hab 76.16637 9.567451 57.09847 95.23428
#> 
#> 
#>   ==================================== 
#> 
#> ***** Optimized Scale of Effect -- Distance *****
#> Scale distances (map units) at 90%, two definitions.
#> 
#> 2D, radius of the circle for the ideal untruncated kernel (selected):
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37
#> 
#> 1D, along a line through the point (comparison):
#>       Mean  2.5%  97.5%
#> hab 125.28 93.92 156.65
#> 
#> 
#>   ==================================== 
#> 
#>  *****     Fitted Model Summary     *****
#> 
#> 
#> Call:
#> stats::glm(formula = y ~ hab, family = poisson, data = data)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   0.4807     0.1021   4.707 2.51e-06 ***
#> hab           0.8691     0.0855  10.166  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 186.585  on 74  degrees of freedom
#> Residual deviance:  70.185  on 73  degrees of freedom
#> AIC: 232.29
#> 
#> Number of Fisher Scoring iterations: 5
#> 
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent. 
#>  Consider increasing max_D in `kernel_prep` to >= 250.56 to ensure accurate estimation of scale.
#> 
#> 
#> 
#> ***** Follow-up recommendation *****
#> Recommended `max_D`: 250.56
#> Recommended starting `sigma` values:
#>      hab 
#> 76.16637 
#> 
#> Action:
#> Refit from the optimized parameters and address the flagged diagnostics. 
```

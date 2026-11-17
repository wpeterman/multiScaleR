# Prepare spatial covariates for native unmarked prediction

Calculate the optimized spatial predictors at new sites and insert them
into an explicitly keyed data frame for native `unmarked` prediction.
The training means and standard deviations are reused, so the new sites
are on the fitted model's covariate scale. Pass the result to
[`unmarked::predict()`](https://rdrr.io/r/stats/predict.html), or give
the same inputs to
[`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
for a reusable snapshot.

## Usage

``` r
prepare_unmarked_data(
  object,
  newdata,
  kernel_inputs,
  join_by,
  include_fixed = FALSE
)
```

## Arguments

- object:

  A fitted `multiScaleR` object whose `opt_mod` is an `unmarkedFit`,
  with its fitted scales, training standardization, and retained
  `kernel_inputs`.

- newdata:

  A nonempty data frame of prediction rows, in the layout that native
  `unmarked` prediction needs for your chosen component. Include the
  site key columns and every non-spatial variable, including factors and
  offset inputs. Existing optimized spatial columns are replaced. Supply
  fixed-radius predictors on their fitted scale unless
  `include_fixed = TRUE`. Other columns and the row order are left
  unchanged.

- kernel_inputs:

  A `multiScaleR_data` object from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  for the new sites. Use the fitted kernel, source and metric
  definitions, buffer extent (`max_D`), raster resolution, and the same
  binned or exact evaluation as the fit. Landscape and surface metrics
  need cell-level inputs, while lean binned inputs are enough for kernel
  means.

- join_by:

  A data frame with one unique key row per prepared site, in
  `kernel_inputs` point order. Its column names must occur in `newdata`,
  and every `newdata` row must match a key. Repeated keys in `newdata`
  repeat that site's spatial values, which lets you write site-by-year
  or site-by-survey scenarios without changing their order. Required.

- include_fixed:

  Logical. The default `FALSE` leaves fixed-radius predictors as native
  inputs that you supply. Set `TRUE` to calculate modeled fixed-radius
  metrics at the new sites too, after verifying that their retained
  preparation values and standardization match the fitted static site
  covariates. Missing or inconsistent training information raises an
  error.

## Value

A data frame with the same rows, row names, and non-spatial columns as
`newdata`. The calculated spatial columns hold the fitted kernel
transforms in standard-deviation units of the training sample, the scale
the native model was fitted on. They are model inputs, not predictions.
Do not standardize them again.

## What the function does and does not do

This function prepares data only. Native prediction determines the
component's meaning, factor handling, offsets, missing-data behavior,
and confidence limits, and neither this function nor
[`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md)
adds uncertainty in the selected scales.
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md)
lists what each `unmarked` component predicts.

Spatial predictors are static site covariates, as the `unmarked` fitting
adapter requires. Supply the native row layout for detection and annual
components yourself. The function does not infer observation times,
expand sites into years, modify an `unmarkedFrame`, or borrow
non-spatial values from the fitting data.

## Match the spatial inputs to the fit

Distances inside `kernel_inputs` are converted with its own `unit_conv`,
and the reported scales remain in map units. Metric definitions and the
choice of binned or exact kernel evaluation must match the fit.
Preparing bins at new sites can shift their representative distances
slightly, as in other binned prediction routes. Source names are
checked, but raster contents are not fingerprinted, so prepare the sites
from the intended source rasters and coordinate system. When the
internal fitted parameters are retained, they are reused directly, which
avoids rounding from a conversion to reported units and back. Stored
parameters that disagree with the reported scales or shapes raise an
error.

The calculated columns use the training sample's mean and standard
deviation. Prediction values are never clamped or standardized within
the new sample. Missing spatial values stay missing for the native
prediction method, while infinite values or invalid training standard
deviations raise errors.

## Recompute fixed-radius predictors

Fixed-radius predictors are recalculated only with
`include_fixed = TRUE`. Their normalization comes from the original
preparation object, and their fitted site values must match that object
after explicit training-site alignment. This rejects manually rescaled
or stale fixed predictors, changed radii, ambiguous training keys, and
missing training buffers rather than guessing their fitted scale. Keep
the original inputs if you need this verification.

## Distance-removal models

For
[`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
fits, only active component formulas identify spatial predictors.
Availability is inactive with one primary period, and the distance
formula is inactive for a uniform key. Optimized static site names must
not also occur in the native yearly or observation covariate tables, so
use distinct names for distinct inputs. Dot formulas use the original
component's available training-column names. GDR models with active
formula offsets or TMB zero inflation are rejected, as
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
explains.

## See also

[`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md),
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md),
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md),
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md),
and
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md)
for fitting and mapping `unmarked` models.

## Examples

``` r
## A stored N-mixture fit used in vignette("v02-user-guide"). Abundance
## depends on the bin1 and cont2 layers of a simulated landscape.
nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
                            package = "multiScaleR"))
layers <- subset(sim_rast(user_seed = 999, dim = 250), c("bin1", "cont2"))

## Three new sites, prepared with the fitted kernel, max_D, covariate
## definitions, and binning
new_pts <- terra::vect(cbind(x = c(800, 1250, 1700),
                             y = c(900, 1250, 1600)))
new_inputs <- kernel_prep(new_pts, layers, max_D = nmix$max_D,
                          kernel = nmix$kernel_inputs$kernel,
                          scale_vars = nmix$kernel_inputs$scale_vars,
                          bin = !is.null(nmix$kernel_inputs$binned),
                          verbose = FALSE)

## Keys follow the prepared point order. Site B appears twice, for
## example in two survey years, and reuses the same spatial values.
keys <- data.frame(site = c("A", "B", "C"))
scenarios <- data.frame(site = c("A", "B", "B", "C"))
nd <- prepare_unmarked_data(nmix, scenarios, new_inputs, join_by = keys)
nd  # bin1 and cont2 in training standard-deviation units
#>   site       bin1     cont2
#> 1    A  0.8463096 0.9121894
#> 2    B  1.6378884 0.6597968
#> 3    B  1.6378884 0.6597968
#> 4    C -1.1986195 0.6290868

## Expected abundance at each row, with native 95% confidence limits
unmarked::predict(nmix$opt_mod, type = "state", newdata = nd)
#>   Predicted        SE     lower     upper
#> 1 2.2541203 0.3245970 1.6998190 2.9891759
#> 2 4.6758934 0.6758364 3.5223707 6.2071771
#> 3 4.6758934 0.6758364 3.5223707 6.2071771
#> 4 0.7324141 0.1144437 0.5392024 0.9948591

## The same rows and component as a reusable snapshot
prepared <- prepare_prediction(nmix, scenarios, new_inputs,
                               component = "state", join_by = keys)
predict(prepared)
#>   Predicted        SE     lower     upper
#> 1 2.2541203 0.3245970 1.6998190 2.9891759
#> 2 4.6758934 0.6758364 3.5223707 6.2071771
#> 3 4.6758934 0.6758364 3.5223707 6.2071771
#> 4 0.7324141 0.1144437 0.5392024 0.9948591
```

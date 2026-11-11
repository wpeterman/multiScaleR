# Plot Marginal Effects from a Fitted Model

Generates marginal effect plots with 95% confidence intervals for each
covariate in the optimized model stored within a `multiScaleR` object.
Each panel sweeps one covariate across its observed range while holding
all other numeric predictors at their sample mean. For `unmarked`
models, other categorical predictors are held at their first observed
level.

## Usage

``` r
plot_marginal_effects(
  x,
  ylab = "Estimated response",
  length.out = 100,
  type = "state",
  link = FALSE
)
```

## Arguments

- x:

  A fitted `multiScaleR` object returned by
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
  Must contain `opt_mod` (the optimized fitted model) and `scl_params`
  (a list with `mean` and `sd` elements for each covariate, used to
  back-transform scaled covariates to their original units for the
  x-axis).

- ylab:

  Character scalar. Y-axis label for all marginal effect panels.
  Default: `"Estimated response"`.

- length.out:

  Positive integer. Number of equally spaced points at which to evaluate
  the marginal effect curve for each covariate. Default: `100`.

- type:

  Character. Prediction type for `unmarked` models. One of `"state"`
  (default) or `"lambda"`. Ignored for non-unmarked models.

- link:

  Logical. If `TRUE`, predictions are passed through the model's inverse
  link function before plotting (i.e., plotted on the response scale).
  For `glm` objects this is applied automatically. Set to `TRUE`
  manually if predicted values appear to be on the link scale. Default:
  `FALSE`.

## Value

A named list of `ggplot2` objects, one per covariate. Plots are printed
as a side effect and the list is returned invisibly. Each plot shows:

- The predicted response as a solid line for numeric predictors or
  points for categorical predictors, with x-axis values back-transformed
  to the original (unscaled) units when applicable.

- A shaded ribbon or categorical error bars for the 95% confidence
  interval (when available from the model's `predict` method).

- A subtitle noting any polynomial (`I(x^2)`) or interaction terms that
  involve the covariate.

## Details

Marginal effects are computed by constructing a `newdata` data frame
where the focal covariate sweeps from its minimum to maximum observed
value and all other predictors are held at their sample mean. For
kernel-scaled covariates, whose scaled mean is 0, this is equivalent to
holding them at their average landscape value. Variables used only in
model offsets are held at their sample mean for prediction but do not
receive marginal effect panels.

For `glm` and most GLM-family models, predictions are transformed via
the inverse link function automatically (e.g.,
[`exp()`](https://rdrr.io/r/base/Log.html) for Poisson log-link models).
Confidence intervals are constructed from `predict(..., se.fit = TRUE)`
using \\\pm 1.96 \times SE\\.

For `unmarked` models, `predict(mod, newdata, type = type)` is used, and
lower/upper bounds come from the `lower` and `upper` columns of the
returned data frame. Predictors are recovered from the selected
component's formula, so factor contrasts such as `Year2` are not
mistaken for data columns. Factors receive one point and confidence
interval per observed level. Site keys not used as fixed predictors and
detection-only predictors do not receive state-effect panels. Numeric
predictors not optimized by
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
retain their original units and are held at their sample mean in other
panels.

For `unmarked` random-effects models, prediction uses `re.form = NA`,
which sets random effects to zero. These curves are conditional
predictions at the stated covariate values, not averages over years,
sites, or the random-effect distribution. Confidence intervals are
conditional on the fitted spatial scale. Explicit nondefault factor
contrasts in ordinary `pcount` (Poisson or negative-binomial) and `occu`
fits use
[`unmarked::linearComb()`](https://ecoverseR.github.io/unmarked/reference/linearComb-methods.html)
and
[`unmarked::backTransform()`](https://ecoverseR.github.io/unmarked/reference/backTransform-methods.html)
to preserve the fitted design and link. Other unmarked classes with such
contrasts stop if their prediction method cannot safely retain that
design.

For models with interaction terms, a message is printed explaining that
each panel represents a conditional slice (at the mean of the
interacting variable) rather than the full interaction surface.

`HLfit` (spaMM) and `zeroinfl` (pscl) models are handled with
class-specific prediction calls. When `x$opt_mod` is a wrapper object
around another fitted model, `plot_marginal_effects()` uses the nested
analysis model for predictor discovery, data recovery, and prediction
when possible.

## Examples

``` r
# \donttest{
data('pts')
data('count_data')
hab <- terra::rast(system.file('extdata', 'hab.tif', package = 'multiScaleR'))

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 250,
                             kernel = 'gaussian')
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...

mod <- glm(y ~ hab, family = poisson, data = count_data)
opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=250.56 to ensure accurate estimation of scale.
#> 
#> 

## Default marginal effect plot on the response scale
plot_marginal_effects(opt)


## Custom y-axis label
plot_marginal_effects(opt, ylab = "Predicted abundance")


## Finer evaluation grid
plot_marginal_effects(opt, length.out = 200)

# }
```

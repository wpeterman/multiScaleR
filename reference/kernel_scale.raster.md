# Apply kernel smoothing to raster layers

Applies a kernel smoothing function to one or more raster layers,
producing a spatially weighted mean (or landscape metric) at the scale
identified by
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
Primarily used to generate prediction rasters for spatial model
projection.

## Usage

``` r
kernel_scale.raster(
  raster_stack,
  sigma = NULL,
  multiScaleR = NULL,
  shape = NULL,
  kernel = c("gaussian", "exp", "expow", "fixed"),
  scale_vars = NULL,
  pct_wt = 0.975,
  fft = TRUE,
  scale_center = FALSE,
  clamp = FALSE,
  pct_mx = 0,
  na.rm = TRUE,
  verbose = TRUE,
  offset_values = NULL,
  ...
)
```

## Arguments

- raster_stack:

  A `SpatRaster` object containing the source raster layer(s) to be
  smoothed. Layer names must match the covariate names used during model
  fitting (or the `source` names in `scale_vars`).

- sigma:

  Numeric vector of kernel scale parameter(s) in the same units as the
  raster projection (e.g., meters). One value per raster layer. Ignored
  when `multiScaleR` is provided; sigma values are extracted from the
  fitted object automatically.

- multiScaleR:

  A fitted object of class `"multiScaleR"` (from
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md))
  or class `"multiScaleR_data"` (from
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)).
  When provided, sigma, shape, kernel, and scale_vars are all extracted
  automatically. Default: `NULL`.

- shape:

  Numeric vector of shape parameters for the exponential power kernel
  (`kernel = "expow"`), one per raster layer. Ignored when `multiScaleR`
  is provided. Default: `NULL`.

- kernel:

  Character. Kernel function used for smoothing. One of `"gaussian"`
  (default), `"exp"`, `"fixed"`, or `"expow"`. Ignored when
  `multiScaleR` is provided.

- scale_vars:

  Optional variable specifications created with
  [`msr_vars`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md).
  Use when projecting landscape metrics or explicitly defined covariates
  without passing a fitted `multiScaleR` object. When `multiScaleR` is
  provided, `scale_vars` is extracted automatically.

- pct_wt:

  Numeric between 0 and 1 (exclusive). Cumulative kernel weight along a
  line through the focal cell used to determine the smoothing window
  radius. A larger value (e.g., `0.99`) includes more of the kernel tail
  but increases computation time. This is the historical one-dimensional
  cutoff, not the proportion of total weight within the planar window.
  Default: `0.975`.

- fft:

  Logical. If `TRUE` (default), smoothing is performed via Fast Fourier
  Transform (FFT) convolution, which is substantially faster for large
  rasters and wide kernels. Some edge effects may occur at raster
  boundaries. Set to `FALSE` to use
  [`terra::focal`](https://rspatial.github.io/terra/reference/focal.html),
  which avoids edge effects but is slower.

- scale_center:

  Logical. If `TRUE`, the smoothed raster values are centered and scaled
  using the mean and standard deviation from the model fitting data
  (extracted from `multiScaleR$scl_params`). Required when using the
  output with
  [`terra::predict`](https://rspatial.github.io/terra/reference/predict.html)
  and a fitted model that was trained on scaled covariates. Requires
  `multiScaleR` to be provided. Default: `FALSE`.

- clamp:

  Logical. If `TRUE`, scaled raster values are clamped to the observed
  covariate range from the model fitting data, preventing extrapolation
  beyond the training range. Only active when `scale_center = TRUE`.
  Default: `FALSE`.

- pct_mx:

  Numeric between -0.99 and 0.99. When `clamp = TRUE`, expands (`> 0`)
  or contracts (`< 0`) the clamping range relative to the observed
  min/max by this proportion. Default: `0` (exact training range).

- na.rm:

  Logical. If `TRUE` (default), `NA` cells are excluded from the
  kernel-weighted mean calculation. If `FALSE`, any window containing a
  `NA` cell will produce a `NA` output.

- verbose:

  Logical. Print layer-level progress messages. Default: `TRUE`.

- offset_values:

  Named numeric vector of values for variables inside model offsets,
  such as `c(n_weeks = 4)` for `offset(log(n_weeks))`. These are raw
  input values, not transformed offsets. When omitted, a log offset uses
  1 (zero on the log scale) and a direct offset uses 0. Set this
  argument for a meaningful projection scenario.

- ...:

  Reserved for deprecated arguments. Currently only `scale_opt`
  (deprecated alias for `multiScaleR`) is handled.

## Value

A `SpatRaster` with one layer per covariate defined in `scale_vars` (or
per raster layer when `scale_vars` is not used). Layer names match the
covariate names from the model. When a fitted `multiScaleR` object is
provided and the model contains site-level predictors (i.e., predictors
not derived from raster layers), constant ("dummy") raster layers filled
with zeros are appended to make the output compatible with
[`terra::predict`](https://rspatial.github.io/terra/reference/predict.html).
Offset variables are appended under their input names so transformed
offset expressions can be evaluated by
[`terra::predict`](https://rspatial.github.io/terra/reference/predict.html).

## Details

**Typical usage** after running
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md):


    opt_hab <- kernel_scale.raster(raster_stack = hab, multiScaleR = opt)
    plot(c(hab, opt_hab))

    ## Scale and center for prediction
    opt_hab_sc <- kernel_scale.raster(raster_stack = hab,
                                      multiScaleR = opt,
                                      scale_center = TRUE)
    preds <- terra::predict(opt_hab_sc, opt$opt_mod, type = "response")

**FFT vs. focal smoothing**

The FFT convolution (`fft = TRUE`) is substantially faster for large
rasters or wide kernels. It produces minor edge artifacts near raster
boundaries, typically within one kernel-width of the edge. For analyses
focused on interior areas this is usually negligible. Use `fft = FALSE`
for exact focal smoothing when edge accuracy matters.

**Dummy rasters for site-level covariates**

When a fitted `multiScaleR` object is supplied, the function inspects
the model frame to find any predictors that are not represented by
raster layers (e.g., survey effort, habitat type measured in the field).
These cannot be projected spatially and are therefore assigned constant
zero rasters, which correspond to the reference value for centered and
scaled covariates. These dummy layers are required for
[`terra::predict`](https://rspatial.github.io/terra/reference/predict.html)
to work but do not represent real spatial variation. Replace them
manually for non-zero projection scenarios. Categorical predictors are
skipped with a warning. Offset terms use their underlying input
variables as layer names. For `offset(log(n_weeks))`, the default layer
is `n_weeks = 1`, representing one week of effort. Set
`offset_values = c(n_weeks = 4)` to predict for four weeks. Direct
offsets default to zero.

## Examples

``` r
## Not Run
r1 <- rast(matrix(rnorm(25^2),
                  nrow = 25))

r1_s <- kernel_scale.raster(r1,
                            sigma = 4,
                            kernel = 'gaussian')
#> 
#> Smoothing spatRaster 1 of 1: lyr.1 at sigma = 4
plot(c(r1, r1_s))

```

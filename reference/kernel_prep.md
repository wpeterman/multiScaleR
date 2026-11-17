# Prepare Kernel Inputs for Scale Optimization

Extracts the raster cells around each survey point, records the distance
from the point to each cell, and computes a first set of kernel-weighted
covariates. Fit your starting model to those covariates, then pass the
model and the returned object to
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
This is the first step of every `multiScaleR` workflow, shown in
[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md).

## Usage

``` r
kernel_prep(
  pts,
  raster_stack,
  max_D,
  kernel = c("gaussian", "exp", "expow", "fixed"),
  scale_vars = NULL,
  sigma = NULL,
  shape = NULL,
  projected = TRUE,
  progress = FALSE,
  verbose = TRUE,
  bin = TRUE,
  nbins = 256L,
  store_cell_data = TRUE
)
```

## Arguments

- pts:

  Survey point locations as a `SpatVector` (terra) or `sf` object, in a
  projected coordinate system whose units match `max_D`. Row names are
  kept as point IDs.

- raster_stack:

  A `SpatRaster` (terra) with one or more layers, on the same projection
  as `pts`. Without `scale_vars`, each layer becomes a covariate with
  the layer's name, so the names must match the model formula. Layer
  names cannot be `"x"`, which is reserved for internal coordinate data.
  The raster must cover a buffer of `max_D` around every point.

- max_D:

  Positive numeric scalar. Radius, in map units, of the circle around
  each point within which raster cells are extracted. It is also the
  upper bound of the sigma search, and kernels are cut off at this
  distance. Set it above the largest plausible scale of effect, at least
  about 3.3 times the largest plausible sigma for a Gaussian kernel. If
  an estimate comes back near `max_D`, increase it and run
  `kernel_prep()` again. See the section on choosing `max_D` in
  [multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md).

- kernel:

  Character. Distance-weighting function used for every covariate:
  `"gaussian"` (default), `"exp"` (negative exponential), `"expow"`
  (exponential power, which also estimates a shape), or `"fixed"` (equal
  weight within a radius). The weight functions are given in
  [multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md).
  To compare kernels, run `kernel_prep()` once per kernel.

- scale_vars:

  Optional covariate specifications created with
  [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md).
  Omit it to make each raster layer one kernel-weighted mean. Supply it
  to define derived covariates, such as a kernel mean and a landscape or
  surface metric from the same layer.

- sigma:

  Optional numeric vector of starting sigma values, in map units, one
  per covariate whose scale is estimated (one per raster layer when
  `scale_vars` is omitted). Default: `NULL`, which starts every sigma at
  `max_D / 2`. The default is recommended. These values set the
  covariates in `kernel_dat` and are not the starting values of
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  which chooses its own.

- shape:

  Optional numeric vector of starting shape values for
  `kernel = "expow"`, with the same length as `sigma`. Default: `NULL`,
  which starts every shape at 2 (the Gaussian form). Ignored for other
  kernels.

- projected:

  Logical. Whether `pts` and `raster_stack` use projected (planar)
  coordinates. Only projected coordinates are supported, so `FALSE` is
  an error. Default: `TRUE`.

- progress:

  Logical. Show progress bars during extraction and weighting. Default:
  `FALSE`.

- verbose:

  Logical. Print status messages. Default: `TRUE`.

- bin:

  Logical. If `TRUE` (default), summarize each kernel mean's cell values
  into distance bins.
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  then evaluates kernel means from the bins instead of every cell, which
  is much faster and no longer depends on `max_D`, at the cost of a
  negligible approximation. The analytic backends and
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  require it for kernel means. Landscape and surface metrics are never
  binned. `FALSE` forces exact cell-by-cell evaluation.

- nbins:

  Integer, at least 2. Number of equal-width distance bins when
  `bin = TRUE`. More bins increase accuracy and memory. The default
  (`256`) reproduces exact optimization results to several significant
  figures. Ignored when `bin = FALSE`.

- store_cell_data:

  Logical. If `TRUE` (default), the per-cell distances (`d_list`) and
  values (`raw_cov`) are kept. If `FALSE`, they are dropped after
  binning, often shrinking the object by one to two orders of magnitude.
  This lean mode needs `bin = TRUE` and only kernel means. With
  landscape or surface metrics, or `bin = FALSE`, the cell data are kept
  and a message says so. The bins suffice for optimization,
  [`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md),
  and the final refit.

## Value

An object of class `multiScaleR_data`, a list. Its full structure is
described in
[multiScaleR_data](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_data-class.md).
The parts you use directly:

- `kernel_dat`: data frame with one row per point and one column per
  covariate. Values are centered and scaled (mean 0, standard
  deviation 1) at the starting sigma. Combine it with your response data
  to fit the starting model.

- `scl_params`: list with the `mean` and `sd` used for that scaling, in
  the covariates' own units.

- `kernel`, `max_D`, `min_D`, and `scale_vars`: the settings passed on
  to
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  and
  [`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md).

Printing the object reports the number of points and covariates, the
kernel, the number of cells extracted around the first point ("Number of
elements"), whether binning is on, the minimum distance (`min_D`, about
one cell width and the lower bound of the sigma search), the maximum
distance (`max_D`), and the unit conversion (distances are stored
divided by `max_D`).

## Prepare the points and rasters

Points and rasters must share a defined, projected CRS. Every distance,
including `max_D` and any `sigma` you supply, is in the linear unit of
that projection, usually meters. Landscape and surface metrics also need
square cells. Row names of `pts` are kept throughout, so model data can
be joined back to the points.

Each point's buffer is the circle of radius `max_D`, and every cell it
intersects is extracted, including partly covered boundary cells, which
keep their full kernel weight. Missing cells drop out of each weighted
mean, and the remaining weights renormalize. Extraction runs in batches
of 25 points to limit the peak memory of temporary tables. Point order,
returned values, and distances are not affected.

## Manage memory

A larger `max_D` extracts more cells per point. Check the size of the
object and the memory needed for parallel optimization with
[`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md).
When every covariate is a kernel mean, `store_cell_data = FALSE` keeps
only the bins.

## Related preparation functions

To assign each observation to a different annual raster while estimating
one pooled landscape relationship, use
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).
To prepare inputs for `unmarked` models, see
[`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md).

## See also

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md),
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md),
[multiScaleR_data](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_data-class.md),
[`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md).
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md)
explains each argument on a worked example.

## Examples

``` r
## Survey points and a binary habitat raster bundled with the package
data("pts")
data("count_data")
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 250,          # map units (meters)
                             kernel = "gaussian",
                             verbose = FALSE)
# 75 points, one covariate, cells within 250 m of each point
kernel_inputs
#> 
#> There are 75 observations at 1 spatial covariate(s): 
#> hab
#> 
#> The specified kernel is:
#> gaussian
#> 
#> Number of elements: 
#> 2052
#> Distance binning:
#> enabled (256 distance bins per covariate)
#> Minimum Distance:
#> 10
#> Maximum Distance:
#> 250
#> Unit Conversion:
#> 250

# Standardized kernel-weighted habitat at the starting sigma (125 m)
head(kernel_inputs$kernel_dat)
#>          hab
#> 1 -0.4206979
#> 2  1.7520402
#> 3  0.5073754
#> 4 -0.8959421
#> 5 -0.1734359
#> 6  0.8108489

# Fit the starting model that multiScale_optim() will optimize
df <- data.frame(count_data["y"], kernel_inputs$kernel_dat)
mod <- glm(y ~ hab, family = poisson, data = df)

## Lean mode keeps only the distance bins
lean_inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                           store_cell_data = FALSE, verbose = FALSE)
# Object sizes in megabytes
c(full = as.numeric(object.size(kernel_inputs)),
  lean = as.numeric(object.size(lean_inputs))) / 1e6
#>     full     lean 
#> 2.547640 0.328192 
```

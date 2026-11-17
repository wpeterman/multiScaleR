# Prepared inputs for scale estimation

[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
and
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
return an object of class `"multiScaleR_data"`, a list that holds
everything the fitting functions need about the survey points and the
raster cells around them. Pass it as `kernel_inputs` to
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
and related functions. This page describes every element.

## Elements

- `kernel_dat`: data frame of the initial kernel-weighted covariate
  values, centered and scaled to unit variance, with one row per point
  and one column per covariate. Row names match the row names of `pts`
  or are sequential integers. Fit your starting model to these columns.

- `d_list`: named list, one element per point, of distances from the
  point to every raster cell within the buffer. `NULL` in lean mode.

- `raw_cov`: named list, one element per point, of matrices holding the
  raw raster cell values within the buffer, aligned with `d_list`. Each
  is stored as a dense or sparse matrix, whichever is smaller. `NULL` in
  lean mode.

- `kernel`: the kernel name (`"gaussian"`, `"exp"`, `"expow"`, or
  `"fixed"`).

- `sigma`: initial sigma values on the internal scale (map units divided
  by `unit_conv`).

- `shape`: initial shape values for `"expow"`, otherwise `NULL`.

- `min_D`: approximate minimum distance between cells, used as the lower
  bound for sigma during the search, in map units.

- `max_D`: the `max_D` you supplied, the upper bound for sigma, in map
  units.

- `n_covs`: number of covariates whose scale is estimated.

- `unit_conv`: internal distance scaling factor, equal to `max_D`.

- `scale_vars`: data frame of class `"multiScaleR_vars"` describing
  every covariate specification (see
  [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)).

- `resolution`: raster cell resolution in map units.

- `cell_size`: cell width and height in map units. Cross-validation uses
  them to keep intersecting cells whose centers lie beyond `max_D`.

- `n_cols`: number of raster columns, used by adjacency-based landscape
  metrics.

- `binned`: when `bin = TRUE`, a list of distance-binned summaries
  (representative bin distances and value and count matrices for each
  covariate) that speed up the search. `NULL` when `bin = FALSE` or when
  no covariate is a kernel-weighted mean.

- `cell_data_stored`: `FALSE` only in lean mode
  (`store_cell_data = FALSE`), when `d_list` and `raw_cov` were dropped.

- `scl_params`: named list with `mean` and `sd`, the centering and
  scaling values applied to `kernel_dat`.

- `raster_group`: from
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
  only, the raster group of each point.

Lean mode keeps only the binned summaries. It works only when every
covariate is a kernel-weighted mean, because landscape and surface
metrics need the cell-level data. Printing the object reports the number
of points and covariates, the kernel, the binning, and the distance
bounds.

## See also

[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md),
[`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md),
[multiScaleR](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md),
and
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md).

## Examples

``` r
data(pts)
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian", verbose = FALSE)
inputs                    # short summary of the prepared inputs
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
head(inputs$kernel_dat)   # standardized covariate values at the start
#>          hab
#> 1 -0.4206979
#> 2  1.7520402
#> 3  0.5073754
#> 4 -0.8959421
#> 5 -0.1734359
#> 6  0.8108489
inputs$max_D              # upper bound for sigma, in map units
#> [1] 250
```

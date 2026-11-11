# Prepare Kernel Inputs from Group-Specific Rasters

Assign each observation to one raster map, then prepare one pooled data
object for
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
Use this when land cover changes across years but the model estimates
shared landscape effects and scales.

## Usage

``` r
kernel_prep_by_group(
  pts,
  raster_stacks,
  group,
  max_D,
  kernel = c("gaussian", "exp", "expow", "fixed"),
  scale_vars = NULL,
  sigma = NULL,
  shape = NULL,
  bin = TRUE,
  nbins = 256L,
  store_cell_data = TRUE,
  progress = FALSE,
  verbose = TRUE
)
```

## Arguments

- pts:

  Projected point locations as an `sf` or `SpatVector` object. Point
  order matches `group` and the returned covariate rows.

- raster_stacks:

  Named list of `SpatRaster` maps. Each map has the same extent, grid,
  CRS, layer names, and class coding. Its name identifies a year or
  other observation group.

- group:

  Character or factor vector, one value per point, identifying the
  corresponding name in `raster_stacks`. Every point uses only its
  assigned map, including cells affected by changes outside the point's
  immediate location.

- max_D:

  Maximum extraction radius in the projected map units.

- kernel:

  Kernel weighting function; see
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

- scale_vars:

  Optional
  [`msr_vars`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  specification. By default, each raster layer becomes a kernel-weighted
  mean covariate.

- sigma:

  Optional starting scales in map units. Defaults to `max_D / 2` for
  each optimized covariate.

- shape:

  Optional starting shapes for the exponential power kernel. Defaults to
  2 for each optimized covariate.

- bin:

  Logical. Build pooled distance bins for kernel covariates. Default:
  `TRUE`. All observations use one set of bin definitions.

- nbins:

  Number of distance bins when `bin = TRUE`; default 256.

- store_cell_data:

  Logical. Retain per-cell values and distances. Default: `TRUE`. Set
  `FALSE` for a lean kernel-only binned object; metrics requiring cell
  detail retain it regardless.

- progress:

  Logical. Display extraction progress. Default `FALSE`.

- verbose:

  Logical. Print preparation messages. Default `TRUE`.

## Value

A `multiScaleR_data` list containing:

- `kernel_dat`:

  Initial covariate data frame, with one row per point. Columns are
  centered and divided by their pooled standard deviations. Use this in
  the starting model.

- `d_list`, `raw_cov`:

  Named per-point distance vectors (divided by `unit_conv`) and
  raster-value matrices. Both are `NULL` in lean mode.

- `binned`:

  Pooled distance-bin summaries for kernel covariates, or `NULL` when
  `bin = FALSE`. Each covariate has finite-value sums and counts by
  point and distance bin.

- `raster_group`:

  Named character vector recording the map assigned to each point, in
  input order.

- `scl_params`:

  Pooled means and standard deviations used to create `kernel_dat`; the
  fitted model replaces these with scaling parameters at the optimized
  scales.

- `kernel`, `sigma`, `shape`:

  Kernel name, initial scales divided by `unit_conv`, and optional
  initial shape values.

- `min_D`, `max_D`, `unit_conv`:

  Lower scale bound, maximum extraction radius, and internal distance
  scaling factor. These are in map units except the internally scaled
  `sigma`.

- `n_covs`, `scale_vars`:

  Number of optimized covariates and the validated covariate
  definitions.

- `resolution`, `n_cols`:

  Raster cell size and column count, used to calculate landscape
  configuration metrics.

- `cell_data_stored`:

  Whether per-cell values and distances are retained.

Pass the full object to
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

## Details

Raster maps must be aligned because landscape configuration metrics use
cell adjacency and the resulting model assumes common map units.
Group-specific maps may contain different missing cells. For kernel
covariates, pooled bins track both finite-cell value sums and counts.
The fitted model has common slopes and scales unless its formula
explicitly includes group interactions. Choose the response likelihood
and any group effect according to the sampling design. A pooled model
does not establish that effects are stable across years.

## See also

[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)

## Examples

``` r
library(terra)
r1 <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 200,
           ymin = 0, ymax = 200, crs = "EPSG:26915")
values(r1) <- seq_len(ncell(r1)) / ncell(r1)
names(r1) <- "habitat"
r2 <- r1
r2[seq_len(30)] <- 0
pts <- vect(cbind(c(50, 70, 90, 110, 130, 150),
                  c(60, 90, 120, 140, 100, 70)),
            type = "points", crs = crs(r1))
prep <- kernel_prep_by_group(pts, list(y1 = r1, y2 = r2),
                             group = rep(c("y1", "y2"), 3),
                             max_D = 25, verbose = FALSE)
head(prep$kernel_dat)
#>       habitat
#> 1  1.14516205
#> 2  0.17317085
#> 3 -0.79882036
#> 4 -1.43564218
#> 5 -0.06144772
#> 6  0.97757736
```

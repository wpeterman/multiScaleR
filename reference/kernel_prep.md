# Kernel Scale Preparation

Prepares the data inputs required for multiscale kernel optimization.
Extracts raster values within a buffer around each point, computes
pairwise distances, and calculates initial kernel-weighted covariate
values. The result is passed directly to
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

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

  Spatial point locations as a `SpatVector` (terra) or `sf` object.
  Points must be in a projected coordinate system that shares units with
  `max_D`.

- raster_stack:

  One or more raster layers of class `SpatRaster` (terra). Layer names
  must match the variable names in the model formula used with
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  unless `scale_vars` is provided to define derived covariates. Layer
  names cannot be `"x"`, which is reserved for internal coordinate data.

- max_D:

  Positive numeric. The maximum radius (in the same units as the
  projection of `pts` and `raster_stack`) around each point to sample
  from the raster. Should be at least 2–3 times the largest expected
  sigma value. Raster layers must extend far enough beyond the points to
  fully cover buffers of this radius.

- kernel:

  Kernel function to be used for weighting raster values by distance.
  One of:

  `"gaussian"`

  :   (Default) Gaussian (normal) decay: weights fall off as a normal
      distribution with standard deviation `sigma`.

  `"exp"`

  :   Negative exponential decay: weights fall off as `exp(-d / sigma)`.

  `"fixed"`

  :   Fixed-radius (step) buffer: all cells within `sigma` receive equal
      weight; cells beyond receive zero.

  `"expow"`

  :   Exponential power kernel parameterized by both `sigma` (scale) and
      `shape` (shape). Requires `shape` to be specified.

- scale_vars:

  Optional variable specifications created with
  [`msr_vars`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md).
  When omitted, each raster layer becomes one kernel-weighted covariate
  with the same name, preserving historical behavior. Provide
  `scale_vars` to define derived model covariates, for example,
  combining a kernel-weighted mean with landscape composition or edge
  metrics from the same source layer. See
  [`kernel_var`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  and
  [`landscape_var`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md).

- sigma:

  Optional numeric vector of initial sigma values for optimization.
  Length must equal the number of covariates with `optimize = TRUE` in
  `scale_vars` (or the number of raster layers when `scale_vars` is not
  provided). Values must be in the same units as the projection.
  Default: `NULL`; initial values are generated automatically as
  `max_D / 2`, which is recommended.

- shape:

  Optional numeric vector of initial shape values when
  `kernel = "expow"`. Length requirements same as `sigma`. Default:
  `NULL`; starting values of 2 are generated automatically.

- projected:

  Logical. Whether `pts` and `raster_stack` are in a projected (planar)
  coordinate system. Currently only projected coordinates are supported.
  Default: `TRUE`.

- progress:

  Logical. Print progress bars to the console during extraction and
  distance calculations. Default: `FALSE`.

- verbose:

  Logical. Print status messages during preparation. Default: `TRUE`.

- bin:

  Logical. If `TRUE` (default), precompute distance-binned summaries of
  the kernel-type covariates so that
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  can evaluate kernel-weighted means in `O(nbins)` per point instead of
  iterating every buffer cell on each optimizer evaluation. This makes
  optimization substantially faster (and independent of `max_D`) at the
  cost of a negligible binning approximation. Landscape-metric
  covariates are never binned and always use the exact per-cell path.

- nbins:

  Positive integer. Number of equal-width distance bins used when
  `bin = TRUE`. More bins increase accuracy and memory; the default
  (`256`) reproduces exact optimization results to several significant
  figures. Ignored when `bin = FALSE`.

- store_cell_data:

  Logical. If `TRUE` (default), the full per-cell `d_list`/`raw_cov`
  data are retained. If `FALSE` and every covariate is kernel-type (no
  landscape metrics) and `bin = TRUE`, the cell-level data are dropped
  after binning, dramatically reducing the object's memory footprint
  (often by 1–2 orders of magnitude). The binned summaries are
  sufficient to drive optimization,
  [`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md),
  and the final refit. Ignored (with a message) when landscape metrics
  are present or `bin = FALSE`.

## Value

A list of class `"multiScaleR_data"` containing:

- `kernel_dat`:

  Data frame of scaled (mean-centered, unit-variance) initial
  kernel-weighted covariate values, one row per point and one column per
  covariate. Row names match `pts` row names or sequential integers. Use
  this directly to fit the initial model passed to
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- `d_list`:

  Named list (one element per point) of numeric distance vectors from
  the point to every raster cell within the buffer.

- `raw_cov`:

  Named list (one element per point) of matrices containing raw raster
  cell values within the buffer, aligned with `d_list`. Each point is
  stored as a dense or sparse matrix, whichever is smaller. `NULL` in
  lean mode.

- `kernel`:

  Character string identifying the kernel used.

- `sigma`:

  Numeric vector of initial sigma values on the internal (scaled)
  parameter space.

- `shape`:

  Numeric vector of initial shape values, or `NULL`.

- `min_D`:

  Numeric. Approximate minimum inter-cell distance, used as the lower
  bound for sigma during optimization.

- `max_D`:

  Numeric. The `max_D` value supplied by the user.

- `n_covs`:

  Integer. Number of covariates with scale being optimized.

- `unit_conv`:

  Numeric. Internal distance scaling factor (equals `max_D`).

- `scale_vars`:

  Data frame of class `"multiScaleR_vars"` describing all covariate
  specifications.

- `resolution`:

  Numeric. Raster cell resolution.

- `n_cols`:

  Integer. Number of raster columns (used for adjacency landscape
  metrics).

- `scl_params`:

  Named list with elements `mean` and `sd` the centering and scaling
  parameters applied to `kernel_dat`. Stored for use by
  [`kernel_scale.raster`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
  when `scale_center = TRUE`.

- `binned`:

  When `bin = TRUE`, a list of precomputed distance-binned summaries
  (representative bin distances and per-covariate value/count matrices)
  used to accelerate optimization. `NULL` when `bin = FALSE`.

- `cell_data_stored`:

  Logical flag indicating whether the per-cell `d_list`/`raw_cov` data
  were retained. `FALSE` only in lean mode (`store_cell_data = FALSE`).

## Details

Point locations and raster layers must share a defined CRS and both must
be projected. All units (including `max_D` and any user-supplied
`sigma`) must be in the same linear unit as the projection (typically
meters or feet).

The `max_D` buffer should be large enough to encompass the plausible
range of the true scale of effect. A common rule of thumb is to set
`max_D` to at least 2–3 times the largest expected sigma. After running
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
the `max_distance` diagnostic will warn if the estimated scale
approaches `max_D`.

Initial `sigma` values do not need to be precise. The optimizer will
refine them. Provide explicit starting values only if the default
(`max_D / 2`) leads to convergence problems.

Row names from `pts` are preserved throughout the returned object so
that downstream model data frames can be joined back to the original
point order. Raster extraction is processed in batches of 25 points to
limit the peak memory used by temporary extraction tables. Point order,
returned values, and distance calculations are unchanged. To assign each
observation to a different annual raster while estimating one pooled
landscape relationship, use
[`kernel_prep_by_group`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).

## Examples

``` r
library(terra)
pts <- vect(cbind(c(3,5,7),
                  c(7,5,3)))

mat_list <- list(r1 = rast(matrix(rnorm(100),
                                  nrow = 10)),
                 r2 = rast(matrix(rnorm(100),
                                  nrow = 10)))
rast_stack <- rast(mat_list)
kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = rast_stack,
                             max_D = 2,
                             kernel = 'gaussian',
                             sigma = NULL)
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...
```

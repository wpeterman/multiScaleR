# Create Spatially Blocked Cross-Validation Folds

Lay a square grid over projected points and assign whole grid blocks to
assessment folds, so that the model must predict parts of the landscape
it was not fitted to. Optional sampling groups join the blocks they
touch into indivisible units. A distance buffer can exclude nearby
training units without dropping assessment rows.

## Usage

``` r
spatial_cv_folds(
  coords,
  block_size,
  k = 5,
  ids = NULL,
  group = NULL,
  buffer = 0,
  origin = c(0, 0),
  crs = NULL
)
```

## Arguments

- coords:

  Projected `sf` POINTs, a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  of points, or a numeric two-column coordinate matrix or data frame.
  Row order must match the data. Meaningful row names must match `ids`.

- block_size:

  Positive width of each grid cell, in the units of the coordinate CRS
  (usually meters). Larger blocks put more distance between held-out
  points and their nearest training points, which moves the task closer
  to predicting a new area.

- k:

  Number of folds, at least two. The default is five. At least `k`
  indivisible spatial units must remain after joining groups.

- ids:

  Unique character row IDs. Defaults to the row names of `coords`, or to
  the row numbers as character strings. Use the same IDs as the modeling
  data.

- group:

  Optional sampling-group labels aligned with `ids`, such as site or
  animal IDs. Blocks joined by any shared group are kept together,
  including through chains of shared groups. This can reduce the number
  of available units.

- buffer:

  Nonnegative training-to-assessment separation, in coordinate CRS
  units. The default 0 adds no exclusions. With a positive buffer, an
  entire training unit is excluded from a fold if any of its points lies
  within, or exactly at, this distance from an assessment point of that
  fold. Assessment rows are never dropped. A buffer that leaves a fold
  with no training rows is an error.

- origin:

  Numeric grid origin of length two, in coordinate CRS units. Defaults
  to `c(0, 0)`. Cells are half-open intervals of width `block_size`,
  indexed by `floor((coordinate - origin) / block_size)`.

- crs:

  Projected CRS declaration for coordinate tables, in any form accepted
  by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).
  For spatial objects, an optional declaration must match their existing
  CRS. Coordinates are never transformed or relabeled here, and
  longitude and latitude are rejected.

## Value

A named list of class `multiScaleR_folds`, structured as described in
[`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md).
Each fold holds `analysis` and `assessment` row IDs. Every row is
assessed exactly once, although a buffer can leave a row out of some
training sets. The attribute `fold_info` holds the `method`
(`"spatial"`), `k`, `block_size`, `buffer`, `origin`, the `crs` (as WKT)
and its distance `units`, an `assignment` data frame with each row's
`row_id`, coordinates `x` and `y`, grid `block`, joined `unit`,
assessment `fold`, and `group` when supplied, and `excluded`, the
training IDs removed by the buffer in each fold.

## Details

Units are assigned largest first to balance row counts, with stable
tie-breaking and no random draws, so the same call always gives the same
folds. Several disjoint blocks can share a fold, so a fold need not be
one contiguous region. Buffers measure the separation between sampled
points, not the distance to block edges. Because whole units are
excluded, the gap can be larger than requested, and a large buffer can
leave few training points. Distance calculations use nearest-point
searches rather than a dense training-by-assessment distance matrix.

Choose the block size, groups, buffer, and prediction target from the
sampling design and the ecological question. A projected CRS and a
distance buffer alone do not establish independence or validate
calibration.
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
maps a set of spatial folds and shows how a buffer changes the training
sets.

## See also

[`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md),
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)

## Examples

``` r
# Six points in three 50 m blocks (crs = 32617 declares meter units)
xy <- cbind(x = rep(c(0, 100, 200), each = 2), y = rep(c(0, 10), 3))
spatial_cv_folds(xy, block_size = 50, k = 3, crs = 32617)
#> Cross-validation folds: spatial 
#>   fold analysis assessment excluded
#>  fold1        4          2        0
#>  fold2        4          2        0
#>  fold3        4          2        0

# The bundled survey points, in 1,000 m blocks. A 320 m buffer removes
# training blocks that have a point near any held-out point.
data("surv_pts")
pts <- crds(vect(surv_pts))
buffered <- spatial_cv_folds(pts, block_size = 1000, k = 3, buffer = 320,
                             crs = 32617)
buffered   # 'excluded' counts rows used for neither fitting nor scoring
#> Cross-validation folds: spatial 
#>   fold analysis assessment excluded
#>  fold1       22         34       44
#>  fold2       11         33       56
#>  fold3       13         33       54
head(attr(buffered, "fold_info")$assignment)   # block and fold of each row
#>   row_id        x         y block unit  fold
#> 1      1 3552.415 2930.8536    15   15 fold2
#> 2      2 3241.634 1376.9512    14   14 fold1
#> 3      3 2620.073 1998.5122    10   10 fold3
#> 4      4 2620.073  755.3902     9    9 fold3
#> 5      5 3552.415 1998.5122    14   14 fold1
#> 6      6 1376.951 1998.5122     6    6 fold1
```

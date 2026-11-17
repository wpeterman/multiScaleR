# Bind Projected Coordinates to a Fitted Response Field

Check the coordinates used to fit a spatial response field and package
them, with their CRS, into a small descriptor for raster prediction. A
response field is a spatially correlated random effect on the response,
separate from the kernel scales of effect. To map such a model, each
raster cell needs its own field value, and the descriptor tells
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
how to compute it from cell coordinates.

## Usage

``` r
response_field_geometry(
  object,
  coords,
  crs,
  batch_size = 128L,
  max_memory_mb = 512
)
```

## Arguments

- object:

  A `multiScaleR` result from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  with the default `backend = "refit"`, whose fitted model is one of:

  - a `glmmTMB` model with one shared exponential, Gaussian or Matern
    field, such as `exp(pos + 0 | grp)`, and no other random effects.
    The family must be Gaussian, Poisson, NB2 or Gamma with a log link,
    or binomial with a logit link.

  - a distance-removal model with a fixed-range abundance field from
    [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).

- coords:

  Numeric matrix or data frame with two finite projected coordinate
  columns (x, y), one row per fitting row in the original row order.
  Repeated locations repeat their coordinates. For `glmmTMB`, use the
  coordinates passed to
  [`glmmTMB::numFactor()`](https://rdrr.io/pkg/glmmTMB/man/numFactor.html),
  in the same units. For a distance-removal field, use exactly the
  `coords` passed to
  [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).
  Those are block locations, which can be coarser than the survey points
  used for the kernel covariates.

- crs:

  The projected CRS of `coords`, as an EPSG code (for example `32617`)
  or any input accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).
  It must have two horizontal axes and linear units that match `coords`.
  Geographic, geocentric, vertical and compound CRSs are rejected.
  Prediction rasters must use an equivalent CRS. Nothing is transformed
  or rescaled, so choose the projected CRS before fitting.

- batch_size:

  Positive integer, default 128. The largest number of raster cells
  passed to each native prediction call. Smaller batches use less memory
  but repeat the native setup more often.

- max_memory_mb:

  Positive number, default 512. A conservative budget, in MiB, for the
  dense arrays of each prediction batch. A batch whose estimate exceeds
  the budget stops before native prediction. Lower `batch_size` if that
  happens. This is a planning estimate, not an operating-system memory
  cap. For a distance-removal field, the budget set in
  [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)
  also applies.

## Value

A `multiScaleR_field_geometry` object, a small list that stores no
model, coordinates or predictions. It contains:

- `schema_version`, and for distance-removal fields `provider`
  (`"native_GDR_fixed_range"`), `range` and `range_estimated` (`FALSE`).

- `crs` (WKT text) and `units`, the linear units of the CRS.

- For `glmmTMB` fields, `position` and `group`, the names of the
  position and group columns, with the fitted `group_value` and
  `group_levels`.

- `training_locations` (distinct field locations) and `training_rows`.

- `provider_version`, the `glmmTMB` or `TMB` version.

- `batch_size` and `max_memory_mb`, the prediction controls.

- `reference_sha256`, `binding_sha256` and `definition_sha256`, SHA-256
  fingerprints of the coordinates, of the fitted context and of the
  descriptor itself. Prediction stops if any of them no longer matches,
  so recreate the descriptor after refitting or changing its inputs.

## Map a fitted field

Pass the descriptor to
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
through `raster_args = list(field_geometry = geometry)`. Prepare the
layers with
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
using `scale_center = TRUE` and `match_buffer = TRUE`. Each cell then
supplies its own projected coordinates to the field, and prediction runs
on one core, in batches of `batch_size` cells. No covariance among cells
is computed or stored. Several rows in `raster_args$const` give one band
(or fit and SE pair) per scenario row, in row order.

## glmmTMB fields

Every raster cell is a new field position, so set
`native_args = list(allow.new.levels = TRUE)` in
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md).
Without it, prediction stops before the first batch. All cells use the
single fitted field group, and a scenario that names another group
stops. Supply ordinary site and offset scenarios as usual. A position or
group scenario passed to
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
is replaced by the cell coordinates during prediction. Do not substitute
a constant coordinate. On the response scale, predictions are the
expected count for Poisson and NB2 models, the mean for Gaussian and
Gamma models, and the per-trial probability for binomial models,
including binomial models fitted to success and failure counts.

## Distance-removal fields

Use `type = "lambda"` for density. Each cell gets the field value
predicted at its center from the fitted block values, with the range
held fixed. `allow.new.levels` is not used. See
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)
for the full workflow, including scale-of-effect estimation with the
field.

## What the standard errors include

Standard errors are conditional on the selected kernel scales. They
include the fitted coefficient and field uncertainty, and for
distance-removal fields the extra field variance between blocks. They
add no observation noise and no kernel-scale uncertainty. `re.form = NA`
gives a zero-field reference, not an average over the field. Checking
the coordinates does not establish that the field is well estimated or
that its intervals are calibrated.

## See also

[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md),
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md),
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
[`prediction_contract()`](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md),
[`vignette("v11-spatial-fields", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v11-spatial-fields.md)

## Examples

``` r
# \donttest{
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  # Give the bundled landscape (in meters) a projected CRS
  land <- terra::rast(system.file("extdata", "landscape.tif",
                                  package = "multiScaleR"))[["land2"]]
  terra::crs(land) <- "EPSG:32617"
  data("surv_pts")
  pts <- terra::vect(surv_pts)
  terra::crs(pts) <- "EPSG:32617"
  pts <- pts[seq(1, 100, by = 2)]   # 50 survey points
  xy <- terra::crds(pts)

  # Counts depend on habitat at a 400 m scale plus a correlated field
  set.seed(3)
  truth <- kernel_scale.raster(land, sigma = 400, kernel = "gaussian",
                               verbose = FALSE)
  x2 <- as.numeric(scale(terra::extract(truth, pts, ID = FALSE)[, 1]))
  field <- drop(t(chol(0.25 * exp(-as.matrix(stats::dist(xy)) / 800))) %*%
                  stats::rnorm(nrow(xy)))
  inputs <- kernel_prep(pts = pts, raster_stack = land, max_D = 1200,
                        kernel = "gaussian", verbose = FALSE)
  df <- data.frame(counts = stats::rpois(nrow(xy), exp(1 + 0.6 * x2 + field)),
                   inputs$kernel_dat)
  df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])   # field positions (m)
  df$grp <- factor("all")                          # one shared field
  m0 <- glmmTMB::glmmTMB(counts ~ land2 + exp(pos + 0 | grp),
                         family = poisson(), data = df,
                         start = list(theta = c(log(0.5), log(500))))
  fit <- multiScale_optim(m0, inputs, n_cores = 1, verbose = FALSE)

  # Bind the fitting coordinates and their CRS to the fitted field
  geometry <- response_field_geometry(fit, coords = xy, crs = 32617)
  print(geometry$training_locations)   # 50 distinct field locations

  # Expected counts with conditional SEs for a 10 x 10 cell window
  layers <- kernel_scale.raster(land, multiScaleR = fit,
    scale_center = TRUE, match_buffer = TRUE, clamp = TRUE, verbose = FALSE,
    site_values = list(pos = as.character(df$pos[1]), grp = "all"))
  window <- terra::crop(layers, terra::ext(1500, 1700, 1500, 1700))
  attr(window, "multiScaleR_projection") <-
    attr(layers, "multiScaleR_projection")
  map <- predict(fit, window, se.fit = TRUE, uncertainty = "conditional",
                 native_args = list(allow.new.levels = TRUE),
                 raster_args = list(field_geometry = geometry))
  print(map)   # bands: fit (expected count) and se.fit
}
#> [1] 50
#> class       : SpatRaster
#> size        : 10, 10, 2  (nrow, ncol, nlyr)
#> resolution  : 20, 20  (x, y)
#> extent      : 1500, 1700, 1500, 1700  (xmin, xmax, ymin, ymax)
#> coord. ref. : WGS 84 / UTM zone 17N (EPSG:32617)
#> source(s)   : memory
#> names       :      fit,   se.fit
#> min values  : 2.103838, 0.379276
#> max values  :  3.97472, 0.585109
# }
```

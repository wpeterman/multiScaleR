# Group-specific raster preparation

Recorded September 23, 2026; updated September 25, 2026.

The public `kernel_prep_by_group()` function now implements this workflow. It
extracts each observation from its assigned map, creates pooled covariate
scaling parameters and distance bins, and returns a standard preparation object
for `multiScale_optim()`. See `vignettes/year-specific-rasters.Rmd` and
`inst/examples/year_specific_rasters.R` for the supported workflow. The manual
combination below records the approach used before the public function.

## Use case

Observations collected across years need landscape covariates from their own
year's raster, while one pooled model estimates shared landscape coefficients
and scales. A year term can allow baseline response differences. This also
applies to other observation-level groups that select different raster maps.

The function assigns each year-2 observation its year-2 raster, including nests
outside management polygons whose buffers may intersect changed cells.

## Binned data are compatible

The current binned representation stores one row per observation in each
covariate's `vsum` and `csum` matrices. These contain finite-cell value sums and
counts by distance bin. Both matrices must follow the observation's raster year.
Changing only `raw_cov` or `kernel_dat` leaves stale binned inputs in use.

For annual maps on an identical grid, prepared at ALL the same point locations
in the same order, with identical `max_D`, kernel settings, and `nbins`, distances
and bin definitions are shared. After verifying this compatibility, replace the
year-2 rows of BOTH binned matrices for every kernel covariate:

```r
# p2024 and p2025 were prepared with bin = TRUE at the same ordered points.
# The maps must have identical geometry and covariate definitions.
stopifnot(identical(p2024$binned$bin_rep, p2025$binned$bin_rep),
          identical(p2024$binned$point_ids, p2025$binned$point_ids),
          identical(p2024$binned$nbins, p2025$binned$nbins),
          identical(p2024$binned$covariates, p2025$binned$covariates),
          identical(p2024$binned$sources, p2025$binned$sources),
          identical(p2024$scale_vars, p2025$scale_vars),
          identical(p2024$unit_conv, p2025$unit_conv))
pooled <- p2024
for (nm in pooled$binned$covariates) {
  pooled$binned$vsum[[nm]][is2025, ] <-
    p2025$binned$vsum[[nm]][is2025, , drop = FALSE]
  pooled$binned$csum[[nm]][is2025, ] <-
    p2025$binned$csum[[nm]][is2025, , drop = FALSE]
}
if (isTRUE(pooled$cell_data_stored)) {
  pooled$raw_cov[is2025] <- p2025$raw_cov[is2025]
  pooled$d_list[is2025] <- p2025$d_list[is2025]
}
# Also rebuild kernel_dat and scl_params using pooled standardization,
# following the unscale/select/rescale steps in the worked example.
```

This can support lean objects (`store_cell_data = FALSE`) when every covariate
is kernel-type. Landscape and surface metrics still require their cell data;
binning only accelerates the kernel covariates in mixed models.

Do not assume that matching `nbins` and `max_D` alone makes separately prepared
objects compatible. The current builder derives bin edges from the maximum
extracted cell distance and representative distances from all prepared points.
Preparing only each year's subset of points can produce different bins even
on the same raster grid. The object does not currently retain the bin edges.
Equal representative distances alone are not a general proof of equal bins.

With cell data retained, a safer general strategy is to combine the selected
cell-level inputs first, then build one pooled set of bins. The current helper,
`.msr_build_kernel_bins()`, is internal, so this is a design direction rather
than a recommended public API. Lean objects with incompatible bin definitions
cannot simply be merged; their original inputs must be prepared again.

## Implemented interface and future extensions

Implemented interface: `kernel_prep_by_group(pts, raster_stacks, group, ...)`,
where `raster_stacks` is a named list of annual SpatRaster objects and `group`
maps each observation to a list name. It extracts only the assigned raster,
supports pooled bins and lean kernel-only objects, and preserves group labels.
Future extensions and checks include:

- Decide whether maps on different grids can support particular metrics and
  what resampling, if any, would preserve interpretation. The current function
  requires one grid and matching layer names.
- Retain explicit bin edges in future preparation objects if users need to
  merge objects prepared separately. The grouped function creates pooled bins
  directly and does not require such a merge.
- Add dedicated prediction helpers for assigning fitted annual surfaces and
  categorical year values. The vignette currently shows explicit scenarios.
- Extend tests to more than two years, varied map units and CRS definitions,
  and additional metric and response-model families.

## Further validation to consider

- Compare pooled binned covariates against an independently assembled annual
  reference at multiple scales and shapes, then compare model likelihoods and
  fitted scales. Quantify binning error against the exact cell-level path.
- Include annual NA changes: both finite-value sums and counts must change.
- Test interleaved years, shuffled IDs, repeated coordinates, missing groups,
  more than two years, and complete-case filtering.
- Cover full, lean, and mixed preparations; reject incompatible lean bins.
- Verify that unchanged annual maps reproduce the ordinary single-map result.
- Test cases where an annual subset has constant covariates but the pooled
  dataset has variation; annual standardization should not block preparation.
- Verify final refits, scale profiling, and predictions that select the correct
  annual map and use fitted pooled scaling parameters.

## Current validation

The original unbinned example passed 10 focused checks. In an additional
check on the current checkout, annual preparations with `bin = TRUE` and
`store_cell_data = FALSE` were combined as above. The combined binned object
matched bins rebuilt from the pooled cell-level inputs. The lean pooled fit
used all 200 observations and estimated sigma = 88.09247 meters, compared with
88.11296 meters for the unbinned example. These simulated results verify this
aligned-grid example, not all candidate workflows or binning accuracy generally.

The worked example now uses `kernel_prep_by_group()` and bins by default. Its
focused test compares assigned annual covariates with independent exact
extractions at three scales and checks the final pooled model refit.

The manual bin-row replacement below remains specific to compatible,
previously prepared objects; use the public function for new analyses.

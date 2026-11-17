# Example landscape raster

Three simulated landscape layers stored as `landscape.tif` in the
package's `extdata` folder, used with
[surv_pts](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md)
and
[landscape_counts](https://wpeterman.github.io/multiScaleR/reference/landscape_counts.md)
in the vignettes.

## Format

A `SpatRaster` with three layers of 200 by 200 cells of 20 map units (an
extent of 0 to 4000 in each direction) and no coordinate reference
system. Values range from 0 to 1:

- `land1`: binary surface with low spatial autocorrelation.

- `land2`: continuous surface with low spatial autocorrelation.

- `land3`: continuous surface with high spatial autocorrelation.

## Source

Simulated for this package.

## See also

[surv_pts](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md),
[landscape_counts](https://wpeterman.github.io/multiScaleR/reference/landscape_counts.md).

## Examples

``` r
land <- terra::rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))
land
#> class       : SpatRaster
#> size        : 200, 200, 3  (nrow, ncol, nlyr)
#> resolution  : 20, 20  (x, y)
#> extent      : 0, 4000, 0, 4000  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source      : landscape.tif
#> names       : land1, land2, land3
#> min values  :     0,     0,     0
#> max values  :     1,     1,     1
terra::plot(land)
```

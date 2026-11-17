# Example habitat raster

A binary habitat raster stored as `hab.tif` in the package's `extdata`
folder, used with
[pts](https://wpeterman.github.io/multiScaleR/reference/pts.md) and
[count_data](https://wpeterman.github.io/multiScaleR/reference/count_data.md).

## Format

A `SpatRaster` with one layer, `hab`, of 400 by 400 cells of 10 map
units (an extent of 0 to 4000 in each direction) and no coordinate
reference system. Cells are 1 for habitat and 0 otherwise.

## Source

Simulated for this package.

## See also

[pts](https://wpeterman.github.io/multiScaleR/reference/pts.md),
[count_data](https://wpeterman.github.io/multiScaleR/reference/count_data.md).

## Examples

``` r
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
hab
#> class       : SpatRaster
#> size        : 400, 400, 1  (nrow, ncol, nlyr)
#> resolution  : 10, 10  (x, y)
#> extent      : 0, 4000, 0, 4000  (xmin, xmax, ymin, ymax)
#> coord. ref. : 
#> source      : hab.tif
#> name        : hab
#> min value   :   0
#> max value   :   1
terra::plot(hab)
```

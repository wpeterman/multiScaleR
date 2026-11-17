# Survey points for the vignette example

Locations of the 100 survey points that go with
[landscape_counts](https://wpeterman.github.io/multiScaleR/reference/landscape_counts.md)
and the `landscape.tif` raster.

## Usage

``` r
data(surv_pts)
```

## Format

An `sf` point object with 100 features, no coordinate reference system,
and coordinates in the map units of `landscape.tif`. Its one attribute
column, `obs`, is a numeric column carried over from the simulation; the
vignettes use `landscape_counts$counts` as the response.

## Source

Simulated for this package.

## See also

[landscape_counts](https://wpeterman.github.io/multiScaleR/reference/landscape_counts.md),
[landscape](https://wpeterman.github.io/multiScaleR/reference/landscape.md),
[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md).

## Examples

``` r
data(surv_pts)
land <- terra::rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))
terra::plot(land$land1)
terra::points(terra::vect(surv_pts), pch = 19)
```

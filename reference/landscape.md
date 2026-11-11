# Simulated raster

Raster data for use with vignette example

## Format

\`landscape_rast\`

A SpatRaster object with three surfaces:

- land1:

  –\> A binary landscape surface with low autocorrelation

- land2:

  –\> A continuous landscape surface with low autocorrelation

- land3:

  –\> A continuous landscape surface with high autocorrelation

## Examples

``` r
land_rast <- terra::rast(system.file("extdata",
                         "landscape.tif", package = 'multiScaleR'))
```

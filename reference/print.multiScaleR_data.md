# Print prepared multiScaleR inputs

Prints a short report of a
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
object: the number of points and covariates, the kernel, the binning,
and the distance bounds.

## Usage

``` r
# S3 method for class 'multiScaleR_data'
print(x, ...)
```

## Arguments

- x:

  A `multiScaleR_data` object.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## See also

[multiScaleR_data](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_data-class.md),
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

## Examples

``` r
data(pts)
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian", verbose = FALSE)
inputs
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
```

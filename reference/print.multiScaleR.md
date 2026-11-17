# Print a multiScaleR fit

Prints a fit from
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md):
the kernel, the estimated scales with their standard errors, the scale
distances, the fitted model, any warnings, and the follow-up
recommendation. Use
[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md)
for intervals.

## Usage

``` r
# S3 method for class 'multiScaleR'
print(x, dimension = c("2d", "1d"), prob = 0.9, ...)
```

## Arguments

- x:

  A `multiScaleR` object.

- dimension:

  Character. `"2d"` (default) prints the ideal planar radius first;
  `"1d"` prints the original line-based distance first.

- prob:

  Numeric between 0 and 1 (exclusive). Proportion of cumulative kernel
  weight used for the distances. Default: `0.9`.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## See also

[multiScaleR](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md),
[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md).

## Examples

``` r
data(pts)
data(count_data)
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian", verbose = FALSE)
fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
                        inputs, verbose = FALSE)
fit                      # same as print(fit)
#> 
#> Call:
#> multiScale_optim(fitted_mod = glm(y ~ hab, family = poisson, 
#>     data = count_data), kernel_inputs = inputs, verbose = FALSE)
#> 
#> 
#> Kernel used: 
#> gaussian
#> 
#> ***** Optimized Scale of Effect *****
#> 
#>         Mean       SE
#> hab 76.16637 9.567451
#> 
#> 
#> ***** Optimized Scale of Effect -- Distance *****
#> Scale distances (map units) at 90%, two definitions.
#> 
#> 2D, radius of the circle for the ideal untruncated kernel (selected):
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37
#> 
#> 1D, along a line through the point (comparison):
#>       Mean  2.5%  97.5%
#> hab 125.28 93.92 156.65
#> 
#>   ==================================== 
#> 
#>  ***** Fitted Model *****
#> 
#> Call:  glm(formula = y ~ hab, family = poisson, data = data)
#> 
#> Coefficients:
#> (Intercept)          hab  
#>      0.4807       0.8691  
#> 
#> Degrees of Freedom: 74 Total (i.e. Null);  73 Residual
#> Null Deviance:       186.6 
#> Residual Deviance: 70.18     AIC: 232.3
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent. 
#>  Consider increasing max_D in `kernel_prep` to >= 250.56 to ensure accurate estimation of scale.
#> 
#> 
#> 
#> ***** Follow-up recommendation *****
#> Recommended `max_D`: 250.56
#> Recommended starting `sigma` values:
#>      hab 
#> 76.16637 
#> 
#> Action:
#> Refit from the optimized parameters and address the flagged diagnostics. 
print(fit, prob = 0.5)   # distances holding half of each kernel's weight
#> 
#> Call:
#> multiScale_optim(fitted_mod = glm(y ~ hab, family = poisson, 
#>     data = count_data), kernel_inputs = inputs, verbose = FALSE)
#> 
#> 
#> Kernel used: 
#> gaussian
#> 
#> ***** Optimized Scale of Effect *****
#> 
#>         Mean       SE
#> hab 76.16637 9.567451
#> 
#> 
#> ***** Optimized Scale of Effect -- Distance *****
#> Scale distances (map units) at 50%, two definitions.
#> 
#> 2D, radius of the circle for the ideal untruncated kernel (selected):
#>      Mean  2.5%  97.5%
#> hab 89.68 67.23 112.13
#> 
#> 1D, along a line through the point (comparison):
#>      Mean  2.5% 97.5%
#> hab 51.37 38.51 64.23
#> 
#>   ==================================== 
#> 
#>  ***** Fitted Model *****
#> 
#> Call:  glm(formula = y ~ hab, family = poisson, data = data)
#> 
#> Coefficients:
#> (Intercept)          hab  
#>      0.4807       0.8691  
#> 
#> Degrees of Freedom: 74 Total (i.e. Null);  73 Residual
#> Null Deviance:       186.6 
#> Residual Deviance: 70.18     AIC: 232.3
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent. 
#>  Consider increasing max_D in `kernel_prep` to >= 250.56 to ensure accurate estimation of scale.
#> 
#> 
#> 
#> ***** Follow-up recommendation *****
#> Recommended `max_D`: 250.56
#> Recommended starting `sigma` values:
#>      hab 
#> 76.16637 
#> 
#> Action:
#> Refit from the optimized parameters and address the flagged diagnostics. 
```

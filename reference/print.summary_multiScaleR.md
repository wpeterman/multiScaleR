# Print a summary of a multiScaleR fit

Prints the object returned by
[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md):
the estimated scales with their intervals, the scale distances under
both definitions, the fitted model, any warnings, and the follow-up
recommendation.

## Usage

``` r
# S3 method for class 'summary_multiScaleR'
print(x, dimension = NULL, ...)
```

## Arguments

- x:

  A `summary_multiScaleR` object from
  [`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md).

- dimension:

  Character, `"2d"` or `"1d"`. Selects which distance table is printed
  first. The default, `NULL`, uses the definition stored in `x`;
  summaries saved before this choice existed print the 1D table first.

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## See also

[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md),
[`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
for the two distance definitions.

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
s <- summary(fit)
print(s, dimension = "1d")   # line-based distances first
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
#>         Mean       SE     2.5%    97.5%
#> hab 76.16637 9.567451 57.09847 95.23428
#> 
#> 
#>   ==================================== 
#> 
#> ***** Optimized Scale of Effect -- Distance *****
#> Scale distances (map units) at 90%, two definitions.
#> 
#> 1D, along a line through the point (selected):
#>       Mean  2.5%  97.5%
#> hab 125.28 93.92 156.65
#> 
#> 2D, radius of the circle for the ideal untruncated kernel (comparison):
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37
#> 
#> 
#>   ==================================== 
#> 
#>  *****     Fitted Model Summary     *****
#> 
#> 
#> Call:
#> glm(formula = y ~ hab, family = poisson, data = data)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   0.4807     0.1021   4.707 2.51e-06 ***
#> hab           0.8691     0.0855  10.166  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 186.585  on 74  degrees of freedom
#> Residual deviance:  70.185  on 73  degrees of freedom
#> AIC: 232.29
#> 
#> Number of Fisher Scoring iterations: 5
#> 
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

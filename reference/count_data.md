# Simulated counts for the hab example

Simulated counts at the 75 points in
[pts](https://wpeterman.github.io/multiScaleR/reference/pts.md), for
learning to estimate a scale of effect with the `hab.tif` raster.

## Usage

``` r
data(count_data)
```

## Format

A data frame with 75 rows, one per point in
[pts](https://wpeterman.github.io/multiScaleR/reference/pts.md), and 2
columns:

- `y`: integer count at each point.

- `hab`: kernel-weighted mean of the binary `hab.tif` raster around each
  point, centered and scaled to unit variance.

## Source

Simulated for this package.

## Details

Counts were simulated from a Poisson distribution with an intercept of
0.5, a `hab` coefficient of 0.75 (log scale), and a Gaussian scale of
effect (sigma) of 75 map units. A fit that recovers a sigma near 75 and
a coefficient near 0.75 is working as intended.

## See also

[pts](https://wpeterman.github.io/multiScaleR/reference/pts.md),
[hab](https://wpeterman.github.io/multiScaleR/reference/hab.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

## Examples

``` r
data(count_data)
head(count_data)
#>   y        hab
#> 1 1 -0.7550406
#> 2 6  1.5678284
#> 3 0  0.3289698
#> 4 0 -0.6590483
#> 5 0 -0.5295792
#> 6 4  1.1883843
# Starting model: hab will be replaced by values at each candidate scale
glm(y ~ hab, family = poisson, data = count_data)
#> 
#> Call:  glm(formula = y ~ hab, family = poisson, data = count_data)
#> 
#> Coefficients:
#> (Intercept)          hab  
#>      0.4802       0.8705  
#> 
#> Degrees of Freedom: 74 Total (i.e. Null);  73 Residual
#> Null Deviance:       186.6 
#> Residual Deviance: 70.24     AIC: 232.3
```

# Simulated counts for the vignette example

Simulated counts at the 100 points in
[surv_pts](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md),
used as the running example in the vignettes with the `landscape.tif`
raster.

## Usage

``` r
data(landscape_counts)
```

## Format

A data frame with 100 rows, one per point in
[surv_pts](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md),
and 2 columns:

- `counts`: integer count at each point.

- `site`: a habitat variable measured at the point itself, which has no
  scale of effect.

## Source

Simulated for this package.

## Details

Counts were simulated from a Poisson distribution with an intercept of
0.25 and log-scale coefficients of -0.5 for `land1`, 0.3 for `site`, and
0.7 for `land2`. The true Gaussian scales of effect (sigma) are 250 map
units for `land1` and 500 for `land2`. The `land1` and `land2`
covariates are not stored here; build them from
[landscape](https://wpeterman.github.io/multiScaleR/reference/landscape.md)
with
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

## See also

[surv_pts](https://wpeterman.github.io/multiScaleR/reference/surv_pts.md),
[landscape](https://wpeterman.github.io/multiScaleR/reference/landscape.md),
[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md).

## Examples

``` r
data(landscape_counts)
summary(landscape_counts)
#>      counts        site        
#>  Min.   : 0   Min.   :-1.7129  
#>  1st Qu.: 0   1st Qu.:-0.6480  
#>  Median : 1   Median :-0.1697  
#>  Mean   : 2   Mean   : 0.0000  
#>  3rd Qu.: 3   3rd Qu.: 0.5189  
#>  Max.   :10   Max.   : 2.9984  
```

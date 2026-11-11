# CI Function

Internal function for calculating confidence intervals

## Usage

``` r
ci_func_r(x, df, min_D = NULL, names = NULL, as_dist = FALSE, ...)
```

## Arguments

- x:

  Estimates

- df:

  Degrees of Freedom

- min_D:

  minimum distance, Default: NULL

- names:

  Names of fitted variables, Default: NULL

- as_dist:

  Calculate distance confidence interval; Default = FALSE

- ...:

  Not used

## Value

data frame with calculated confidence intervals

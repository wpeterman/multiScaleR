# Scale Function

Scaling function to be applied to rasters

## Usage

``` r
scale_type(
  d,
  kernel = c("gaussian", "exp", "expow", "fixed"),
  sigma,
  shape = NULL,
  r_stack.df = NULL,
  output = NULL
)
```

## Arguments

- d:

  Vector of distances

- kernel:

  Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow';
  Default: 'gaussian')

- sigma:

  Scaling parameter

- shape:

  Shape parameter if using exponential power kernel

- r_stack.df:

  Dataframe values extracted from rasters

- output:

  If NULL, a vector of weights is returned, otherwise a weighted raster
  values are returned, Default: NULL

## Value

A vector of weights or vector of weighted raster values

## Details

DETAILS

## Examples

``` r
### TO BE COMPLETED ###
```

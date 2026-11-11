# Plot Sigma Profile

Plots the profiled log-likelihood or AICc across sigma values for each
spatial covariate. The optimized sigma is marked with a vertical red
line.

## Usage

``` r
# S3 method for class 'sigma_profile'
plot(x, ...)
```

## Arguments

- x:

  An object of class `sigma_profile` from
  [`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md).

- ...:

  Additional arguments (not currently used).

## Value

A list of `ggplot` objects (one per covariate), returned invisibly.
Plots are printed as a side effect.

## See also

[`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)

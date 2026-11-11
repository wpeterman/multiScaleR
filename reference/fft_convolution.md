# Fast Fourier Transformation

Kernel smoothing with fft

## Usage

``` r
fft_convolution(x, kernel, fun = "mean", na.rm = TRUE)
```

## Arguments

- x:

  A matrix of the raster layer

- kernel:

  A weight matrix for the smoothing

- fun:

  Either "mean" or "sum"

- na.rm:

  Logical; should NA values be removed from mean calculation?

## Value

A matrix

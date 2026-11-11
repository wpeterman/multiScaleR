# Print method for multiScaleR

Print method for objects of class `multiScaleR`.

## Usage

``` r
# S3 method for class 'multiScaleR'
print(x, dimension = c("2d", "1d"), prob = 0.9, ...)
```

## Arguments

- x:

  A `multiScaleR` object

- dimension:

  Character. `"2d"` (default) prints the ideal planar radius first;
  `"1d"` prints the original line-based distance first.

- prob:

  Numeric between 0 and 1 (exclusive). Cumulative kernel-weight
  proportion used for distances. Default: `0.9`.

- ...:

  Ignored

## Value

Invisibly returns the input `multiScaleR` object

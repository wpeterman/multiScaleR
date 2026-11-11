# Print method for summary_multiScaleR

Print method for objects of class `summary_multiScaleR`.

## Usage

``` r
# S3 method for class 'summary_multiScaleR'
print(x, dimension = NULL, ...)
```

## Arguments

- x:

  A `summary_multiScaleR` object

- dimension:

  Character. `"2d"` or `"1d"`. Selects which distance table is printed
  first. By default, use the definition stored in `x`; summaries saved
  before this switch existed retain 1D reporting.

- ...:

  Ignored

## Value

Invisibly returns the input `summary_multiScaleR` object

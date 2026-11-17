# Create Grouped Cross-Validation Folds

Keep all rows from a sampling group in the same assessment fold, so that
a held-out row is never predicted with help from its own group. Whole
groups are assigned to folds to balance the row counts, largest groups
first.

## Usage

``` r
grouped_cv_folds(group, k = 5, ids = NULL)

# S3 method for class 'multiScaleR_folds'
print(x, ...)
```

## Arguments

- group:

  Vector of sampling-group labels, one per row, such as site, animal,
  year, or choice-set IDs. Missing or empty labels are rejected.

- k:

  Number of folds, at least two and no greater than the number of
  groups. The default is five.

- ids:

  Unique character row IDs, in the order of `group`. Use the row names
  of the modeling data. The default uses the names of `group`, or the
  row numbers as character strings when `group` is unnamed. Named groups
  must align exactly with supplied IDs.

- x:

  A `multiScaleR_folds` object, for
  [`print()`](https://rdrr.io/r/base/print.html).

- ...:

  Not used. Reserved for future print controls.

## Value

A named list of class `multiScaleR_folds` with one element per fold
(`fold1`, `fold2`, and so on), ready for
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
and the other cross-validation functions. Each element holds two
character vectors of row IDs, `analysis` (rows used to fit) and
`assessment` (rows held out). Every row is assessed exactly once. The
attribute `fold_info` holds the `method` (`"grouped"`), `k`, an
`assignment` data frame with the `row_id`, `group`, `unit`, and
assessment `fold` of every row, and `excluded`, the IDs left out of each
fold's training rows (none for grouped folds).
[`print()`](https://rdrr.io/r/base/print.html) shows the number of
analysis, assessment, and excluded rows in each fold, and returns `x`
invisibly.

## Details

Assignment is deterministic and does not use the random-number
generator, so the same call always gives the same folds. Groups of equal
size are ordered by their labels. Balancing row counts does not stratify
the responses. Choose the sampling groups and the prediction question
before you look at any scores. To keep groups together within spatial
blocks, pass them as `group` to
[`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md).
Building folds does not decide which models can be scored. For matched
choice sets, score the folds with
[`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md).
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
compares grouped, spatial, and random fold designs.

## See also

[`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md),
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)

## Examples

``` r
# Twelve records from six sites, two visits each
ids <- paste0("row", 1:12)
folds <- grouped_cv_folds(rep(letters[1:6], each = 2), k = 3, ids = ids)
folds                    # rows fitted and held out in each fold
#> Cross-validation folds: grouped 
#>   fold analysis assessment excluded
#>  fold1        8          4        0
#>  fold2        8          4        0
#>  fold3        8          4        0
folds$fold1$assessment   # both visits to a site are held out together
#> [1] "row1" "row2" "row7" "row8"
```

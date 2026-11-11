# Kernel scaling function

Function for internal use with optim

## Usage

``` r
kernel_scale_fn(
  par,
  d_list,
  cov_df,
  kernel,
  fitted_mod,
  join_by = NULL,
  mod_return = NULL,
  opt_context = NULL,
  cov_w = NULL
)
```

## Arguments

- par:

  list of parameters

- d_list:

  List of distance vectors

- cov_df:

  List of data frames with values extracted from rasters

- kernel:

  Kernel used

- fitted_mod:

  fitted model object

- join_by:

  Optional data frame of site keys used to map prepared spatial
  covariates to the unmarked frame without changing its row order.

- mod_return:

  Default: NULL

- opt_context:

  Cached optimization context created internally by
  \`multiScale_optim()\`

- cov_w:

  Optional precomputed (unscaled) kernel-weighted covariate matrix. When
  supplied, the per-covariate weighting step is skipped and this matrix
  is scaled and refit directly. Used by \`profile_sigma()\` to reuse the
  columns of covariates that are held at their optimum. Default: NULL
  (compute internally).

## Value

Either estimated parameters or the fitted model using provided
parameters

## Details

For internal use

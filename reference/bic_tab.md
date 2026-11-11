# multiScaleR model selection

Creates a model selection table ranked by BIC for a list of fitted
models. Mixed lists containing both `multiScaleR` objects and plain
model objects are supported.

## Usage

``` r
bic_tab(mod_list,
        mod_names = NULL,
        verbose = FALSE,
        ...)
```

## Arguments

- mod_list:

  List containing fitted `multiScaleR` objects (from
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md))
  and/or plain fitted model objects. All models must have been fit to
  the same number of observations.

- mod_names:

  Optional character vector of model names. Length must equal
  `length(mod_list)`. By default, names are constructed as
  `[kernel]right-hand-side-formula` for `multiScaleR` objects.

- verbose:

  Logical. If `TRUE`, the table is also printed to the console. Default:
  `FALSE`.

- ...:

  Additional arguments (not used).

## Value

A data frame of class `"bictab"` (from the `AICcmodavg` package) sorted
by ascending BIC, containing:

- `Modnames`:

  Model name.

- `K`:

  Number of fitted likelihood parameters, including ancillary parameters
  such as negative-binomial dispersion, plus optimized scale and shape
  parameters for `multiScaleR` models.

- `BIC`:

  Bayesian information criterion value.

- `Delta_BIC`:

  Difference from the top model.

- `BICWt`:

  BIC weight.

- `Cum.Wt`:

  Cumulative BIC weight.

- `LL`:

  Log-likelihood.

## Details

The table is constructed using
[`bictabCustom`](https://rdrr.io/pkg/AICcmodavg/man/bictabCustom.html)
from the `AICcmodavg` package. BIC penalizes model complexity more
heavily than AIC as sample size grows, and is often preferred when the
goal is identifying the single best-supported model rather than model
averaging. Sigma (and shape for `"expow"`) is counted in K. Fitted
parameter counts and log-likelihood validation follow
[`aic_tab`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
including all estimated components of unmarked models.

## Author

Bill Peterman

## Examples

``` r
## Simulate data
set.seed(555)

points <- vect(cbind(c(5,7,9,11,13),
                     c(13,11,9,7,5)))

mat_list <- list(r1 = rast(matrix(rnorm(20^2),
                                  nrow = 20)),
                 r2 = rast(matrix(rnorm(20^2),
                                  nrow = 20)))
rast_stack <- rast(mat_list)
kernel_inputs <- kernel_prep(pts = points,
                             raster_stack = rast_stack,
                             max_D = 5,
                             kernel = 'gaussian',
                             sigma = NULL)
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...
## Example response data
y <- rnorm(5)

## Create data frame with raster variables
dat <- data.frame(y = y,
                  kernel_inputs$kernel_dat)
mod1 <- glm(y ~ r1,
            data = dat)
mod2 <- glm(y ~ r2,
            data = dat)
mod3 <- glm(y ~ r1 + r2,
            data = dat)

## NOTE: This code is only for demonstration
## Optimization results will have no meaning

opt_mod1 <- multiScale_optim(fitted_mod = mod1,
                             kernel_inputs = kernel_inputs,
                             par = NULL,
                             n_cores = NULL)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=11.9 to ensure accurate estimation of scale.
#> 
#> 
#>  WARNING!!!
#>  The standard error of one or more `sigma` estimates is >= 50% of the estimated mean value.
#>  Carefully assess whether or not this variable is meaningful in your analysis and interpret with caution.
#> 
#> 

opt_mod2 <- multiScale_optim(fitted_mod = mod2,
                             kernel_inputs = kernel_inputs,
                             par = NULL,
                             n_cores = NULL)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=16.44 to ensure accurate estimation of scale.
#> 
#> 
#>  WARNING!!!
#>  The standard error of one or more `sigma` estimates is >= 50% of the estimated mean value.
#>  Carefully assess whether or not this variable is meaningful in your analysis and interpret with caution.
#> 
#> 
opt_mod3 <- multiScale_optim(fitted_mod = mod3,
                             kernel_inputs = kernel_inputs,
                             par = NULL,
                             n_cores = NULL)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=16.44 to ensure accurate estimation of scale.
#> 
#> 

 ## BIC table
 mod_list <- list(opt_mod1, opt_mod2, opt_mod3)

 bic_tab(mod_list = mod_list)
#> 
#> Model selection based on BIC:
#> 
#>                     K   BIC Delta_BIC BICWt Cum.Wt    LL
#> [gaussian]~ r1 + r2 6  3.26      0.00  0.98   0.98  3.20
#> [gaussian]~ r2      4 11.56      8.30  0.02   0.99 -2.56
#> [gaussian]~ r1      4 13.41     10.15  0.01   1.00 -3.49
#> 

 ## BIC table with specified names
 bic_tab(mod_list = mod_list,
         mod_names = c('mod1', 'mod2', 'mod3'))
#> 
#> Model selection based on BIC:
#> 
#>      K   BIC Delta_BIC BICWt Cum.Wt    LL
#> mod3 6  3.26      0.00  0.98   0.98  3.20
#> mod2 4 11.56      8.30  0.02   0.99 -2.56
#> mod1 4 13.41     10.15  0.01   1.00 -3.49
#> 
```

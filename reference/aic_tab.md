# multiScaleR model selection

Creates a model selection table ranked by AIC or AICc for a list of
fitted models. Mixed lists containing both `multiScaleR` objects and
plain model objects (e.g., `glm`) are supported, allowing comparison of
scale-optimized models against fixed-scale alternatives.

## Usage

``` r
aic_tab(mod_list,
        AICc = TRUE,
        mod_names = NULL,
        verbose = FALSE,
        ...)
```

## Arguments

- mod_list:

  List containing fitted `multiScaleR` objects (from
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md))
  and/or plain fitted model objects (e.g., `glm`, `lm`). All models must
  have been fit to the same set of observations; an error is raised if
  sample sizes differ.

- AICc:

  Logical. If `TRUE` (default), the second-order corrected AIC (AICc) is
  used. Set to `FALSE` for standard AIC. AICc is recommended when the
  ratio of observations to parameters is small (n / K \< 40 is a common
  rule of thumb).

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

A data frame of class `"aictab"` (from the `AICcmodavg` package) sorted
by ascending AIC(c), containing:

- `Modnames`:

  Model name (from `mod_names` or auto-generated).

- `K`:

  Number of fitted likelihood parameters, including ancillary parameters
  such as negative-binomial dispersion, plus optimized scale and shape
  parameters for `multiScaleR` models.

- `AICc` or `AIC`:

  Information criterion value.

- `Delta_AICc` or `Delta_AIC`:

  Difference from the top model.

- `ModelLik`:

  Relative likelihood (`exp(-0.5 * Delta)`).

- `AICcWt` or `AICWt`:

  Akaike weight.

- `LL`:

  Log-likelihood.

- `Cum.Wt`:

  Cumulative Akaike weight.

## Details

The table is constructed using
[`aictabCustom`](https://rdrr.io/pkg/AICcmodavg/man/aictabCustom.html)
from the `AICcmodavg` package. For `multiScaleR` objects, the optimized
(refitted) model stored in `opt_mod` is used for log-likelihood and
parameter count. The sigma parameter (and shape for `"expow"`) is added
to K automatically. For plain model objects in a mixed list, sigma is
not added to K. The fitted model's log-likelihood degrees of freedom are
used when available, so negative-binomial dispersion contributes one
parameter. For `unmarkedFit` models, native AIC and log-likelihood
metadata determine the count, including detection coefficients, factor
contrasts, mixture parameters, and estimated random-effect variances.
Optimizer and TMB free-parameter metadata are checked for agreement when
present. Latent random effects and the integration limit `K` are not
counted. Missing or inconsistent parameter metadata and non-finite
log-likelihoods raise an error instead of guessing from predictor names.
AICc tables raise an error if any model has `n <= K + 1`; use
`AICc = FALSE` for an AIC table in that case.

## Author

Bill Peterman

## Examples

``` r
## Bundled counts simulated with a habitat scale of 75 map units
data("pts")
data("count_data")
hab <- terra::rast(system.file("extdata", "hab.tif", package = "multiScaleR"))
prepared <- kernel_prep(pts, hab, max_D = 300, sigma = 75, verbose = FALSE)
dat <- data.frame(y = count_data$y, prepared$kernel_dat)
fixed_model <- glm(y ~ hab, family = poisson(), data = dat)
null_model <- glm(y ~ 1, family = poisson(), data = dat)
spatial_model <- multiScale_optim(fixed_model, prepared, par = 75 / 300,
                                  n_cores = 1, verbose = FALSE)
#> 
#> 
#> Optimization complete
models <- list(spatial = spatial_model, fixed = fixed_model, null = null_model)

## The spatial model pays for one extra estimated sigma parameter.
aic_tab(models, AICc = FALSE, mod_names = names(models))
#> 
#> Model selection based on AIC:
#> 
#>         K    AIC Delta_AIC AICWt Cum.Wt      LL
#> fixed   2 232.34      0.00  0.73   0.73 -114.17
#> spatial 3 234.33      1.99  0.27   1.00 -114.17
#> null    1 346.69    114.35  0.00   1.00 -172.34
#> 
aic_tab(models, AICc = TRUE, mod_names = names(models))
#> 
#> Model selection based on AICc:
#> 
#>         K   AICc Delta_AICc AICcWt Cum.Wt      LL
#> fixed   2 232.51       0.00   0.75   0.75 -114.17
#> spatial 3 234.67       2.16   0.25   1.00 -114.17
#> null    1 346.74     114.23   0.00   1.00 -172.34
#> 
```

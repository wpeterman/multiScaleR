# Check Optional Software for Bayesian Scale Models

Check whether the R packages that a
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
computation method needs can be loaded and meet the package's minimum
versions. Nothing is installed, no model is fitted, and neither the
compiler nor sampling performance is tested.

## Usage

``` r
check_bayes_setup(
  method = c("sample", "laplace", "marginal", "prepare", "slice"),
  correction = c("none", "importance"),
  cores = 1L
)

# S3 method for class 'multiScaleR_bayes_setup'
print(x, ...)
```

## Arguments

- method:

  Computation method to check, as in
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).
  `"sample"` (default, NUTS) needs RTMB, tmbstan, rstan, and posterior.
  `"laplace"` and `"marginal"` need RTMB. `"slice"` and `"prepare"` need
  no optional packages to run in one process. Slice chain diagnostics
  need posterior (version 1.6.0 or later).

- correction:

  `"none"` (default) or `"importance"`. Importance correction applies
  only to `method = "marginal"` and also needs posterior.

- cores:

  Requested number of worker processes; default 1. A value greater than
  one also checks processx, which parallel fits need. It has no effect
  for `"prepare"`, which always uses one process.

- x:

  A `multiScaleR_bayes_setup` object from `check_bayes_setup()`.

- ...:

  Not used. Supplying any further argument is an error.

## Value

An object of class `multiScaleR_bayes_setup`, a list with these
elements.

- `method`, `correction`: the settings checked.

- `ready`: `TRUE` when every listed package loads and meets its minimum
  version.

- `packages`: a data frame with one row per package and columns
  `package`, `minimum` (the minimum version, or `NA` when any version
  works), installed `version`, `loadable`, `meets_version`, `passed`,
  and a diagnostic `message` that says what to install or update.

- `notes`: what the check does and does not cover.

Printing the object lists the packages and the notes, and returns it
invisibly. A passed check does not certify a working compiler,
successful sampling, convergence, or adequate memory. Fit a small model
next.

## Details

Preparation (`method = "prepare"`) validates your data and formulas
before you install any sampling software. Installing the NUTS packages
from source can require a C++ compiler, which on Windows comes with
Rtools.
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
covers installation, and
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
compares the methods.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)

## Examples

``` r
check_bayes_setup(method = "prepare")  # always ready: nothing to install
#> Bayesian package setup: prepare 
#> Package checks passed.
#> Preparation needs no optional Bayesian packages.
#> Checks R packages only; does not install software or run a model.
#> Compiler configuration, memory and sampling diagnostics require a separate small fit. 

# Packages for slice sampling on two workers. `ready` is FALSE when
# processx is missing. The message column says what to install.
setup <- check_bayes_setup(method = "slice", cores = 2)
setup$ready
#> [1] TRUE
setup$packages[c("package", "loadable", "passed")]
#>    package loadable passed
#> 1 processx     TRUE   TRUE
```

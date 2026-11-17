# multiScaleR

README last updated: 2026-11-14

| CI | Development | CRAN | License |
|----|----|----|----|
| [![R-CMD-check](https://github.com/wpeterman/multiScaleR/actions/workflows/GitHub.yaml/badge.svg)](https://github.com/wpeterman/multiScaleR/actions/workflows/GitHub.yaml) | [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable) | [![CRAN status](https://www.r-pkg.org/badges/version/multiScaleR)](https://cran.r-project.org/package=multiScaleR) | [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) |
| [![codecov](https://codecov.io/gh/wpeterman/multiScaleR/graph/badge.svg)](https://app.codecov.io/gh/wpeterman/multiScaleR) | [![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) | [![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/multiScaleR)](https://cran.r-project.org/package=multiScaleR) | [![DOI](https://img.shields.io/badge/DOI-10.1007/s10980--025--02267--x-yellow.svg)](https://doi.org/10.1007/s10980-025-02267-x) |
| [![r-universe checks](https://wpeterman.r-universe.dev/multiScaleR/badges/checks)](https://wpeterman.r-universe.dev/multiScaleR) | [![r-universe version](https://wpeterman.r-universe.dev/multiScaleR/badges/version)](https://wpeterman.r-universe.dev/multiScaleR) |  |  |

## Overview

The *scale of effect* is the distance over which a landscape variable
influences an ecological response. multiScaleR estimates it directly. It
weights the raster cells around each survey point with a distance kernel
and chooses the kernel scale (sigma) that maximizes the likelihood of
your regression model. You fit your model once at starting values, and
multiScaleR searches over scales, reports each one with its uncertainty,
and helps you compare, check, and map the result.

multiScaleR works with:

- linear, generalized linear, additive, and mixed models;
- occupancy and abundance models from `unmarked`, including
  distance-removal models;
- step-selection models fitted with
  [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html)
  or `amt`;
- Bayesian models in which a scale of effect can vary among years or
  with a covariate.

It also computes kernel-weighted landscape and surface metrics,
cross-validates models with spatial or grouped folds, maps predictions
across a landscape, and simulates data with known scales of effect so
you can test a design before collecting data.

## Installation

Install the CRAN release with:

``` r

install.packages("multiScaleR")
```

Install the latest build from R-universe with:

``` r

install.packages(
  "multiScaleR",
  repos = c("https://wpeterman.r-universe.dev", "https://cloud.r-project.org")
)
```

R-universe rebuilds multiScaleR from the `master` branch on GitHub
whenever that branch changes. This version can therefore be newer than
the CRAN release. R-universe provides precompiled Windows and macOS
binaries with the vignettes already built, so this route does not
require Rtools.

If installing from GitHub on Windows, you may need Rtools to compile the
package from source. Choose the Rtools release that matches your version
of R and follow the installation instructions at
<https://cran.r-project.org/bin/windows/Rtools/>.

Then execute the following commands in R:

``` r

# Install 'remotes' package, if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

remotes::install_github("wpeterman/multiScaleR", 
                        build_vignettes = TRUE) # Install & build vignette

library(multiScaleR) # Load the package
```

## A first example

The package includes simulated counts at 75 points, generated with a
habitat effect whose true scale of effect (sigma) is 75 map units.

``` r

library(multiScaleR)
data(pts)          # survey points
data(count_data)   # counts at each point
hab <- terra::rast(system.file("extdata", "hab.tif", package = "multiScaleR"))

# 1. Weight the habitat raster around each point with a Gaussian kernel
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian")

# 2. Fit your model once, at the starting scale
start_model <- glm(y ~ hab, family = poisson, data = count_data)

# 3. Estimate the scale of effect
fit <- multiScale_optim(start_model, inputs)
summary(fit)   # sigma, its interval, and the distance holding 90% of the weight
plot(fit)      # the fitted kernel weight curve
```

The estimated sigma is close to the simulated value of 75, and
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
converts it into the distance that contains 90% of the kernel’s weight,
which is the number to report as the scale of effect.

## Learn more

The vignettes work through simulated data and explain how to read each
result. Start with the quick start and the user guide, then turn to the
topics you need. Open one in R with
[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md),
or read them all on the [package
website](https://wpeterman.github.io/multiScaleR/articles/).

| Vignette | Name in [`vignette()`](https://rdrr.io/r/utils/vignette.html) |
|:---|:---|
| [01. Quick start](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.html) | `v01-quickstart` |
| [02. User guide](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.html) | `v02-user-guide` |
| [03. Mapping predictions and limiting extrapolation](https://wpeterman.github.io/multiScaleR/articles/v03-map-predictions.html) | `v03-map-predictions` |
| [04. Landscape metrics as covariates](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.html) | `v04-landscape-metrics` |
| [05. Surface texture metrics as covariates](https://wpeterman.github.io/multiScaleR/articles/v05-surface-metrics.html) | `v05-surface-metrics` |
| [06. Using year-specific rasters in one model](https://wpeterman.github.io/multiScaleR/articles/v06-year-specific-rasters.html) | `v06-year-specific-rasters` |
| [07. Choose a fitting method for your model](https://wpeterman.github.io/multiScaleR/articles/v07-fitting-methods.html) | `v07-fitting-methods` |
| [08. Cross-validate scale-of-effect models](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.html) | `v08-cross-validation` |
| [09. What your estimates and intervals mean](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.html) | `v09-estimates-and-intervals` |
| [10. Estimate scales of effect in step-selection models](https://wpeterman.github.io/multiScaleR/articles/v10-step-selection.html) | `v10-step-selection` |
| [11. Account for spatial autocorrelation with a spatial field](https://wpeterman.github.io/multiScaleR/articles/v11-spatial-fields.html) | `v11-spatial-fields` |
| [12. Fit distance-removal models with scales of effect](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.html) | `v12-distance-removal` |
| [13. Estimate scales of effect with a Bayesian model](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.html) | `v13-bayesian-scales` |
| [14. Predict new sites and years from a Bayesian model](https://wpeterman.github.io/multiScaleR/articles/v14-bayesian-prediction.html) | `v14-bayesian-prediction` |
| [15. Choose and check a Bayesian computation method](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.html) | `v15-bayesian-computation` |

Reference topics in the help system give the details that the vignettes
leave out:

| Help topic | What it covers |
|:---|:---|
| [`?multiScaleR`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-package.md) | Package overview and the main functions, step by step |
| [`?"multiScaleR-kernels"`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md) | Kernel families, sigma, hard radii, and effective distance |
| [`?"multiScaleR-fitting-methods"`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md) | The fitting backends and the uncertainty each reports |
| [`?"multiScaleR-prediction"`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md) | What each prediction represents and how uncertain it is |
| [`?"multiScaleR-glossary"`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md) | Plain-language definitions of the terms used |

`NEWS.md` lists the changes in each version, including the switch to
planar (2D) kernel distances in version 0.7.4.

## Experimental tools

[`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md),
[`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md),
the count, binary, and positive families of
[`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md),
[`diagnose_bayes_laplace()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_laplace.md),
and the `radius_grid` argument of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
are experimental. Their arguments and results may change in a later
release without a deprecation period.

## Citation

Run `citation("multiScaleR")` to cite the paper that describes the
method and the package version you used.

![multiScaleR package logo](reference/figures/logo.png)

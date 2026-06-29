## Test environments

- local Windows 11, R 4.6.0

## R CMD check results

Local validation for this release candidate completed with:

- 0 errors
- 0 warnings
- 0 notes

Command used:

```r
rcmdcheck::rcmdcheck(args = c("--as-cran", "--no-manual"),
                     check_dir = "checks",
                     error_on = "never")
```

## Changes from Previous Submission

This update is the 0.7.0 release of multiScaleR. Major changes include:

- self-contained vignettes for CRAN rebuilds
- new categorical landscape metrics, including edge, adjacency,
  information-theory, and class-level metrics
- new continuous surface texture metrics, including weighted surface
  heterogeneity metrics
- faster `profile_sigma()` evaluation for models with multiple expensive
  covariates, with a reproducible benchmark included under `tools/benchmarks/`
- improved diagnostics for landscape metrics that are undefined or uninformative
  during scale optimization
- expanded unit tests and vignette documentation for landscape and surface
  metric workflows

### Addressed CRAN feedback

CRAN checks for version 0.6.13 reported an ERROR on
`r-devel-linux-x86_64-fedora-gcc` while rebuilding `multiScaleR_Guide.Rmd`.
The failure occurred because the vignette setup attempted to download
precomputed `.RData` objects from GitHub during the CRAN rebuild.

Version 0.7.0 removes remote vignette downloads. The guide now relies on package
data, lightweight live examples, and precomputed objects bundled in
`inst/extdata`, so vignette rebuilds do not depend on network access.

## Downstream Dependencies

There are currently no downstream dependencies on CRAN.

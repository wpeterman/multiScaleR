## Test environments

-   local Windows 11, R 4.4.3
-   GitHub Actions (macOS, Windows, Ubuntu)

## R CMD check results

0 errors | 0 warnings | 1 note

The local note was:

- unable to verify current time

## Changes from Previous Submission

This update includes changes made since version 0.5.0:

- expanded unit test coverage
- improved error and diagnostic messaging
- added optional profile-likelihood scale intervals
- added sigma profiling and plot methods
- made profile-likelihood confidence intervals opt-in
- fixed complete-case alignment in multiscale optimization inputs
- added structured optimization diagnostics and a `diagnostics()` accessor
- added an optional custom refit hook for model classes that cannot use default model updates
- fixed singular-Hessian fallback SE values to remain numeric
- excluded missing raster cells from sparse kernel weighted averages
- added linear and user-specified sigma grids to `profile_sigma()`
- preserved original sparse kernel dot-product behavior for complete raster layers
- preserved point row identities across `kernel_prep()` outputs used during optimization
- ensured PSOCK workers use the same multiScaleR code as the main R session
- fixed complete-case row alignment when fitted model frames retain original row names
- updated vignettes and documentation

### Addressed CRAN feedback:

-   NA

## Downstream Dependencies

There are currently no downstream dependencies on CRAN.

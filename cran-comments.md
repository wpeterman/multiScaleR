## Test environments

-   local Windows 11, R 4.5.3
-   GitHub Actions (macOS, Windows, Ubuntu)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes from Previous Submission

This update includes changes made since version 0.5.0:

- added a Windows R-devel compatibility shim for Rcpp header compilation
- forced C++17 compilation for compatibility with CRAN Windows R-devel checks
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
- fixed PSOCK optimization for unqualified model calls such as `glm.nb()` after `library(MASS)`
- updated vignettes and documentation

### Addressed CRAN feedback:

-   Version 0.6.11 failed the CRAN incoming Windows R-devel pretest during
    package installation because the R-devel toolchain compiled the package
    with C++20 and the installed Rcpp headers failed before package code was
    compiled. Version 0.6.12 explicitly sets `CXX_STD = CXX17` in both
    `src/Makevars` and `src/Makevars.win`.
-   Version 0.6.12 still failed the CRAN incoming Windows R-devel pretest
    during package installation because the Windows R-devel headers available
    to win-builder did not declare `R_NamespaceRegistry`, while the installed
    Rcpp headers still referenced it. Version 0.6.13 adds a small
    Windows/R-devel-only compatibility header that declares the symbol before
    Rcpp headers are included. Debian R-devel checks for 0.6.12 otherwise
    passed with only the recent-update incoming note.

## Downstream Dependencies

There are currently no downstream dependencies on CRAN.

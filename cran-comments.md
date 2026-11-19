## Submission

This updates multiScaleR from 0.7.2, the current CRAN version, to 1.0.1.
Versions 0.7.3, 0.7.4, and 1.0.0 were released on GitHub and R-universe but
not published on CRAN, so this submission includes their changes. NEWS.md
lists them in full.

Main changes since 0.7.2:

* Bayesian estimation of scales of effect (`multiScale_bayes()`). Its slice
  sampler needs no Stan toolchain.
* New opt-in fitting backends for `multiScale_optim()`, cross-validation that
  re-estimates the scales in every fold, and prediction at new points and
  rasters.
* Reorganized documentation: fifteen vignettes and expanded help pages.
* Some results change for existing scripts. NEWS.md lists each case under
  "Changes that can affect existing scripts".

1.0.1 fixes a failure in the R-universe Windows builds of 1.0.0. Total RAM is
now read from the operating system in compiled code (src/msr_os_ram.c), and
the earlier PowerShell, /proc/meminfo, and sysctl queries run only as a
fallback.

## Test environments

* Local: Windows 11 x64, R 4.6.2, `R CMD check --as-cran` on the submitted
  tarball, with OMP_THREAD_LIMIT=2. Examples, tests, vignette rebuilds, and
  the PDF and HTML manuals passed. A second run with pandoc on the PATH
  checked README.md, NEWS.md, and their URLs.
* GitHub Actions on the same package code (only the version, NEWS.md, and
  README differ): macOS (R release), Windows (R release and R-devel), and
  Ubuntu (R release and R-devel).

## R CMD check results

0 errors | 0 warnings | 0 notes

The installed size is about 10 MB: R code 2.8 MB, vignettes 2.8 MB, example
data in inst/extdata 1.9 MB, and compiled code 1.4 MB. The example data are
the rasters and fitted objects that the examples and vignettes use, so they
run without network access.

## Reverse dependencies

There are no reverse dependencies on CRAN.

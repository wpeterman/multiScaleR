## Resubmission

This is a resubmission of 1.0.1. The incoming checks of the first upload
reported "Overall checktime 24 min > 10 min" on r-devel-windows-x86_64. The
tests (332 s) and the vignette rebuild (340 s) took 11 of those minutes. In
this version:

* The 77 slowest test files call `skip_on_cran()`. Locally, the tests now take
  75 s in CRAN mode instead of 224 s.
* Vignettes 2 to 15 open with a chunk that calls `knitr::knit_exit()` when
  R CMD check rebuilds them without NOT_CRAN=true, so the check neither loads
  the package nor runs their code. The vignettes in the tarball are built in
  full, and our GitHub Actions checks still run all of them. The quick-start
  vignette still runs in full during the check. Locally, the vignette rebuild
  step now takes 45 s instead of 249 s.
* The slow examples were already in `\donttest{}`. The remaining examples run
  in about 15 s.

Most of the remaining check time goes to installation, the package load
checks, and the R code checks of a large package (about 47,000 lines of R).

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
  tarball, with OMP_THREAD_LIMIT=2 and pandoc on the PATH. Examples
  (including `\donttest{}`), tests, vignette rebuilds, and the PDF and HTML
  manuals passed.
* GitHub Actions on the same R and C code: macOS (R release), Windows
  (R release and R-devel), and Ubuntu (R release and R-devel). These runs set
  NOT_CRAN=true, so they run every test and vignette.

## R CMD check results

0 errors | 0 warnings | 0 notes

The installed size is about 10 MB: R code 2.8 MB, vignettes 2.8 MB, example
data in inst/extdata 1.9 MB, and compiled code 1.4 MB. The example data are
the rasters and fitted objects that the examples and vignettes use, so they
run without network access.

## Reverse dependencies

There are no reverse dependencies on CRAN.

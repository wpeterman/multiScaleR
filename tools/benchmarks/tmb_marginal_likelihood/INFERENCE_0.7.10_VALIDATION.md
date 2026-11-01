# Inference integration 0.7.10

Validated November 1, 2026 on inference-dev, parent `3ab4d2d86aa6042dbe9e7858157f2f102a673243`.

The release integrates native formula designs and expanded coefficient counts,
supported weighted-SQ likelihood prerequisites, private spatial-only helpers,
and step-selection point, raster, and snapshot contracts. Fixed training
references, priors, diagnostics, starts, cold refits, zero and limiting cases,
curvature checks, and distinct uncertainty targets are preserved. Private
spatial or marginal fitting admission has not been promoted.

The normal vignette-enabled `R CMD check --as-cran --no-manual --no-tests`
returned Status OK for 0.7.10. Its archive SHA256 is
`587232438e1193c60519e1a4ebad79f152f75455efd6952b51904d579c2c6cb4`.
The separate full installed suite passed 19,326 assertions with no failures,
the same 12 warning identities, and 10 development-only skips. That suite ran
on the exact full56 archive
`ce109ef15c52794025774fd38ebe0a96d1a2421861b306e5cb80ae3e8a3b7c46`.
Only DESCRIPTION version and NEWS changed afterward. Runtime, tests, generated
help, and vignette files were independently compared across the 432-file source
set. The release check rebuilt vignettes normally; the unchanged runtime suite
was not replayed solely for the version bump. Timing comparisons are ineligible
under foreign R concurrency. All owned R exited; unrelated RStudio work was
untouched.

Qualified fixed-scale normalization outside the checked unmarked curves remains
a known gap in this release. Isolated candidate 105 passed 855 focused assertions
but is not included here. Its shared metadata and Cox repair still needs broader
adapter conformance and fresh full validation. Expired raw curve inputs, broader
amt movement updates, and private spatial reference qualification remain explicit
gates. See `unified_development_20261028/DEVELOPMENT.md` and the local continuation
records. Universal formula or model support has not been established.

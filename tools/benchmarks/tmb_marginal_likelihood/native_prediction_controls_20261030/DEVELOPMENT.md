# Native prediction controls, version 0.7.8

Zero-random-effect prediction recognizes formula syntax independently of its
environment. Random-only grouping/slope inputs may be omitted where the native
provider permits it; fixed variables and offsets remain checked. Point, snapshot
and raster paths retain native results. Zeroing effects is not integration.

GAM snapshots bind native term selection and smoothing-covariance controls.
Component outputs, unavailable native SE errors, warnings, and compact RDS
snapshots are preserved. Smoothing uncertainty remains distinct from selected
kernel-scale uncertainty. Provider-only contract metadata remains supported.

Validation before the required metadata-only version bump: package check
Status: OK with normal vignette/index rebuilding; full exact-archive suite
6,904 passed, zero failed, 12 unchanged warnings, 10 development-only skips.
The final focused run passed 294 with eight development-load namespace warnings
that did not recur in the installed suite. All 355 source hashes were verified.
Only DESCRIPTION version/NEWS metadata changed afterward; validation.json
records the exact final identities and qualifier. Prior failures are retained
locally. Priors, diagnostics, numerical gates and unrelated work are preserved.

Next priorities: selected random-effect and selected/excluded GAM term row
requirements, using direct native oracles and retaining strict fixed/offset checks.

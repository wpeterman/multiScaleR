# Functional probability precision

Validated October 28, 2026. The standard package check and complete installed
tests passed in separate phases, preserving test coverage within resource caps.

The shared prediction-contract batch closed successfully before these changes.
Its exact archive, numerical parity checks, resource observations and package
validation remain authoritative for its frozen source. Its live verifier will
correctly flag this later deliberate source drift.

The new diagnose_bayes_probabilities report addresses a documented integration
gap: passing coefficient and derived-scale diagnostics does not establish the
precision of a scientific probability. Users provide named fitted-scale
thresholds and an optional absolute probability MCSE target. Original chain
order and the fitted scale design, random-effect map and unit conversion define
the binary events. Mean ESS, MCSE, Rhat and per-chain probabilities describe
those events. The report retains the original diagnostics without a new overall
pass flag. A target comparison never overrides existing gates.

Thresholds concern kernel scale parameters, not transformed kernel-weight radii.
Constant indicators, too few chains/iterations and nonfinite derived scales
remain explicit limitations. Compact/approximate fits and private Gaussian
projection blocks report unavailable; importance candidates are never MCMC
chains. The result retains summaries only, with no sampler or chain arrays.

unit02 passes 145 focused assertions, including independent ordinary-array
reconstruction and posterior-package indicator MCSE/ESS/Rhat references, original
chain order, named thresholds and fitted/crossed designs, map-unit equivalence,
scalar/indexed original columns, constant/stuck/nonfinite/undersized chains,
malformed layouts and requests, raw failures, and object/RNG immutability.
compatibility01 passes 290 assertions across the existing derived diagnostics,
Bayesian scales, log summaries, prediction contracts and strict export tests.
Both runs have no errors, failures, warnings or skips. unit01 retains five
failures: a matrix-shaped target rejection edge fixed in code, and four test
fixtures recentered to have variable indicators. No gate was weakened.

documentation01 generated the new Rd topic and only the intended extra export,
checked help syntax, and rendered the computation vignette (64,012 bytes) using
compact stored examples. Visual browser review remains policy-blocked from the
preceding batch; execution/help validity are not a claim of visual inspection.

snapshot01-input.json pins 311 files. Its audit found only four expected existing
file changes and three additions, all unrelated source intact, the original
derived-diagnostics code unchanged, and all closed prediction-contract artifacts
unchanged. A fresh source archive was built with vignettes, followed by a
standard check without nested tests and the complete unchanged installed test
tree with NOT_CRAN=true. This split preserves parallel-test coverage under the
global interpreter limit. Do not add focused counts to complete-suite totals.

The complete gate passes: 5,016 assertions, zero failures, 12 existing test
warnings and 10 development-only skips. The standard vignette-enabled check
reports zero errors, warnings or notes. All 311 source files match, there are
no additional runtime sources, archive SHA256 is
E8BC0137D2073AB6C1CFA0E01146415C9EDFBFAC89DD529D7C63BA054DE1154F,
and the unchanged extracted archive test tree was executed with NOT_CRAN=true.
Prior closed artifacts and historical pinned evidence remain unchanged; index
and HEAD are unchanged. All 47 resource observations are known/foreign-R-free:
maximum four actual R interpreters, minimum 16.08 GiB available RAM and 51.24 GiB
free disk. No package source, test or numerical gate changed during validation.
final-source-evidence.json freezes this exact scope before later development.

No calibration
study, sampler launch, fitting speed claim or public spatial promotion is part
of this batch. No staging, commit, version bump, merge or push.

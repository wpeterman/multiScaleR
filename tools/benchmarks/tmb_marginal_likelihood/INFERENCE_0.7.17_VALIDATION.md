# Version 0.7.17 validation and handoff

Branch: `inference-dev`. This release publishes the accepted prediction and
training-only CV metadata work based on commit `457afeb1`. The November 6 user
request explicitly authorizes the patch version bump, commit and push. There
is no master merge or further numerical development in this closeout.

## Published behavior

- Supported nested scaling, polynomial and spline calls and formula offsets
  retain training definitions across native, analytic, Bayesian and private
  predictions. Missing or ambiguous metadata produces an explicit error.
- GLM, NB2/Gamma, weighted-SQ and complete-set choice CV use final analysis-fold
  predictors. Assessment rows never supply learned transformation constants.
- Ordinary LM/GLM/NB-GLM and qualified time-fixed Cox predictions and curves use
  original fitted terms for dot expansion, exclusions, factors and offsets.
  Unknown subclasses retain their provider paths.
- Ordinary unmarked occupancy and pcount predictions use checked original
  state/detection transformations, factor contrasts and native covariance.
  Point rows, prepared outputs, raster chunks and curves share the native target,
  including overall ZIP abundance. Detection curves use observation covariates;
  `lambda` aliases the state component in ordinary unmarked plotting.
- Private fixed-state likelihood snapshots require verified training inputs and
  explicit realized effects. They provide conditional arithmetic, without
  admitting new fitting methods, SEs or latent-distribution integration.

Full model parameter counts, priors, fixed reference normalization, diagnostics,
coefficient SEs and existing scientific/admission guards are preserved. Native
confidence limits still condition on selected scales. Setting effects to zero
does not integrate over them. Bayesian uncertainty remains separately defined.

## Validation evidence

The accepted 484-file source was validated at version **0.7.16**, then applied
without source drift. Release 0.7.17 changes only DESCRIPTION/NEWS metadata
after that gate, plus development handoff documents outside the package build.
No R validation was launched after the previous 06:00 November 6 logoff cutoff.

| Gate | Result |
|---|---|
| Normal vignette-enabled package check, `--no-manual --no-tests` | OK |
| Complete installed suite from the exact checked archive | 37,858 passes, 0 failures, 0 warnings, 10 existing helper skips |
| Focused consistency suite | 2,593 passes, no failures, errors, warnings or skips |
| Source, checked archive and installed dependency identity | Independently verified |
| Original scientific markers | All 14 preserved |
| Original skip identities and installed test ledgers | Preserved, with reviewed regression additions |
| `--as-cran` incoming checks | Older external metadata timeouts remain unresolved |

Accepted archive SHA256:
`99b936e13264f662547d8ab7ac251f54ed787e8973bf453ab1acc711f13ab765`.

Detailed evidence is retained locally in `utarget_20261106/`, particularly
`freeze08-source01.json`, `full08-verified-state01.json`, `static-review13.json`
and `FINAL14.json`. The accepted archive and the earlier failed `udot_20261106`
archive/logs are retained locally. Bulk evidence, build trees and intermediate
objects are not part of this release commit. This is a normal package-validation
result, not a CRAN incoming or scientific-calibration claim. The version-only
bump has static verification, not a new version-0.7.17 package check.

## Start here on resume

Read [the next compatibility plan](continuation_20261030/NEXT_SUPPORTED_COMPATIBILITY_20261106.md)
and [the development history](unified_development_20261028/DEVELOPMENT.md).

1. Inspect branch/source identity, other agents' work, global R processes and
   available resources before starting a new bounded validation producer. Reload
   the updated package in a fresh user session; existing sessions were untouched.
2. Audit component-specific unmarked prediction preparation first. Currently
   state and detection contexts are both prepared before the requested type is
   evaluated. Qualify native target dependencies before narrowing this work,
   retaining ZIP derived-target requirements and conditional random-design
   guards. Reuse the existing deterministic fixtures. This finding is source
   inspection only, not a validated extension of support.
3. Then qualify original random-design metadata for conditional learned
   unmarked slopes. Preserve full model K and original likelihood/diagnostics
   even when a requested prediction uses only one component.
4. Keep latent-model CV pending a proper held-out joint density and explicit
   site/observation folds, truncation and new-group policies. Component means or
   SEs do not establish a CV scoring distribution. Additional subclasses,
   weighted metrics and private no-year spatial methods retain their own gates.
5. For the next coherent source edit, regenerate affected help, run focused
   regressions, then a fresh vignette-enabled check and complete installed suite
   from its exact archive. Prior check coverage becomes historical after edits.

Keep one validation pipeline, bounded workers and one thread per interpreter;
inspect foreign processes and leave them untouched. Retain compact atomic
records, unique inputs and all failed/interrupted attempts. Delete only exactly
identified, verified regenerable artifacts. No assertion-count target is needed:
consolidate duplication only after a coverage audit retains the scientific
cases, numerical oracles and failure behavior.

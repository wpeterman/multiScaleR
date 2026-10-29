# Resumable stopping point: shared Gaussian layer and spatial-scale prototype

Finalized October 28, 2026, before the 3:30 PM Eastern logoff. Read REPORT.md,
PROTOTYPE.md, PLAN.md, final-source-audit.json and final-process-audit.json first.

## Safe stopped state

At 14:42 Eastern all owned numerical jobs/controllers had exited. The local
15:25 deadline watchdog was stopped only after exact PID, creation-time, command
and parent identity checks; its absence was confirmed. Fresh inventory contained
no R interpreter. No foreign process or automation was stopped or modified.
Available RAM was 18.66 GiB and disk 53.43 GiB. No queued job is intended to resume.
Do not restart any retained manager, check, pilot or promotion blindly.

No staging, commit, push, merge or version bump. Branch remains
codex/bayesian-hierarchical-scales, HEAD87bf8779fe571921f2fd5d3e45c0f9e478530022,
version0.7.4. The index is empty. Preserve all concurrent/unrelated edits, old
studies, original failures, source snapshots and producer records. No subagents
or cross-chat messages were used. Previously policy-blocked DLL cleanup was not
retried. Do not delete another owner's checked source or stop their jobs.

## Completed implementation

New private modules R/latent_gaussian.R and R/spatial_scale_prototype.R provide
proper fixed-covariance Gaussian projection blocks and one spatial radial-basis
field on log kernel scale, plus optional existing annual effects. New SD prior
is explicit. Range/covariance are fixed; first prototype fixed scale mean is
intercept-only. Public formulas/interfaces/priors/diagnostic thresholds stay
unchanged. Shared scale/response components branch only for new Gaussian blocks;
legacy index arithmetic is identical. Unsupported scalar/legacy prediction
routes reject blocks. A separate private class and spatial prediction adapter
preserve frozen training projection and coherent new annual groups.

Bayesian joint Laplace uses existing numerical gates. The normalized likelihood
has joint/latent-Laplace objective adapters only, not a qualified ML fitter,
confidence intervals or limiting-model witnesses. No public prototype API,
random-slope syntax, temporal covariance estimation, estimated spatial range,
new response family or generic factorization is claimed.

## Retained validation

- mathematical-gates.json:111assertions pass,0failure/error/warning/skip.
  Includes90new block/spatial/integration assertions and21strict API assertions.
- legacy-parity.json:8original cases, ordinary components identical;
  ADvalue/gradient/Hessian differences0.
- focused-complete.json:1639existing assertions/27files pass; this preceded the
  final private constructor guard, which is covered by current mathematical and
  full package checks. Do not sum overlapping test counts.
- prototype-check-complete.json:StatusOK,0check errors/warnings/notes; installed
  tests4763pass/0fail/12warnings/10development-helper skips; examples/donttest and
  rebuilt vignettes pass. Reused a completed concurrent same-package check only
  after strict source auditing, not an additional duplicate check.
- check-reuse-audit.json/reused-* logs and manifest:273shared files byte-identical.
  All model/native/tests/help/source-vignettes/runtime inputs match.32external-only
  build-excluded files and21own-only generated/baseline auxiliaries are disclosed.
  No archive-byte identity claim. The actual checked candidate is recorded there.
- final-source-audit.json:all273shared files still match live source and frozen
  numerical engine; no missing/changed files. All pilot source identities stable.

Frozen numerical engine:
C:/Users/peterman.73/.codex/scratch/msr-spatial-20261028-source

Earlier own baseline check failed strict export inventory for the concurrently
added prepare_unmarked_data(). Only that exact expected name was added; strict
set equality is preserved. Failed baseline and prior operational exits remain
retained and must not be relabeled as passes. See OPERATIONS.md.

## Scientific outcomes: no pilots may repeat automatically

PILOT_PROTOCOL.md fixed null and smooth spatial cases:24sites,3years,72Gaussian
observations,rank3field,range350,population scale220,annualSD.15,fieldSD0/.45,
explicit half-normal field prior scale.4. Both use the same inputs for Laplace
and NUTS; they are synthetic fixture illustrations, not calibration.

Original NUTS:4serial chains8000/2000,.99/depth12. Original raw/derived diagnostic
gates and functionalCDF MCSE<=.005 remain unchanged. Neither reference qualifies.

- pilot-null.json:Laplace numericalpass,.72sec; NUTS277.95sec,18divergences,
  derivedpass,2median-CDFprecision failures.204/228CDFflags,2unavailable.
- pilot-smooth.json:Laplace numericalpass,.54sec; NUTS423.80sec,0divergences,
  raw/derivedpass,but1field-SDmedian-CDFprecision failure.194/228flags,1unavailable.
- pilot-null-affine.json:mean-Hessian audit stopped before sampling; min eigenvalue
  -1.414376. No repaired covariance or affine sampler fit.
- pilot-null-dense.json:one separately declared dense-metric null diagnostic,
  exact original input hash,newseed202610292,same8000/2000controls/priors/gates.
  NUTS241.22sec,23divergences,derivedpass,2CDFprecision failures.205/228flags,
  2unavailable. It did not resolve reliability. Original Laplace output reused.

All completed outputs/failures retained; no further sampler configurations/fits
are authorized by these frozen protocols. No sampler/rawdrawarray/ADtape saved.
final-pilot-summary.json provides compact metrics and independent-repeat median
changes. Median site-year log-interval width ratios(Laplace/NUTS)3.29null,4.11smooth.
Maps can look similar while SD uncertainty differs materially. No calibration,
false-positive, nonstationarity detection or validated NUTS-equivalence claim.

Report figures read compact summaries only, with fixture factor-order mapping
corrected explicitly; both PNGs were visually inspected. They label references
provisional. report-figures.R performs no fitting.

## Next authorized development requires a new bounded design

First qualify a spatial reference and investigate null SD/latent geometry under
unchanged target/gates. Then develop conditional/block integration with latent
reconstruction, zero-SD limits and honest uncertainty targets. A likelihood fitter
requires cold score/stationarity, identification and limiting-witness gates.
Do not pass a coupled spatial field into independent-year scalar factorization.
Broader null/weak/confounded/spatial-holdout calibration precedes public API and
estimated-covariance extensions. Existing exact-response marginalization research
is not assumed to be more efficient just because it reduces dimension.

Preserve original annual-study decisions including sbc123/187 limiting-witness
and prior018 competing rejections. No archived fit needs to run to resume source
work. Recheck live source identity/resources before any future test or fit;
source drift requires qualifying the old check scope, not overwriting changes.

The study is compact(about1MiB). Own source/check trees contain regenerable
compiled/check products(about31.4/28.0/31.5MiB) but also retained provenance.
No cleanup is needed for current headroom. Preserve manifests/scripts/results;
classify exact owned unloaded build products before any future deletion.

# NB2/Gamma CV completed locally

Applied and validated on inference-dev at e54cff6d, version 0.7.15.
No staging, commit, version bump, merge or push. All existing dirty work remains.
All 472 curated root source hashes match the accepted final source manifest;
only 18 reviewed source/help/test/vignette/NEWS files changed in this batch.

## Supported behavior

- Ordinary and spatial NB2/Gamma log-link GLM CV refit within each training set.
  Theta and Gamma density dispersion are frozen for held-out scoring. Gamma and
  fixed-theta NB2 require native refit; estimated NB2 also supports analytic
  profile/joint engines. All existing admission/start/cold/curvature checks remain.
- NB2/Poisson count comparisons and Gamma/Gaussian continuous mean/training-IC
  comparisons preserve likelihood-measure guards. Selected folds retain their
  actual family, nuisance parameters and controls. Compact/full AIC/AICc/BIC agree.
- Estimated nuisance parameters add one to native K; fixed ones add zero.
  Regression rank, optimized scales/shapes/radii and native positive-weight nobs
  retain their original conventions. No sum-frequency substitution for nobs.
- Estimated Gamma density uses native AIC's deviance/sum-frequency scale,
  separately from native Pearson coefficient SEs. Known Gamma dispersion uses
  a scoped normalized density before fitting, avoiding the inconsistent native
  dispersion-field-only AIC convention. Existing fits/global methods are unchanged.
- Seventeen deterministic spatial cases cover all four kernels through admitted
  backends, Gaussian weighted-SQ, fixed-radius SE/count rules, original exact
  Poisson-limit admission and assessment perturbation. Radius Hessian SEs remain
  absent; coefficient covariance and optimized-radius K remain. Ordinary native
  fixtures also verify factors/products, polynomial/dot/offset/saved-scale terms,
  rank aliases, zero frequencies, response support and retained native failures.
- Priors, normalized targets, reference transforms, private structural guards
  and original diagnostics are unchanged. Scores remain plug-in quantities,
  without parameter/selection integration, observation intervals or calibration.

## Validation

R 4.6.1 on Windows 11, one thread per interpreter, no foreign jobs altered.
Fresh normal vignette-enabled final check: **Status: OK**.
COMPLETE tests installed from that exact archive: **34,928 passes,
0 failures/warnings, 10 unchanged helper skips**. Every old per-file
assertion count and skip identity is preserved; the only final increment is 67
reviewed radius/exponential assertions. All 14 scientific markers remain.
Focus19 passed 1,894 with no failures, errors, warnings or skips. Generated help and worked
NB/Gamma CV-selection examples built successfully; UTF-8 and original R chunks
outside the deliberately expanded CV vignette were preserved.

Final archive SHA256: 62eaf30c0ebc6f7ed2a07fe1ee9c3e042e92d65d502ec9dac04517bab2ef112a
Durable source archive: accepted_multiScaleR_0.7.15.tar.gz.
See full20-verified-state01.json, applied16.json and final26.json for exact
source/archive/installed-library/test identities. Historical full14 also passed
Status: OK,34,861 assertions and the same10 skips before the final test additions.
This is software qualification, not a recovery, calibration, speed or new
cross-platform/CRAN-incoming claim. Concurrent foreign-R timing is ineligible;
the older external incoming-metadata timeout gate remains separate.

## Preserved attempts

| Phase | Pass | Fail | Error | Warning | Skip |
|:--|--:|--:|--:|--:|--:|
| focus05 | 21 | 3 | 4 | 0 | 0 |
| focus06 | 1078 | 11 | 2 | 1 | 0 |
| focus07 | 1490 | 0 | 0 | 0 | 0 |
| focus10 | 274 | 7 | 1 | 0 | 0 |
| focus12 | 3534 | 0 | 0 | 0 | 0 |
| focus19 | 1894 | 0 | 0 | 0 | 0 |

Probe02 retained all eight native references, including the fixed-Gamma AIC
trap; its rounded JSON is descriptive, while tests provide exact density/SE/K
oracles. Focus05 exposed an oracle weight-symbol environment error. Focus06
exposed MASS's previous IRLS family size versus final native theta/twologlik.
Focus10 exposed incorrect test expectations about valid estimated-scale SEs,
elapsed-time equality and error text. Corrections preserve native scientific
behavior; every failed log/input/result remains. One native default weighted/
rich NB alternation-limit failure stays an expected failed/unrankable fold.
No seed search, automatic retry, tolerance/prior/diagnostic relaxation or
private fitting promotion occurred.

## Resources and continuation

No owned numerical work remains. Cleanup removed only 36
hash-verified regenerable object files (16,911,715 logical bytes)
inside this batch's sealed temporary workspaces. All sources, DLLs, archives,
unique inputs and failures remain. Observed free disk: 72.63GiB.
No full chains, tapes, candidate arrays or fit-object checkpoints were saved.

Read RESUME.md and NEXT_CV_COMPATIBILITY.md before further development. Future
edits require fresh identities/validation; never replay a closed producer.
No commit is authorized by this checkpoint. A later explicit commit request
requires the repository's version/NEWS and validation rules.

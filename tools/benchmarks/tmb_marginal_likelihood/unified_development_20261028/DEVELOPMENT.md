# Current validated continuation: multiScaleR 0.7.6

Use **inference-dev**. Read `../fixed_metric_integration_20261030/DEVELOPMENT.md` and its adjacent `validation.json`
for the current shared prediction, weighted-SQ, fixed-metric, and legacy-input
contracts. The Windows check completed with zero errors/warnings and one elapsed-example
timing NOTE under concurrent work; exact-archive tests passed
6799 assertions, zero failures, 12 existing warnings,
and 10 development-only benchmark skips. Priors, diagnostics, starts/cold/limiting
gates, and uncertainty distinctions remain intact. No master merge is implied.

The earlier consolidation below is historical source/branch context; its byte
manifest must not be substituted for current-source coverage or used to reset
later work. Detailed local study records and unrelated edits remain preserved.

## Historical consolidation

# Unified scale-inference development

Consolidated October 28, 2026. This is the continuation index for the marginal
likelihood, integrated likelihood, and Bayesian tracks. Read this document and
run `verify_checkpoint.py` before relying on the retained validation. Historical
study reports and failures remain authoritative for their own frozen sources.

## Branch decision and scope

Use **codex/inference-dev** in the existing multiScaleR checkout. It was created
from `87bf8779fe571921f2fd5d3e45c0f9e478530022`, the common HEAD of master,
codex/bayesian-hierarchical-scales, and codex/performance-joint-scale-assessment.
All three tracks already occupied this working tree. No code merge or source
replacement was needed. The old branch references are preserved.

The development implementation is still largely uncommitted or untracked.
Creating this branch does not store those changes in Git history. Do not switch
or reset the checkout assuming another branch contains a saved copy. The source
manifest in checkpoint.json records the dependency-complete current candidate;
it is an integrity record, not a backup. No staging, commit, version bump, push,
or merge to master is authorized by this consolidation.

Do not merge the old dev/development branches wholesale. They have diverged
histories and are not the location of the current inference implementation.
Other worktrees and unrelated edits remain untouched.

There is no present source-level reason to keep the three tracks on different
branches. Share preparation, conditional response densities, kernel arithmetic,
Gaussian blocks, and explicit capability/target descriptions. Keep experimental
methods private until their own gates pass. Separate future worktrees are useful
only for actual concurrent edits or a frozen release candidate; they are not a
reason to maintain three competing versions of the shared model.

## Necessary memory from each track

| Track | What is implemented | What must remain qualified |
|---|---|---|
| Marginal likelihood | Normalized response density; neutral model preparation; target metadata and capability checks; private joint latent-Laplace and accurate independent-group integration; scalar cache; shared proper Gaussian blocks and a private spatial-scale basis prototype | Factorized integration requires independent scalar components. Spatial likelihood adapters are objectives, not a qualified ML fitter. Population-parameter uncertainty and likelihood intervals remain incomplete. |
| Integrated likelihood | Specialized profile GLM/NB2, joint NB2/mixed and conditional-logistic routes; scale inference and prediction; compact native prediction snapshots; explicit unmarked preparation; fixed-radius prediction and opt-in raster support matching | Native conditioning, random-effect, detection, offset, trial, site, and uncertainty contracts belong to each backend. Snapshots reuse fixed outputs. Neither snapshot speed nor small-map timing establishes fitting speed or large-map performance. |
| Bayesian | Public hierarchical scale model with NUTS, Laplace and corrected marginal routes; setup/prior checks, joint posterior summaries, new-observation/year prediction, simulation, plots and worked vignettes | Importance acceptance is model dependent. NUTS diagnostics alone do not establish probability or quantile precision. The coefficient-integrated posterior lacks reconstruction for public prediction. Spatial references remain unqualified. |

The common model is not one inferential target. For population parameters theta
and latent effects b, the hierarchical likelihood integrates
`p(y | b, theta) p(b | theta)` over b. The posterior also includes declared proper
population priors and propagates uncertainty in theta. A Laplace approximation
is an integration method that can serve either target. Existing profiled and
conditional likelihoods retain their own conditioning. Do not substitute a
posterior callback for a likelihood objective or expose posterior objectives as
likelihoods for AIC.

### Current implementation map

| Layer | Source to read first |
|---|---|
| Shared specification and response density | R/scale_model_prepare.R; R/scale_response_likelihood.R; R/scale_model_capabilities.R; R/scale_inference_metadata.R |
| Hierarchical likelihood | R/marginal_scale_model.R; R/marginal_scale_likelihood.R; R/marginal_scale_fit.R |
| Independent scalar integration | R/factorized_scale_*.R, including limits, score, cache, optimizer and posterior summaries |
| Existing integrated/profile routes | R/multiScale_optim.R; R/scale_inference*.R; R/profile_glm_*.R; R/joint_nb2*.R; R/joint_mixed*.R; R/conditional_logit.R |
| Posterior computation and diagnostics | R/multiScale_bayes.R; R/bayes_scale_likelihood.R; R/bayes_scale_laplace.R; R/bayes_scale_marginal.R; R/bayes_scale_importance.R; R/bayes_scale_diagnostics.R |
| Prediction contracts | R/bayes_scale_prediction.R; R/bayes_scale_predict_scales.R; R/prepare_prediction.R; R/prepare_native_prediction.R; R/prepare_unmarked_data.R; R/predict_profile_glm.R |
| Proper blocks and private spatial model | R/latent_gaussian.R; R/spatial_scale_prototype.R |

### Scientific decisions that must survive consolidation

* Accurate annual integration accepts sbc-111 and the exact zero-SD sbc-057
  case. Retain prior_shift-018's competing-optimum rejection and sbc-123/187's
  limiting-witness rejections. A better ineligible objective is not an accepted
  optimum. Accurate integration does not cure weak identification.
* **Supersession:** the October 27 weak-case transport median discrepancy was
  investigated by four much longer matched-start, independent-seed references
  in consolidation_20261028. They passed sampled, derived, and fixed-CDF
  precision gates, and the discrepancy did not persist. Do not describe that
  particular discrepancy as still unresolved. This is selected-case agreement,
  not general transport calibration or public readiness.
* Exact Gaussian coefficient integration passed density/gradient/Hessian audits,
  but both fixed sampling pilots failed the probability-precision target and
  had worse useful information per second. Dimension reduction alone is not a
  performance argument. Original joint coefficient reconstruction is missing.
* The completed predictive study used 56 matched inputs: original NUTS qualified
  on 54; the marginal approximation on 10, all Gaussian. Original marginal
  outcomes were 10 eligible, 39 rejected, and seven execution failures. Preserve
  these original denominators separately from repair rechecks. All eight
  repaired Poisson cases compute but still fail unchanged importance gates;
  one repaired comparison has a substantial predictive-tail discrepancy.
* The spatial prototype has a frozen proper basis/covariance and fixed range,
  a field on one selected kernel with an intercept-only fixed scale mean, optional IID annual
  effects, and an explicit field-SD prior. It is private. Ecological kernel
  distance, spatial smoothing range, and log-scale field amplitude differ.
  Laplace interval widths differed materially from provisional NUTS outputs.
  Original null/smooth and the declared dense null reference did not all pass
  the strict reference gates. An affine attempt stopped on an indefinite
  Hessian. No covariance repair, automatic pilot replay, or calibrated
  nonstationarity claim follows from those runs.
* Phase63's equivalent transforms on 2,500 existing cells took 9.95 versus
  0.79 seconds; repeated native predictions took 3.21 versus 0.22 seconds
  including preparation. These are separate single-pair workloads. Cold
  startup was slower for the snapshot workflow. `include_fixed = FALSE` and
  `match_buffer = FALSE` remain the defaults.

### Validation status and superseded handoffs

The final phase63 attempt 2 completed a vignette-enabled package check with
**0 errors, 0 warnings, 0 notes**. Installed tests reported **4,763 passing
assertions, 12 test-level warnings, and 10 development-helper skips**. Test-level
warnings are disclosed even though package-check status is clean. The spatial
track reused that same check after a strict source audit; it is not a second
independent check. Do not sum overlapping focused/full assertion counts.

This consolidation independently compared all **305 files** in the checked
source manifest with the live checkout and found no differences before branch
creation. verify_checkpoint.py repeats that read-only comparison, checks for
additional runtime sources, and pins the evidence and continuation files. No R
fit, benchmark, documentation regeneration, or package check was repeated because
this batch changes no package source, help, examples, or vignette. These new
handoff files live under an already build-excluded research directory.

Older notes saying the full combined check is pending, or reporting the failed
export-inventory attempt as the last gate, are superseded by phase63 attempt 2
and latent_gaussian_20261028/check-reuse-audit.json. Preserve those failures as
history. A later source change invalidates reuse for the changed candidate;
validate that candidate rather than extending the claim by branch name.

## Next feature batches, in priority order

1. **Finish shared prediction contracts and adapters.** Neutral preparation,
   seven-route capability checks, and inference metadata already exist. Extend
   the remaining prediction description with target population, known/new effect
   policy, conditioning, integrated/omitted uncertainty, and original-coordinate
   readiness. Reuse current adapters and preserve public formals, classes and
   original diagnostics. Cover known/new rows, complete joint draws, offsets,
   trials, training transformations, and old saved fits. Metadata must never
   grant numerical eligibility or manufacture an interval. Do not add a unified
   convenience front end or automatic method fallback in this batch.
2. **Expose functional precision and qualify the spatial reference.** Add a
   bounded diagnostic interface for requested scale/probability thresholds,
   returning chain-based MCSE/ESS or explicit unavailability. Preserve existing
   sampled/derived gates. Declare a new bounded spatial reference design before
   fitting, then investigate near-zero field-SD geometry without changing the
   posterior, repairing Hessians, or selecting passing seeds. This is the first
   scientific blocker for the flagship spatial extension; the older transport
   discrepancy is not the blocker now.
3. **Develop boundary-aware likelihood uncertainty.** Establish population-SD
   uncertainty and predictive intervals with a declared bootstrap/profile or
   test-inversion design, retaining exact zero-SD faces, all starts, cold scores,
   competing modes, and limiting witnesses. A spatial ML fitter also needs those
   gates before an accepted-fit interface. Label plug-in conditional predictions
   honestly until population-parameter uncertainty is propagated. Do not use
   NUTS as exact truth for unpenalized ML or attach an unvalidated Wilks cutoff.
4. **Calibrate the spatial model before broadening it.** With a qualified
   reference, use predeclared constant-scale, smooth alternative, weak-signal,
   response-only spatial-confounding, and spatial-holdout cases. Assess false
   detection, scale/SD recovery, intervals, prior/basis sensitivity, failure rates
   and computational cost. A positive equal-tail SD interval is not a test of
   nonstationarity. Separate scientific datasets from deterministic fixtures.
5. **Add further Gaussian blocks in small batches.** Random slopes first, then
   explicit temporal covariance with time-gap and SD semantics; estimated spatial
   range follows identifiability studies. Coupled effects require joint or newly
   validated block/sequential integration. Generalize scalar factorization only
   after proving independent one-dimensional components and tail bounds. Add
   response families through density/data/derivative/prediction/simulation
   adapters with their own audits. Formula syntax does not imply support.

Spatial workflow source/CRS/unit and coverage diagnostics can proceed as a small
independent batch on the same branch. A larger-map benchmark should use an
appropriate existing map in an uncontended window, with equivalence checked
before accepting timing or peak-memory gains. Do not rerun accepted small-map
pairs without a new validation question.

Exact coefficient reconstruction is a later opt-in posterior batch, justified
only if its cost and precision warrant it. Draw one full conditional coefficient
vector per retained latent/hyperparameter vector and share it across prediction
rows. Do not substitute modes or independently resample coefficients per row.

## Assertions required to continue

| Contract | Required assertion and existing test anchors |
|---|---|
| Candidate identity | Branch/HEAD/version, unchanged checked-source and evidence identities, additional-source inventory, and preserved index/working files. `verify_checkpoint.py`; checkpoint.json. Source drift is a reason to revalidate, not overwrite edits. |
| Shared statistical model | Normalized ordinary/AD values, scores and Hessians agree at physically matched vectors. Offsets, trials, NB2 size and Gaussian SD retain their meanings. Tests: test-scale-response-likelihood.R; test-scale-model-preparation.R; test-latent-gaussian.R. |
| Target and capability separation | Proper population priors only for posterior targets; normalized latent densities retained in both. Metadata is descriptive. Unsupported combinations fail before numerical work; no fallback after rejection. Tests: test-scale-inference-metadata.R; test-scale-model-capabilities.R. |
| Factorization and boundaries | Only independent scalar blocks factorize; crossed response/site, spatial or temporal coupling is rejected. Preserve exact zero SD, integration tails/refinement, cold score/Hessian, competing and limiting witnesses, and all six retained scientific decisions. Tests: test-factorized-scale-likelihood.R; test-factorized-limits.R; test-factorized-score.R; test-factorized-refinement.R; test-marginal-scale-likelihood.R. |
| Sampling and correction | Preserve every prior and gate: Rhat, bulk/tail/derived ESS, divergences, depth, E-BFMI, importance population/ESS/tail diagnostics and functional precision. Invalid AD importance fallback is value-only; it cannot validate invalid derivatives or curvature. Constant indicators report nonestimable MCSE. Tests: test-bayesian-derived-diagnostics.R; test-bayesian-target-evaluation.R; test-bayesian-importance.R. |
| Prediction and saved fits | Freeze training normalization, contrasts, bases, distance convention, row/site IDs and parameter order. Reject missing scale covariates and unknown fixed levels. Unused future responses cannot affect predictions. Preserve joint vectors/weights and one shared effect per new group per replicate, including row reordering/splitting and reload/RNG contracts. Tests: test-bayesian-prediction.R; test-prepared-prediction.R; test-prepared-unmarked.R; test-unmarked-fixed-prediction.R; test-scale-inference-metadata.R. |
| Native and spatial predictions | Preserve native component/re.form/detection/offset/trial controls and conditional-SE restrictions. Snapshots invalidate changed inputs/controls and are fixed-output reuse. Matched raster support uses intersecting-cell buffers, finite-cell normalization and center-cell NA masks. Never silently cap a fixed radius above max_D or call relative risk/selection a probability. Tests: test-unmarked-spatial-prediction.R; test-projection-buffer.R; test-spatial-scale-prototype.R; test-spatial-scale-integration.R. |
| Uncertainty and scientific claims | Keep posterior expected response, posterior predictive distribution, conditional empirical-Bayes summaries, native conditional SE, joint frequentist uncertainty and plug-in intervals distinct. Preserve failures in planned denominators; use dataset-level calibration uncertainty and matched inputs/priors/seeds. Numerical eligibility does not establish coverage or pooling superiority. |
| Package completeness | Include dependency-complete R/native registration, Rcpp exports, tests/fixtures, generated help, examples, runtime example data and vignettes. A successful full dirty-tree check does not validate a selected subset for a future commit. |

## Operational continuation

Read the current process inventory before launching anything. Historical PIDs,
cutoffs and one-night automation instructions are closed-study evidence, not live
authorization to restart producers or controllers. No work is queued by this
consolidation. Do not message or resume the three old chats automatically.

For new numerical work, use the conservative common limits: at most two owned
fit workers, at most four actual R interpreters globally, one numerical thread
each, available RAM greater than 12 GiB and free disk greater than 5 GiB.
Recognize rsession variants and Rserve. Unknown inventory closes admission.
Full checks and timing/peak-memory pairs require a genuinely foreign-R-free,
uncontended window. Never stop another owner's process to create that window.
Use compact atomic checkpoints, and retain rejected/failed/pending outcomes.
Avoid serialized full chains, candidate populations, samplers and AD tapes.
Do not retry historical policy-blocked cleanup by another mechanism.

Existing study scripts pin R 4.6.1 explicitly; the general local-environment
skill records R 4.6.0 as its default. Preserve each evidence run's recorded
toolchain and verify the selected executable and optional glmmTMB/TMB binary
compatibility for a new run. Use a usable Pandoc path for vignette checks and
clear inherited incompatible locale variables as described by the Windows
validation skill. No toolchain change is made by this handoff.

For an implementation batch, run the verifier first, freeze a new source record,
then run the relevant focused tests above and update roxygen/help/examples if
behavior changes. The final dependency-complete candidate should run:

```r
devtools::document()
devtools::test()
devtools::check(args = "--no-manual")
```

Use a clean source snapshot without .git if OneDrive paths interfere. Keep
vignettes enabled for the final check when they change. A future authorized
commit must receive a single 0.0.1 version bump and matching NEWS entry after
validation. Do not stage the entire dirty checkout mechanically.

## Evidence index and precedence

Paths below are relative to the repository. Follow the latest final checkpoint,
not a historical intermediate instruction embedded in an older RESUME file.

* Common exact check: tools/benchmarks/joint_scale_assessment/phase63/RESUME.md,
  REPORT.md, check-release-2/checked-source-md5.csv, check-summary.txt and
  installed-tests.log.
* Shared block/spatial closeout: tools/benchmarks/tmb_marginal_likelihood/
  latent_gaussian_20261028/RESUME.md, REPORT.md, PROTOTYPE.md,
  check-reuse-audit.json, mathematical-gates.json and final-pilot-summary.json.
* Posterior prediction calibration: tools/benchmarks/bayesian_scales/
  predictive_20261028/RESUME.md, REPORT.md, INTEGRATION_PLAN.md and
  ARTIFACT_MANIFEST.csv. Its read-only verify-manifest.py audits the larger
  closed study. Do not run its writers or cleanup scripts after closure.
* Shared density and transport correction: tools/benchmarks/
  tmb_marginal_likelihood/consolidation_20261028/RESUME.md, REPORT.md and
  INTEGRATION_STRATEGY.md. Later preparation/capabilities and the combined
  package check supersede its earlier pending architecture/check priorities.
* Neutral preparation and capability implementation: tools/benchmarks/
  tmb_marginal_likelihood/preparation_20261028 and capabilities_20261028.
* Origins: Codex chats Marginal Likelihood (01a170db-7999-7c91-b92b-7344f026541c),
  Integrated likelihood (01a1595a-b988-7b31-aaa2-4259513b41d0), and
  Assess Bayesian or Laplace Flex (01a166f8-e2c8-7640-9483-57f3226de9c0).
  This consolidation read their latest turns; it did not send messages to them.

# Shared Gaussian layer and spatially varying scale prototype

October 28, 2026. Development batch prepared before the 3:30 PM Eastern logoff.
No staging, commit, push, merge or version change. The existing dirty development
checkout and concurrent work are preserved.

## Decision

The shared latent-Gaussian representation and its private spatial-scale consumer
are implemented and covered by mathematical, regression and package checks.
Existing public interfaces and numerical diagnostics are preserved. This is a
substantive model extension: ecological kernel scale can now vary smoothly by
location while retaining independent annual effects.

The spatial prototype is **not ready for a public inference or nonstationarity
detection interface**. The two illustrations expose substantial disagreement
between joint Laplace and NUTS uncertainty. Neither original NUTS run passes every
reference gate, and one predeclared dense-metric null repeat also fails. All
failures are retained. Successful optimization and positive curvature establish
numerical eligibility, not approximation accuracy or interval calibration.

## Implemented model and integration boundary

`R/latent_gaussian.R` represents a proper finite Gaussian coefficient block with
fixed projection matrix A and fixed positive-definite coefficient covariance C.
With C = LL', its observation effect is tau A L z, z ~ N(0,I). An independent
centered density includes the determinant and scale normalizers. Singular,
indefinite and ill-conditioned covariance matrices are rejected without repair.
Scale and response components share this representation. Existing indexed IID
effects retain their original arithmetic.

`R/spatial_scale_prototype.R` adds one fixed-range radial-basis field to a selected
predictor's log kernel scale. In schematic form,

    log s_i = log(unit conversion) + gamma_0
              + tau_year z_year[i] + tau_field b(location_i)' L z_field.

This models a varying ecological distance scale, rather than a varying landscape
response coefficient. The field smoothing range is distinct from kernel scale.
Coordinates, knots and range are supplied explicitly. Coordinates must be in a
known projected system; numeric columns do not carry verifiable CRS metadata.
The first prototype uses an intercept-only fixed scale mean with optional
existing annual intercepts. Ordinary public fixed scale formulas remain intact.

The basis is projected off the training intercept and normalized to unit mean
pointwise prior variance over training rows. Consequently tau_field is the RMS
pointwise prior log-scale SD on those rows, not the realized surface SD. Training
centering and normalization are frozen at prediction. Shared coefficients induce
covariance across prediction locations; this finite basis adds no independent
new-location innovation. Extrapolation follows the selected basis.

The posterior target retains existing proper priors and requires an explicit
half-normal field-SD prior. The likelihood target includes normalized Gaussian
latent densities but no population priors. Separate target checks prevent
silently treating a posterior as an unpenalized likelihood.

Available private operations are preparation, joint posterior Laplace fitting,
joint/latent-Laplace AD objectives, and dedicated spatial scale prediction. The
likelihood adapter is **not** a qualified maximum-likelihood fitter, confidence
interval method or limiting-model witness. The separate prototype class and
explicit guards prevent accidental use by existing scalar integrators and public
prediction methods. Public formals, defaults and diagnostic thresholds are
unchanged. No additional package dependency was introduced.

The block representation can support future random-slope projection matrices and
proper temporal/spatial covariance structures. This batch does not implement
their public syntax, covariance estimation, intrinsic fields, additional response
families or general graph factorization. A shared spatial field couples sites and
years, so the independent-year scalar factorization cannot simply be reused.
Joint or block integration must respect that dependence.

## Validation and exact source scope

The focused mathematical run passed 111 assertions: 90 new block/prototype/
integration assertions and 21 existing simulation/selection assertions. It had
zero failures, errors, warnings or skips. It checks normalizers and Jacobians,
proper covariance handling, independent derivative references, frozen prediction,
coherent new annual effects, target/prior separation and unsupported routes.

An exact limiting test sets the landscape response coefficient to zero. The
response then contains no scale-field information. Latent integration returns
the correct normalizer and a flat field-SD likelihood, rather than manufacturing
curvature. Eight predecessor cases, covering both kernels and all four existing
response densities, have identical ordinary components and zero differences in
AD values, gradients and Hessians. A separate focused regression run passed
1,639 assertions in 27 existing files, with no failures, warnings or skips.
These counts overlap package testing and should not be summed as independent
evidence.

The completed full package check has **0 check errors, 0 check warnings and
0 notes**; examples, donttest examples and rebuilt vignettes passed. Its installed
tests report **4,763 passes, 0 failures, 12 test warnings and 10 skips**. The test
warnings concern existing CRS fixtures and surface/refit tests. All skips concern
development benchmark helpers or policy files excluded from the built package.
The skipped benchmark-policy helpers are not counted as validated by this package
check. Separate focused Bayesian/factorized checks are reported above.

This full check was performed by concurrent same-package work. Rather than start
a duplicate, the batch adopted its actual checked source after a byte audit:
273 shared files match, including all model R/native source, tests, help,
source vignettes, runtime inputs, DESCRIPTION, NAMESPACE and build exclusions.
The 32 external-only files match existing build exclusions. Our 21 additional
files are 12 prebuilt vignette outputs and 9 unchanged baseline benchmark
auxiliaries. There is **no archive-byte identity claim**. The extra leading
validation library contains glmmTMB; the Bayesian dependencies use the same
private library. Copied logs, manifest and the exact reuse audit are retained.

The immutable numerical engine is
`C:/Users/peterman.73/.codex/scratch/msr-spatial-20261028-source`.
`final-source-audit.json` confirms all 273 shared files still match the live
checkout and frozen engine, with no missing or changed files at the audit.
Every completed pilot reports unchanged source identity.

The earlier own baseline check is retained as a failure: its strict export list
omitted the concurrently documented/exported `prepare_unmarked_data()`. Only
that exact expected name was added; the strict set-equality test remains. Its
4,605 passes do not turn the failed baseline into a pass. A previous scheduling
interruption and a later unavailable resource-inventory exit are also retained.
No foreign process was stopped. See OPERATIONS.md and CHECK_REUSE_PROTOCOL.md.

## Fixed same-input spatial illustrations

PILOT_PROTOCOL.md froze two cases before fitting: 24 sites, 3 years, 72 Gaussian
observations, a rank-3 field, range 350, population kernel scale 220, annual
log-scale SD 0.15, and spatial SD either 0 or 0.45. The explicit new field-SD
prior scale is 0.4. Inputs use a synthetic binned fixture, rather than realistic
overlapping raster neighborhoods. Laplace and NUTS receive exactly the same input
within each case. No archived annual study or scientific fit was rerun.

Laplace uses six starts and 2,000 draws. NUTS uses four serial chains of 8,000
iterations with 2,000 warmup, adapt_delta 0.99 and depth 12. Original gates remain
Rhat <= 1.01, bulk/tail ESS >= 400, zero divergences/depth hits and E-BFMI >= 0.3.
All 72 site-year log scales, both SDs and two regional contrasts also receive
derived diagnostics. Functional CDF MCSE must be <= 0.005 at lower, median and
upper thresholds. Constant indicators are unavailable, not zero-error estimates.

| Case | Joint Laplace seconds | NUTS seconds | Divergences | Derived gates | CDF precision failures | Reference qualified |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| Constant spatial field | 0.72 | 277.95 | 18 | Pass | 2 | No |
| Smooth spatial field | 0.54 | 423.80 | 0 | Pass | 1 | No |
| Constant field, dense-metric repeat | Original fit reused | 241.22 | 23 | Pass | 2 | No |

Times are the numerical fitting/sampling sections in a loaded process; setup and
compilation are excluded. End-to-end case times were 285.30, 430.66 and 248.55
seconds. Numerical speed alone does not make Laplace an equivalent substitute.
The dense repeat verifies the original null input SHA256 exactly, retains all
original thresholds and uses a new seed. It does not repair the failure.

The original constant-field precision failures are the annual-SD median CDF
(MCSE 0.005171) and spatial-SD median CDF (0.005240). The smooth case fails the
spatial-SD median CDF (0.005428). The dense repeat fails the corresponding two
null CDF checks (0.005288 and 0.005314). Small excesses remain failures.

The original Laplace comparisons flag 204/228 constant-field CDF checks and
194/228 smooth-field checks, with 2 and 1 unavailable comparisons respectively.
The dense repeat flags 205/228 at the original frozen thresholds, with 2
unavailable. Flags use the predeclared three-combined-MCSE rule. These correlated
pointwise checks are approximation diagnostics, not independent hypothesis tests
or a coverage estimate. The reference qualifications above remain essential.

The medians of site-year 95% log-scale interval width ratios, Laplace/NUTS, are
3.29 and 4.11 for the original two cases. In the constant-field illustration,
NUTS estimates spatial SD median 0.099 with interval [0.006, 0.407], while joint
Laplace gives 0.445 [0.147, 1.404]. In the smooth illustration, the corresponding
values are 0.300 [0.153, 0.693] and 0.531 [0.222, 1.252], for truth 0.45.
This is a material uncertainty discrepancy even though maps of medians look
similar. A Gaussian approximation in noncentered log-SD coordinates need not
capture the SD marginal distributions; positive local curvature does not certify
that approximation. The precise geometry causing the NUTS divergences remains
unresolved.

Descriptive log-scale median RMSE is 0.120 versus 0.136 (NUTS versus Laplace) in
the constant-field case and 0.063 versus 0.071 in the smooth case. These are two
illustrations with unqualified references, not evidence of average recovery,
calibration, false-positive control or superiority to constant-scale models.

An additional cold affine-metric audit stopped before sampling because the
posterior-mean Hessian had minimum eigenvalue -1.414376. No covariance was repaired
and no affine sampler was launched. This rejects that particular preconditioner;
it does not prove multimodality or invalidate positive curvature at the MAP.
The final dense-metric diagnostic was declared separately and exhausted the
bounded follow-up. No further sampler configurations were attempted.

Figures were generated only from retained compact summaries, then inspected:

![First-year scale medians, with provisional reference labels](<C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/multiScaleR/tools/benchmarks/tmb_marginal_likelihood/latent_gaussian_20261028/spatial-pilot-surfaces.png>)

![SD medians and intervals, with provisional reference labels](<C:/Users/peterman.73/OneDrive - The Ohio State University/R/Packages/multiScaleR/tools/benchmarks/tmb_marginal_likelihood/latent_gaussian_20261028/spatial-pilot-sd-intervals.png>)

Site coordinates were reordered explicitly to account for the fixture's
lexicographic factor levels. Shared map colors and annotated values represent
distance medians. The smooth-looking null map illustrates why visual patterns
alone cannot establish nonstationarity. A positive equal-tail field-SD interval
under a continuous positive prior does not reject a zero-field null.

## When pooling helps, and when it cannot

Annual and spatial pooling offer a model for borrowing information among
exchangeable years or nearby locations. It is useful when the landscape response
is identifiable, neighborhoods contain scale-sensitive contrast, and sparse
years/locations share a scientifically defensible structure. Bayesian inference
also propagates uncertainty in the shrinkage SDs through scale maps and regional
contrasts. These are model capabilities, not a new calibrated performance claim
from the two pilots.

Pooling cannot supply missing scale information when the landscape coefficient
is zero or the neighborhood response is insensitive to scale. Priors and the
basis can still generate smooth maps in weak data. Response-only spatial
structure, omitted variables, spatial sampling bias, fixed range choice and
limited extracted neighborhoods may confound interpretation. A prior-sensitive
surface is not proof that ecological scale varies geographically. Spatially
blocked prediction and null/weak/confounded simulations are required.

Full Bayes is flexible but is not the only mathematical route for this model.
The shared normalized likelihood also admits joint or block Laplace and other
multivariate integration. The present independent-year scalar algorithm is not
valid for the coupled spatial field. Whether a useful frequentist approximation
is faster and sufficiently accurate remains an open development question.

## Next priorities and stopping boundary

1. Obtain a qualified spatial reference under the unchanged target and gates.
   Diagnose the null SD/latent geometry with a declared design and an independent
   reference; do not extend chains or change metrics opportunistically. Existing
   exact-response integration research is not assumed to improve information per
   second, and should not be promoted merely because dimension decreases.
2. Develop conditional/block marginal integration of the shared field and SDs,
   including coefficient reconstruction for prediction. Test conditional values,
   normalizers, derivatives, zero-SD limits and resulting marginal uncertainty
   before treating it as an improved approximation.
3. For likelihood fitting, add explicit zero-variance faces, cold stationarity
   and limiting-model witnesses, identification checks and a justified uncertainty
   target. The current adapter supplies none of those fit-level guarantees.
4. Freeze broader spatial-null, weak-effect, response-confounding and blocked-
   holdout experiments before exposing a public constructor. Only then expand
   random slopes, temporal covariance or estimated spatial range with independent
   covariance/normalizer tests. General factorization must describe conditional
   dependence, rather than assume annual independence.

All prior annual-study failures and protected decisions remain untouched,
including sbc-123/sbc-187 limiting witnesses and prior-018 competing rejection.
No original diagnostic cutoff or prior was weakened to obtain a pass.

## Disk, ownership and reproducibility

At the report inventory, available RAM was 18.72 GiB and disk 53.43 GiB. Numerical
jobs used one thread and ran serially through an exclusive owned lock, with
known-inventory admission gates. No raw samplers, draw arrays or AD tapes were
serialized. Original seeds, input hashes, warnings, compact summaries and failed
outcomes are retained. No repeated numerical design overwrote an original fit.

The study directory was about 0.78 MiB before this report. The own baseline source,
baseline check and prototype engine trees occupied approximately 31.42, 28.00 and
31.48 MiB including regenerable compilation/check products. Their source-only
freezes were about 8.3 MiB each. No duplicate full prototype check was created.

Retain the source manifests, producer scripts, protocols, original failures,
mathematical/scientific results, check logs and reports. Compiled objects/DLLs,
installed-check trees and generated vignette outputs are regenerable cleanup
candidates after exact ownership and no-loader verification, but they were not
deleted in this batch. Other owners' check sources and jobs remain untouched.
Previously policy-blocked DLL deletions were not retried. Current free space
does not justify risking the provenance retained here.

`RESUME.md` records the final process/source audit and safe next steps. No queue
or sampler should be restarted automatically. The local deadline watchdog is
identity-stopped once all owned jobs are complete; no heartbeat was created or
modified by this batch.

Final process audit at 2:42 PM Eastern confirmed no owned numerical job or R
interpreter remained. The watcher was identity-verified, stopped and confirmed
absent. Available RAM was 18.66 GiB and disk remained 53.43 GiB. HEAD and version
remain unchanged, and the Git index is empty. `change-inventory.json` identifies
this batch's source/help/test files without attributing the checkout's broader
dirty state to this work.

For revalidation after future source changes, use the configured Windows R and
Pandoc environment, refresh resource/process admission, and run:

```r
devtools::document()
devtools::test()
devtools::check(args = "--no-manual")
```

These are future validation commands, not a request to repeat the completed
frozen check or pilots. Preserve the existing check evidence when source changes.

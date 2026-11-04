# Inference development contracts

Keep native refitting, analytic likelihood development, and Bayesian development
on `inference-dev`. Share preparation and prediction contracts while retaining
each method's likelihood, prior, diagnostics, and uncertainty target. A common
interface does not establish scientific calibration or authorize an experimental
backend for another model class.

## Fixed-kernel radius uncertainty

With `kernel = "fixed"`, each estimated kernel scale is a cutoff radius. Raster
support changes as cells enter or leave that radius. Kernel means and weighted
surface metrics retain their existing formulas and training normalization.

The native refit omits radius Hessian SEs and Wald limits, following the same
contract as landscape and unweighted surface radii. The raw optimizer Hessian,
point estimates, convergence results, and native coefficient covariance remain
available. Coefficient SEs remain conditional on the selected spatial parameters.
Requested profile endpoints retain their method and boundary information;
profiling failures remain explicit.

Fitted summaries, distance tables, and kernel plots show the actual radius in
both dimensional definitions. Old saved fits cannot regain Wald limits through
legacy finite SEs. Standalone `kernel_dist(sigma = ..., kernel = "fixed")` and
`plot_kernel()` continue to describe uniform-kernel weight quantiles.

Count every estimated radius once in the likelihood penalty. Keep the native
model's coefficient, dispersion, covariance, or effective-degree-of-freedom
count. Preset radii and inapplicable shapes add no estimated spatial parameter.
Unavailable radius uncertainty does not make an estimated radius a known input.

The regression contract is in
[test-fixed-kernel-radius.R](../../../../tests/testthat/test-fixed-kernel-radius.R).
It covers old objects, mean and weighted-SQ fits, coefficient covariance,
profile endpoints and failures, plots, and AIC/BIC penalties. Smooth-kernel
curvature failures retain their existing diagnostics and missing-SE behavior.

## Prediction and native model behavior

Preserve source names, metric definitions, physical scales, row/site alignment,
and the original reference normalization. Reuse retained transformations and
factor/logical coding rather than estimating them from prediction rows.

Native GAM/BAM predictions retain their supported output types, term selections,
matrix geometry, covariance decisions, and warnings. Single-predictor models
retain their inherited response-curve behavior. Multi-component GAMs require an
explicit component where applicable and retain their distinct link/output
conventions. Conditional native SEs do not include ecological-scale uncertainty,
observation noise, or an entire Bayesian posterior.

The native binomial/Poisson BAM regression fixture covers ordinary/discrete fits
over Gaussian, exponential, exponential-power, and fixed kernels. Keep native
trial initialization, convergence, penalized score, positive information, and
covariance assertions. Model-selection tests retain native family labels,
observation counts, and spatial-parameter penalties.

Gamma and `mgcv::nb()` BAM checks extend points, numerical snapshots, maps,
raw-reference curves, ordinary workers, and model selection. Estimated-theta
cases and Gamma cases cover all four kernels; fixed-theta cases use Gaussian
kernels. Retain the native convergence, score, information, and covariance
gates. Tighter fitting controls used by declared software fixtures do not change
package fitting defaults or establish ecological-scale calibration.

For verified native NB BAM unit-scale bookkeeping discrepancies, selection
tables use the normalized count density and actual effective-degree-of-freedom
count. Estimated theta counts once; fixed theta and the known unit scale do
not add estimated parameters. Record native and used values in the table's
`likelihood_adjustments` attribute. Other model counts and all fitted objects
retain their existing definitions. Unknown discrepancies fail explicitly.

The uncommitted weighted-NB batch covers known-scale legacy `mgcv::negbin`
with explicit `scale = 1`, plus fixed/estimated-theta extended NB with zero
and fractional case weights. These eight Gaussian-kernel fixtures retain the
complete original reference rows, native covariance and observation counts.
They do not reinterpret the legacy default estimated working scale. Large-weight
rounding verification checks canonical native AIC arithmetic with explicit theta
and rescaled density at the unchanged tolerance, reporting the original weighted
target and actual count. Inconsistent getter/AIC state remains an error.

Effect-curve verification keeps its complete native basis comparison and
1e-10 tolerance. Large exact finite ordinary numeric matches can avoid temporary
comparison copies. Every difference, custom class or non-native comparison method
retains original comparison behavior; no native prediction batch is split.
Verified matrices are released before link prediction. Further blockwise basis
construction still requires separate equivalence and failure tests.

Numeric 0/1 Bernoulli responses use binomial/logit, with success probability on
the response scale. Conditional GLM likelihood comparisons require the same
selected predictors, rows, design, offset, encoding and nuisance conventions.
Integrated random effects and Bayesian population priors retain their distinct
targets. Logical FALSE/TRUE scalar outcomes and validated logical count columns now
use the same zero/one likelihood target in profile and shared preparation.
Native data/model-frame types, trial totals and priors are retained. Selection
signatures compare supported logical outcomes as zero/one copies, preserving
original fitted likelihoods, covariance and parameter penalties. Changed
outcome/trial/frequency-weight targets still fail explicitly. Compare
native covariance within the same response representation; equal likelihoods
across encodings do not replace native IRLS covariance. Additional links
and weighted-proportion admissions remain separate, explicit audits.

Pcount model comparisons recover actual used visits from native component
formulas, checking row layout and stored omitted-site metadata. Different
missing visits are incompatible despite equal raw counts/site totals. Native
marginal likelihoods and free variance/mixture penalties remain unchanged;
realized latent effects are excluded. Site-row AICc/BIC conventions do not
certify repeated-group small-sample assumptions, fit quality or truncation
stability. Default siteCovs=NULL frame counts and identities are covered by full83:
checked native numSites/sampleSize, original response IDs or original site
positions, and validated omitted-site indices. Null abundance/detection-only
models do not bypass pcount actual-visit checks. This shared ID/count repair
does not qualify other unmarked classes' formula-specific visit masks.

Stage1 cross_validate() refits GLM coefficients and scales on explicit training
rows only. Representative distance bins and spatial/fixed prediction
transforms are learned from training buffers and frozen for assessment.
The predeclared physical buffer domain is unchanged. Whole planned
assessment occurrences, count constants, frequency/trial distinctions and
failed folds remain visible; any failed fold makes pooled scores NA.
Gaussian identity, Poisson log and binomial logit scores are plug-in
targets, not new uncertainty/calibration evidence. Existing backend
starts/cold/curvature/admission gates apply. Raw cell inputs are required;
new fold generators, additional model classes and predictive integration
need separate qualification. CV selection performance needs outer
validation or an independent test set. Full fold models are opt-in.

## Private inference gates

Keep existing analytic/Bayesian kernel and family admissions. The fixed-kernel
reporting correction does not admit a discontinuous radius to a smooth analytic
or Bayesian fitting backend. Preserve every prior and normalized target.

Private weighted-SQ and spatial-scale work retains its structural admission
guard, conservative AD-memory prerequisites, starts, cold checks, zero/limiting
cases, and curvature requirements. The spatial-only preparation has no annual
latent coordinates; structural software checks do not establish field-SD or
scale-surface recovery. Retain null failures and approximation qualifications.

## Continued validation and refinement

Optimized discrete BAM points, snapshots, and maps have regression coverage.
For further family or control extensions, compare the actual native provider on the exact
same batch of rows. An ordinary prediction target cannot substitute for an
optimized one. Prove raster instrumentation preserves the unmodified public
route, and retain warning identities, geometry, and missing-value positions.
Build any new audit from the accepted release source and commit identity.

Logical/numeric binomial selection compatibility is covered by full81;
mixed pcount actual-visit and default-frame contracts are covered by
full82/full83. Next review native used-visit semantics of other unmarked
classes with a declared zero-fit check, then trial/weight encoding against
the native GLM provider. Keep frequency weights and grouped-binomial
constants distinct.
Derive native row-batch semantics before replacing full basis matrices with
blocks. Keep private fitting admissions and every scientific numerical gate
unchanged.

Use focused regression tests and generated documentation, followed by a normal
vignette-enabled package check and the complete installed suite from that exact
archive. Preserve old scientific gates and all failed/interrupted attempts.
Assertion counts record coverage; they are not development targets. Remove
redundancy only after demonstrating that unique behavior and failure coverage
remain intact. Never prune a failing test to obtain a passing check.

Inspect live source, actual R processes, memory, and disk before numerical work.
Use one numerical thread, at most two owned fit workers and three owned R
interpreters including drivers. Leave foreign work untouched. Keep compact
atomic records, avoid full fit/chain/tape/array checkpoints, and omit build-ignored
legacy example fits from temporary source copies. Publish only the reviewed
package source, tests, documentation, and concise development guidance.

## CV publication and next fold contracts

Version 0.7.14 publishes explicit-fold canonical GLM CV with training-only
scale/coefficient fitting and frozen training preprocessing. Future group or
spatial generators must return the same disjoint analysis/assessment ID lists.
CV score comparisons must match planned targets, splits, responses, trials
and weights and retain failures. Predictive scores are not AIC weights or
additional parameter penalties. Further model-class scoring requires separate
predictive-target qualification; all existing private admission guards remain.

## Qualified CV folds and fixed map geometry

Grouped/spatial generators return explicit whole-unit row-ID splits. Spatial
buffers exclude complete training units and use projected point distances. CV
refits scales, coefficients and learned transforms within training rows. Native
grouped rasters select fixed input maps; their row labels are validated and
subset with their raw buffers. Neither map labels nor pixel geometry add fitted
parameters, annual effects or spatially varying scales.

CV bin domains use declared extraction radius plus half the fixed cell diagonal
when known. The original distance tolerance stays unchanged, all intersecting
cells remain, and no held-out maximum selects the domain. Legacy unknown grid
geometry keeps its strict old gate. Freeze geometry with training bin references.

cv_tab matches splits/outcomes/trials/weights and ranks compatible held-out
losses without penalties or Akaike weights. select_cv uses only retained native
training fits and their full K/likelihood guards. Preserve exact ties/failures,
selected row/map records, candidate settings and the inherited-controls label.
Point scoring does not integrate parameter/selection uncertainty, add intervals
or establish calibration. Matched-set, mixed, GAM, unmarked and Bayesian scoring
need separate target qualification; never substitute zero random effects for
integration or posterior objectives for maximized likelihoods.

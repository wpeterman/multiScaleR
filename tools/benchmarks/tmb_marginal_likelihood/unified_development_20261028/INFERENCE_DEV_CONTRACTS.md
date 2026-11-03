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

The next prepared native prediction audit concerns optimized discrete BAM
points, snapshots, and maps. Compare the actual native provider on the exact
same batch of rows. An ordinary prediction target cannot substitute for an
optimized one. Prove raster instrumentation preserves the unmodified public
route, and retain warning identities, geometry, and missing-value positions.
Build any new audit from the accepted release source and commit identity.

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

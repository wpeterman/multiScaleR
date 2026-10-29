# Bounded Gaussian response-SD profile trace

Declared before new implementation or numerical fitting. This is the next
small marginal-likelihood uncertainty primitive on codex/inference-dev. It is
private, produces a diagnostic trace rather than a confidence interval, and
does not change any existing fitter, prior, public interface or numerical gate.

## Scope and target

Use the existing normalized unpenalized marginal-likelihood specification:
Gaussian response, one kernel with a shared intercept-only scale, one indexed
Gaussian response random intercept, no scale random effects or private Gaussian
projection blocks. Reject posterior priors, other families, coupled/spatial
blocks and unsupported structures before constructing an objective. This first
trace concerns the response-group SD in response/linear-predictor units, not an
annual kernel log-scale SD or the spatial field amplitude. Those remain later
extensions. The target is full ML, not REML or a posterior callback.

Proposed private function .msr_factorized_response_sd_profile(spec, sd,
control=list(), initial_parameters=list()). The caller supplies distinct finite
nonnegative natural-SD values and must include exact zero. Preserve requested
order. Positive SDs are fixed, not optimized; zero uses .msr_marginal_reduce to
omit the response latent/SD block exactly. No epsilon floor or repaired Hessian.
Optimize the remaining population/nuisance parameters using the existing
factorized objective, retaining the full normalized Gaussian latent integral.

## Conditional numerical contract

Each grid point has all declared deterministic starts (default three) and every
supplied warm start, with failures preserved. Keep names/indices when assembling
full vectors. At positive SD, remove only the fixed log_response_sd coordinate
from nuisance optimization, and always restore its fixed value in fn/gr/he calls.
Re-evaluate the accepted vector before extracting parameters; rejected trials
must not contaminate TMB state. Do not return AD objects/tapes in the trace.

Retain original tolerances: nuisance/cold score <=1e-4; optimizing/cold/ordinary
objective agreement <=1e-6; gradient agreement <=1e-4; unmodified positive-definite
nuisance Hessian with condition <=1e12; existing Newton polishing certificate.
Record the full gradient and curvature separately. A constrained nuisance
optimum can have nonzero score in the fixed-SD direction and is not a global
MLE. Never label a trace point globally eligible or reuse it as an accepted fit.
The original full-parameter fitter's gates remain completely unchanged.

Preserve lower ineligible candidates and competing eligible nuisance optima
using the existing objective/curvature-distance rules. Every requested point
remains in the output, including rejected/unavailable points. No interpolation,
discarding failed starts or selecting a favorable grid subset.

Required broad/nearest kernel witnesses must hold the SAME fixed response SD.
Do not call the existing unrestricted limit audit unchanged: it reoptimizes that
SD and would compare a different conditional target. Build a small conditional
witness adapter with all three starts, cold/ordinary value checks and retained
curvature diagnostics. A feasible audited limiting value can reject finite-kernel
identification even when its Hessian is singular; it is an exclusion witness,
never an accepted finite fit. Preserve every failure and require both witnesses
to complete. Do not regularize curvature or loosen the existing 1e-5 witness rule.

Reuse the existing Gaussian singleton/unbounded-response identification guards.
Return a distinct private trace class, per-point attempts/parameters/diagnostics
and numeric curve. Any relative objective must say "relative to the best
requested valid point", not a global likelihood ratio. No confidence level,
Wilks cutoff, p-value, interval endpoints, predictive interval or calibration
claim is exposed. Boundary-aware confidence procedures require a later declared
bootstrap/test-inversion design. All existing zero-SD/competing/limiting decisions
and prior studies remain untouched.

## Declared deterministic software fixture and resource bound

Use bayesian_scale_fixture(n_site=6,n_year=5,seed=713403) only to construct fixed
binned landscapes/covariates. Gaussian formula count~habitat+temp+(1|year)+
offset(log(effort)), shared habitat scale ~1, reference/kernel scale200 map units.
For a deterministic informative response use beta(.5,.25,2), response SD.35,
fixed year z(-1,-.3,0,.4,.9), residual SD.2 and a permuted fixed normal-quantile
residual vector. This is a software fixture, not a calibration dataset/panel.
The declared trace grid is c(0,.2,.5), default three starts plus one supplied
true-parameter warm start. Retain all attempts. No sampler, new scientific
replicate, seed search or hyperparameter-prior study is declared.

Bound each focused numerical validation phase to600seconds and use one owned
worker, one thread, at most four actual R interpreters globally, known identity,
RAM>12GiB/disk>5GiB. Unknown/foreign inventory blocks admission; leave foreign
work untouched. Full package checks require an R-free uncontended window.
No long launch after04:00Eastern; stop launches at first wake>=04:15, checkpoint
by04:35 and delete the one-night heartbeat by04:45. No staging/commit/version/
merge/push. Preserve source/evidence identities and every failed attempt.

## Necessary assertions before promotion

* At each fixed SD, independently evaluate the normalized Gaussian covariance
  density with Sigma=sigma^2 I + sd^2 Z Z' using base-R Cholesky/solve, including
  exact zero and its small-positive limiting evaluation at the same nuisance
  vector. Compare objective, nuisance finite-difference score and the analytic
  beta-Hessian block X' Sigma^-1 X. Preserve offsets/training normalization.
* Verify fixed SD and coordinate maps in every candidate and limiting witness;
  retain full-direction scores rather than force them to zero. Check all starts,
  candidate rejections, failed/competing optima and incomplete limiting audits.
* Reject unsupported targets/structures and identification failures early; never
  attach a posterior prior, silently remove groups or manufacture a likelihood CI.
* Confirm existing fitter/helper source identities and public formals unchanged.
  Add focused tests and generated internal Rd. Keep NAMESPACE unchanged. Then
  build a fresh dependency-complete source and run the complete package gate.

If this bounded primitive cannot be validated, retain the implementation/failures
as a private unresolved prototype and checkpoint exact remaining gates; do not
promote it, weaken tests or claim uncertainty coverage.

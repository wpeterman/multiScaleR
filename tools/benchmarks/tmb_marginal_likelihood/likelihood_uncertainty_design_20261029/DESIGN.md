# Boundary-aware likelihood uncertainty: next development design

Design only, October29,2026. No bootstrap, simulation cohort, sampler, interval
implementation or numerical producer is declared by this document. Finish the
annual-SD candidate's current gate first. Keep the package source frozen afterward;
use the remaining night for static review and a safe checkpoint. Calibration needs
a separately declared precision/resource budget and coverage study. A few valid
profile points cannot supply that evidence within the remaining window.

## Target and first scope

Retain normalized unpenalized full marginal likelihood. Proper latent Gaussian
densities belong to the likelihood; population-parameter priors do not. Bayesian
method="marginal" keeps its proper population priors and is not this target.
Native/profile likelihoods retain their own conditioning and are not substituted
for the private hierarchical likelihood.

Begin with either of two already scoped Gaussian models, separately:

* One response random-intercept SD and a shared intercept-only kernel scale,
  with no annual scale group.
* One independent kernel log-scale group SD, with no response random effect.

Do not expand simultaneously to multiple SDs, crossed response groups, spatial
fields/range, non-Gaussian response integration or temporal dependence. Response
SD uses linear-predictor units; scale-group SD uses log kernel-scale deviations.
The ecological kernel distance, its covariance range in a future spatial model,
and either SD are different quantities.

## The required statistic is not the present curve

For an SD tau>=0 and nuisance parameters psi, let

    P(tau) = min_psi -log L(psi,tau)
    M      = min_(psi,tau>=0) -log L(psi,tau)
    T(tau) = 2 * (P(tau) - M).

The current private curves report objective differences from the best requested
valid grid point. That value is not M. A constrained point can retain a nonzero
fixed-SD score and an indefinite full Hessian while its nuisance block is valid.
It is not an eligible unrestricted fit or a likelihood confidence interval.

Before forming T, obtain a qualified unrestricted full-parameter candidate with
the original exact-zero faces, starts, cold scores, curvature and limiting gates.
The profile solution supplies a feasible full-model embedding and an additional
warm start, never a replacement for existing starts. Assert M<=P(tau). Retain
negative discrepancies with both ordinary objective evaluations and audit errors;
do not silently clamp a negative statistic or conceal a lower feasible candidate.
An unresolved unrestricted optimum makes the statistic unavailable.

Retain derivative coordinates. At a positive SD, the stored full fixed-direction
score uses the log-SD coordinate named in the parameter map; it is not a derivative
with respect to natural tau. At exact zero that coordinate is absent. A full
fixed-parameter Hessian entry is also not the reoptimized profile curvature.
Under regularity and a unique nuisance optimum, use eta=log(tau): profile
curvature in eta involves the Schur complement
H_eta,eta - H_eta,psi * H_psi,psi^-1 * H_psi,eta. Natural-SD derivatives additionally
require P_tau=P_eta/tau and P_tau,tau=(P_eta,eta-P_eta)/tau^2. Boundaries, competing
optima or failed numerical gates invalidate a routine curvature interpretation.
No Wald interval or cutoff follows from retaining these matrices. Posterior
curvature, with its priors and population-coordinate Jacobians, cannot replace
unpenalized likelihood curvature.

Exact zero remains a model reduction, not log(epsilon). At an accepted zero face,
use the existing reduced-model eligibility rule rather than inventing positive
curvature in an absent coordinate. Interior unrestricted candidates keep all
original full-parameter gates. Profile candidates keep the conditional nuisance
gates and full-direction records. Neither rule can replace the other.

## Preserve the two limiting contracts

The response-SD profile must hold the same response SD in broad/nearest kernel
witnesses, because that covariance remains active. The annual-SD profile loses
its annual coordinate in those exact kernel limits: the normalized latent prior
integrates to one for every tau. Record that structural omission, not a fitted
zero SD or an artificial zero gradient/Hessian.

All original and supplied starts, lower ineligible candidates, competing optima,
cold/ordinary values and gradients, unmodified curvature and every failed limiting
audit remain. A feasible singular-curvature limit can exclude finite-kernel
identification but cannot become an accepted finite fit. An unavailable or rejected
observed-data fit is not repaired by choosing a confidence cutoff.

## Candidate calibration procedure to validate later

At each tested null tau0, condition on the actual extracted landscape summaries,
fixed covariates, row/group IDs, offsets and the frozen training normalization.
Fit the nuisance parameters under that null. A proposed parametric bootstrap
would then generate Gaussian responses from those fitted parameters, refitting
both the null and unrestricted models to each same generated response vector.
This is a plug-in composite-null procedure, not a proven coverage result.

For response groups, draw one Gaussian group effect shared across that group's
rows, then observation noise from the response model. For annual scales, draw
one standard-normal scale effect per group, form the corresponding log kernel
scale, compute means from the retained binned landscapes and add response noise.
At tau0=0 omit the effect exactly. Do not simulate population parameters from
Bayesian priors, resample groups independently by row, re-extract new rasters,
change support/standardization, or turn a Gaussian response into a count model.

For each replicate retain the null fit, unrestricted fit and T under unchanged
gates. Predeclare tie handling, finite-replicate test calibration and the nominal
test level before launching. Boundary mass, weak identification and nuisance
estimation require validation; no generic chi-square or mixture cutoff is attached
by this design. NUTS agreement is not frequentist coverage truth for unpenalized ML.

Monte Carlo precision must be specified for the calibrated test or critical
quantile, in addition to numerical fitting accuracy. For an independently
simulated exceedance probability q, its binomial sampling error is approximately
sqrt(q*(1-q)/B). Precision of a critical quantile is a separate requirement. The
MCMC ESS/MCSE functions do not turn these bootstrap records into original chains.
Zero observed exceedances do not prove q=0 or zero numerical error. Keep such
outcomes explicit and require a declared finite-sample precision rule rather
than substituting a plug-in zero MCSE.
No replicate count, seed set, cohort or throughput claim is specified here.

## Invert tests without manufacturing endpoints

Confidence-set inversion requires a validated acceptance rule over tau, including
exact zero, and numerical uncertainty from both fitting and calibration. Preserve
all requested grid points and refine only under a declared rule. Do not assume
the set is connected or that there are exactly two crossings. A failed grid point
is unresolved, not a value to interpolate across. Acceptance at the largest tested
SD does not establish a finite upper endpoint or prove an unbounded endpoint.
Tail/parameter-range reasoning is required before either conclusion.

An interval for a population SD would not by itself provide joint prediction
uncertainty for response coefficients, residual dispersion, scales and new groups.
It cannot be relabeled as posterior draws, latent conditional quantiles or native
conditional standard errors. Existing prediction contracts remain unchanged.

## Failure denominators and evidence

Predeclare original input/seed identities and run slots before a producer exists.
Keep execution failures, numerical rejections, unbounded/limiting cases and
completed fits separate. Never replace a failed replicate with another seed,
discard a bad start, condition the primary denominator on successful diagnostics,
or rerun only favorable inputs. Any repair policy requires its own declared
method label and retains the original failure and paired outcomes.

Coverage studies must span declared nuisance settings and identification regimes,
including true zero, small positive SD, weak landscape response and finite-support
limitations. Report interval availability and failures alongside coverage among
available results; that conditional coverage is not unconditional calibration.
The policy for unresolved observed/bootstrap fits must be justified before public
promotion, not inferred after seeing the results.

Store compact atomic records with source/runtime/target identities, seed/input
keys, all start outcomes, selected parameters, ordinary/cold/coordinate errors,
curvature, zero-face and limiting decisions, statistic/calibration values and
explicit status. Save no sampler, AD tape, full chain or large candidate array.
Use unique checkpoint names and resume pending slots only; preserve all completed
and failed records. Verify live owned-process identities before any stop.

## Assertions before an implementation can be promoted

| Contract | Required evidence |
|---|---|
| Target | Normalized Gaussian densities/latent measure; no population priors or REML adjustment; native/Bayesian targets unchanged. |
| Profile and unrestricted fits | Same physical vectors and names, exact zero, all starts, feasible embedding, P>=M, separate nuisance/full stationarity and unmodified curvature. |
| Numerical integration | Cold/ordinary/AD values and scores, original tail/partition and large-SD coordinate gates, independent Gaussian covariance or scalar-integral references. |
| Limits and identification | Same-SD response limits versus SD-independent annual limits, all witnesses/failures and no accepted singular finite fit. |
| Generating mechanism | Shared group draws, correct response noise, fixed offsets/design/normalization/landscape support and no accidental population-prior sampling. |
| Calibration | Predeclared null/nuisance cohorts, seed/slot identities, precision for the reported test/quantile and complete failure denominators. |
| Set inversion | Exact-zero inclusion, multiple/unresolved crossings, open/unbounded limits and retained numerical/calibration uncertainty. |
| Compatibility | Existing response/annual profile behavior, priors, diagnostics, public formals, prediction dependence and uncertainty distinctions preserved; fresh focused/full code validation. |

Read R/factorized_scale_profile.R, R/factorized_scale_sd_profile.R,
R/factorized_profile_attempt.R and the original factorized fitter/limit/score/
refinement sources, plus their permanent tests. Use the current CONTINUE.md and
exact final-source-evidence records for source-qualified validation. Do not replay
closed scientific producers, refresh old baselines or mistake metadata for gates.

The one-night authority still ends04:45 Eastern. No long launches after04:00;
first wake>=04:15 stops launches, safe checkpoint by04:35 and heartbeat deletion
by04:45. This design declares no additional numerical study. The already declared
annual build/check/installed-test pipeline may complete under its existing guards.

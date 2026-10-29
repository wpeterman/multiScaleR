# Bounded annual kernel-scale SD profile trace

Declared before implementation or fitting. Start only after the shared inference
guide is closed. This is the next small likelihood-uncertainty primitive, shared
with the already validated private Gaussian response-SD trace. It supplies no
confidence cutoff, accepted global fit, posterior prior or calibration study.
No public export, existing interface change, stage, commit, bump, merge or push.

## Target and structural scope

One Gaussian response, one kernel predictor, an intercept-only fixed kernel
scale and one independent indexed scale-group intercept/SD. No response random
effects, additional scale groups, Gaussian projection/spatial blocks, temporal
correlation or population priors. The group can represent years, which remain
exchangeable; the SD is spread in log kernel scale, not response SD or a distance
radius. Use normalized full marginal likelihood, not REML or a posterior target.

Proposed private .msr_factorized_scale_sd_profile(spec, sd, control=list(),
initial_parameters=list()). Natural SD values must be distinct finite >=0,
include exact zero and preserve requested order. Positive log_scale_sd is fixed
in every objective/gradient/Hessian and nuisance audit. Zero removes the annual
latent/SD block exactly with .msr_marginal_reduce(spec,"scale"), without a floor.
Integrate independent scale effects with the existing accurate factorized
objective and its original integration/tail gates. Do not introduce new priors,
clipping, integration tolerances, Hessian repair or a fallback after rejection.

## Shared numerical contract and compatibility

Prefer a small common internal fixed-SD problem/attempt helper for the annual
and existing response trace. Preserve every existing response helper's formals,
output fields, values, warning/failure records and RNG behavior. Avoid refactoring
unrelated fitters or prediction. The response trace can delegate through wrappers
with its original interface. Its existing217assertion module must still pass.
At fixed named vectors, compare old frozen response helper fn/gr/he/parameter
extraction with the new wrapper exactly, using predecessor source only as a
read-only function reference. Do not invoke any closed study driver or replay a
scientific producer. Reuse frozen response trace records and independent covariance
tests as regression evidence. Any parity failure blocks promotion.

Each annual SD point retains all default3starts plus every supplied warm start,
including duplicates, warnings and failures. Retain nuisance/cold scores<=1e-4,
optimizing/cold/ordinary objective agreement<=1e-6, gradient agreement<=1e-4,
unmodified positive nuisance Hessian with condition<=1e12 and existing Newton
polishing. Re-evaluate the accepted full named vector before extracting blocks.
Keep full fixed-direction scores/Hessians even when their curvature differs from
the nuisance block. Never force the fixed SD score to zero or call a constrained
point a globally eligible fit. Preserve lower ineligible candidates (1e-5 rule)
and competing nuisance optima (existing2loglik/curvature-distance rule).

## The kernel-limit distinction must be proved

Source inspection of .msr_factorized_limit_spec establishes that broad/nearest
kernel limits remove the annual scale block and fixed scale coordinate. Unlike
the response-SD trace, these limits make annual SD identically irrelevant: the
normalized latent prior integrates to one at every requested SD. Do not leave a
spurious fixed log_scale_sd coordinate or silently reoptimize it. Explicitly
record its structural omission and why the same conditional target no longer
depends on that SD. No artificial zero-score/positive-curvature certificate.

Both limits still need all3starts, cold/ordinary audits and retained unmodified
curvature/failures. A feasible audited singular-curvature limit may exclude finite
kernel identification at the original1e-5 tolerance, but is never an accepted
finite fit. Incomplete limiting audits remain unavailable. Preserve every grid
point and all original zero/limiting/competing decisions. Relative objectives
refer only to the best requested valid point, never a global likelihood ratio.

## One deterministic software fixture

Use the same binned-landscape helper settings as the preceding response-profile
fixture: bayesian_scale_fixture(n_site=6,n_year=5,seed=713403), Gaussian kernel,
reference/kernel scale200 map units. Formula count~habitat+temp+offset(log(effort)),
scale formula habitat~1+(1|year), no response group. Set beta(.5,.25,2), internal
scale_beta log(.2), annual log-scale SD .35, annual z(-1,-.3,0,.4,.9), residual
SD .2 and the same permuted fixed normal-quantile residual vector used by the
response fixture. This is one declared deterministic software case, not a new
scientific simulation/calibration panel, sampler configuration or seed search.
Grid c(0,.2,.5), default3starts plus the supplied generating population vector.

Independently evaluate the normalized product of one-dimensional Gaussian-group
integrals using ordinary base-R integration and the existing response/kernel
components at held-fixed population vectors. Check objective and nuisance
finite-difference derivatives, exact-zero agreement with direct Gaussian density,
small-positive continuity at the same nuisance vector and fixed SD in every
candidate. Prove broad/nearest values are SD-independent after exact reduction.
Test unsupported structures/priors, names/dimensions, lower/competing candidates,
failed/incomplete limits and absence of intervals/global-fit eligibility. Do not
change the fixture, seeds or gates to obtain passing/available grid points.

## Validation and bounds

Freeze the predecessor source/evidence and compact old versions before edits.
All closed manifests remain immutable. Focused phases max600seconds, one owned
worker/thread, known identities, global<=4actual R (<=2owned fit workers),
RAM>12GiB/disk>5GiB. Unknown inventory blocks launch; leave foreign jobs untouched.
No raw chains, samplers, AD tapes or candidate arrays in checkpoints. Keep all
failed attempts. If accurate quadrature or parity cannot pass these fixed gates,
retain the private prototype as unresolved; do not weaken tests or claim coverage.

Generate internal Rd, run response/annual factorization and limiting compatibility,
then a fresh exact complete candidate check including unchanged archive tests
in an R-free window. This code change cannot borrow the guide's unchanged-code
coverage. No long launch after04:00 Eastern; first wake>=04:15 stops launches,
checkpoint by04:35 and delete heartbeat by04:45. After cutoff checkpoint/delete
only. Start no work whose bounded validation cannot fit the remaining window.

# Private Gaussian response-group SD trace

Implementation and focused validation are complete; documentation and a fresh
full package gate remain pending at 01:22 Eastern on October 29, 2026.

The trace holds natural response-group SD fixed while optimizing population
nuisance parameters in the normalized unpenalized full Gaussian marginal
likelihood. Exact zero removes the latent/SD block. Both kernel-limit witnesses
hold that same SD. Every default and supplied start, rejection, cold/ordinary
audit and unmodified curvature diagnostic remains explicit. Only nuisance
curvature is required at a constrained point; the full fixed-direction score
and Hessian remain recorded. Global-fit eligibility and interval availability
are always false. Existing fitters, priors, public interfaces and gates are
unchanged.

The single declared deterministic software fixture has three requested points
and four starts per point. All 12 nuisance starts pass. Both required limiting
models retain three starts at each point, and no limiting witness rejects the
selected finite kernel in this fixture. Objectives are 6.339741560234627 at
SD 0, -1.809600483787065 at .2 and -.277298925700718 at .5. Relative objectives
refer only to the best requested valid point. They do not supply confidence
endpoints, a Wilks cutoff or evidence of coverage.

unit04 passes 217 assertions with zero failures, errors, warnings or skips.
Independent base-R covariance likelihoods verify every audited objective;
finite differences verify all nuisance directions, and X' Sigma^-1 X verifies
the beta Hessian. Tests retain fixed SD in every candidate/witness, nonzero
fixed-direction scores, exact-zero/small-positive agreement, offsets, unchanged
input/RNG, early unsupported and unbounded/singleton rejection, lower ineligible
and competing candidates, incomplete audits and limiting exclusions.

Failures are retained: unit01 had two fixture assertions caused by floating-point
offset subtraction; unit02 failed custom-class JSON serialization. Neither
changed a scientific gate. unit03 passed 214 assertions. unit04 uses a base-R
test observer for external compact checkpointing, so no package dependency is
introduced. The numeric checkpoint contains no AD objects or raw chains.

The previous clean check and 5,044 passing installed assertions belong to the
closed probability_precision_20261029 candidate. They do not validate this new
source. This trace remains private and has no uncertainty promotion pending
the fresh candidate check and a later boundary-aware interval design.

At 01:25 Eastern, documentation01 is complete: generated internal Rd passes
checkRd and NAMESPACE is byte-identical. compatibility01 passes 429 assertions
across Gaussian response, controls, limits, refinement, posterior reconstruction,
capabilities and simulation/selection, with zero failures/errors/warnings/skips.
A fresh 314-file snapshot has exactly one changed existing file (NEWS) and three
new files (private R helper, test and generated internal help). All prior fitter,
prior, public-interface and numerical-gate sources are byte-identical, and all
four closed predecessor evidence manifests pass. The build rebuilt vignettes;
standard as-cran check with examples/vignettes and separate exact archive tests
are running under the global interpreter cap. Full validation remains pending.

## Final validation, October 29 at 01:39 Eastern

The fresh exact 314-file candidate passes the standard as-cran check, including
examples, donttest examples and rebuilt vignettes. The unchanged exact archive
tests ran separately against that check-installed package with NOT_CRAN=true:
5,261 passing assertions, zero failures, 12 existing warnings and 10 development
skips. Focused counts (217 new, 429 compatibility) are not additional to that
total. Archive SHA256 is
D0F606A9A7EE7E5BC62F872D7D1BE0F52E89816611FF0EF2DE69A441BC8FF31F.
Source and archive-test identities remain exact, branch/HEAD/version/index are
unchanged, and all four closed predecessor/historic evidence audits pass.

The full candidate gate has44resource samples, max3actual R interpreters,
minimum13.70GiB RAM and46.24GiB free disk; no unknown/foreign process occurred.
Closeout confirms no R process remains. Compact numeric traces preserve all12
nuisance attempts and18limiting attempts, without AD objects, chains or samplers.
This is validated private software support for a conditional profile trace.
It supplies no accepted global fit, likelihood CI, boundary calibration or
public spatial method. A later boundary-aware interval design remains necessary.
The next source preflight is separately declared and does not inherit this check
for any changed source. Frozen reports/manifests are never refreshed or replayed.

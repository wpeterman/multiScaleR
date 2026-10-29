# Bounded independent spatial reference qualification

Declared October 28, 2026, before new sampling. This is a new selected-case
reference study, not a replay of the closed latent_gaussian_20261028 producers.
No numerical launch is allowed until the functional_precision_20261028 package
gate and this study's deterministic adapter tests pass. Keep the model private.

## Question and fixed target

Can longer, more conservative original-coordinate NUTS runs qualify a reference
for each of the two existing spatial illustrations? Do independent seeds agree
at fixed probability thresholds? Where do any remaining divergences occur?
The earlier null, smooth and dense-metric runs remain failures under their own
frozen protocols. Their diagnostics and denominators will not be replaced.

Use the exact existing null and smooth data, binned landscapes, truth, knots,
fixed covariance/range, priors and full normalized posterior. Regenerate inputs
only through the documented deterministic fixture, then require each original
pilot's SHA256 input hash before sampling. Use the newly checked installed
package, with its archive and 311-file source identity recorded. No source,
population prior, Jacobian, latent density, likelihood, parameterization or
diagnostic cutoff changes. A hash mismatch is a failed admission, not permission
to fit an approximately matched dataset.

The fixed design has 24 sites on a 4 by 6 projected grid and three years, giving
72 Gaussian observations. Gaussian kernel; fixed population scale 220 map units;
annual log-scale SD 0.15; rank-three spatial log-scale field with SD zero (null)
or 0.45 (smooth); knots (200,200), (800,200), (500,800); fixed smoothing range
350; field-SD half-normal prior scale 0.4. All remaining existing priors and
response parameters follow the original pilot. The ecological scale, smoothing
range and field amplitude retain distinct meanings.

## Finite run list and stopping rule

Exactly four planned fits, in this order, with one owned fit worker:

| Run | Dataset | Sampler seed |
|---|---|---:|
| null-a | Original null | 202610301 |
| null-b | Original null | 202610302 |
| smooth-a | Original smooth | 202610303 |
| smooth-b | Original smooth | 202610304 |

Each uses four serial chains, 16,000 iterations with 4,000 warmup, diagonal
Euclidean metric, adapt_delta 0.999 and max_treedepth 14. Preserve the original
deterministic starts: prepared unconstrained vector plus
0.05*sin(parameter_index + chain_index). Use full posterior sampling with
laplace=FALSE, one core and one numerical thread. Do not rerun Laplace; retain
its already recorded compact summaries and CDF values as qualified historical
approximation outputs. No Hessian repair, affine transformation, seed selection,
automatic extension, metric search, extra fits or discarded first attempt.

The per-fit wall limit is 2,700 seconds and the cumulative fit-phase limit is
10,800 seconds, including failures. The external guard samples every 15 seconds
and stops verified owned processes at the first over-limit observation. A
limit/resource/deadline interruption keeps
the last atomic started/chain-summary/output records and remains interrupted.
Unstarted fits remain pending in the planned denominator of four. Do not shorten
chains to manufacture a completed result. All four outcomes must be reported.

## Qualification and summaries

Use the package's unchanged raw diagnostics: Rhat <=1.01, bulk/tail ESS >=400,
zero divergences and depth hits, and E-BFMI >=0.3. Apply the same Rhat/ESS limits
separately to all 72 fitted site-year log scales, both SDs and the east-west and
north-south log-scale contrasts. Reconstruct the field using its frozen
projection and Cholesky factor. Never coerce the private object into the public
Bayesian class or apply independent-year factorization to the coupled field.

Primary functional events use all 228 thresholds frozen in each original
pilot's functional table. Report chain probabilities, binary-indicator mean
ESS, Rhat and MCSE. Every MCSE must be finite and <=0.005; constant indicators
have unavailable precision. Also report each new run's 2.5%, 50% and 97.5%
quantiles, quantile MCSE, and CDF precision at those quantiles, with the same
0.005 CDF criterion. No diagnostic is inferred from permuted/thinned draws.

A run qualifies only if raw, derived, fixed-threshold and own-quantile CDF
criteria all pass and exact input/source identity is retained. A dataset's
reference additionally requires both planned independent runs to qualify, and
their fixed-threshold probabilities to agree within three combined MCSEs for
every event. This is a conservative numerical agreement screen, not a collection
of independent hypothesis tests or a calibration result. All unavailable and
disagreeing events remain explicit. Quantile MCSE has the quantity's own units;
binary-event MCSE uses the probability axis. Neither is posterior uncertainty.

Compare retained Laplace CDFs only at those identical fixed thresholds, using
the original independent Laplace simulation MCSE and the new chain MCSE. Keep
the three-combined-MCSE flag as descriptive approximation evidence, conditional
on the new reference qualification. Do not pool old unqualified chains with
new chains, or claim four new datasets from four seeds.

The old JSON serialized thresholds to 12 decimal digits. The two new runs use
the exact same parsed stored values, so their agreement screen shares exact
thresholds. The retained Laplace CDF values were computed before serialization;
without old draws they cannot be recomputed at the rounded values. Label that
comparison historical and retain this precision qualification. Do not claim
bitwise re-evaluation of the original Laplace CDFs.

For geometry, save compact per-chain energy/E-BFMI, acceptance, step size,
treedepth/leapfrog and SD summaries, and divergent versus nondivergent SD/latent
norm summaries. These locate associations; they cannot prove a cause of a
divergence. Save no sampler, original/derived chain arrays, importance candidates
or AD tapes. Retain errors, warnings, exact identities, timestamps and all gates.

## Operational admission and limits

At most two owned fit workers are allowed globally, but this protocol uses one.
At most four actual R interpreters globally, including rsession/Rserve, one
thread each; available RAM >12 GiB and disk >5 GiB. Unknown inventory blocks
launch. Require an R-free window before the first fit and each later fit;
leave foreign work untouched. Validate source/input identities before sampling
and source identities after each run. Monitor actual R.dll hosts rather than
counting verified forwarding executables as interpreters.

Do not start long phases at or after 04:00 Eastern October 29. At the first wake
at/after 04:15 stop launches and finish a safe checkpoint by 04:35. Stop only
identity-verified owned unfinished workers if necessary; no owned numerical work
past 04:45. Delete the one-night heartbeat by then. No staging, commit, version
bump, merge or push. Public spatial inference, likelihood uncertainty and broad
null/weak/confounded/holdout calibration remain separate unmet gates.

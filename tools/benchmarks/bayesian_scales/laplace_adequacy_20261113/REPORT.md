# Fitted scale-probability checks for uncorrected Laplace

This study evaluates an opt-in assessment of fitted scale probabilities.
The original approximation remains the fit. The assessment estimates the
same posterior CDF from a fresh importance proposal, with a separate
outcome for a discrepancy, a difference within the requested tolerance,
or a result that remains unresolved.

## Method

The proposal uses a heavy-tailed version of the fitted joint or marginal
geometry and a 20% defensive latent-prior component. Response coefficients
are sampled conditionally on the outer parameters. Their conditional is
Gaussian for Gaussian responses. Other responses use a normalized Student t
proposal fitted to the conditional geometry. Each vector receives the ratio
of the exact unconstrained posterior density to the actual proposal density.
The calculation never drops an invalid ratio. The original approximation
draws are independent of the fresh reference simulation.

The reference uses the package's existing importance overlap checks. Event
classification additionally requires an indicator ESS of at least 400 and
an absolute difference between independent halves of no more than the smaller
of 0.05 and the requested tolerance. For each event, a numerical interval for
the difference subtracts the weighted reference CDF and its three-MCSE band
from the approximate CDF's Wilson interval using a multiplier of three.
The procedure declares a discrepancy when the entire interval is outside
the tolerance on the same side. It declares within-tolerance when the whole
interval lies inside the tolerance. An original fit that failed its numerical
diagnostics always leaves the assessment unresolved.

These numerical bounds are conditional on the sampled proposal. Neither
finite-variance approximations for weighted estimates nor diagnostic
thresholds can establish that the proposal visited every important mode.
Events can be highly dependent within a dataset. Report counts alongside
the number of independent datasets and interpret the numbers descriptively.

## Frozen simulation

Five independent datasets are generated in each of seven scenarios. Their
declared priors and synthetic response geometry span strong and weak scale
information, narrow and broad scale priors, annual variation close to zero,
an informative annual design, a repeated-site response effect, and two NB2
designs. Gaussian and exponential kernels are represented. The simulator
computes weighted covariates and responses independently of package fitting.
The exact settings and seeds are in `design.R`. They are fixed generating
configurations, with independent simulated datasets within each scenario.

The fit plan includes 35 NUTS references, 35 uncorrected joint Laplace fits,
and 25 uncorrected marginal Laplace fits. Marginal Laplace needs a scale
random effect, so ten shared-scale cases are outside that method's scope.
NUTS runs four chains with 4,000 iterations per chain, including half warmup,
with target acceptance 0.99 and maximum tree depth 12. Eligibility requires
the package's raw and derived diagnostics. Approximations retain 8,000
draws. The assessment uses 65,536 new candidates and a tolerance of 0.05
at 150, 250, and 350 map units. No tuning follows the fit outcomes.

NUTS probabilities use their retained-chain Monte Carlo standard errors.
For the independent comparator, the same interval construction substitutes
NUTS for the importance reference. An ineligible fit, undefined reference
precision, or ambiguous difference stays unresolved. Evaluation reports
the full three-way comparison rather than treating an unresolved decision
as agreement.

## Reproduce

Run the two fitting shards from the package root under R 4.6.1 with the
development DLL available. The runner requires RTMB, tmbstan, rstan, posterior
and digest. The archived configuration records the package source hashes,
locale, package versions, seed plan and candidate budget.

```text
Rscript tools/benchmarks/bayesian_scales/laplace_adequacy_20261113/run.R fit 1
Rscript tools/benchmarks/bayesian_scales/laplace_adequacy_20261113/run.R fit 2
Rscript tools/benchmarks/bayesian_scales/laplace_adequacy_20261113/run.R assess
```

Each attempt gets one compact RDS record. Rerunning a shard skips completed
records only when the source fingerprint matches. The report preserves
warnings, failed fits, rejected fits, and unresolved events. Neither
posterior draws nor differentiation tapes are archived.

## Results

All 95 planned attempts returned. The original numerical diagnostics
accepted 35/35 NUTS references, 31/35 joint approximations, and 17/25
marginal approximations. The importance calculation returned finite
reference summaries for all 60 approximations. Its global checks passed
for 20/35 joint and 16/25 marginal references, including some with
rejected original fits.

| Method | NUTS comparator | Assessment discrepancy | Assessment unresolved | Assessment within tolerance |
|---|---|---:|---:|---:|
| Joint | Discrepancy | 54 | 81 | 0 |
| Joint | Unresolved | 3 | 126 | 11 |
| Joint | Within tolerance | 0 | 28 | 27 |
| Marginal | Discrepancy | 21 | 46 | 0 |
| Marginal | Unresolved | 1 | 140 | 39 |
| Marginal | Within tolerance | 0 | 27 | 26 |

The assessment detected 54 of 135 NUTS-supported discrepancies for joint
Laplace and 21 of 67 for marginal Laplace, across 12 and 5 independent
datasets respectively. No NUTS-supported discrepancy was reported within
tolerance, and no NUTS-supported within-tolerance event was flagged as a
discrepancy. Four flagged events and 50 within-tolerance decisions had
unresolved NUTS comparators. The assessment remained unresolved on 448 of
630 matched events. These are counts for the frozen panel, without treating
different thresholds or annual scale estimates within one dataset as
independent replicates.

Results differed by geometry. The joint assessment flagged 11 of 15
events in the broad-prior shared Gaussian panel and 30 of 60 in annual
NB2. The marginal assessment resolved 56 of 60 events in the informative
annual Gaussian panel. Both references failed overlap checks in the
Gaussian repeated-site panel. In the repeated NB2 panel every event
remained unresolved, even though a single marginal reference passed
global overlap checks. In those cases, a full posterior fit remained
necessary for the probability comparison.

The results support an optional screen for selected fitted scale
probabilities. They do not qualify uncorrected Laplace as a replacement
for full posterior sampling or establish broad accuracy for the screen.
The source and results preserve all unresolved events. A separate Gaussian
test integrates the coefficients analytically and checks the weighted
scale CDF against numerical quadrature.

Both fitting shards and assessment exited successfully. The independent
`verify.R` audit confirmed the frozen hashes, all 95 fit identities,
and 630 paired events before the roxygen-only documentation cleanup.
The review verification matches every comparison's key, score, precision
and decision to its original approximation and NUTS records. The tested
non-roxygen lines of `R/bayes_laplace_diagnostic.R` have
SHA-256 `a9308cbbe610956c5bc7cacfc90628c04bcf609242536687f0142df4191e080b`
when joined with LF and no terminal newline. The prose-only change leaves
that fingerprint identical. The archived `frozen-diagnostic.R` preserves the original full-file SHA-256.
For a hash-identical replay, start from commit `8f7a26a8e`. Copy this study
directory and `R/bayes_scale_probabilities.R` into that checkout. Use
`frozen-diagnostic.R` as `R/bayes_laplace_diagnostic.R`, and register the
`diagnose_bayes_laplace` export in `NAMESPACE`. The audit can verify the
frozen source in the review worktree without replacing the current package
files. The `frozen-source/` directory retains original versions of the
files touched by the other inference updates.

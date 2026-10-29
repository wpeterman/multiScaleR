# Annual-scale inference and consolidation

Updated October 28, 2026. Numerical studies, the final frozen package check,
five-file live integration and source/evidence audit are complete. The final
operational audit found that the original process guard missed three existing
RStudio interpreters. The check passed, but its required foreign-interpreter-free
window was not established. That operational gate remains unmet.
Nothing has been staged, committed, pushed or merged. Version 0.7.4
and HEAD `87bf8779fe571921f2fd5d3e45c0f9e478530022` are unchanged by this work.
The deadline is 05:30 America/New_York. Expensive launches stop by 04:15.

## Assessment

The useful gains concern reliability and a maintainable model definition. The
earlier weak-case median discrepancy does not persist in the longer independent
comparison. A common response-density layer passes its mathematical and focused
tests, and all six original scientific decisions pass. The final package check
also passes. Further small cache speedups have lower priority.

The new exact Gaussian coefficient integration is mathematically sound on its
independent checks, but its first sampling implementation is a poor performance
choice on this difficult input. Both pilots miss the fixed precision target.
It stays experimental. Reducing parameter dimension alone is insufficient.

## Long NUTS comparison on the disputed input

All methods use exactly the archived weak sparse-year2 input4, with 48 observations
in five years (3/5/8/12/20), kernel coefficient 0.4, unchanged complete priors and
fixed full-scaffold predictor standardization. Four physically matched chain
starts are reused. Two independent seeds each run original noncentered and exact
transported coordinates, with four serial chains,128,000 iterations and 16,000 warmup iterations,
adapt_delta = 0.99 and depth 12. No fit was retried or selected after its outcome.

The screen checks agreement between chains (Rhat), effective sample size (ESS)
and Monte Carlo standard error (MCSE), which measures uncertainty from finite
simulation. BFMI screens exploration of the sampler's energy distribution.
Every fit passes the original sampled and derived checks: Rhat <= 1.01, bulk/tail
ESS >= 400, zero divergences, zero depth hits and BFMI >= 0.3. Both primary probability
targets also pass their separate Rhat/ESS requirements and MCSE <= 0.005. These
targets are probabilities that year1 log distance is <= 4.85 and year4 is <= 5.3.
Constant outer-tail indicators remain nonestimable rather than falsely passing.

| Fit | Whole-fit seconds | Sampler/derived checks | Primary probability precision |
|---|---:|---|---|
| Noncentered A | 1441.25 | Pass | Pass |
| Transport A | 652.87 | Pass | Pass |
| Noncentered B | 1429.66 | Pass | Pass |
| Transport B | 676.87 | Pass | Pass |

All four predeclared independent comparisons pass. None of the 72 annual/SD
endpoint and median comparisons exceeds three combined estimated MCSE. Every
mean and primary probability comparison also passes. The largest standardized
mean difference is 2.072. The largest primary probability difference is 1.008.
Quartile comparisons and all original flags remain in the compact evidence.

The first long pair's year4 medians are 5.508410 and5.511379 log units, compared
with the earlier short-run medians 5.52273 and5.10250. The corresponding long-run
probabilities at 5.3 are 0.463900/0.462435, with MCSE 0.003822/0.004326. The
discrepancy does not persist under the longer, matched-start, independent-seed
protocol. This supports a finite-run explanation, but the experiment does not
isolate run length, warmup, initialization or quantile-MCSE estimation as the
specific cause. Original short-run traces were not retained, so their compact
diagnostics cannot resolve those mechanisms retrospectively. The original
failure remains. Agreement across these runs cannot prove discovery of every mode.

The noncentered runs took about 2.16 times as long as the transported runs. Both
methods have identical iteration counts and posterior targets, but machine load
was shared. This is a selected-case timing result, not a general speed guarantee
or interval-calibration study. Public coordinate/draw integration remains
separate work. Evidence: long-reference-analysis.json/rds and four long-reference
fit records. No samplers or raw draw arrays were saved.

![Agreement and uncertainty on the selected difficult input](long-reference-comparison-v3.png)

Top points are cumulative probabilities. Bottom colored points are posterior
medians, and their intervals are marginal 95% credible intervals. Top error
bars measure simulation precision. The unpooled
likelihood contributes points only. The broad population-SD interval illustrates
remaining uncertainty with five years. This selected input does not establish
coverage, and a positive SD credible interval is not a test against zero.

## Exact Gaussian coefficient integration

Normal response-coefficient priors permit exact integration of those coefficients
for Gaussian responses. Annual effects, SDs, dispersion and the original priors
retain their meanings. Eight independent dense observation-covariance checks
compare values, gradients and Hessians. A scalar coefficient integral and an
unsupported-family rejection also pass. The scalar difference is 0 at recorded
precision. This mathematical audit took 0.6 seconds.

Two fixed pilots use the same disputed data and remaining physical starts,
32,000 iterations and 8,000 warmup iterations, four serial chains, unchanged settings, and 900 seconds
per-fit limits. Both pass sampled/derived diagnostics. None of 162 prespecified
endpoint/median comparisons flags disagreement, and means/probabilities agree
within their estimated simulation error.

| Pilot | Seconds | Worst bulk/tail ESS | Fixed probability precision |
|---|---:|---|---|
| Beta A | 771.53 | 4758/24926 | Fail |
| Beta B | 596.91 | 5266/26768 | Fail |

The overall pilot gate is FALSE because both miss MCSE <= 0.005. Year4 probability
MCSE is 0.009011/0.008574, despite ordinary convergence passing. Runtime is
comparable with the much longer transport references. It is not a preferred
sampler for this input. The prototype does not reconstruct coefficient draws
and is not a public prediction-capable fit object.

A deterministic 12,000-evaluation audit at the same four physical starts took
0.63 seconds and created no fit. Across 2,000 distinct perturbed gradient calls,
noncentering took 0.09 seconds, transport 0.06 and exact-beta integration 0.16.
Clock resolution and shared load limit these short timings, but they identify
extra evaluation cost rather than an automatic benefit from fewer parameters.

| Method | Year4 probability ESS per whole-fit second |
|---|---:|
| Noncentered long reference | 11.69 |
| Transport long reference | 20.12 |
| Exact-beta pilot | 4.71 |

These observed rates include warmup and diagnostic costs, with different chain
lengths reported above. They are descriptive efficiency estimates for one
functional, not a controlled universal ranking. No further pilot or tuning is
planned. The exact integral remains a possible building block for other models.

## What is being consolidated

The same hierarchical model can support likelihood and Bayesian inference. If
a denotes annual effects and theta the population parameters, the normalized
hierarchical likelihood is L(theta)=integral p(y|a,theta)p(a|theta)da. A Bayesian
posterior adds proper population priors, p(theta|y) proportional to L(theta)pi(theta).
Both models pool. Laplace describes how an integral is approximated. The
inferential target depends on which distributions enter that integral.

| Implementation | Target and role | Present limit |
|---|---|---|
| Existing optimized/profiled likelihood and mixed backends | Data likelihood, with specialized profile or random-effect integration | Identification and confidence-inference gates still matter |
| Bayesian Laplace | Local Gaussian approximation to the proper posterior | Skewness and competing modes can limit interval accuracy |
| Bayesian nested/marginal approximation | Integrates selected population uncertainty, with correction diagnostics | Annual-tail calibration and precision remain incomplete |
| Accurate factorized annual likelihood (RTMB prototype) | Higher-accuracy independent annual integrals, with unpenalized population parameters | Factorization restrictions and separate population-SD uncertainty inference |
| Full Bayesian NUTS | Joint proper posterior, including population-SD uncertainty | Computational cost and exploration require explicit checks |

The candidate now shares normalized Gaussian, Poisson, binomial and NB2 response
densities across Bayesian, Laplace-likelihood and factorized-likelihood adapters.
It preserves likelihood constants, prior penalties, latent densities, exact
Gaussian covariance integration and each kernel's existing arithmetic. Ordinary
NB2 computes its existing size parameter directly, avoiding predictor cancellation.
No public interface or default changes in this five-file batch.

Forty predecessor/current value, gradient and Hessian comparisons pass at the
original 1e-10 tolerance, with maximum difference 7.11e-15. Response attributes
are checked separately and unchanged. All 1,290 focused assertions pass, including
the concurrent prediction tests, with zero failures/errors. All six fixed
scientific decisions also pass. The final current-source package check finished
at 01:54 Eastern with Status: OK, zero package-check errors, warnings or notes,
successful examples and rebuilt vignettes. Installed tests report 4,225 passing
assertions, zero failures, 13 test-level warnings and 10 development-only skips.
These test warnings are distinct from package-check warnings: nine concern
missing test-polygon coordinate systems, three concern fitted binomial
probabilities near zero or one, and one reports the installed glmmTMB/TMB build
version mismatch (1.9.21 versus 1.9.25). Dependencies were not changed during
concurrent work. This successful check must be read with the process-isolation
limitation described below.

Exact predecessor and protected-source identities passed before the five-file
promotion at 02:03 Eastern. All 1,316 live focused assertions then passed, with
zero failures or errors. The final evidence audit verifies 277 live code, help,
test, example, vignette and data files against the frozen candidate, all six
retained scientific outputs and the unchanged frozen producers. No scientific
fits were repeated. One descriptive NEWS addition from concurrent prediction
work differs from the checked snapshot, with both exact hashes retained in
concurrent-news-delta.json. It was preserved. The whole live tree is therefore
not byte-identical to the checked snapshot.

| Scientific case | Required decision retained | Seconds |
|---|---|---:|
| gaussian-positive | Accepted | 203.11 |
| sbc-111 | Accepted | 241.53 |
| sbc-057 | Accepted, exact zero-SD model available | 94.19 |
| prior_shift-018 | Rejected, competing optima | 290.32 |
| sbc-123 | Rejected, limiting-model witness | 123.58 |
| sbc-187 | Rejected, limiting-model witness | 109.32 |

These are numerical regression decisions. Accepting sbc-111 in this gate does
not resolve its separate marginal-posterior lower-tail calibration question.

The initial new test failed only because legacy binomial values retain row names
while its base-dbinom reference did not. That failed record, test, logs and entire
snapshot remain intact. The successor removes incidental names from the value
comparison and adds exact predecessor/current attribute checks. No prior,
mathematical target or diagnostic cutoff changed, and no scientific fit was
dispatched by the failed queue.

The combined successor snapshot includes independent prediction development.
That work has its own owner, source freezes, example and 16-dataset predictive
protocol. It is preserved and not claimed as this study's implementation. Its
new-group prediction behavior is not permission to invent future fixed-factor
coefficients or temporal autocorrelation. Full checking waits for foreign
interpreters to avoid disrupting that work. Promotion requires exact live
predecessor and final protected-source/help identities and copies only the five
response-layer files. Concurrent edits will not be overwritten to force a pass.

Five prediction/output or documentation files changed again after the combined
snapshot. Its waiting-only checker and dispatcher were identity-stopped before
any package build or test started. All completed scientific results are retained.
One final current-source check ran after the original guard reported a free
interpreter window. The final loader audit subsequently found that this report
omitted RStudio's rsession-utf8.exe hosts.
Scientific evidence was reused only after every numerical R/src file was byte-identical,
with three explicit output-file exceptions and an identical parsed fitting entry.
Current public tests, examples and vignettes ran afresh. This avoids both
stale-source promotion and repeated scientific fits.

## When Bayesian pooling helps, and when it does not

Bayesian inference is useful when years plausibly share a scale distribution,
annual measurements are weak, and uncertainty in the population SD matters.
Proper priors can express defensible information and support a coherent posterior
for annual scales and new groups. Hierarchical likelihood can also pool, and
may be faster when its likelihood is well identified and appropriate uncertainty
inference is available. Computational convergence cannot establish identification.

The closed October27 results provide mixed evidence, which must stay visible.
Transport improved eligibility on ten stronger fixed inputs to 10/10, compared
with 4/10 original runs. All 50 annual truths were inside its nominal 95-percent
intervals, but ten independent datasets do not establish calibration. At zero
SD, mean annual-spread/unpooled ratios were 0.408, versus 0.993 at positive SD. This
shows adaptive pooling in those examples. The frequentist zero-SD point RMSE
was 0.00665 versus Bayesian 0.01210, so pooling does not guarantee lower point error.

Four weak unequal-year illustrations favored unpooled point estimation overall.
Their mean sparse/dense log-scale RMSE was 0.1160/0.0844 for unpooled ML and
0.1502/0.1011 for Bayesian posterior mean log scales. The prior median matched
the generating population location, which does not establish robustness to the
other prior assumptions. No prior sensitivity result is claimed in this
unchanged-prior study.

On selected weak input4, greater integration accuracy alone did little for the
hierarchical likelihood's point error. The archived accurate fit took 55.70
seconds with the same six original starts. Accurate modes differed from Laplace by
at most 0.00883 log units, with log-likelihood gain 0.12747. Dense-year RMSE was
0.49447 for accurate modes, 0.45441 for accurate conditional means and 0.49559 for
Laplace modes. The original Bayesian mean's dense-year RMSE was 0.16695. This
contrasting selected case does not override the four-input average. Point
estimands differ, and neither comparison supplies a general ranking.

Proper priors can make inference well defined even when an unpenalized likelihood
has no useful finite optimum. This follows analytically for the implemented
families, rather than from a new simulation. Discrete likelihood masses are at
most1. For n conditionally independent Gaussian observations, the likelihood is
bounded by (2*pi)^(-n/2)*sigma^(-n). The implemented log-sigma normal prior has
finite negative moments: E[sigma^(-n)]=exp(-n*m+n^2*s^2/2), for log-sigma~N(m,s^2).
Together with the other normalized proper priors, this bounds the evidence.
It does not create data information, guarantee finite natural-scale moments,
or make sampling reliable. Log-scale summaries and finite posterior quantiles
remain distinct from potentially unstable arithmetic moments.

## Flexibility and access

Formula-based hierarchical scale estimation is available through specifications
such as `scales = list(habitat = ~1 + (1 | year))`. The nonlinear kernel still requires
numerical estimation. There is no general closed-form annual-scale estimator.
Exact Gaussian coefficient/response-effect integrals remove particular blocks,
not all annual or population uncertainty.

The joint Bayesian representation supports more implemented combinations than
independent annual factorization, including supported multiple kernels and
crossed/additive grouping structures. The scalar factorized prototype requires
one kernel and at most one scale grouping term. Its current response-effect
integration is Gaussian and excludes response effects that couple years. Joint inference
is needed when those independence conditions fail. Temporal covariance, random
slopes, custom families and other unsupported models still need new code and
validation. The factorized likelihood remains a private prototype, with a stable
exported interface reserved for a later integration batch. No model here has an
exclusive mathematical claim to Bayes.

## Development line and remaining work

INTEGRATION_STRATEGY.md recommends one future development line rooted in current
master ancestry, with reviewed batches for the shared response layer, neutral
model preparation, explicit targets/integration/uncertainty metadata and separate
research improvements. Current master and this branch share HEAD, and much
inference code is uncommitted. Old dev/development have substantial divergent history and
should not be merged wholesale. No actual branch merge is performed here.

Remaining scientific limits include broader independent calibration,
boundary-aware population-SD confidence inference, prior sensitivity and public
coordinate/draw reconstruction. Existing sbc123/187 limiting-witness rejections
and prior018 competing-optimum rejection remain mandatory. A positive Bayesian
SD credible interval is not a test against exact zero under a continuous prior.
Posterior credible, confidence and repeated-latent prediction-error intervals
must retain their different meanings.

## Resources, cleanup and reproducibility

Workers use one numerical thread, at most four global numerical leaves, known
process inventories, >12GiB available RAM and >5GiB disk, with default 16 GiB
admission. Exact creation-time/command/ancestry records identify owned processes.
Owned heavy work is serialized. No raw sampler, draw array, AD tape or duplicate
simulation design is written. Compact checkpoints are atomic.

The original 113.26 MiB failed-test snapshot and 113.97 MiB combined successor are
retained as provenance. Closed studies and every original failure are preserved.
No deletion succeeded in this continuation. Previously policy-blocked DLL
folders were not retried. The final new owned DLL folder contains one 1,488,384
byte copy identical to retained src/multiScaleR.dll. Its producer exited, and
the exact active loader inventory points to an installed library copy rather
than this temporary copy. It is a verified regenerable cleanup candidate.
Automatic approval review rejected its deletion with only "blocked by policy."
The folder is retained, and no alternate deletion method was attempted.

| Artifact | Final classification |
|---|---|
| Temp/RtmpqmZsrr/pkgload84644af271cb | Verified new owned duplicate DLL, 1.42 MiB; deletion policy-blocked, retained |
| Temp/RtmpO8L8GI/pkgload7ce45781267a and Temp/RtmpQ7MDBG/pkgloadebd42ea0746 | Earlier policy-blocked folders remain untouched; no retry |
| Failed, combined and final source snapshots; closed studies; compact fit/check evidence | Retain as provenance |
| Compiler objects, build archive and installed check tree | Regenerable in principle; retained pending an exact separate cleanup review |
| Sampler objects, raw draw arrays, serialized AD tapes or duplicate pilot designs | None written by this study |

The process correction recognizes rsession variants and Rserve hosts, counts
them conservatively, and throws on unavailable process or RAM inventory. Its
12 name cases, unknown-inventory failure test and real foreign-interpreter
blocking check pass without launching R. Three RStudio sessions created on
October26 remain untouched. Their omission means the earlier global-interpreter
counts and foreign-R-free claims cannot be certified retrospectively. Numerical
results and source identities still pass their own gates. No claim of complete
operational compliance is made. A later isolated integration check requires an
actually free interpreter window; the completed check is not blindly restarted.

Final audit readers are also retained. The first stopped at the NEWS difference.
The second correctly pinned that exact descriptive exception but assumed the
beta summary held only its two pilot rows; it actually includes four reference
rows too. The final reader verifies all six identities before checking the two
pilots separately. These were read-only audit failures, with no fitting retries,
prior changes or relaxed numerical diagnostics. Final evidence is in
final-audit-evidence-v3.json, process-guard-validation.json and the operational
audit. Resource headroom was 21.8 GiB RAM and 61.7 GiB disk at 02:15 Eastern.
Producer snapshots, study inputs/results, source, reports and fingerprints remain
retained. There is no further prespecified fit to run productively in this window.
Owned workers and controllers have exited. The exact-identity-verified watchdog
was stopped, and this study's heartbeat was deleted. Foreign RStudio sessions
remain untouched. The source/check evidence is ready for review, with the
operational isolation gate explicitly incomplete.

Primary evidence is in PLAN.md, RESUME.md, integration-inventory.json,
concurrent-source-delta.json, long-reference-analysis.json/rds,
beta-math.rds, beta-pilot-analysis.json/rds and combined-* validation records.
Earlier comparisons remain in the closed factorized_20261027 report. RESUME.md
records the stopping point and next development priorities before the 05:30 cutoff.

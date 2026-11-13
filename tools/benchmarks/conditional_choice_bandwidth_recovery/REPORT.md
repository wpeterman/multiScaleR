# Conditional-choice bandwidth recovery

## Design and numerical comparison

This study simulated 384 independent datasets with one used and three
available alternatives per choice set. Gaussian and exponential kernel
designs crossed two training sample sizes and four habitat signals.
Every dataset had 160 independent held-out choice sets. A generating
bandwidth exists for the three nonzero habitat effects. See the
[protocol](PROTOCOL.md) for the simulation design and fixed seeds.

The public conditional-logit backend used its existing six starts.
The independent reference evaluated the matched-set likelihood in base R
and reestimated every coefficient at each bandwidth. It used a fixed
121-point log grid, additional refinement inside sampled local minima
and every adjacent bracket where the independent bandwidth score crossed
from negative to positive, and both bounds. It provides a checked
numerical profile over the declared domain. The maximum absolute projected
bandwidth score fell from 0.0218 in the initial reference to
`3.71e-6` after refining adjacent score brackets. The package fit,
assessment probabilities and native diagnostics apart from run time
reproduced exactly across both executions.

## Recovery and assessment predictions

All 384 public fits returned an eligible score. All 384 independent
likelihood profiles returned. Each public fit retained four parameters,
including the bandwidth, and counted one event per training choice set.
The held-out probabilities independently reproduced to within
`3.8e-15` absolute probability. This checks normalization and
prediction at the returned bandwidth even if it differs from the
best independent profile point.

| Habitat signal | Estimates within 25% of truth | Estimates at a search bound | Profiles supporting both bounds | Mean fitted minus no-habitat held-out loss |
|---|---:|---:|---:|---:|
| Strong interior | 78/96 | 0/96 | 0/96 | -0.325 |
| Weak interior | 22/96 | 27/96 | 47/96 | -0.0151 |
| Strong near lower bound | 88/96 | 9/96 | 0/96 | -0.321 |
| Exact null | Undefined | 46/96 | 93/96 | +0.00482 |

The last column is natural-log loss per held-out choice set, with a
training-fitted no-habitat model as the comparator. Each signal row
averages its four kernel by sample-size cells. Supporting both bounds
means the sampled profile at each bound lies within 1.92 log-likelihood
units of the best independent point. This descriptive threshold does not
supply interval coverage. An exact habitat null has no generating
bandwidth because the habitat term vanishes.

For strong interior effects, the median fitted-to-generating bandwidth
ratios were 1.03 and 0.98 for Gaussian kernels with 100 and 450 training
sets, and 1.07 and 1.05 for the exponential kernel. Increasing training
sets from 100 to 450 reduced weak-signal boundary estimates from
20/48 to 7/48 and increased recovery within 25% from 8/48 to 14/48.
The full [recovery table](results/recovery-combined.csv) separates kernel,
sample size and signal and gives percentile ranges. The
[prediction table](results/prediction-combined.csv) gives held-out losses
and simulation-bootstrap intervals for paired contrasts with the model
at the known bandwidth, the oracle choice probabilities, the
no-habitat model, and the independent profiled fit.

## A weakly identified fit can have several eligible solutions

In 383 of 384 cases, the best independently checked objective was within
0.0001 total nats of the public result. The exception was a Gaussian
kernel with 450 training sets and an exact habitat null (`choice_047`).
All six default starts passed the existing native checks. Four reached
a local optimum near 148 units and two reached one near 323 units. The
public fit selected 147.78. The independent profile had a better solution
at the upper limit of 800 units by 0.01689 total training nats.

A separately recorded six-start diagnostic included both search limits.
It selected 800 units and agreed with the independent objective at the
displayed precision. Held-out loss was 1.272263 nats per set, compared
with 1.272848 for the primary fit. The training improvement was small
and the generating habitat effect was zero, so neither bandwidth
represents a recovered ecological scale. The original starts and primary
fit remain in the main cohort. No package defaults or uncertainty checks
were changed.

## Interpretation

These simulations isolate one ecological bandwidth under a correctly
specified response model. They establish that the package evaluates the
one-case likelihood and predicts complete assessment sets correctly for
the tested data. They also show how a small habitat effect can produce
a broad likelihood profile. More training choice sets do not resolve
every weak signal. Boundary estimates can be useful descriptions of the
fitted score surface, but they are limited by the supplied bandwidth
range. Assessment loss describes prediction for a matched set and cannot
establish identification of its ecological bandwidth.

The designs use independent sets with synthetic cell values and two
movement covariates. They do not address animal-level transfer, serial
dependence, availability-proposal misspecification, landscape metric
radii, or profiles involving multiple bandwidths. The study reports
finite-sample recovery and conditional prediction under these simulation
designs. The profile span is descriptive, and scale confidence intervals
remain unavailable for this backend.

## Reproduce

Run from the package root with the development package and survival
available. Use two disjoint parts or run both commands sequentially.
Preserve the raw records before rerunning.

```sh
Rscript tools/benchmarks/conditional_choice_bandwidth_recovery/run.R /path/to/results/part1 1 2
Rscript tools/benchmarks/conditional_choice_bandwidth_recovery/run.R /path/to/results/part2 2 2
Rscript tools/benchmarks/conditional_choice_bandwidth_recovery/summarize.R /path/to/results
Rscript tools/benchmarks/conditional_choice_bandwidth_recovery/run.R /path/to/results/boundary-starts 1 1 384 /path/to/results/cases-combined.csv
Rscript tools/benchmarks/conditional_choice_bandwidth_recovery/summarize-boundary.R /path/to/results
```

The initial-grid script and output are also preserved. The raw archive
contains complete held-out probabilities, all independent profile points,
refinement candidates and starts for both executions and the diagnostic.
Unpack it inside `results` to reproduce summaries against the archived inputs. To reproduce the
paired replay, run `run-initial-grid.R` into `initial_grid_part1` and
`initial_grid_part2` with the same shard arguments. Then run
`verify-replay.R /path/to/results`. The original executions saved every
case before a trailing reader error when their source file was edited.
The replay uses the final script and independent refinement.

- [Archived raw profiles, predictions and starts](results/raw-profiles-and-predictions.zip)
- [Raw archive index](results/raw-archive-index.csv)
- [Case-level results](results/cases-combined.csv)
- [All native eligibility records](results/starts-combined.csv)
- [Recovery and profile support](results/recovery-combined.csv)
- [Held-out choice log loss](results/prediction-combined.csv)
- [Independent likelihood and prediction agreement](results/agreement-combined.csv)
- [Reference refinement comparison](results/independent-reference-refinement.csv)
- [Boundary-inclusive diagnostic comparison](results/boundary-start-comparison.csv)

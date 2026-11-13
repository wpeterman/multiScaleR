# Fitted spatial response-field predictive assessment

## Main result

The public field-CV function scored 341 of 360 prespecified simulations.
All 19 failed cases stopped at the initial native training fit. All 360
known-parameter calculations returned. The fitted joint score was worse on
average than the paired known-parameter score in all nine family-design
cells. The empirical 95% field-bootstrap intervals for those mean excess
losses stayed above zero.

Each simulation generated a new shared field, independent covariates and
complete repeated-location observations. The study compares interpolation,
a sparse stronger field and buffered transfer. Each cell starts with 40
independent simulations. The paired results below use the same assessment
outcomes and score target at the fitted and known parameter values. See
[the protocol](PROTOCOL.md) for the design and eligibility rules.

| Family | Field design | Fitted scored | Mean excess joint loss per row, 95% interval | Fitted / known mean squared standardized residual |
|---|---|---:|---:|---:|
| Gaussian | Interpolation | 39/40 | 0.027 [0.015, 0.039] | 1.06 / 0.96 |
| Gaussian | Sparse stronger field | 38/40 | 0.058 [0.023, 0.096] | 1.03 / 0.93 |
| Gaussian | Buffered transfer | 38/40 | 0.060 [0.034, 0.089] | 1.46 / 0.99 |
| Poisson | Interpolation | 40/40 | 0.021 [0.007, 0.037] | 1.05 / 0.97 |
| Poisson | Sparse stronger field | 39/40 | 0.060 [0.025, 0.099] | 1.20 / 0.85 |
| Poisson | Buffered transfer | 39/40 | 0.040 [0.017, 0.070] | 2.68 / 1.33 |
| NB2 | Interpolation | 36/40 | 0.015 [0.002, 0.028] | 1.02 / 0.98 |
| NB2 | Sparse stronger field | 34/40 | 0.119 [0.042, 0.226] | 1.82 / 1.38 |
| NB2 | Buffered transfer | 38/40 | 0.095 [0.043, 0.160] | 1.88 / 0.68 |

Loss is in natural log units. A lower score is better. The squared
standardized residual is the squared prediction error divided by the
predicted observation variance. A calibrated predictive mean and variance
give an expected value of one. The statistic can also rise when the mean
is wrong. Average first within each independent field to account for shared
locations and visits. The table's known-parameter column uses the scored
fitted subset. The [calibration table](results/calibration-combined.csv)
also includes all 40 generating-parameter predictions per cell.

## Native fitting and parameter recovery

Gaussian, Poisson and NB2 scored 115/120, 118/120 and 108/120 simulations,
respectively. Three Gaussian failures came from disagreement between eligible
starts, and two came from a start failing the original eligibility checks.
The two Poisson and 12 NB2 failures had no eligible start. Of the 24 starts in failed NB2 simulations, 23 lacked positive fixed
covariance and 22 had a nonpositive native Hessian or joint precision.
The returned non-Gaussian fits include 21 with one eligible start and
21 with two eligible but disagreeing starts. All starts remain available for review. The study
used the existing fit-selection and reporting rules.

Among eligible fits in the sparse stronger design, median fitted-to-generating
range ratios were 0.48 for Gaussian, 0.25 for Poisson and 0.37 for NB2.
The range fell below one third or above three times its generating value
in 12/38 Gaussian, 21/39 Poisson and 16/34 NB2 scored fits. The NB2 size
had a median of 1.20 times the generating value, but its 90th percentile
reached 25.2 million times that value. That tail includes near-Poisson
variance behavior even though the data were generated under NB2. The
[recovery table](results/recovery-combined.csv) retains signed coefficient
errors, field SD ratios, observation-dispersion ratios and percentile ranges.

For the exact Gaussian predictive comparator, fitted nominal 90% marginal
interval coverage was 88.8% in interpolation, 90.0% in sparse stronger
fields and 84.4% in buffered transfer. Paired known-parameter coverage
was 89.7%, 88.9% and 86.6%, respectively. These percentages average
coverage by field. The buffered Gaussian fitted coverage interval
from resampling fields was 80.4% to 88.4%. Its known-parameter analogue
was 82.9% to 90.1%. The lower fitted coverage occurs even when the
prediction distribution calculation is exact at the specified parameters.

## Independent count reference

The prespecified subset covered the first two simulations per count family
and design. Independent elliptical slice and likelihood-power bridge
estimates resolved all 12 generating and all 11 eligible fitted score
references. One planned fitted NB2 reference had no eligible training fit.
No attempted reference failed. The largest fitted joint loss difference
was 0.00303 nats per observation, with reference MCSE 0.00728 in that case.
The largest generating-parameter difference was 0.01445 nats per
observation, with MCSE 0.00388. Thus small approximation differences at
fitted parameters coexist with larger fitted-minus-generating loss in
the 360-case panel.

All 556 available count rows met the prespecified mean, variance and
randomized-PIT reference precision criteria. The largest relative
discrepancies at fitted parameters were 1.92% for the mean and 2.82% for
the observation variance. For the randomized PIT, the fraction between
0.05 and 0.95 averaged 87.4% over six Poisson fields and 89.9% over five
NB2 fields at fitted parameters. The generating-parameter counterparts
averaged 88.4% and 89.3%, respectively, over six fields per family.
These are illustrative counts from a small cohort. The primary count
results assess moments and joint scores, and leave tail and interval
calibration unresolved.

## What the comparison measures

For Gaussian outcomes the generating-parameter reference uses independent
dense Gaussian conditioning. For counts the main cohort uses the same
Laplace score and tilted-moment calculation at fitted and generating
parameters. This paired comparison measures the practical consequence
of estimating parameters under these designs while holding the numerical
method fixed. The separate sampler cohort below checks the approximation
at both parameter vectors. Neither comparison isolates a single
parameter's effect on prediction.

These outcomes support the public workflow as a way to test a specified
model on complete held-out locations. They do not establish general
predictive calibration or stable estimation of the range and dispersion.
The cases use the correct response and covariance family. They do not
cover covariance misspecification, estimation of ecological scales, or
new geographic study populations. Failures remain part of the practical
availability denominator, and all score differences refer to the paired
scored subset.

## Reproduce

From the package root, make sure the development package and matched
native providers are available. Run each command with a fresh output
directory. Two single-threaded workers were used for the disjoint primary
parts and then for the disjoint sampler parts. The session records include
package versions and source fingerprints. The first sampler attempt in each
shard completed its first generating reference, then stopped on a driver
column-name error before reading a fitted status. Those partial files remain
under `results/reference1_initial_interrupted` and
`results/reference2_initial_interrupted`. The corrected runner repeated the
same prespecified cases and seeds.

```sh
Rscript tools/benchmarks/spatial_field_cv_fitted/run.R /path/to/results/part1 1 2
Rscript tools/benchmarks/spatial_field_cv_fitted/run.R /path/to/results/part2 2 2
Rscript tools/benchmarks/spatial_field_cv_fitted/summarize.R /path/to/results
Rscript tools/benchmarks/spatial_field_cv_fitted/reference.R /path/to/results /path/to/results/reference1 1 2
Rscript tools/benchmarks/spatial_field_cv_fitted/reference.R /path/to/results /path/to/results/reference2 2 2
Rscript tools/benchmarks/spatial_field_cv_fitted/reference-summary.R /path/to/results
```

- [All simulated fields and outcomes](results/scenarios-combined.csv)
- [Native attempts](results/starts-combined.csv)
- [Complete native audit packets](results/native-audits.zip)
- [Audit packet hashes](results/native-audit-index.csv)
- [Predictive rows](results/predictions-combined.csv)
- [Full availability denominators](results/availability-combined.csv)
- [Paired joint score summary](results/scores-combined.csv)
- [Moment and interval summary](results/calibration-combined.csv)
- [Parameter recovery](results/recovery-combined.csv)

- [Independent count reference scores and failures](results/reference-scores-combined.csv)
- [Independent count reference moments and randomized PIT](results/reference-moments-combined.csv)
- [Count reference score summary](results/reference-score-summary.csv)
- [Count reference moment summary](results/reference-moment-summary.csv)

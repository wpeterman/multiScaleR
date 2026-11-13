# Replicated mixed-model CV numerical assessment

The assessment ran all 224 planned folds through `cross_validate_mixed()` with
the default native optimizer. It crossed seven model configurations, four data
scenarios, four seeds, and two prediction targets. Each dataset had 18 groups
and six visits per group. New-group folds held out four entire groups.
Existing-group folds held out the last two visits per group. In the sparse
scenario, half the groups contributed only one analysis visit to the
existing-group fit. Weak and strong intercept SDs were 0.12 and 1.15 on the
linear predictor scale. Sparse data used a lower baseline response, and the
extreme scenario perturbed held-out responses. Grouped binomial rows had 4,
12, or 25 trials. All models included a training-standardized predictor and a
log exposure offset.

## Independent references

The Gaussian reference formed each group's dense observation covariance and
computed its conditional mean and joint Normal density with direct linear
solves. For GLMMs, the reference integrated native R response densities over
the fitted Gaussian intercept using composite eight-point Gauss–Legendre
quadrature. It compared a fixed range of 16 prior SDs on either side of zero
with panels 0.125 SD wide against a range of 20 SDs with panels 0.0625 SD wide.
The reference conditioned on analysis responses for existing groups and used
the fitted prior directly for new groups. It evaluated the predictive mean
directly, and computed the joint negative log score as the analysis log
integral minus the combined analysis and assessment log integral. Refinement
had to agree within 1e-8 after scaling by 1 + the absolute reference value.
The largest observed scaled refinement discrepancy was 2.18e-14.

## Availability and error

214 of 224 folds scored. The native optimizer reported singular convergence
or a nonpositive Hessian in the ten failed folds. All ten remain in
[fold-results.csv](fold-results.csv). Their 152 planned group assessments
remain in [group-results.csv](group-results.csv) with the reference marked
as unrun. References were available for all 2,312 groups in scored folds.

| Model | Folds scored / planned | Available group references |
| --- | ---: | ---: |
| Gaussian random intercept | 32 / 32 | 352 |
| Gaussian random intercept and slope | 24 / 32 | 236 |
| Poisson random intercept | 32 / 32 | 352 |
| Bernoulli random intercept | 32 / 32 | 352 |
| Grouped binomial random intercept | 32 / 32 | 352 |
| NB2 random intercept | 30 / 32 | 316 |
| Gamma random intercept | 32 / 32 | 352 |

The eight Gaussian random-slope failures included four weak, two sparse, one
strong, and one extreme scenario. The two NB2 failures occurred in sparse
existing-group folds. New-group availability was 110 of 112 folds.
Existing-group availability was 104 of 112 folds.

| Error against reference | Comparisons | 95th percentile | 99th percentile | Maximum |
| --- | ---: | ---: | ---: | ---: |
| Absolute joint negative log score | 2,312 | 7.11e-15 | 2.27e-13 | 3.29e-10 |
| Absolute response mean | 6,384 | 3.55e-15 | 1.78e-14 | 1.18e-9 |
| Scaled response mean | 6,384 | 1.05e-15 | 3.76e-15 | 1.60e-10 |

The largest score difference was in sparse Poisson fold `fold129`, group 14.
The largest absolute and scaled mean differences were in extreme NB2 fold
`fold195`, group 8, row `g08v06`. The scaled error divides the absolute
difference by 1 + the absolute reference value. These comparisons cover the
held-out calculations at the same fitted training parameters. Training
parameter recovery, selection performance, and predictive calibration
require separate assessments.

## Reproduce the run

Run both scripts from the package source root with a compatible glmmTMB and
TMB installation. The recorded run used R 4.6.1, glmmTMB 1.1.14, and TMB
1.9.21. [Run metadata](run-metadata.txt) records the source fingerprints.

```sh
Rscript --vanilla tools/benchmarks/mixed_cv_replicated/run.R
Rscript --vanilla tools/benchmarks/mixed_cv_replicated/summarize.R
```

The [plan](plan.csv) lists every intended fold. The
[availability table](availability.csv) shows counts by model, scenario, and
target. The [error table](error-summary.csv) gives quantiles by model.
[Row results](row-results.csv) preserve the scored response means and their
references.

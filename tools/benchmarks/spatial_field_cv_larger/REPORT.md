# Numerical assessment of larger non-Gaussian spatial fields

## Study target

The reference samples the latent posterior at frozen model parameters. It
uses separately evaluated normalized response densities and an independent
elliptical slice sampler. Adjacent likelihood-power bridges give the joint
assessment loss. The analysis posterior gives response means and variances
after integration over new-location innovations.

The 28 larger-field cases contained eight or 16 training nodes and four or
six assessment nodes. They included one or four observations per node, two
field SDs, two covariance functions, and weak and strong nearest-neighbor
correlation. Four difficult small-field quadrature cases checked the
reference. See the [protocol](PROTOCOL.md) for fixed settings and the
reference precision criteria.

## Primary results

All 28 package calculations returned. The independent reference resolved 27
joint scores, 330 of 332 response means, and 330 of 332 observation variances.
All temperature chains used independent numerical Gaussian proposals. The
maximum split rank-normal R-hat was 1.0095, and the lowest adjacent bridge
overlap was 0.0754.

Among resolved joint references, the median, 95th percentile and maximum
absolute differences were 0.00112, 0.00763 and 0.04515 nats per assessment
observation. The median reference Monte Carlo standard error was 0.00219 nats
per observation. Two differences exceeded twice their Monte Carlo standard
error. The larger difference came from a sparse NB2 field with Gaussian
correlation, field SD 1.8, and nearest-neighbor correlation 0.25
(`larger_016`). Its reference and package block losses were 14.7520 and
14.9326 across four rows, a difference of 0.04515 nats per row with MCSE
0.00800. The other case (`larger_008`) was an NB2 field with 16 training
nodes, six assessment nodes and stronger dependence. Its difference was
0.00774 nats per observation with MCSE 0.00356.

The sparse Poisson case `larger_015` missed the primary joint precision
criterion. Its apparent difference was 0.03970 nats per observation with
MCSE 0.01879. It is excluded from the resolved quantiles. One Poisson
response mean and one NB2 mean missed their precision criterion. Two
NB2 variance references also missed theirs. The raw records include those
references and their MCSEs.

| Response | Joint scores resolved | Joint loss error per observation, p95 / max | Mean error, p95 / max | Relative variance error, p95 / max |
|---|---:|---:|---:|---:|
| Poisson | 11/12 | 0.00685 / 0.00725 | 2.27% / 5.82% | 3.97% / 8.77% |
| NB2 | 12/12 | 0.02457 / 0.04515 | 2.69% / 3.82% | 3.22% / 5.99% |
| Bernoulli | 2/2 | 0.00110 / 0.00111 | 2.21 / 2.81 percentage points | 6.15% / 8.29% |
| Gamma | 2/2 | 0.00542 / 0.00562 | 2.12% / 3.00% | 5.52% / 8.21% |

Joint errors are absolute nats. Log-link mean errors are relative. These
differences contain approximation error and finite reference Monte Carlo
error. The Bernoulli and Gamma results have two independent scenarios
each. Their row-level tails describe these two scenarios.

The worst resolved Poisson mean had reference 0.35005, package 0.32968
and reference MCSE 0.00958. The maximum resolved Poisson variance
difference was 0.17476 versus 0.15943, with MCSE 0.00830. The independent
sampler saw counts as high as 59 in analysis and 43 in assessment for
the two sparse count cases with the largest apparent joint discrepancies.

## Reference checks

The small-field bridge scores agreed with independent quadrature to within
1.4 reference Monte Carlo standard errors in all four control cases.
Absolute block differences were 0.000694 nats for Poisson, 0.00977 nats
for NB2, 0.01551 nats for Bernoulli and 0.01025 nats for Gamma. The NB2
control mean and variance were too noisy for the primary moment precision
criteria. The separate normalized Gaussian bridge check differed from
its analytic log marginal likelihood by 0.00462 nats.

## Independently seeded precision follow-up

The longer Poisson reference for `larger_015` met the joint precision
criterion. The package score exceeded that reference by 0.01862 nats per
assessment observation with MCSE 0.01148. The independent NB2 reference
for `larger_016` also passed and reproduced the direction of the discrepancy.
Its difference was 0.02384 nats per observation with MCSE 0.00729,
compared with 0.04515 and MCSE 0.00800 in the primary run. The smaller
replicated magnitude makes the original maximum an uncertain estimate of
the exact approximation error. Both independent references indicate that
the package loss is higher in this case.

The longer analysis-posterior chains resolved all previously unresolved
larger-field moments. The Poisson mean MCSE was 2.04% of the mean. The
NB2 mean MCSE was 1.16%, and the two NB2 variance MCSEs were 3.84% and
3.61%. These follow-up results remain separate from the primary summaries.

## Interpretation

The study gives supporting evidence for the approximate joint score and
tilted moments at 12 and 22 field locations. It also identifies a sparse
NB2 case whose score discrepancy exceeds the reference precision by
several standard errors. Differences on that scale can matter when
ranking closely competing models. Moment accuracy requires a separate
assessment from the joint score.

The reference used known parameters and observed outcomes from a bounded
set of simulated fields. These comparisons do not establish field recovery,
uncertainty in fitted parameters, calibration, or accuracy under a new
study design. Stronger dispersion, more extreme correlations, and larger
or differently configured assessment fields remain outside this study.

## Replay after the response-moment arithmetic fix

The final package evaluation reproduced all 28 larger-field outcomes
without resampling the independent references. The maximum absolute score
change was 4.62e-14 nats. Changes in the response means and variances,
divided by the larger of one and the previous absolute value, were below
4.45e-15. Regression tests separately cover near-certain Bernoulli outcomes
and representable log-link variance when a raw second moment overflows.

## Reproduce

Run each command from the package root. The scripts use the development
package and the earlier quadrature records in this repository.

```sh
Rscript tools/benchmarks/spatial_field_cv_larger/run.R /path/to/fresh-results
Rscript tools/benchmarks/spatial_field_cv_larger/summarize.R /path/to/fresh-results
Rscript tools/benchmarks/spatial_field_cv_larger/check-bridge.R
Rscript tools/benchmarks/spatial_field_cv_larger/followup.R /path/to/fresh-results
Rscript tools/benchmarks/spatial_field_cv_larger/replay-package.R
```

- [Design and seeds](results/design.csv)
- [Joint scores, failures and reference diagnostics](results/scenarios.csv)
- [Assessment-row moments and Monte Carlo precision](results/predictions.csv)
- [Temperature-specific diagnostics](results/temperatures.csv)
- [Grouped tail summaries](results/error-summary.csv)
- [Quadrature controls for joint scores](results/score-controls.csv)
- [Quadrature controls for response moments](results/moment-controls.csv)
- [Environment and source fingerprints](results/session.txt)
- [Independent follow-up joint scores](results/followup-scores.csv)
- [Independent follow-up moment precision](results/followup-moments.csv)
- [Independent follow-up temperature diagnostics](results/followup-temperatures.csv)

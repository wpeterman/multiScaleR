# Numerical accuracy of spatial response-field CV

## Joint held-out block scoring

All 1,152 planned non-Gaussian scenarios returned a finite Laplace score and a
converged independent quadrature reference. This includes 288 cases for each
family. Every zero-field case (192) and dense Gaussian reference check (288)
also completed. The maximum absolute joint loss difference was
`7.11e-15` nats at zero field variance and `1.42e-14` nats for the Gaussian
identity. The largest mean discrepancy in those checks was `7.55e-15`.

Over the non-Gaussian factorial design, the median, 95th percentile, 99th
percentile and maximum absolute joint score errors were 0.00240, 0.0304,
0.0485 and 0.0594 nats per block. Dividing each block error by its assessment
row count gave 0.000415, 0.0177, 0.0314 and 0.0483 nats per observation,
respectively. All 117 cases above 0.01 nats per observation had one visit per
location. With eight visits per location, the maximum was 0.00352 nats per
observation. These empirical quantiles weight the simulated design evenly and
do not estimate the frequency of errors in field data.

| Response | 95th percentile joint error per observation | Worst joint error per observation | 95th percentile response mean error | Worst response mean error |
|---|---:|---:|---:|---:|
| Poisson | 0.0131 | 0.0309 | 8.86% | 74.53% |
| Bernoulli | 0.0242 | 0.0483 | 1.96 percentage points | 4.42 percentage points |
| NB2 | 0.0223 | 0.0397 | 10.45% | 50.30% |
| Gamma | 0.0154 | 0.0175 | 2.31% | 15.79% |

Errors in the first two columns are absolute nats. Count and Gamma response
mean errors are absolute relative to the reference mean. Bernoulli mean
errors are absolute probability errors. Joint scores retain dependence
across locations and visits.

The largest absolute total block error was 0.0594 nats for a two-row NB2
assessment (`field_0539`). The worst per-observation error was 0.0483 nats
for one Bernoulli observation (`field_0586`, field SD 2 and correlation
0.20). One Gamma reference emitted an intermediate numerical warning
(`field_1008`). Its final integral met both checks. Across all quadrature
cases the maximum final refinement and proposal checks remained below the
predeclared `5e-6` tolerance. Four non-Gaussian cases used 80 nodes per
dimension.

## Response means need their own accuracy assessment

| Response and case | Reference mean | Gaussian-posterior mean | Absolute response difference | Absolute joint score error |
|---|---:|---:|---:|---:|
| Poisson `field_0357` | 0.78746 | 1.37436 | 0.58690 | 0.00495 nats |
| NB2 `field_0167` | 2.10753 | 3.16765 | 1.06012 | 0.03726 nats |
| Gamma `field_0792` | 0.009275 | 0.007810 | 0.001465 | 0.00934 nats |
| Bernoulli `field_0142` | 0.37078 | 0.41498 | 0.04420 | 0.00413 nats across 16 rows |

The Poisson, NB2 and Gamma rows give the highest relative mean errors for
their respective families. All three used field SD 2, correlation 0.85 and
one visit per location. The Bernoulli row gives the largest absolute
probability error. It used field SD 2, correlation 0.85, and eight visits at
each of two assessment locations.

The largest absolute errors across response rows were 2.35119 Poisson
counts (15.54343 versus 17.89462), 5.30230 NB2 counts (50.58472 versus
45.28242), and 18.06384 Gamma units (1425.56213 versus 1407.49829).
Their relative errors were 15.13%, 10.48%, and 1.27%, respectively. The
raw tables retain both error scales and the joint score for each case.

## Tilted log-link response moments

For a positive integer moment order, completing the square in the prior
moves the analysis integral to a shifted latent normal distribution. The
same field Cholesky factor evaluates the shifted likelihood. The first
follow-up used one additional Laplace solve per assessment node and compared
its mean with the independent quadrature mean. On the 5,184 log-link rows,
the worst relative errors at field SD 2 fell from 74.53% to 3.66% for
Poisson, from 50.30% to 4.16% for NB2, and from 15.79% to 0.335% for Gamma.
This calculation leaves the joint scoring likelihood ratio untouched.
The second-moment follow-up compared first moments, second moments and total
predictive variances for all 864 log-link scenarios and 5,184 assessment
rows. All quadrature references converged and all warm-started tilted
evaluations returned positive predictive variances. At field SD 2, the
largest relative second-moment errors fell from 405.3% to 4.65% for Poisson,
from 291.5% to 4.82% for NB2, and from 35.77% to 0.696% for Gamma.

| Family at field SD 2 | Gaussian approximation variance error, 95th percentile / max | Tilted variance error, 95th percentile / max |
|---|---:|---:|
| Poisson | 33.87% / 403.10% | 0.599% / 3.79% |
| NB2 | 48.17% / 295.91% | 0.805% / 4.51% |
| Gamma | 14.06% / 38.37% | 0.307% / 0.726% |

The executable Poisson vignette has 50 training nodes and six assessment
nodes. In 30 alternating-order timings at its fitted parameters, the
two moment orders added a median of 0.02 seconds. The benchmark times the
scoring stage after native fitting. The increment consists of two latent
solves for each distinct assessment location, reusing the training
covariance factor and warming each solve from the preceding mode in shifted
coordinates. Those timings provide evidence for the vignette configuration. Repeating
the full moment study through the package implementation reproduced the
prototype moment table and summary byte for byte and again added a median
of 0.02 seconds on the vignette.

## Scope

The results support explicit approximate labeling of the Laplace joint
block score for the four response families and use of tilted Laplace
moments for the three log-link families. Bernoulli predictions retain the
Gaussian-posterior probability integral and had a worst absolute
probability difference of 4.42 percentage points in this design.
Differences on the scale of the observed numerical errors can alter
rankings between closely competing models.

These are known-parameter comparisons with up to three field locations
and up to 16 assessment observations. They leave larger fields, extreme
correlation, different observation dispersions, nonstationary fields,
parameter estimation or selection uncertainty, and predictive calibration
untested. Separate package fixtures verify native fitting, the held-out
likelihood calculation, leakage, field covariance, and ecological scale
refitting for their declared datasets. Those checks do not measure the
frequency of convergence failures on this simulation design.

## Reproduce the study

Run from the package root with the development package and glmmTMB loaded
through the usual R environment. See [the pre-run design](PROTOCOL.md) and
[recorded environment and package source fingerprints](results/session.txt).

```sh
Rscript tools/benchmarks/spatial_field_cv_quadrature/run.R /path/to/fresh-results
Rscript tools/benchmarks/spatial_field_cv_quadrature/summarize.R /path/to/fresh-results
Rscript tools/benchmarks/spatial_field_cv_quadrature/tilted-means.R /path/to/fresh-results
Rscript tools/benchmarks/spatial_field_cv_quadrature/tilted-moments.R /path/to/fresh-results
```

- [Prespecified scenarios and seeds](results/design.csv)
- [All score outcomes and quadrature diagnostics](results/scenarios.csv)
- [All response mean comparisons](results/predictions.csv)
- [Tail error summaries by family and design](results/error-summary.csv)
- [First-moment tilt comparisons](results/tilted-means.csv)
- [First-moment tilt summary](results/tilted-summary.csv)
- [First and second moments and observation variance](results/tilted-moments.csv)
- [Tilted-moment tail errors](results/tilted-moments-summary.csv)
- [Fitted-vignette prototype timing pairs](results/tilted-runtime.csv)
- [Fitted-vignette implementation timing pairs](results/tilted-runtime-implementation.csv)
- [Tilted implementation environment and source fingerprints](results/session-tilted.txt)

The first execution recorded all outcomes, then reported a parse error when
the runner was refactored before its source file was closed. A clean full
execution of the finalized runner reproduced the design, outcomes,
predictions and environment records byte for byte. All analyses use those
complete records.

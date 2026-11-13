# Bayesian approximation assessment

## Scope and decision

The replicated assessment covered Gaussian, Poisson, grouped binomial and NB2
responses with one annual scale-group SD. Each family used an annual model
with a Gaussian kernel and a repeated-site model with an exponential kernel.
The repeated-site design included a scale trend and a response random
intercept. Both importance-corrected proposals agreed closely with NUTS on
the prespecified fitted-scale probabilities when both fits were eligible.
Availability and predictive-score precision varied by family and design.

Both proposals passed the availability criterion in seven panels. In the
repeated-site NB2 panel, only the conditional coefficient proposal passed that
criterion. Annual Poisson had 46 eligible NUTS references, so the matched
comparison describes that subset. These results apply to the simulated
designs, priors and computation budgets.

## Gaussian annual qualification cohort

| Method | Eligible fits | One-sided 95% lower bound for eligibility | Median seconds per fit |
|:--|--:|--:|--:|
| NUTS | 48 of 48 | 0.939 | 9.5 |
| Joint importance proposal | 48 of 48 | 0.939 | 7.1 |
| Conditional coefficient proposal | 48 of 48 | 0.939 | 27.0 |

The 144 fits used 48 independent prior-predictive datasets. The
[fit and calibration tables](results/annual/) preserve the outcomes. No
fits were rejected or failed. Runtime describes this Windows run and is affected by
the configured priors, design, budgets, and concurrent work.

Each fitted year contributed two probability targets at 150 and 350 map units.
Among 288 possible comparisons for each proposal, 287 had estimable combined
Monte Carlo error. The median absolute differences from NUTS were 0.00376
for the joint proposal and 0.00336 for conditional coefficients. The largest
differences were 0.01981 and 0.02274, respectively. Neither proposal had a
comparison above the declared material-difference threshold.

| Quantity | NUTS 90% coverage | Joint proposal 90% coverage | Conditional coefficients 90% coverage |
|:--|--:|--:|--:|
| Habitat coefficient | 41 of 48 | 41 of 48 | 40 of 48 |
| First-year log distance | 44 of 48 | 44 of 48 | 44 of 48 |
| Between-year scale SD | 43 of 48 | 44 of 48 | 44 of 48 |

All exact 95% coverage intervals contained 0.90. The largest empirical CDF
deviation from uniformity was 0.166 and no monitored curve exceeded its
descriptive 0.223 simultaneous band. Finite posterior Monte Carlo error
also contributes to these empirical CDFs.

Joint future response scores were averaged over coherent posterior vectors.
The median log-loss differences from NUTS were 0.00038 for the joint proposal
and 0.00079 for conditional coefficients. The largest absolute differences
were 0.375 and 0.367, both on dataset `gaussian-annual-33`. The planned study
retains these scores.

A separate diagnostic integrated the new year effect by Gauss-Hermite
quadrature. From 128 to 256 nodes, all four log scores changed by less than
0.001. At 256 nodes the original NUTS fit scored 8.765, the joint proposal
8.740, and conditional coefficients 8.696. A longer independent NUTS fit
scored 8.705. Relative to that run, the proposal differences were 0.035 and
-0.009. The original and independent NUTS results still differ by 0.060
after quadrature. This illustrates finite posterior simulation precision in
a rare held-out block. The diagnostic retains its separate seed and
does not replace the original paired comparison.

## Replication across response families and repeated-site designs

All eight panels completed with 48 independent datasets each and three
methods per dataset. The full archive contains 1,152 fits. Five were
rejected by their existing numerical checks and no execution failed.

| Family and design | Eligible NUTS | Eligible joint | Eligible conditional coefficients | Largest paired probability difference |
|:--|--:|--:|--:|--:|
| Gaussian annual | 48 | 48 | 48 | 0.023 |
| Poisson annual | 46 | 48 | 48 | 0.030 |
| Grouped binomial annual | 48 | 48 | 48 | 0.025 |
| NB2 annual | 48 | 48 | 48 | 0.021 |
| Gaussian repeated site | 48 | 48 | 47 | 0.024 |
| Poisson repeated site | 48 | 48 | 48 | 0.021 |
| Grouped binomial repeated site | 48 | 48 | 48 | 0.033 |
| NB2 repeated site | 48 | 46 | 48 | 0.025 |

The exact one-sided 95% eligibility lower bound is 0.939 for 48 of 48,
0.905 for 47 of 48, and 0.875 for 46 of 48. Repeated-site NB2 joint
importance fell short of the prespecified 0.90 availability bound.
Annual Poisson NUTS also fell short, with 46 paired datasets remaining.

None of the 5,282 estimable fitted-scale comparisons passed the threshold for
a material discrepancy. There were 94 comparisons without estimable combined
Monte Carlo error. The 72 coverage summaries use eligible fits and three
quantities per panel and method. For repeated-site Poisson, the first-year
scale was inside the nominal 90% interval in 38 of 48 NUTS fits, with
an exact 95% interval of 0.650 to 0.895 for coverage. Both proposals covered
39 of 48. The other 71 exact intervals contained 0.90. No monitored
empirical CDF exceeded its descriptive simultaneous band.

The largest original joint predictive-score difference across the eight
panels was 1.918 nats on one three-observation repeated-site Poisson block.
The absolute median paired score difference was below 0.003 nats for both
proposals in every panel. The one augmentation per posterior vector used in
the primary analysis supplies a finite simulation estimate. Interpret
predictive-score tails separately from fitted-scale probabilities.

In the repeated-site Poisson case, independent Gaussian quadrature integrated
the shared new annual and site effects. Increasing the order from 48 to 72
nodes on each axis changed the block loss by less than 0.001 for every fit.
At order 72, the original NUTS fit scored 15.030, the joint proposal 15.210,
the conditional proposal 15.029, and the independently seeded longer NUTS
fit 14.858. The primary difference of 1.918 nats shrank after integration.
The gap between the two NUTS fits and the remaining differences among
methods show that this rare block is sensitive to posterior simulation as
well as the new-effect integration budget. The [diagnostic records](results/replicated/poisson-tail/)
retain all estimates.

The [replicated tables and fit archive](results/replicated/) preserve all
outcomes, numerical diagnostics, and derived probability records. The
independent verification script reconstructs eligibility, coverage and
the paired discrepancy counts from the archived tables. Scope includes the
tested priors, kernels, group structure and importance budgets.

## Eight-dataset pilot

The [pilot](results/pilot/) used a disjoint seed series with one dataset for each combination
of response family and design. All 24 fits passed their own numerical screens.
For the 112 fitted-scale probability comparisons across both proposals, the
largest absolute difference from NUTS was 0.0165. The largest joint log-loss
difference was 0.0314. These results establish that the simulator, fitting
routines, and assessments ran across the tested cases. A single dataset per
cell gives limited information on calibration and availability.

The pilot's compact records were retained before regenerating its summary
from the frozen pilot script. The source and configuration fingerprint is
recorded with each output directory. The assessment script produces coverage
and paired probability tables and rank-CDF figures.

# Fitted Bayesian CV assessment

This study uses 16 independent simulated datasets, two training splits, and
two fitted candidates. It covers Gaussian and negative binomial (NB2)
responses under common and annual bandwidth generators. The [protocol](PROTOCOL.md)
describes the frozen inputs, refits, numerical references, and oracle repair.

## Fit availability and replay

All 64 training refits passed the package's raw and derived NUTS diagnostics
and produced finite assessment scores without captured warnings. Each cell
has four independent datasets, and there are no missing or rejected fold
attempts. The independent calculation checked the response order, point IDs,
bin distances, reference center and spread, and new-group flags in every
fit. Using the original score seed, it replayed all package scores to within
1.11 × 10^-13 nats per block. The fitted quadrature refinement changed observed
scores by at most 2.61 × 10^-9 and stress scores by at most 1.78 × 10^-7.

The original independent oracle returned 96 of 128 finite joint scores
and 384 of 512 finite marginal probabilities. Every original missing oracle
value came from a whole-year fold under the annual generator. The separately
fingerprinted numerical repair returned all 128 scores and 512 probabilities.
The repaired results agreed with formerly finite results to within
7.11 × 10^-14 nats and 3.22 × 10^-15 in probability. The maximum oracle
quadrature refinement difference was 6.42 × 10^-6 nats for observed outcomes
and 0.00161 for stress outcomes. These oracle results use the repaired
values. The original outputs remain in the frozen archive. All fitted-model
scores and replay comparisons are from the original production run.

## Score precision

The fitted-score reference integrates over the new effects conditional on
the same 6,000 posterior vectors used by the package. The budget experiment
therefore measures new-group integration error. It does not include Monte
Carlo error in fitting the posterior. The scores refer to one whole
eight-visit block. Divide by eight to get the per-visit scale in `cv_tab()`.
The tables include diagnostics-eligible blocks only and four scoring
seeds per latent block.

| Response | Outcomes | New draws per posterior vector | RMSE | 95th percentile absolute error | Max absolute error |
|:--|:--|--:|--:|--:|--:|
| Gaussian | Observed | 1 | 0.0261 | 0.0510 | 0.0980 |
| Gaussian | Observed | 8 | 0.0078 | 0.0156 | 0.0327 |
| Gaussian | Observed | 32 | 0.0038 | 0.0079 | 0.0100 |
| Gaussian | Shifted | 1 | 0.4376 | 0.9363 | 1.5404 |
| Gaussian | Shifted | 8 | 0.2040 | 0.3946 | 1.1673 |
| Gaussian | Shifted | 32 | 0.1067 | 0.2626 | 0.4704 |
| NB2 | Observed | 1 | 0.0146 | 0.0298 | 0.0617 |
| NB2 | Observed | 8 | 0.0054 | 0.0110 | 0.0206 |
| NB2 | Observed | 32 | 0.0025 | 0.0052 | 0.0095 |
| NB2 | Shifted | 1 | 0.3396 | 0.4579 | 2.6602 |
| NB2 | Shifted | 8 | 0.0841 | 0.1867 | 0.3497 |
| NB2 | Shifted | 32 | 0.0368 | 0.0668 | 0.1727 |

Each table row uses 24 assessment blocks from eight independent datasets.
The 96 scored streams within a row share data and posterior vectors. A
Gaussian shifted score differed by more than 0.1 nats per visit in 7 of 96
streams at one draw, and 1 of 96 at eight draws. The NB2 counts were 2 and
0, respectively. At 32 draws, neither family exceeded that descriptive
band, though 19 Gaussian and 3 NB2 streams still missed the reference by
more than 0.1 nats per block. The worst NB2 stream was a whole-year annual
candidate for dataset `nbinom2-02`. Its shifted score exceeded the reference
by 2.66 nats per block, and the density contribution count was 1.86 among
6,000 equally weighted posterior vectors. The minimum contribution counts
were 1.00 Gaussian and 1.06 NB2 under the shift, versus 1,044 and 417 for
observed outcomes. Increasing the augmentation budget only addresses the
integration over new effects. It cannot add posterior mass in an
undersampled parameter region.

## Predictive comparison and marginal coverage

The table reports the annual candidate's score minus the common candidate's
score using the same assessment visits and fitted quadrature reference.
Negative values favor the annual candidate. The SD is computed across four
independent paired datasets in each row. The descriptive paired t intervals
in `summary-paired_summary.csv` all span zero.

| Response | Generating bandwidth | Fold | Mean difference (nats per visit) | Dataset SD | Annual lower in |
|:--|:--|:--|--:|--:|--:|
| Gaussian | Common | Whole year | -0.0148 | 0.0495 | 2 of 4 |
| Gaussian | Common | Seen years, new site | -0.0046 | 0.0081 | 2 of 4 |
| Gaussian | Annual | Whole year | +0.0059 | 0.0131 | 1 of 4 |
| Gaussian | Annual | Seen years, new site | -0.0276 | 0.0412 | 2 of 4 |
| NB2 | Common | Whole year | +0.0005 | 0.0069 | 1 of 4 |
| NB2 | Common | Seen years, new site | +0.0040 | 0.0062 | 2 of 4 |
| NB2 | Annual | Whole year | -0.0072 | 0.0447 | 2 of 4 |
| NB2 | Annual | Seen years, new site | +0.0151 | 0.0637 | 2 of 4 |

The score at the package seed reversed three of the 32 dataset-by-fold
comparisons relative to the same-fit reference. All three were new-site
predictions in seen years, with reference gaps of 0.0021 to 0.0047 nats per
visit. Of 128 matched seed comparisons at each budget, the signs differed
from the reference 11 times at 1 draw, 4 times at 8 draws, and 5 times at 32
draws. The scores can remain practically tied even as individual errors shrink.

Against the known-parameter oracle, the fitted-model mean observed score
gaps are tabulated by response family, generating bandwidth, candidate and
fold in `summary-oracle_gaps.csv`, along with SDs and extrema. The oracle
knows the generating global parameters and training-year bandwidths, which
helps separate the limited fit information from the new-effect integration.
Mean central 90% marginal coverage across the individual cells was 87.5% to
96.9% for Gaussian fitted models and 84.4% to 93.8% for NB2 fitted models.
The corresponding oracle ranges were 90.6% to 93.8% and 87.5% to 90.6%.
Each cell has 32 visits from four independent datasets, so these are small
descriptive checks. Marginal coverage does not measure joint coverage.

## Interpretation

The Gaussian and NB2 panels examine joint forecasting of one new year at
known sites and one new site in known years. The former tests annual
bandwidth extrapolation. The latter tests whether the fitted annual
bandwidths aid prediction across sites, after integrating a new response
effect. The two questions can give different score differences. Choose
the assessment target that matches the intended use.

Start with the default augmentation budget for ordinary outcomes. For a
close comparison, or if held-out responses are jointly rare, increase
`new_group_draws` for both candidates with the same seed and fit settings
and check whether the score difference is stable at the required precision.
A second seed checks the repeatability of both fitting and augmentation.
The fixed-fit stress cases isolate one source of error, and they do not
estimate how often such extreme blocks occur. This study uses one new
latent effect per block, Gaussian kernels, and NUTS with equal-weight draws.
The [earlier integration benchmark](../cv_new_group_integration/REPORT.md)
addresses multiple new effects and weighted posterior mixtures. Neither
study justifies a universal draw budget or a claim of calibration across
other data-generating regimes.

## Archive verification

`verify.R frozen` verified all 32 fit-call records, 64 fold attempts,
128 joint references, 1,184 augmentation scores and 512 marginal
probabilities against the archived records. No package methods, thresholds
or defaults were changed in this assessment.

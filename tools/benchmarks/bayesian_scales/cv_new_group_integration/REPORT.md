# Joint new-group prediction in Bayesian cross-validation

This assessment isolates Monte Carlo integration over unseen **response**
random intercepts in `cross_validate_bayes()`. The posterior vectors and
weights are fixed, and the scores are conditional on those finite mixtures.
Fitting error, scientific predictive calibration, and unseen ecological-scale
effects require separate evidence.

## Design

Each assessment set contains eight observations, including repeated visits
within two unseen sites and observations at two known sites. The cases use one
joint eight-row block or two four-row blocks. Gaussian cases vary the response
group SD and observation SD. Count cases use either sparse outcomes or one
unseen site with 36, 44, and 39 events, against a baseline expected count of
about 0.4. This last case is a deliberate stress test.

The main run fixes 96 synthetic posterior vectors and evaluates 1, 4, 16, and
64 new-effect draws per vector with 32 independent seeds per case and budget.
Two cases use uniform weights. Three assign extra mass to one vector to
represent concentrated importance weights. A second run uses 4,000 vectors,
the nominal count retained by the default four chains and 2,000 iterations
with half discarded as warmup. That run uses 16 seeds and one new-effect
draw per vector for four selected cases. A third run uses 4,000 vectors and
16 seeds for three diffuse weighted mixtures. It spreads 14% of the total
weight across the first 32 vectors, then distributes the rest uniformly. The
resulting weight ESS is 1,166 and the maximum weight is 0.00459. These meet
the existing importance fit's raw ESS and maximum-weight cutoffs. The
concentrated synthetic mixtures exceed the maximum-weight cutoff and serve
as stress tests. All tiers condition on fixed posterior vectors, isolating
new-effect simulation variability from posterior sampling variability.

## References and scoring

The Gaussian reference integrates independent unseen intercepts using
a conditional multivariate normal covariance per posterior vector. It adds
the residual variance to the diagonal and the group variance to every entry
for observations sharing an unseen group. The Poisson reference integrates each unseen group over
its standard-normal effect with adaptive quadrature, centered and scaled
around the integrand mode. Both references multiply all observation densities
for a block within each vector, then take their posterior-weighted average.
Conditional expected responses use the Gaussian linear mean or the Poisson
lognormal mean. The script independently evaluates the kernel, training
standardization, and known-site response effects from its synthetic inputs.

Fold and block log scores are divided by their observation counts. Signed
error is computed minus reference negative log score, so a positive error
means the Monte Carlo score overstated loss. RMSE, empirical 95th percentile
of absolute error, maximum absolute error, and the fraction of absolute errors
above 0.1 are reported. The 0.1 band is descriptive and represents about a
10% multiplicative difference in the geometric mean density per observation.

## Results

All 752 scoring calls produced finite estimates. The Poisson
reference changed by at most 7.1 × 10⁻¹⁴ when quadrature tolerances were
tightened. The reference also reproduced each univariate Gaussian density
within 1 × 10⁻¹¹. A replay of the main tier reproduced all 11,776
records exactly, excluding elapsed time.

### Budget sensitivity with 96 vectors

| Case | New draws per vector | Weight ESS | Mean error | Error SD | 95th percentile absolute error | Maximum absolute error | Fraction above 0.1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Gaussian, extreme, weighted | 1 | 13.84 | 9.6003 | 5.8967 | 18.1576 | 23.5094 | 93.8% |
| Gaussian, extreme, weighted | 4 | 13.84 | 4.5271 | 3.0565 | 9.8763 | 13.1263 | 96.9% |
| Gaussian, extreme, weighted | 16 | 13.84 | 1.4612 | 1.4822 | 3.9795 | 5.2027 | 93.8% |
| Gaussian, extreme, weighted | 64 | 13.84 | 0.4375 | 0.3428 | 1.0150 | 1.2820 | 90.6% |
| Gaussian, large group SD | 1 | 96.00 | 0.0198 | 0.0768 | 0.1388 | 0.2417 | 21.9% |
| Gaussian, large group SD | 4 | 96.00 | -0.0003 | 0.0358 | 0.0716 | 0.0952 | 0.0% |
| Gaussian, large group SD | 16 | 96.00 | 0.0034 | 0.0152 | 0.0280 | 0.0330 | 0.0% |
| Gaussian, large group SD | 64 | 96.00 | -0.0018 | 0.0096 | 0.0191 | 0.0234 | 0.0% |
| Gaussian, small group SD | 1 | 96.00 | 0.0004 | 0.0060 | 0.0119 | 0.0132 | 0.0% |
| Gaussian, small group SD | 4 | 96.00 | 0.0009 | 0.0032 | 0.0058 | 0.0079 | 0.0% |
| Gaussian, small group SD | 16 | 96.00 | 0.0003 | 0.0015 | 0.0027 | 0.0047 | 0.0% |
| Gaussian, small group SD | 64 | 96.00 | 0.0004 | 0.0005 | 0.0011 | 0.0012 | 0.0% |
| Poisson, extreme, weighted | 1 | 13.84 | 6.2783 | 5.6343 | 17.5536 | 19.7627 | 96.9% |
| Poisson, extreme, weighted | 4 | 13.84 | 1.0886 | 1.3329 | 3.1322 | 5.0989 | 84.4% |
| Poisson, extreme, weighted | 16 | 13.84 | 0.1660 | 0.3075 | 0.8215 | 1.1557 | 50.0% |
| Poisson, extreme, weighted | 64 | 13.84 | -0.0004 | 0.0567 | 0.1252 | 0.1348 | 9.4% |
| Poisson, sparse, weighted | 1 | 30.60 | 0.0005 | 0.0147 | 0.0250 | 0.0326 | 0.0% |
| Poisson, sparse, weighted | 4 | 30.60 | 0.0002 | 0.0078 | 0.0143 | 0.0179 | 0.0% |
| Poisson, sparse, weighted | 16 | 30.60 | 0.0007 | 0.0036 | 0.0064 | 0.0093 | 0.0% |
| Poisson, sparse, weighted | 64 | 30.60 | 0.0001 | 0.0018 | 0.0034 | 0.0042 | 0.0% |

### One augmentation with 4,000 vectors

| Case | Weight ESS | Mean error | Error SD | 95th percentile absolute error | Maximum absolute error | Fraction above 0.1 |
|---|---:|---:|---:|---:|---:|---:|
| Gaussian, extreme, weighted | 15.94 | 0.4884 | 0.4653 | 1.2567 | 1.2841 | 87.5% |
| Gaussian, large group SD | 4000.00 | 0.0053 | 0.0090 | 0.0195 | 0.0211 | 0.0% |
| Gaussian, small group SD | 4000.00 | 0.0002 | 0.0007 | 0.0015 | 0.0018 | 0.0% |
| Poisson, extreme, weighted | 15.94 | 0.0161 | 0.1064 | 0.1827 | 0.3556 | 18.8% |

### One augmentation with 4,000 diffuse weighted vectors

| Case | Mean error | Error SD | 95th percentile absolute error | Maximum absolute error | Fraction above 0.1 |
|---|---:|---:|---:|---:|---:|
| Gaussian, extreme | 0.5418 | 0.5968 | 1.6863 | 1.9888 | 81.2% |
| Gaussian, large group SD | 0.0038 | 0.0130 | 0.0185 | 0.0190 | 0.0% |
| Poisson, extreme | 0.0602 | 0.1288 | 0.2671 | 0.3498 | 37.5% |

The largest fold error was 23.509 nats per observation in the weighted
extreme Gaussian case at one draw per vector and seed 351117. The exact
score was 1.450 and the estimated score was 24.959. The largest marginal
row error was 61.692 nats at seed 371105 in the weighted extreme Poisson
case. The complete archive retains the individual scores and seeds.

## What the default supports

One augmentation per posterior vector is a low-cost starting budget.
With 4,000 uniformly weighted vectors, the two Gaussian cases with
regular outcomes had maximum observed fold errors of 0.00181 and
0.0211 nats per observation. The sparse count case was also stable in
the 96-vector panel at a weight ESS of 30.6. The 4,000-vector concentrated
extreme cases had a weight ESS of about 16 and showed large errors. Diffuse
weighting improved precision for the large-SD Gaussian case with regular
outcomes. The extreme Gaussian and Poisson cases still exceeded 0.1 nats per
observation in 13 of 16 and 6 of 16 repeats, respectively, at a weight ESS
of 1,166. Weight ESS and a finite score alone are insufficient for a close
model comparison.

The average of conditional densities is unbiased for the fixed mixture
under independent augmentation. Taking its negative logarithm can
introduce upward bias. Increasing the budget substantially reduced error
in the two extreme scenarios, but 64 draws did not reliably bring the
extreme Gaussian error below 0.1. For rare joint outcomes or concentrated
weights, rerun with the same seed and fitting settings and increase
new_group_draws until the scores and candidate ranking stabilize to the
precision needed for the comparison. If they do not stabilize at a
practical budget, treat the comparison as unresolved. The budget needed
for stable scores varies with outcome rarity and the posterior weights.

## Reproduce

From the package root, load a compatible compiled development DLL and run:

```text
Rscript tools/benchmarks/bayesian_scales/cv_new_group_integration/run.R
Rscript tools/benchmarks/bayesian_scales/cv_new_group_integration/run.R tools/benchmarks/bayesian_scales/cv_new_group_integration/default_vectors 16 4000 1 gaussian_small,gaussian_large,gaussian_extreme_weighted,poisson_extreme_weighted
Rscript tools/benchmarks/bayesian_scales/cv_new_group_integration/run.R tools/benchmarks/bayesian_scales/cv_new_group_integration/diffuse_vectors 16 4000 1 gaussian_large_diffuse_weighted,gaussian_extreme_diffuse_weighted,poisson_extreme_diffuse_weighted
```

Each output directory contains `all-results.csv.gz` with every estimate,
reference, seed, and outcome status, `summary.csv` with the grouped error
statistics, and `session-info.txt` with the software environment.

The first two archived tiers predate the redundant maximum-weight column
in the runner. Their weights are fixed by the same scenario definitions.

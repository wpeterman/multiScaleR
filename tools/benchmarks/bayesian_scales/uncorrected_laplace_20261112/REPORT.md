# Uncorrected Laplace approximations for Bayesian scale models

This bounded assessment complements the replicated importance-correction
study and the existing annual NB2 vignette comparison. It compares uncorrected
joint Laplace and uncorrected marginal Laplace with NUTS on the same newly
generated training and future observations, using the same model, priors,
offsets and ecological covariates. Agreement is measured conditional on the
tested settings. A numerical diagnostic pass indicates an eligible fit.

## Design

The twelve independent prior-predictive datasets cover three scenarios with
four replicates each.

| Scenario | Training design | Response |
|---|---|---|
| Annual | Six sites by three years, one partially pooled annual scale and Gaussian kernel | Gaussian |
| Repeated | Unbalanced panel of nine sites and four years, fixed scale trend, annual scale random intercept, site response random intercept and exponential kernel | Gaussian |
| Repeated | Same unbalanced design | NB2 |

The data generator is adapted from the sibling qualification study. It draws
the true parameter vector from the declared priors and independently
calculates kernel values, training normalization and response simulation.
The assessment set includes a known site-year and two observations sharing
a new site and new year. All methods use the same three observations. Seeds
910020 to 910240 distinguish the twelve replicates from the sibling study.

NUTS uses four chains, 3,000 iterations per chain, half warmup, target
acceptance 0.99 and maximum tree depth 12. Each approximation uses 6,000
posterior draws and its native optimizer controls. Marginal Laplace uses
the default uncorrected grid. Eligibility requires the native diagnostics
for the method. NUTS additionally requires passing derived scale diagnostics.
Warnings, failed fits and rejected fits remain in the archive.

## Outcomes

- Training landscape coefficient, first-year log distance, between-year
  scale SD and a new-year log distance. Retain 90% intervals, median,
  CDF at simulated truth and coverage for each dataset.
- Predeclared first-to-last annual probabilities of log distance corresponding
  to 150 and 350 map units. Compare posterior estimates on matched cases.
  A numerical difference is flagged when it exceeds both 0.05 and five
  combined Monte Carlo standard errors. Approximate-fit MCSE measures
  simulation variability conditional on its approximation. Undefined MCSE
  remains missing.
- Joint future negative log score over the three observations using coherent
  joint vectors and one augmentation per retained vector. These scores are
  descriptive. An isolated run has no separate estimate of augmentation
  precision.

Report the number of eligible, rejected and failed fits for each scenario
and method. Compare approximation agreement on eligible pairs and show
coverage denominators. Four replicates per cell can expose substantial
approximation errors. They have limited power to assess calibration or
rule out practically important errors. Preserve the unsuccessful attempts
when interpreting the eligible pairs.

## Reproduce

Run from the package root with the development DLL available. The
approximations need RTMB. NUTS also needs tmbstan, rstan and posterior.

```text
Rscript tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112/run.R fit 1
Rscript tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112/run.R fit 2
Rscript tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112/run.R assess
```

Both fitting shards record SHA-256 hashes of the scripts and the fitted
package's relevant R sources. The runner writes one compact RDS per attempt.
Those records exclude raw posterior draws and sampler objects. Rerunning
skips records with matching hashes. The assessment CSVs retain fit outcomes,
individual intervals and probabilities, and matched differences.

## Results

All 36 attempted fits returned. NUTS and joint Laplace passed their numerical
checks on all 12 datasets. Marginal Laplace passed on ten. Two fits, Gaussian
annual replicate 1 and Gaussian repeated replicate 3, returned with a warning
about the marginal numerical checks and were retained as rejected. No
attempt was replaced.

| Panel | Method | Eligible pairs | Material CDF differences / assessable | CDF comparisons | Largest absolute CDF difference | Median new-year interval width ratio | Median scale-SD interval width ratio |
|---|---|---:|---:|---:|---:|---:|---:|
| Gaussian annual | Joint | 4 | 21 / 24 | 24 | 0.334 | 1.359 | 1.606 |
| Gaussian annual | Marginal | 3 | 12 / 18 | 18 | 0.138 | 0.959 | 0.886 |
| Gaussian repeated | Joint | 4 | 9 / 31 | 32 | 0.634 | 1.329 | 1.514 |
| Gaussian repeated | Marginal | 3 | 11 / 21 | 24 | 0.167 | 0.892 | 0.928 |
| NB2 repeated | Joint | 4 | 19 / 32 | 32 | 0.175 | 1.156 | 1.583 |
| NB2 repeated | Marginal | 4 | 16 / 32 | 32 | 0.099 | 0.978 | 0.970 |

Interval widths are divided by those from NUTS on the same dataset. Joint
Laplace widened the scale-SD interval in all three panels on the median
ratio, while the marginal method had median ratios closer to one. Close
interval widths still coexisted with material probability differences.
For scale-SD interval coverage, the joint method covered 2/4 cases in each
panel. The NUTS comparison covered 4/4, 3/4 and 2/4 in the same panels.
On the marginal eligible pairs, coverage was 3/3, 2/3 and 2/4 for both
methods. These denominators are small.

Future negative log-score differences varied in direction across datasets.
The largest eligible absolute difference was 1.033 nats per three-observation
block for the joint method and 0.429 for the marginal method. Both
approximations were faster on these fixtures, with median elapsed times
1.9–3.4 seconds per panel and method versus 17.3–64.3 seconds for NUTS
under concurrent workloads.

The approximations are useful for quickly exploring posterior geometry
and sensitivity. Neither qualifies as a general substitute for sampled
posterior intervals, probability estimates or predictions in the designs
tested here. A passed numerical check and agreement on one summary
cannot establish approximation accuracy for another estimand.

`verify.R` checked the recorded source hashes, the 36 fit identities,
176 paired probability rows and 96 paired interval rows. The assessment
script emitted a benign data-frame row-name recycling warning when adding
case labels to multiple rows. The CSV audit verified the attached identities
and recomputed the reported comparison arithmetic.

## Supplemental slice sampling comparison

The package already verifies a shared-scale Gaussian slice posterior against
quadrature and tests a binomial coefficient update against a normalized
reference. This supplemental run checks the full hierarchical slice sampler
against NUTS on Gaussian annual replicates 1 and 2 and NB2 repeated replicates
1 and 2. It reuses those four NUTS records. The slice fits use independent
seeds, four chains of 3,000 iterations, and half warmup. Both samplers must
pass their own raw and derived diagnostics.

All four slice fits passed, with 28 assessable matched scale probabilities.
None exceeded both 0.05 and five combined MCSEs. The largest absolute
difference was 0.030 for Gaussian and 0.0268 for NB2. The runs took
58.4 to 69.0 seconds under concurrent workloads. This supports agreement
on the four tested datasets and thresholds. The package documents slice
sampling as a full-posterior method, and its status stays unchanged.

```text
Rscript tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112/run_slice.R fit
Rscript tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112/run_slice.R assess
```

The `slice_agreement/` directory contains all four compact fit records and the
two assessment CSVs.

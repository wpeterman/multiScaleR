# Conditional-choice bandwidth recovery study

## Prespecified design

Simulate 384 independent one-case conditional-choice datasets. Gaussian and
exponential kernels cross 100 and 450 training choice sets, four habitat
signals and 24 outcome replicates. Each choice set has four alternatives.
The independent assessment sample has 160 complete sets. The simulated
buffers contain one cell per training distance bin. This keeps the
generating kernel calculation and the fitted bin-distance calculation on
the same target. The cell values mix independent alternative-level habitat
features over distance and small independent cell variation. Two movement
covariates and a known offset enter the true and fitted choice model. The habitat coefficient is standardized by the generating
standard deviation at the declared bandwidth, shared across simulations in
a kernel and signal cell.
This is a synthetic matched-set design with independent sets. It does not
simulate trajectories, available-step proposals or animal-level variation.

| Signal | Generating bandwidth | Standardized habitat coefficient |
|---|---:|---:|
| Strong interior | 250 map units | 1.2 |
| Weak interior | 250 map units | 0.25 |
| Strong near lower bound | 75 map units | 1.2 |
| Null habitat | Undefined | 0 |

The fitted bandwidth range is 60 to 800 map units. At the exact null,
the label in the run table marks a convenient generating design value.
It does not define a recoverable parameter.

## Primary calculations

Call public `cross_validate_clogit()` with the conditional-logit backend
and its existing six deterministic starts. Preserve all starts and eligibility
decisions. Compare the returned objective, kernel bandwidth, coefficients and
within-set assessment probabilities with independent calculations. Recalculate
the held-out softmax with the reported coefficients and bandwidth to verify
prediction even when the public optimum differs from the reference optimum. Report
complete primary and scored denominators.

Build a reference using the exact finite-choice-set likelihood and analytic
score in base R. For every point of a fixed 121-point log bandwidth grid,
estimate all three response coefficients to the same numerical tolerance.
Refine every sampled local minimum inside its neighboring grid interval and
compare with both bounds. Check the constrained gradient and retain each
grid evaluation and refinement. The best independently evaluated value is a
numerical reference for this one-dimensional search. It is not a mathematical
proof that no narrower feature exists between grid points.

Report both bandwidth error and selected-boundary frequency. Describe the
grid span within 1.92 log-likelihood units of its best point, with an explicit
count of disconnected supported regions. This span describes flatness and is
not a calibrated confidence interval. For fitted-scale choice scores,
compare held-out loss against three predeclared references using the same
sets: the model at the true bandwidth with coefficients learned on analysis
sets, the known-generating probability, and the training-fitted model without
habitat. At the exact null, the no-habitat fit is the fixed-truth comparator.

If the public optimizer fails to return or the independently profiled
objective is better by more than 0.0001 total nats, retain that primary
result. In a separate diagnostic, add explicit starts at the bandwidth
bounds and four interior values. The diagnostic does not replace primary
estimates, alter numerical eligibility or change package defaults.

## Independent boundary-bracketing refinement

The first execution used interior sampled minima and the two endpoints as
refinement candidates. An intermediate numerical self-check found that a
sampled endpoint can be the best grid point even when the analytic score
points toward a narrow adjacent interior minimum. Preserve that original
execution. Repeat the same frozen design after adding a deterministic
analytic-score root search for every adjacent grid interval whose score
changes from negative to positive. This supplements interior grid minima
and both boundaries. It changes neither the simulated outcomes nor the
package optimizer, starts or gates. Compare the package outputs across
the two executions before using the corrected independent reference.
